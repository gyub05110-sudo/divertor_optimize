"""Matched exploratory searches for direct and L-dependent deposition models."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import differential_evolution
from scipy.stats import qmc
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.broadening import BroadeningConfig, PathLengths
from d3d_optimization.broadened_evaluator import BroadenedEvaluator
from d3d_optimization.broadening_study import save, load_json, capture_provenance, summarize

ROOT = Path(__file__).resolve().parent
HISTORICAL = {
    'IMP': 'input/D3D/optimized/continuous-boundary-20261024/IMP/manifest.json',
    'OMP': 'input/D3D/optimized/continuous-boundary-incumbent-20261024/OMP/manifest.json',
}


def model_key(sigma_mm):
    return 'direct' if sigma_mm == 0 else f'sigma_{sigma_mm:.9g}mm'


def select_warm_start(record, dimension, bound):
    c = record.get('verified_coefficients_m', record.get('coefficients_m'))
    if c is None or np.shape(c) != (dimension,) or not np.all(np.isfinite(c)):
        return None, 'missing_or_invalid_historical_coefficients'
    c = np.asarray(c, float)
    if np.max(abs(c)) > bound:
        return None, 'outside_requested_coefficient_bound'
    return c, 'included'


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mode', choices=['preflight', 'optimize', 'verify'], default='optimize')
    p.add_argument('--source', choices=['IMP', 'OMP', 'both'], default='both')
    p.add_argument('--sigma-mm', type=float, nargs='+', default=[0., .1, .25, .5])
    p.add_argument('--length-ref-m', type=float, default=10.)
    p.add_argument('--exponent', type=float, choices=[.5, 1.], default=1.)
    p.add_argument('--tubes', type=int, default=160)
    p.add_argument('--nodes', type=int, default=81)
    p.add_argument('--step-mm', type=float, default=2.)
    p.add_argument('--wall-step-mm', type=float, default=3.)
    p.add_argument('--anchor-mm', type=float, default=5.)
    p.add_argument('--sobol-power', type=int, default=7)
    p.add_argument('--population', type=int, default=24)
    p.add_argument('--generations', type=int, default=8)
    p.add_argument('--seed', type=int, default=20261027)
    p.add_argument('--bound-mm', type=float, default=100.)
    p.add_argument('--smoke', action='store_true')
    a = p.parse_args(argv)
    numeric = [*a.sigma_mm, a.length_ref_m, a.step_mm, a.wall_step_mm, a.anchor_mm, a.bound_mm]
    if (not all(np.isfinite(v) for v in numeric) or any(v < 0 for v in a.sigma_mm)
            or min(a.length_ref_m,a.step_mm,a.wall_step_mm,a.anchor_mm,a.bound_mm) <= 0
            or a.tubes < 8 or a.nodes < 9 or not 0 <= a.sobol_power <= 20
            or a.population < 5 or a.generations < 0 or not 0 <= a.seed < 2**32):
        p.error('finite physical/numerical controls and valid sample counts required')
    a.sigma_mm = sorted(set([0., *a.sigma_mm]))
    if a.smoke:
        a.tubes, a.nodes, a.sobol_power, a.population, a.generations = 48, 41, 3, 8, 1
    return a


def export_arrays(path, evaluator, result):
    if not result['geometry_valid'] or 'peak_q' not in result:
        return
    if 'deposition' in result:
        dep = result['deposition']
        points, heat = dep['points'], dep['q']
        extra = dict(receiving_wall=dep['wall'], cell_power=dep['cell_power'], areas=dep['areas'],
                     connection_length_m=dep['lengths_m'], sigma_m=dep['sigma_m'])
    else:
        points = result['maps'][0].points
        heat = np.column_stack([m.q for m in result['maps']])
        witness = result.get('boundary_witness')
        if witness and witness['resolved']:
            points = np.vstack([points, witness['point']])
            heat = np.vstack([heat, witness['q']])
        extra = {}
    np.savez_compressed(path, wall=result['wall'], points=points, q_over_P=heat,
                        peak_q=result['peak_q'], lcfs=evaluator.lcfs, **extra)


def search(e, directory, args, warm):
    """Same candidate proposal budget for each physical model; no gradient claims."""
    directory.mkdir(parents=True, exist_ok=True)
    if not e.baseline_admissible:
        return dict(status='baseline_inadmissible', best_coefficients_m=None, evaluated=0)
    start = time.monotonic()
    seen, pool = {}, []
    best = None
    stage = 'initial'
    with (directory/'evaluations.jsonl').open('x', encoding='utf-8') as log:
        def objective(c):
            nonlocal best
            c = np.asarray(c, float)
            key = tuple(np.round(c, 13))
            if key in seen:
                return seen[key]
            result = e.evaluate(c)
            row = dict(case=len(seen)+1, stage=stage, coefficients_m=c.tolist(), **summarize(e, result))
            if 'worst_peak_ratio' in row:
                value = row['worst_peak_ratio'] + result['penalty']
            else:
                value = 1e6 + result['penalty']
            if row['feasible'] and 'worst_peak_ratio' in row:
                score = row['worst_peak_ratio']
                pool.append((score, c.copy()))
                if best is None or score < best[0]:
                    best = score, c.copy()
                    save(directory/'best_checkpoint.json', row)
            row['search_objective'] = float(value)
            seen[key] = float(value)
            log.write(json.dumps(row, allow_nan=False, ensure_ascii=False)+'\n')
            if len(seen) % 16 == 0:
                log.flush()
                print(e.metadata['source'], directory.name, stage, len(seen),
                      'best', None if best is None else round(best[0], 7), flush=True)
            return value
        objective(np.zeros(e.dimension))
        for c in warm:
            objective(c)
        stage = 'sobol'
        samples = 2*qmc.Sobol(e.dimension, scramble=True, seed=args.seed).random_base2(args.sobol_power)-1
        # Include local perturbations of both baseline and historical candidate;
        # short design components make uniformly global samples mostly invalid.
        local = min(.005, .1*min(e.lengths.values()), e.bound)
        for i, sample in enumerate(samples):
            center = warm[i % len(warm)] if warm and i % 2 else np.zeros(e.dimension)
            amplitude = local if i % 4 else min(.02, e.bound)
            objective(np.clip(center + amplitude*sample, -e.bound, e.bound))
        rng = np.random.default_rng(args.seed)
        ranked = sorted(pool, key=lambda item: item[0])
        init = [c for _, c in ranked[:args.population]]
        center = best[1] if best else np.zeros(e.dimension)
        while len(init) < args.population:
            init.append(np.clip(center+rng.uniform(-local, local, e.dimension), -e.bound, e.bound))
        stage = 'differential_evolution'
        opt = differential_evolution(objective, [(-e.bound, e.bound)]*e.dimension,
            init=np.asarray(init), seed=args.seed, maxiter=args.generations,
            tol=0., atol=0., polish=False)
    result = dict(status='candidate_found' if best else 'no_feasible_candidate',
                  best_coefficients_m=None if best is None else best[1].tolist(),
                  best_worst_peak_ratio=None if best is None else best[0],
                  evaluated=len(seen), feasible_evaluations=len(pool),
                  elapsed_seconds=time.monotonic()-start, generations=int(opt.nit),
                  optimizer_message=str(opt.message), global_optimum_proven=False,
                  search_scope='bounded exploratory matched-budget search')
    save(directory/'search.json', result)
    return result


def run_source(args, source):
    folder = args.output.resolve()/source
    if folder.exists() and any(folder.iterdir()):
        raise FileExistsError(f'refusing to overwrite existing study: {folder}')
    warm_path = HISTORICAL[source]
    extra = [warm_path] if (ROOT/warm_path).is_file() else []
    capture_provenance(ROOT, folder, extra)
    start = time.monotonic()
    direct, meta = prepare_evaluator(ROOT, source, count=args.tubes, nodes=args.nodes,
        step=args.step_mm/1000, bound=args.bound_mm/1000, boundary_mode='continuous')
    lengths = PathLengths(direct.field, direct.paths)
    warm = []
    warm_status = 'no_historical_manifest'
    if extra:
        historical = load_json(ROOT/warm_path)
        if any(historical[k] != meta[k] for k in ('shot_sha256','geometry_sha256','surfaces')):
            raise ValueError('historical warm start uses different inputs or design freedoms')
        c, warm_status = select_warm_start(historical, direct.dimension, direct.bound)
        if c is not None:
            warm.append(c)
    metadata = dict(schema='connection-broadening-study.v1', source=source, mode=args.mode,
        smoke_test=args.smoke, metadata=meta, models={}, candidates={'baseline': [0.]*direct.dimension},
        settings={key: str(value) if isinstance(value, Path) else value for key,value in vars(args).items()},
        physical_scope='fixed equilibrium; phenomenological endpoint deposition; no solved transport PDE',
        physical_validation=False, global_optimum_proven=False, validation_status='not_run',
        historical_warm_start_status=warm_status)
    if warm:
        metadata['candidates']['historical_direct'] = warm[0].tolist()
    evaluators = {}
    for sigma in args.sigma_mm:
        key = model_key(sigma)
        config = BroadeningConfig(sigma/1000, length_ref_m=args.length_ref_m,
            exponent=args.exponent, wall_step_m=args.wall_step_mm/1000, anchor_m=args.anchor_mm/1000)
        e = BroadenedEvaluator(direct, meta, config, lengths=lengths)
        evaluators[key] = e
        model_dir = folder/key
        model_dir.mkdir(exist_ok=True)
        baseline = summarize(e, e.baseline)
        metadata['models'][key] = dict(config=asdict(config), baseline=baseline)
        export_arrays(model_dir/'baseline.npz', e, e.baseline)
        print(source, key, 'baseline', 'feasible', baseline['feasible'],
              'missing', baseline.get('missing_power'), 'errors', baseline.get('model_errors'), flush=True)
        if args.mode == 'optimize':
            outcome = search(e, model_dir, args, warm)
            metadata['models'][key]['search'] = outcome
            if outcome['best_coefficients_m'] is not None:
                metadata['candidates'][key+'_best'] = outcome['best_coefficients_m']
        save(folder/'manifest.json', metadata)
    cross = {}
    direct_choice = metadata['candidates'].get('direct_best', metadata['candidates'].get('historical_direct',metadata['candidates']['baseline']))
    direct_selected_wall = direct.evaluate(direct_choice)['wall']
    for model, e in evaluators.items():
        cross[model] = {}
        for name, c in metadata['candidates'].items():
            r = e.baseline if name == 'baseline' else e.evaluate(c)
            row = summarize(e, r)
            if r['geometry_valid']:
                base_wall = e.direct.baseline['wall']
                row['max_wall_displacement_m'] = float(np.max(np.linalg.norm(r['wall']-base_wall, axis=1)))
                row['max_wall_difference_from_direct_selected_m'] = float(np.max(np.linalg.norm(r['wall']-direct_selected_wall,axis=1)))
            cross[model][name] = row
            export_arrays(folder/model/(name+'.npz'), e, r)
    metadata['elapsed_seconds'] = time.monotonic()-start
    metadata['cross_evaluation'] = cross
    save(folder/'manifest.json', metadata)
    print(source, 'study finished', round(metadata['elapsed_seconds'], 2), 'seconds', flush=True)
    return metadata


def main(argv=None):
    args = parse_args(argv)
    if args.mode == 'verify':
        from verification.validate_broadening import validate_study
        validate_study(args.output, args.source)
        return
    sources = ['IMP','OMP'] if args.source == 'both' else [args.source]
    for source in sources:
        run_source(args, source)
    from verification.plot_broadening import write_report
    write_report(args.output)


if __name__ == '__main__':
    main()
