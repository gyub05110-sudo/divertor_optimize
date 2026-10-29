"""Compare plate geometries at a fixed hypothetical operating point."""
import argparse
from dataclasses import asdict
from pathlib import Path
import json
import time
import numpy as np
from scipy.optimize import differential_evolution
from scipy.stats import qmc
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.transport_evaluator import TransportEvaluator,summarize
from d3d_optimization.reference_study import capture,export
from d3d_optimization.broadening_study import save,load_json
from compare_broadening import HISTORICAL,select_warm_start

ROOT=Path(__file__).resolve().parent
MODES=('direct','frozen_transport','thermal')


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=['preflight','optimize','verify'],default='optimize')
    p.add_argument('--source',choices=['IMP','OMP','both'],default='both')
    p.add_argument('--power-mw',type=float,default=1.)
    p.add_argument('--density',type=float,default=1e19)
    p.add_argument('--chi',type=float,default=.3)
    p.add_argument('--tubes',type=int,default=320)
    p.add_argument('--nodes',type=int,default=161)
    p.add_argument('--step-mm',type=float,default=.5)
    p.add_argument('--wall-step-mm',type=float,default=1.)
    p.add_argument('--quadrature',type=int,default=128)
    p.add_argument('--sobol-power',type=int,default=7)
    p.add_argument('--population',type=int,default=24)
    p.add_argument('--generations',type=int,default=8)
    p.add_argument('--seed',type=int,default=20261028)
    p.add_argument('--bound-mm',type=float,default=100.)
    p.add_argument('--smoke',action='store_true')
    a=p.parse_args(argv)
    if (not all(np.isfinite(x) and x>0 for x in [a.power_mw,a.density,a.chi,a.step_mm,a.wall_step_mm,a.bound_mm])
        or a.tubes<8 or a.nodes<9 or a.quadrature<8 or a.population<5 or a.generations<0
        or not 0<=a.sobol_power<=20 or not 0<=a.seed<2**32):
        p.error('finite positive controls and valid counts required')
    if a.smoke:a.tubes,a.nodes,a.sobol_power,a.population,a.generations=48,41,3,8,1
    return a


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


def run_source(args,source):
    folder=args.output.resolve()/source
    if folder.exists() and any(folder.iterdir()):raise FileExistsError(f'existing study: {folder}')
    warm_path=HISTORICAL[source]
    capture(ROOT,folder,[warm_path])
    ref=ReferencePoint(power_W=args.power_mw*1e6,density_m3=args.density,chi_perp_m2_s=args.chi)
    start=time.monotonic()
    direct,meta=prepare_evaluator(ROOT,source,count=args.tubes,nodes=args.nodes,
        step=args.step_mm/1000,bound=args.bound_mm/1000,boundary_mode='continuous')
    historical=load_json(ROOT/warm_path)
    if any(historical[k]!=meta[k] for k in ['shot_sha256','geometry_sha256','surfaces']):
        raise ValueError('warm start geometry/input mismatch')
    warm,status=select_warm_start(historical,direct.dimension,direct.bound)
    warm=[] if warm is None else [warm]
    manifest=dict(schema='reference-transport-study.v1',source=source,reference=asdict(ref),
        mode=args.mode,smoke_test=args.smoke,metadata=meta,models={},
        candidates={'baseline':[0.]*direct.dimension},physical_validation=False,
        global_optimum_proven=False,validation_status='not_run',historical_warm_start_status=status,
        physical_scope='one-way effective two-point closure then endpoint deposition; no feedback, particles, momentum or radiation',
        settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()})
    if warm:manifest['candidates']['historical_direct']=warm[0].tolist()
    evaluators={};cached=None
    for mode in MODES:
        e=TransportEvaluator(direct,meta,ref,mode=mode,wall_step_m=args.wall_step_mm/1000,
                             thermal_paths=cached,quadrature_order=args.quadrature)
        cached=e.thermal_paths;evaluators[mode]=e
        destination=folder/mode;destination.mkdir()
        baseline=summarize(e,e.baseline)
        manifest['models'][mode]=dict(baseline=baseline)
        export(destination/'baseline.npz',e,e.baseline)
        print(source,mode,'baseline feasible',baseline['feasible'],'missing',baseline.get('missing_power'),flush=True)
        if args.mode=='optimize':
            outcome=search(e,destination,args,warm)
            manifest['models'][mode]['search']=outcome
            if outcome['best_coefficients_m'] is not None:
                manifest['candidates'][mode+'_best']=outcome['best_coefficients_m']
        save(folder/'manifest.json',manifest)
    cross={}
    for mode,e in evaluators.items():
        cross[mode]={}
        for name,c in manifest['candidates'].items():
            result=e.baseline if name=='baseline' else e.evaluate(c)
            row=summarize(e,result)
            if result['geometry_valid']:
                row['max_wall_displacement_m']=float(np.max(np.linalg.norm(result['wall']-direct.baseline['wall'],axis=1)))
            cross[mode][name]=row
            export(folder/mode/(name+'.npz'),e,result)
    manifest['cross_evaluation']=cross;manifest['elapsed_seconds']=time.monotonic()-start
    save(folder/'manifest.json',manifest)
    print(source,'reference study saved',round(manifest['elapsed_seconds'],1),'s',flush=True)
    return manifest


def main(argv=None):
    args=parse_args(argv)
    if args.mode=='verify':
        from verification.validate_reference_transport import validate_study
        return validate_study(args.output,args.source)
    for source in ['IMP','OMP'] if args.source=='both' else [args.source]:run_source(args,source)
    from verification.plot_reference_transport import write_report
    write_report(args.output)


if __name__=='__main__':main()


