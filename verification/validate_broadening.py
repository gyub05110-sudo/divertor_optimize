"""Independent first-hit, quadrature and refinement screens for the surrogate."""
import argparse
from dataclasses import replace
from pathlib import Path
import sys
import time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.integrate import quad
from shapely.geometry import LineString, Polygon
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.frozen_mapping import first_hits
from d3d_optimization.broadening import BroadeningConfig, PathLengths
from d3d_optimization.broadened_evaluator import BroadenedEvaluator
from d3d_optimization.broadening_study import load_json, save, verify_provenance, summarize, acceptance, gate_original_screen

LEVELS = [('cubic',320,.002,161,.002), ('cubic',640,.001,321,.001),
          ('cubic',1280,.0005,641,.0005), ('quintic',640,.001,321,.001)]


def independent_lengths(field, paths, hits):
    """Reintegrate B/Bp on cut polyline with midpoint Simpson quadrature."""
    output = np.full(len(paths), np.nan)
    for i, (path, distance, segment) in enumerate(zip(paths, hits.distance, hits.segment)):
        if segment < 0:
            continue
        points = np.asarray(path.coords)
        arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
        cut = np.searchsorted(arc, distance, side='left')
        p = np.vstack([points[:cut], hits.points[i]])
        ds = np.linalg.norm(np.diff(p, axis=0), axis=1)
        mids = (p[:-1]+p[1:])/2
        def metric(v):
            b = field.magnetic(v)
            return np.linalg.norm(b,axis=1)/np.linalg.norm(b[:,:2],axis=1)
        edge, middle = metric(p), metric(mids)
        output[i] = np.sum(ds*(edge[:-1]+4*middle+edge[1:])/6)
    return output


def audit(e, result, baseline_audit=True):
    row = summarize(e, result)
    if not result['geometry_valid'] or 'peak_q' not in result:
        row['audit_passed'] = False
        return row
    full = first_hits(e.paths, result['wall'])
    maps = result.get('direct_maps', result.get('maps'))
    mismatch = int(np.count_nonzero(full.segment != maps[0].segment))
    both = full.segment >= 0
    point_error = float(np.max(np.linalg.norm(full.points[both]-maps[0].points[both],axis=1))) if np.any(both) else 0.
    wall = LineString(result['wall'])
    geometry = bool(wall.is_simple and Polygon(result['wall']).is_valid
                    and not wall.intersects(Polygon(e.lcfs)) and not wall.intersects(e.vessel_solids))
    fixed_error = 0.
    for i in e.fixed_indices:
        ids = np.flatnonzero(result['original_segment'] == i)
        fixed_error = max(fixed_error, float(np.max(abs(result['wall'][ids[0]:ids[0]+2]-e.wall[i:i+2]))))
    source_errors = []
    rsep,zup = e.metadata['source_RZ_m']
    sign = -1 if e.metadata['source'] == 'IMP' else 1
    for source in e.sources:
        def density(x):
            p = [[rsep+sign*x,zup]]
            return float(np.exp(-x/source.width)*2*np.pi*abs(e.field.gradient(p)[0,0])/np.linalg.norm(e.field.magnetic(p)[0]))
        integrated = quad(density,0,source.edges[-1],epsabs=1e-13,epsrel=1e-10)[0]
        source_errors.append(abs(integrated-source.normalization)/integrated)
    extra_ok = True
    if 'deposition' in result:
        dep = result['deposition']
        quadrature_error = float(np.max(abs(dep['quadrature_integrated_power']-dep['ledger']['deposited'])))
        length_independent = independent_lengths(e.field,e.paths,full)
        good = both & np.isfinite(dep['lengths_m'])
        delta = abs(length_independent[good]-dep['lengths_m'][good])/length_independent[good]
        max_length_error = float(np.max(delta)) if np.any(good) else 1.
        drift = np.asarray(dep['diagnostics']['path_flux_drift_xi_m'])
        sigma = np.asarray(dep['sigma_m'])
        drift_ratio = float(np.max(drift[good]/sigma[good]))
        row.update(independent_density_integral=dep['quadrature_integrated_power'].tolist(),
            density_vs_cdf_power_max_error=quadrature_error,
            max_independent_length_relative_error=max_length_error,
            max_path_flux_drift_over_sigma=drift_ratio)
        extra_ok = quadrature_error <= 1e-4 and max_length_error <= .005 and drift_ratio <= .05
    row.update(full_first_hit_segment_mismatches=mismatch, full_first_hit_point_error_m=point_error,
        full_wall_geometry_passed=geometry, fixed_wall_coordinate_error_m=fixed_error,
        source_normalization_relative_errors=source_errors,
        audit_passed=bool(baseline_audit and e.baseline_admissible and extra_ok and mismatch == 0
                          and point_error < 1e-10 and geometry and fixed_error < 1e-12
                          and max(source_errors) < 1e-8))
    return row


def validate_source(folder, levels=LEVELS):
    folder = Path(folder)
    if (folder/'validation.json').exists() or any(folder.glob('*/*_refined.npz')):
        raise FileExistsError(f'refusing to overwrite completed or partial validation: {folder}')
    verify_provenance(ROOT,folder)
    manifest = load_json(folder/'manifest.json')
    if manifest.get('schema') != 'connection-broadening-study.v1':
        raise ValueError('not a broadening study')
    from compare_broadening import export_arrays
    result = dict(schema='connection-broadening-validation.v1', complete=False, levels=[],
                  comparisons={}, anchor_sensitivity={}, exponent_sensitivity={})
    for index, (method,count,step,nodes,wall_step) in enumerate(levels):
        start = time.monotonic()
        d, meta = prepare_evaluator(ROOT, manifest['source'],count=count,step=step,nodes=nodes,
            method=method,bound=manifest['settings']['bound_mm']/1000,boundary_mode='continuous')
        if any(meta[k] != manifest['metadata'][k] for k in ['shot_sha256','geometry_sha256','surfaces']):
            raise ValueError('validation inputs differ from study')
        lengths = PathLengths(d.field,d.paths)
        level = dict(interpolation=method,tubes=count,step_m=step,target_nodes=nodes,wall_step_m=wall_step,models={})
        for model, saved in manifest['models'].items():
            config = replace(BroadeningConfig(**saved['config']), wall_step_m=wall_step)
            e = BroadenedEvaluator(d,meta,config,lengths=lengths)
            base_row = audit(e,e.baseline)
            rows = {'baseline': base_row}
            for name,c in manifest['candidates'].items():
                if name == 'baseline':
                    continue
                evaluated = e.evaluate(c)
                rows[name] = audit(e,evaluated,base_row['audit_passed'])
                if index == 2:
                    export_arrays(folder/model/(name+'_refined.npz'),e,evaluated)
            if index == 2:
                export_arrays(folder/model/'baseline_refined.npz',e,e.baseline)
            level['models'][model] = rows
            print(manifest['source'],method,count,model,'baseline audit',base_row['audit_passed'],
                  'q integral error',base_row.get('density_vs_cdf_power_max_error'),flush=True)
            # Closure-sensitivity checks are separate from numerical convergence.
            if index == 1 and config.sigma_ref_m > 0 and model+'_best' in manifest['candidates']:
                coeff = manifest['candidates'][model+'_best']
                anchor_rows = {}
                for anchor in [.0025,.01]:
                    other = BroadenedEvaluator(d,meta,replace(config,anchor_m=anchor),lengths=lengths)
                    anchor_rows[str(anchor)] = summarize(other,other.evaluate(coeff))
                result['anchor_sensitivity'][model] = anchor_rows
                other = BroadenedEvaluator(d,meta,replace(config,exponent=.5 if config.exponent==1 else 1),lengths=lengths)
                result['exponent_sensitivity'][model] = summarize(other,other.evaluate(coeff))
        level['elapsed_seconds'] = time.monotonic()-start
        result['levels'].append(level)
        save(folder/'validation.json',result)
    for model in manifest['models']:
        result['comparisons'][model] = {}
        for candidate in manifest['candidates']:
            evidence = [level['models'][model][candidate] for level in result['levels']]
            verdict = gate_original_screen(acceptance(evidence),
                manifest['models'][model]['baseline'], manifest['cross_evaluation'][model][candidate])
            if manifest['smoke_test']:
                verdict['improvement_confirmed'] = False
                verdict['smoke_test'] = True
            verdict['anchor_sensitive'] = None
            if candidate == model+'_best' and model in result['anchor_sensitivity']:
                central = evidence[1].get('peak_ratios')
                variants = list(result['anchor_sensitivity'][model].values())
                valid = central is not None and all(r['feasible'] and r['baseline_admissible'] and 'peak_ratios' in r for r in variants)
                spread = max(float(np.max(abs(np.asarray(r['peak_ratios'])-central))) for r in variants) if valid else None
                verdict['max_anchor_ratio_change'] = spread
                verdict['anchor_sensitive'] = not valid or spread > .002
                verdict['robust_to_anchor_choice'] = not verdict['anchor_sensitive']
            result['comparisons'][model][candidate] = verdict
    result['complete'] = True
    result['physical_validation'] = False
    save(folder/'validation.json',result)
    manifest['validation_status'] = 'complete_see_per_candidate_verdicts'
    save(folder/'manifest.json',manifest)
    return result


def validate_study(output, source='both'):
    summaries = {}
    for name in ['IMP','OMP'] if source=='both' else [source]:
        summaries[name] = validate_source(Path(output)/name)['comparisons']
    save(Path(output)/'validation_summary.json',summaries)
    from verification.plot_broadening import write_report
    write_report(output)
    return summaries


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source',choices=['IMP','OMP','both'],default='both')
    args = parser.parse_args()
    validate_study(args.output,args.source)
