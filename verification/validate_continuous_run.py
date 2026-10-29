"""Full-wall and resolution audit for the versioned continuous-boundary search."""
import argparse
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import quad
from shapely.geometry import LineString, Polygon
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.frozen_mapping import first_hits, map_heat
from d3d_optimization.case import sha256
from d3d_optimization.run_provenance import verify as verify_provenance
from optimize_independent_sources import summary, export, save

ROOT=Path(__file__).resolve().parents[1]
LEVELS=[('cubic',320,.002,161),('cubic',640,.001,321),('cubic',1280,.0005,641),
        ('cubic',2560,.00025,1281),('quintic',1280,.0005,641)]


def acceptance(levels):
    valid=len(levels)==5 and all(r['feasible'] and r.get('audit_passed',False) for r in levels)
    if not valid:return dict(numerically_verified_feasible=False,improvement_confirmed=False,
                             verification_status='constraints_or_independent_audit_failed')
    ratios=np.array([r['peak_ratios'] for r in levels])
    refinement=float(np.max(abs(ratios[3]-ratios[2])))
    interpolation=float(np.max(abs(ratios[4]-ratios[2])))
    worst=float(np.max(ratios[3])); observed=max(refinement,interpolation)
    converged=observed<=.002
    nonworsening=bool(np.all(ratios<=1+1e-10))
    return dict(numerically_verified_feasible=True,resolution_screen_passed=bool(converged),
                nonworsening=nonworsening,improvement_confirmed=bool(converged and nonworsening and 1-worst>2*observed+.0001),
                finest_worst_peak_ratio=worst,finest_reduction_percent=100*(1-worst),
                max_fine_grid_ratio_change=refinement,max_interpolant_ratio_change=interpolation,
                measured_uncertainty_indicator=observed,
                verification_status='passed' if converged and nonworsening else 'resolution_or_nonworsening_failed',
                global_optimum_proven=False,experimental_validation=False)


def audit_level(e, meta, coefficients):
    r=e.evaluate(coefficients); row=summary(e,r)
    row.update(tubes=meta['tubes'],step_m=meta['step_m'],target_nodes=meta['target_nodes'],interpolation=meta['interpolation'])
    if not r['geometry_valid']:
        row['audit_passed']=False
        return row,r
    full=first_hits(e.paths,r['wall'])
    maps=[map_heat(e.field,full,s.weights,s.h,r['wall']) for s in e.sources]
    mismatches=int(np.count_nonzero(full.segment!=r['maps'][0].segment))
    point_delta=float(np.nanmax(np.linalg.norm(full.points-r['maps'][0].points,axis=1)))
    q_delta=float(max(np.max(abs(a.q-b.q)) for a,b in zip(maps,r['maps'])))
    fixed_delta=0.
    for i in e.fixed_indices:
        ids=np.flatnonzero(r['original_segment']==i)
        fixed_delta=max(fixed_delta,float(np.max(abs(r['wall'][ids[0]:ids[0]+2]-e.wall[i:i+2]))))
    wall=LineString(r['wall'])
    geometry_ok=bool(wall.is_simple and Polygon(r['wall']).is_valid and not wall.intersects(Polygon(e.lcfs))
                     and not wall.intersects(e.vessel_solids))
    errors=[]; rsep,zup=meta['source_RZ_m']; sign=-1 if meta['source']=='IMP' else 1
    for source in e.sources:
        def density(x):
            p=[[rsep+sign*x,zup]]
            return float(np.exp(-x/source.width)*2*np.pi*abs(e.field.gradient(p)[0,0])/np.linalg.norm(e.field.magnetic(p)[0]))
        integral,_=quad(density,0,source.edges[-1],epsabs=1e-13,epsrel=1e-10,limit=100)
        errors.append(abs(integral-source.normalization)/integral)
    row.update(full_first_hit_segment_mismatches=mismatches,full_first_hit_point_max_delta_m=point_delta,
               independent_q_max_delta=q_delta,fixed_segment_max_coordinate_change_m=fixed_delta,
               full_wall_geometry_passed=geometry_ok,source_normalization_relative_errors=errors,
               source_weight_sums=[float(s.weights.sum()) for s in e.sources],
               candidate_peak_q=r['peak_q'].tolist(),baseline_peak_q=e.baseline_peaks.tolist())
    row['audit_passed']=bool(mismatches==0 and point_delta<1e-10 and q_delta<1e-10 and fixed_delta<1e-12
                            and geometry_ok and max(errors)<1e-8 and r['boundary_witness']['resolved'])
    return row,r


def validate_run(p):
    p=Path(p); manifest=json.loads((p/'manifest.json').read_text(encoding='utf-8'))
    if manifest.get('boundary_mode')!='continuous' or manifest['mode']!='optimize':
        raise ValueError('requires a continuous-boundary optimization')
    provenance=verify_provenance(ROOT,p)
    for name,digest in manifest['implementation_sha256'].items():
        if sha256(ROOT/name)!=digest:raise ValueError(f'implementation changed since search: {name}')
    coefficients=manifest['coefficients_m']; levels=[]
    for method,count,step,nodes in LEVELS:
        start=time.monotonic()
        e,meta=prepare_evaluator(ROOT,manifest['source'],count,step,nodes,method,manifest['bound_m'],
                                 boundary_mode=manifest['boundary_mode'],imp_buffer_m=manifest['imp_buffer_m'])
        if any(meta[k]!=manifest[k] for k in ['shot_sha256','geometry_sha256','policy','surfaces','objective_version']):
            raise ValueError('saved inputs or design policy changed')
        row,r=audit_level(e,meta,coefficients);row['elapsed_seconds']=time.monotonic()-start
        levels.append(row)
        if count==2560 and r['geometry_valid']:
            # Export finest cubic arrays; verification status is assigned after all checks.
            export(e,r,p,manifest)
        save(p/'continuous_validation.json',dict(levels=levels,complete=False))
        print(manifest['source'],method,count,'ratio',row.get('worst_all_wall_peak_ratio'),
              'feasible',row['feasible'],'audit',row['audit_passed'],flush=True)
        del e,r
    verdict=acceptance(levels)
    note=(f"연속 경계·2560 경로 검증: 최대 열유속 감소 {verdict.get('finest_reduction_percent',0):.3f}%. "
          +('수치 개선 확인.' if verdict['improvement_confirmed'] else '유의한 개선 미확인.')
          +f" 검증 상태: {verdict['verification_status']}. 실험 검증·전역 최적해 증명 아님.")
    manifest.update(verdict,validity_note=note,validity_audit_summary=verdict)
    if verdict['verification_status']=='passed':
        manifest['verified_coefficients_m']=coefficients
    else:
        manifest.pop('verified_coefficients_m',None)
    save(p/'manifest.json',manifest)
    viewer=json.loads((p/'viewer_result.json').read_text(encoding='utf-8'))
    viewer.update(numerically_verified_feasible=verdict['numerically_verified_feasible'],
                  validity_note=note,validity_audit_summary=verdict)
    save(p/'viewer_result.json',viewer)
    save(p/'continuous_validation.json',dict(levels=levels,complete=True,summary=verdict,
         validation_script_sha256=sha256(Path(__file__)),input_provenance=provenance,
         physical_scope='fixed equilibrium, direct mapping, no transverse transport'))
    return verdict


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source',choices=['IMP','OMP','both'],default='both')
    args=parser.parse_args()
    summary_path=args.output/'validated_summary.json'
    results=json.loads(summary_path.read_text(encoding='utf-8')) if summary_path.exists() else {}
    for s in ['IMP','OMP'] if args.source=='both' else [args.source]:
        results[s]=validate_run(args.output/s)
    save(args.output/'validated_summary.json',results)
    if not all(r['verification_status']=='passed' for r in results.values()):raise SystemExit(1)


if __name__=='__main__':main()
