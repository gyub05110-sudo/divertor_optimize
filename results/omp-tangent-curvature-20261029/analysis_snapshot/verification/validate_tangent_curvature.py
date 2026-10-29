"""Shared-grid validation of the three fixed tangent/radius candidates."""
import argparse,time
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from d3d_optimization.boundary_evaluator import prepare_boundary,BoundaryEvaluator
from d3d_optimization.tangent_evaluator import constrained_evaluator,ACTIVE
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.broadening_study import save,load_json,digest
from optimize_regional_transport import write_arrays
from optimize_boundary_transport import verify
from verification.validate_boundary_transport import LEVELS,audit,verdict


def geometry_audit(e,c,r):
    check=e.constraints(c); dense={}; u=np.linspace(0,1,50001)
    for n,curve in e.curves(c).items():
        if n not in ACTIVE:continue
        d=curve.point(u); t=np.gradient(d,u,axis=0,edge_order=2); acc=np.gradient(t,u,axis=0,edge_order=2)
        curvature=abs(t[:,0]*acc[:,1]-t[:,1]*acc[:,0])/np.linalg.norm(t,axis=1)**3
        dense[n]=float(np.max(curvature[2:-2]))
    floor=e.curves(c)['floor']; distances=np.linalg.norm(floor.point(u)-e.lip,axis=1)
    opposite=max(float(np.max(abs(r['targets'][n]-e.baseline['targets'][n]))) for n in ['imp_sloped','imp_vertical'])
    residual=max(dense[n]-check['curvature_max_1_m'][n] for n in ACTIVE)
    return dict(**check,dense_finite_difference_max_curvature_1_m=dense,
        dense_vs_analytic_max_excess_1_m=residual,dense_gap_m=float(distances.min()),opposite_error_m=opposite,
        independent_geometry_passed=bool(check['passed'] and residual<.01 and
            distances.min()>=check['gap_m']-1e-10 and distances.min()-check['gap_m']<1e-6 and opposite<1e-12))


def run(folder):
    verify(folder);study=load_json(folder/'study.json'); ref=ReferencePoint(**study['reference'])
    records={mm:dict(complete=False,levels=[],validator_sha256=digest(Path(__file__))) for mm in [20,50,100]}
    manifests={mm:load_json(folder/f'R{mm:03d}'/'manifest.json') for mm in records}
    for mm in list(records):
        if 'coefficients_m' not in manifests[mm]:
            save(folder/f'R{mm:03d}'/'validation.json',dict(complete=True,status='no_candidate_to_validate',
                search_status=manifests[mm]['search'].get('status'),verdict=dict(resolution_screen_passed=False,improvement_confirmed=False)))
            del records[mm]
    if not records:return
    for index,(method,count,step,nodes,wall_step) in enumerate(LEVELS):
        start=time.perf_counter()
        original=prepare_boundary(ROOT,'OMP',ref,count=count,step=step,nodes=nodes,wall_step=wall_step,method=method)
        base=audit(original,original.baseline)
        if index==2:write_arrays(folder,original,original.baseline,'baseline_refined')
        if index==0:
            fine_original=BoundaryEvaluator(original.direct,original.metadata,ref,wall_step/2,thermal_paths=original.thermal_paths)
            fine_base=audit(fine_original,fine_original.baseline)
        for mm,record in records.items():
            e=constrained_evaluator(original,mm/1000,'gpu'); c=np.array(manifests[mm]['coefficients_m']);r=e.evaluate(c)
            selected=audit(e,r);geom=geometry_audit(e,c,r)
            selected['audit_passed']=bool(selected['audit_passed'] and geom['independent_geometry_passed'])
            rows=dict(baseline=base,selected=selected)
            record['levels'].append(dict(method=method,tubes=count,trace_step_m=step,nodes=nodes,
                wall_step_m=wall_step,rows=rows,geometry_audit=geom,elapsed_seconds=time.perf_counter()-start))
            if index==0:
                fine=constrained_evaluator(fine_original,mm/1000,'gpu');fr=fine.evaluate(c);fs=audit(fine,fr)
                record['wall_only_refinement']=dict(wall_step_m=wall_step/2,rows=dict(baseline=fine_base,selected=fs),
                    max_peak_relative_change=max(abs(fine_base['peak_q']/base['peak_q']-1),abs(fs['peak_q']/selected['peak_q']-1)))
            if index==2:write_arrays(folder/f'R{mm:03d}',e,r,'selected_refined')
            save(folder/f'R{mm:03d}'/'validation.json',record)
            print('VALIDATE',mm,method,count,'audit',base['audit_passed'],selected['audit_passed'],
                'peak',selected.get('peak_MW_m2'),'gap',geom['gap_m'],'radius',geom['actual_minimum_radius_m']['floor'],flush=True)
    for mm,record in records.items():
        record['verdict']=verdict(record['levels'],record['wall_only_refinement']);record['complete']=True
        save(folder/f'R{mm:03d}'/'validation.json',record);print('VERDICT',mm,record['verdict'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.output)
