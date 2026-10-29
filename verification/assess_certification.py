"""Assess evidence for experiments/global optimality without changing the design."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.optimize import brentq
from d3d_optimization.case import load_case, separatrix, sha256
import pyarrow.parquet as pq
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.shape import NormalDisplacementBSpline

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'input/D3D/optimized/independent-production'
OUT=RUN/'certification_assessment_20261024'
OUT.mkdir(exist_ok=True)

def save(name,data):
    (OUT/name).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')

def boundaries():
    review=json.loads((ROOT/'baseline_review/baseline_review.json').read_text(encoding='utf-8'))
    shot=ROOT/'input/D3D/d3d/d3d_shot_1a3155370c.parquet'
    f,lcfs,axis=load_case(shot,1400.)
    _,xp,psi=separatrix(f,lcfs,axis)
    a,b=np.asarray(review['wall'])[69:71]
    def residual(t):return float(f.values([a+t*(b-a)])[0]-psi)
    ts=np.linspace(0,1,1001);values=[residual(t) for t in ts]
    roots=[brentq(residual,lo,hi,xtol=1e-14) for lo,hi,u,v in zip(ts[:-1],ts[1:],values[:-1],values[1:]) if u*v<0]
    if len(roots)!=1:raise ValueError(('ambiguous segment 69 crossing',roots))
    strike=a+roots[0]*(b-a)
    selected=next(r for r in review['optimization_footprints'] if r['segment']==69)
    endpoint=min(np.asarray(selected['points']),key=lambda p:np.linalg.norm(p-strike))
    endpoint_constraints=[]
    for c in np.eye(4):
        basis=NormalDisplacementBSpline(c)
        endpoint_constraints.append(dict(displacement=basis.displacement([0,1]).tolist(),slope=basis.first_derivative([0,1]).tolist()))
    fields=pq.ParquetFile(shot).schema_arrow.names
    data=dict(shot_path=str(shot),shot_sha256=sha256(shot),time_ms=1400.,available_columns=fields,
              divertor_surface_heat_measurement_columns=[],experimentally_validated=False,
              xpoint_RZ_m=xp.tolist(),continuous_strike_RZ_m=strike.tolist(),sampled_design_endpoint_RZ_m=endpoint.tolist(),
              endpoint_gap_m=float(np.linalg.norm(endpoint-strike)),flux_residual=float(residual(roots[0])),
              basis_endpoint_constraints=endpoint_constraints,
              implication='Moving the endpoint to the separatrix strike alone still fixes both its position and tangent. A peak reduction is not guaranteed.',
              design_policy_changed=False)
    save('boundary_and_data.json',data)
    print('IMP continuous strike',strike,'gap_m',data['endpoint_gap_m'],flush=True)

def challenge():
    source='OMP'
    m=json.loads((RUN/source/'manifest.json').read_text(encoding='utf-8'))
    c=np.asarray(m['verified_coefficients_m']);neighbor=c.copy();neighbor[14]+=.0001
    checks=[]
    for count,step,nodes in [(320,.002,161),(640,.001,321),(1280,.0005,641),(2560,.00025,1281)]:
        e,meta=prepare_evaluator(ROOT,source,count,step,nodes,bound=m['bound_m'])
        results=[]
        for name,coef in [('saved',c),('neighbor',neighbor)]:
            r=e.evaluate(coef)
            row=dict(name=name,geometry_valid=r['geometry_valid'],feasible=r['feasible'])
            if r['geometry_valid']:
                ratios=[float(max(h.q)/b) for h,b in zip(r['maps'],e.baseline_peaks)]
                row.update(worst_ratio=max(ratios),ratios=ratios,unmapped_power=r['missing'],accounting_error=r['power_balance_error'],
                           shelf_min_angle_deg=r['angles']['shelf']['minimum_angle_deg'])
            results.append(row)
        checks.append(dict(tubes=count,step_m=step,nodes=nodes,results=results,
                           better_feasible_neighbor=bool(results[1]['feasible'] and results[1]['worst_ratio']<results[0]['worst_ratio'])))
        save('optimality_challenge.json',dict(source=source,saved_coefficients_m=c.tolist(),neighbor_coefficients_m=neighbor.tolist(),
             changed_coefficient_index=14,change_m=.0001,checks=checks,rigorous_interval_certificate=False,
             original_results_replaced=False))
        print('OMP',count,[(r['name'],r['feasible'],r.get('worst_ratio')) for r in results],flush=True)

if __name__=='__main__':
    boundaries()
    challenge()
