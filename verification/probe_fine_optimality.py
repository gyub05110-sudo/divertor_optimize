"""Local feasible-neighbor challenge at the audited resolution; not a global solver."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from d3d_optimization.independent_sources import prepare_evaluator

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'input/D3D/optimized/independent-production'
OUT=RUN/'certification_assessment_20261024'
m=json.loads((RUN/'OMP/manifest.json').read_text(encoding='utf-8'))
c=np.asarray(m['verified_coefficients_m'])
e,meta=prepare_evaluator(ROOT,'OMP',2560,.00025,1281,bound=m['bound_m'])

def result(coef):
    r=e.evaluate(coef)
    row=dict(feasible=r['feasible'],geometry_valid=r['geometry_valid'])
    if r['geometry_valid']:
        ratios=[float(max(h.q)/b) for h,b in zip(r['maps'],e.baseline_peaks)]
        row.update(ratios=ratios,worst_ratio=max(ratios),shelf_min_angle_deg=r['angles']['shelf']['minimum_angle_deg'],
                   unmapped_power=r['missing'])
    return row

original=result(c);best=dict(original,coefficients_m=c.tolist());records=[]
for index in [12,13,14,15]:
    for delta in [-.001,-.0001,-.00001,.00001,.0001,.001]:
        x=c.copy();x[index]+=delta
        if max(abs(x))>e.bound:continue
        row=dict(index=index,delta_m=delta,**result(x));records.append(row)
        if row['feasible'] and row['worst_ratio']<best['worst_ratio']:
            best=dict(row,coefficients_m=x.tolist())
    print('tested',index,'best',best['worst_ratio'],flush=True)
report=dict(source='OMP',tubes=2560,step_m=.00025,nodes=1281,original=original,best=best,probes=records,
            rigorous_global_certificate=False,interpretation='Finite feasible-neighbor test only; does not bound unexplored coefficient space.')
(OUT/'fine_optimality_probe.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
print('DONE',original['worst_ratio'],best['worst_ratio'],flush=True)
