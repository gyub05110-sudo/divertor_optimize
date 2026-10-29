"""Rerank a small search shortlist on an independent 640-tube grid."""
import argparse
from pathlib import Path
import numpy as np
from d3d_optimization.boundary_evaluator import prepare_boundary
from d3d_optimization.tangent_evaluator import constrained_evaluator
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.broadening_study import save,load_json,digest
from optimize_boundary_transport import verify
from optimize_regional_transport import write_arrays

ROOT=Path(__file__).resolve().parent


def run(folder,limit=5):
    verify(folder);study=load_json(folder/'study.json');ref=ReferencePoint(**study['reference'])
    fine=prepare_boundary(ROOT,'OMP',ref,count=640,step=.00025,nodes=321,wall_step=.0005)
    search=prepare_boundary(ROOT,'OMP',ref)
    inherited=[]
    for mm in [100,50,20]:
        sub=folder/f'R{mm:03d}';m=load_json(sub/'manifest.json')
        if 'coefficients_m' not in m:continue
        if 'qualification' in m:raise FileExistsError('candidate qualification already recorded')
        e=constrained_evaluator(fine,mm/1000,'gpu');chosen=[];seen=set()
        for row in m['search']['feasible_pool']:
            c=np.asarray(row['coefficients_m']);key=tuple(c)
            if key in seen:continue
            seen.add(key);chosen.append((c,'own_search'))
            if len(chosen)>=limit:break
        for c in inherited:
            if tuple(c) not in seen:chosen.append((c,'stricter_radius_selected'))
        rows=[];best=None
        for index,(c,origin) in enumerate(chosen):
            r=e.evaluate(c);row=dict(rank=index+1,origin=origin,coefficients_m=c.tolist(),**e.describe(r));rows.append(row)
            if r['feasible'] and (best is None or r['peak_q']<best[0]):best=(r['peak_q'],c)
            print('SHORTLIST',mm,index+1,r['feasible'],row.get('peak_MW_m2'),flush=True)
        if best is None:raise RuntimeError('All shortlisted candidates fail refined grid; search results retained')
        c=best[1];coarse=constrained_evaluator(search,mm/1000,'gpu');r=coarse.evaluate(c)
        if not r['feasible']:raise RuntimeError('Chosen refined candidate failed original grid')
        qualification=dict(tubes=640,nodes=321,wall_step_m=.0005,own_shortlist_limit=limit,rows=rows,
            qualifier_sha256=digest(Path(__file__)),selected_peak_MW_m2=best[0]*ref.power_W/1e6)
        save(sub/'qualification.json',qualification)
        m['search_selected_coefficients_m']=m['coefficients_m'];m['search_selected']=m['selected']
        m['coefficients_m']=c.tolist();m['selected']=coarse.describe(r);m['qualification']=qualification
        # Recheck backend parity for the final selected geometry too.
        cpu=constrained_evaluator(search,mm/1000,'cpu');rc=cpu.evaluate(c)
        parity=abs(r['peak_q']/rc['peak_q']-1)
        if parity>1e-10 or not rc['feasible']:raise ValueError('qualified CPU/GPU parity failed')
        m['cpu_gpu_peak_relative_difference']=parity
        save(sub/'manifest.json',m);write_arrays(sub,coarse,r,'selected');inherited.append(c.copy())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--limit',type=int,default=5)
    a=p.parse_args()
    if a.limit<1:p.error('positive shortlist limit required')
    run(a.output,a.limit)
