"""Diagnose the single-branch pre-impact spreading hypothesis using archived real donors.
No real-wall heat flux is computed or ranked with an incomplete first-hit map.
"""
import argparse
from pathlib import Path
import numpy as np
from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.case import load_case
from d3d_optimization.flux_footprint import branch_power

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence',required=True,type=Path)
    p.add_argument('--source-csv',required=True,type=Path)
    p.add_argument('--shot',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    table=np.genfromtxt(a.source_csv,names=True,delimiter=',')
    field,_,_=load_case(a.shot,1400)
    # x increases radially outward on the OMP. Translate the dimensional flux label
    # using the SAME separatrix gradient used to define thermal_sigma_m.
    rsep=table['R_m'][0]-table['x_mid_m'][0]
    z=table['Z_m'][0]
    sep=np.array([[rsep,z]])
    dpsi=float(field.gradient(sep)[0,0])
    def label(r):
        return (field.values(np.c_[r,np.full(len(r),z)])-field.values(sep)[0])/dpsi
    edges=label(rsep+np.r_[table['x_low_m'][0],table['x_high_m']])
    assert np.all(np.diff(edges)>0)
    centers=label(table['R_m'])
    weights=table['power_fraction']
    assert abs(weights.sum()-1)<1e-12
    result={}
    for name,desc in load_json(a.evidence)['inputs'].items():
        f=Path(desc['path'])
        if digest(f)!=desc['sha256']:raise ValueError('archived candidate changed')
        with np.load(f,allow_pickle=False) as raw:
            np.testing.assert_allclose(table['x_mid_m'],raw['source_x_m'],rtol=0,atol=1e-14)
            seg=raw['hit_segment']
            if np.any(seg<0):raise ValueError('missing first hits')
            wall=raw['wall'];t=wall[seg+1]-wall[seg]
            grad=field.gradient(raw['hit_points'])
            g=np.abs((grad*t).sum(axis=1))/(np.linalg.norm(t,axis=1)*np.linalg.norm(grad,axis=1))
            sigma=raw['thermal_sigma_m']
            masses=branch_power(edges[[0,-1]],centers,sigma,weights)
            result[name]=dict(power_g_below_025=float(weights[g<.25].sum()),power_g_below_05=float(weights[g<.5].sum()),
                              below_separatrix=masses['below'],above_traced_outer=masses['above'],
                              inside_label_range=float(masses['cell_power'][0]),file_sha256=desc['sha256'])
    save(a.output,dict(cases=result,scope='alternative pre-impact single displacement; outside label range is untraced, not physical loss',
         source_csv_sha256=digest(a.source_csv),shot_sha256=digest(a.shot)),exclusive=True)
    print(result)

if __name__=='__main__':main()

