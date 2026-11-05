"""Fixed-input deposition comparison on actual DIII-D first-hit branches; all lengths are SI."""
import argparse
from pathlib import Path
from time import perf_counter
from functools import partial
import numpy as np
import shapely
from shapely.geometry import LineString,Polygon,Point
from scipy.optimize import brentq
from d3d_optimization.case import load_case
from d3d_optimization.frozen_mapping import trace_paths,first_hits
from d3d_optimization.broadening_study import load_json,save,digest
from verification.preimpact_routing import Wall,GaussianMass,ring_grid,route_bins,split_routes,scores,complete_bounds

ROOT=Path(__file__).resolve().parents[1]
METHODS=((.002,0.),(.001,0.),(.004,0.),(.002,.5))
CUTS=(-1.285,-1.305)


class Flux:
    def __init__(self,field,rz):
        self.field=field
        self.sep=field.values([rz])[0]
        self.grad=field.gradient([rz])[0,0]

    def values(self,p):return (self.field.values(p)-self.sep)/self.grad
    def gradient(self,p):return self.field.gradient(p)/self.grad


class Section:
    def __init__(self,z,lo,hi,flux):
        self.z,self.lo,self.hi,self.flux=z,lo,hi,flux
        u=flux.values(np.c_[np.linspace(lo,hi,257),np.full(257,z)])
        if not (np.diff(u)>0).all():raise ValueError('launch branch must be monotone')
        self.umin,self.umax=u[0],u[-1]

    def starts(self,u):
        radii=[brentq(lambda r:self.flux.values([[r,self.z]])[0]-v,self.lo,self.hi,xtol=1e-14) for v in u]
        return np.c_[radii,np.full(len(radii),self.z)]


def private_section(z,wall,xp):
    # Outboard side of the negative-label private-flux branch on this horizontal plane.
    flux=wall.flux
    rcrit=brentq(lambda r:flux.gradient([[r,z]])[0,0],xp[0]-.12,xp[0]+.10)
    outer_sep=brentq(lambda r:flux.values([[r,z]])[0],rcrit,xp[0]+.19)
    segments=Polygon(wall.p).intersection(LineString([(flux.field.r[0],z),(flux.field.r[-1],z)]))
    parts=segments.geoms if hasattr(segments,'geoms') else [segments]
    part,=[q for q in parts if q.bounds[0]<outer_sep<q.bounds[2]]
    return Section(z,max(rcrit,part.bounds[0]),outer_sep,flux)


def traced_arcs(labels,section,wall,step):
    ids=[]
    for part in np.array_split(labels,max(1,int(np.ceil(len(labels)/80)))):
        starts=section.starts(part)
        arcs=np.full(len(part),-1)
        inside=shapely.contains_xy(Polygon(wall.p),starts[:,0],starts[:,1])
        if inside.any():
            # The original solver's direction, complete path, all wall intersections, earliest hit.
            paths=trace_paths(wall.flux.field,starts[inside],step=step,max_length=8.)
            hits=first_hits(paths,wall.p)
            ok=hits.segment>=0
            s=wall.hit_s(hits.points[ok],hits.segment[ok])
            arcs[np.flatnonzero(inside)[ok]]=wall.arc_id(s)
        ids.extend(arcs)
    return np.asarray(ids,int)


def ranges(wall,section,low,high,step,dense):
    # Numerical grouping preserves both endpoints and never changes source normalization.
    linear=np.linspace(low,high,33 if dense else 17)
    cuts=np.sort(np.r_[low,high,linear,wall.level[(wall.level>low)&(wall.level<high)]])
    cuts=np.r_[low,cuts[(cuts>low+1e-12)&(cuts<high-1e-12)],high]
    cuts=cuts[np.r_[True,np.diff(cuts)>1e-12]]
    samples=(.007,.30,.68,.993) if dense else (.002,.5,.998)
    trace=partial(traced_arcs,section=section,wall=wall,step=step)
    good,bad=split_routes(cuts,trace,depth=5,fractions=samples)
    accepted=[]
    for lo,hi,j in good:
        v0,v1=sorted(wall.level[j:j+2]);a=max(lo,v0);b=min(hi,v1)
        if b<=a:bad.append((lo,hi));continue
        accepted.append((a,b,j))
        if lo<a:bad.append((lo,a))
        if b<hi:bad.append((b,hi))
    return accepted,bad


def old_power(wall,archive,grid,spacing,phase):
    # Integrate archived native cell powers in physical area onto the identical new ring edges.
    source_area=archive['area_edges'];source_power=archive['cell_power']
    if abs(source_area[-1]-wall.area[-1])>1e-8:raise ValueError('archived material area mismatch')
    n=len(grid['area']);s=wall.s[-1]*((np.arange(n+1)+phase)/n)
    cycles=np.floor((np.arange(n+1)+phase)/n)
    u=wall.area_at(np.clip(s-cycles*wall.s[-1],0,wall.s[-1]))
    sums=cycles*source_power.sum()+np.interp(u,source_area,np.r_[0.,np.cumsum(source_power)])
    return np.diff(sums)


def case_run(folder,case,raw,field,flux,table,omp,xp,quality):
    start=perf_counter();out=folder/(case+'_'+quality);out.mkdir(exist_ok=False)
    wall=Wall(raw['wall'],flux)
    u=flux.values(np.c_[table['R_m'],table['Z_m']]);weights=table['power_fraction']
    np.testing.assert_allclose(raw['source_x_m'],table['x_mid_m'],rtol=0,atol=1e-14)
    np.testing.assert_allclose(raw['target_gradient'],np.linalg.norm(field.gradient(raw['hit_points']),axis=1),rtol=1e-10)
    assert np.all(raw['hit_segment']>=0) and np.isclose(table['power_W'].sum(),1e6,rtol=1e-13)
    mass=GaussianMass(u,raw['thermal_sigma_m'],weights)
    fine=quality=='fine';step=.0005 if fine else .001
    good,bad=ranges(wall,omp,0.,omp.umax,step,fine)
    # Check all independent archived donor rays against the new numerical routing, and record
    # any discrepant interval as unresolved (including power from all displaced donors).
    hit_s=wall.hit_s(raw['hit_points'],raw['hit_segment']);arc=wall.arc_id(hit_s)
    accepted=[];mismatched=[];point_errors=[];matched=0
    for l,h,j in good:
        ids=np.flatnonzero((u>=l)&(u<h))
        if np.any(arc[ids]!=j):
            mismatched.append((l,h));continue
        accepted.append((l,h,j))
        for i in ids:
            fn=lambda s:float(flux.values(wall.points_at([s]))[0]-u[i])
            s=brentq(fn,wall.turn[j],wall.turn[j+1],xtol=5e-13)
            point_errors.append(float(np.linalg.norm(wall.points_at([s])[0]-raw['hit_points'][i])))
        matched+=len(ids)
    bad+=mismatched;good=accepted
    route_cases={'omp':(good,[(-np.inf,0.),*bad,(omp.umax,np.inf)])}
    for z in CUTS:
        sec=private_section(z,wall,xp)
        assert abs(sec.umax)<1e-10
        ng,nb=ranges(wall,sec,sec.umin,0.,step,fine)
        route_cases[f'leg{z:g}']=(good+ng,[(-np.inf,sec.umin),*nb,*bad,(omp.umax,np.inf)])
    report=dict(case=case,quality=quality,step_m=step,donors=len(u),archived_hit_arc_matches=matched,
                mismatch_intervals=mismatched,archived_hit_max_distance_m=max(point_errors,default=None),
                p_below_sep=float(mass(-np.inf,0.)),p_above_range=float(mass(omp.umax,np.inf)),measures={})
    routes={}
    for name,(assigned,unknown) in route_cases.items():
        routes[name]=dict(assigned=assigned,unknown=[[float(l) if np.isfinite(l) else None,float(h) if np.isfinite(h) else None] for l,h in unknown])
        interval_list=sorted([(l,h) for l,h,_ in assigned]+list(unknown))
        assert np.isneginf(interval_list[0][0]) and np.isposinf(interval_list[-1][1])
        assert max(abs(x[1]-y[0]) for x,y in zip(interval_list,interval_list[1:]))<5e-12
        missing=sum(float(mass(l,h)) for l,h in unknown)
        metrics={}
        for width,phase in METHODS:
            key=f'{width*1000:g}mm_{phase:g}'
            grid=ring_grid(wall,width,phase);area=grid['area']
            old=old_power(wall,raw,grid,width,phase)
            p,cap=route_bins(wall,grid,assigned,unknown,mass)
            if abs(old.sum()-1)>3e-10 or abs(p.sum()+missing-1)>3e-9:
                raise ValueError(f'wall ledger mismatch {case} {name} {p.sum()} {missing}')
            m=scores(p,area)
            metrics[key]=dict(old=scores(old,area),routed=m,resolved=float(p.sum()),unknown=missing,
                              qmax_RZ=wall.points_at([grid['s_mid'][np.argmax(p/area)]])[0].tolist(),
                              arbitrary=complete_bounds(p,area,missing),flux_conserving=complete_bounds(p,area,missing,cap))
            np.savez_compressed(out/(name+'_'+key+'.npz'),s=grid['s_mid'],rz=wall.points_at(grid['s_mid']),
                                area=area,original_MW=old,resolved_MW=p,unknown_bin_cap_MW=cap)
        report['measures'][name]=metrics
    save(out/'routes.json',routes,exclusive=True)
    report['seconds']=perf_counter()-start
    save(out/'result.json',report,exclusive=True)
    print(case,quality,'matched',matched,'time',round(report['seconds'],1),flush=True)


def main():
    from d3d_optimization.case import separatrix
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence',type=Path,required=True);p.add_argument('--source-csv',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--case',choices=['baseline','OMP20','OMP50'],required=True)
    p.add_argument('--quality',choices=['coarse','fine'],required=True)
    a=p.parse_args()
    shot=ROOT/'input/D3D/d3d/d3d_shot_1a3155370c.parquet';field,lcfs,axis=load_case(shot,1400)
    sep,xp,ps=separatrix(field,lcfs,axis)
    table=np.genfromtxt(a.source_csv,delimiter=',',names=True)
    rsep=table['R_m'][0]-table['x_mid_m'][0];z=table['Z_m'][0]
    assert abs(sep-rsep)<1e-11 and abs(z-axis[1])<1e-12 and np.all(table['sign_along_B']==-1)
    flux=Flux(field,(rsep,z));omp=Section(z,rsep,rsep+table['x_high_m'][-1],flux)
    entry=load_json(a.evidence)['inputs'][a.case]
    if digest(entry['path'])!=entry['sha256']:raise ValueError('archived source changed')
    with np.load(entry['path']) as raw:case_run(a.output,a.case,raw,field,flux,table,omp,xp,a.quality)

if __name__=='__main__':main()

