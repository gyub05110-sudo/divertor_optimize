"""Exploratory baseline geometry review; no optimization or transport solve."""
from pathlib import Path
import json
import time
import hashlib
import numpy as np
import shapely
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import substring
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from d3d_optimization.case import load_case, separatrix
from d3d_optimization.frozen_mapping import trace_paths, first_hits

ROOT=Path(__file__).resolve().parent
RUN=ROOT/'input/D3D/optimized/d3d-divertor-floor-shelf'
OUT=ROOT/'baseline_review';OUT.mkdir(exist_ok=True)
SHOT=ROOT/'input/D3D/d3d/d3d_shot_1a3155370c.parquet'
manifest=json.loads((RUN/'manifest.json').read_text(encoding='utf-8'))
assert hashlib.sha256(SHOT.read_bytes()).hexdigest()==manifest['shot_sha256']
field,lcfs,axis=load_case(SHOT,1400.)
r_omp,xp,psi_sep=separatrix(field,lcfs,axis)
arrays=np.load(RUN/'case_arrays.npz'); wall=arrays['baseline_wall']
original=np.asarray(json.loads((RUN/'DIIID_geom.json').read_text(encoding='utf-8'))['limiter'])
original=np.vstack([original,original[0]])
assert np.allclose(wall,original,atol=1e-14)
polygon=Polygon(wall);boundary=LineString(wall)
psi_axis=float(field.values([axis])[0]);span=psi_sep-psi_axis
rr=np.linspace(field.r[0]+1e-5,axis[0]-.02,500)
vals=field.values(np.c_[rr,np.full(len(rr),axis[1])])-psi_sep
roots=[brentq(lambda r:float(field.values([[r,axis[1]]])[0]-psi_sep),a,b,xtol=1e-13)
       for a,b,u,v in zip(rr[:-1],rr[1:],vals[:-1],vals[1:]) if u*v<0]
if len(roots)!=1:raise ValueError(('ambiguous IMP crossing',roots))
r_imp=roots[0]

def clip(path,point):
    line=path if isinstance(path,LineString) else LineString(path)
    p=np.asarray(substring(line,0,line.project(Point(point))).coords)
    sampled=p[::max(1,int(len(p)/220))]
    if np.linalg.norm(sampled[-1]-p[-1])>1e-12:sampled=np.vstack([sampled,p[-1]])
    return sampled

def leg(direction,epsilon,spacing):
    start=xp+epsilon*direction
    # Project the local saddle ray onto the exact cubic psi_X contour.
    for _ in range(6):
        grad=field.gradient([start])[0]
        residual=field.values([start])[0]-psi_sep
        start-=residual*grad/np.dot(grad,grad)
    grad=field.gradient([start])[0]; tangent=np.array([-grad[1],grad[0]])
    sign=1. if np.dot(tangent,direction)>0 else -1.
    def rhs(s,p):
        grad=field.gradient([p])[0];t=np.array([-grad[1],grad[0]])
        return sign*t/np.linalg.norm(t)
    def bottom(s,p):return p[1]+1.445
    bottom.terminal=True;bottom.direction=-1
    def left(s,p):return p[0]-field.r[0]-.002
    left.terminal=True;left.direction=-1
    def right(s,p):return field.r[-1]-.002-p[0]
    right.terminal=True;right.direction=-1
    sol=solve_ivp(rhs,(0,3),start,rtol=2e-10,atol=1e-12,method='DOP853',max_step=spacing,
                  events=[bottom,left,right],dense_output=True)
    ss=np.r_[np.arange(0,sol.t[-1],spacing),sol.t[-1]]
    path=np.vstack([xp,sol.sol(ss).T])
    hits=first_hits([path],wall)
    if hits.segment[0]<0:raise ValueError('leg did not hit wall')
    path=clip(path,hits.points[0])
    return dict(points=path.tolist(),strike_RZ_m=hits.points[0].tolist(),
                segment=int(hits.segment[0]),max_psiN_residual=float(np.max(abs((field.values(path)-psi_sep)/span))))

hessian=np.array([[field.spline.ev(xp[1],xp[0],dy=2),field.spline.ev(xp[1],xp[0],dx=1,dy=1)],
                  [field.spline.ev(xp[1],xp[0],dx=1,dy=1),field.spline.ev(xp[1],xp[0],dx=2)]])
eigen,basis=np.linalg.eigh(hessian);assert eigen[0]<0<eigen[1]
directions=[]
for sg in [-1,1]:
    v=np.sqrt(-eigen[0])*basis[:,1]+sg*np.sqrt(eigen[1])*basis[:,0]
    v/=np.linalg.norm(v)
    if v[1]>0:v=-v
    directions.append(v)
legs={}
for v in sorted(directions,key=lambda v:v[0]):
    name='inner' if v[0]<0 else 'outer'
    result=leg(v,.0002,.0005);fine=leg(v,.0001,.00025)
    result['strike_refinement_delta_m']=float(np.linalg.norm(np.array(result['strike_RZ_m'])-fine['strike_RZ_m']))
    assert result['strike_refinement_delta_m']<1e-5 and result['max_psiN_residual']<1e-6
    legs[name]=fine
    legs[name]['strike_refinement_delta_m']=result['strike_refinement_delta_m']
    print('leg',name,json.dumps({k:v for k,v in legs[name].items() if k!='points'}),flush=True)

sources={}; footprints=[]
for name,r_sep,sg in [('IMP',r_imp,-1),('OMP',r_omp,1)]:
    crossing=shapely.intersection(boundary,LineString([[field.r[0],axis[1]],[field.r[-1],axis[1]]]))
    pts=list(crossing.geoms) if hasattr(crossing,'geoms') else [crossing]
    candidates=[p.x for p in pts if (p.x-r_sep)*sg>0]
    wall_r=min(candidates,key=lambda r:abs(r-r_sep))
    limit=min(.1,abs(wall_r-r_sep));count=800
    edges=limit*np.linspace(0,1,count+1)**2;centers=(edges[:-1]+edges[1:])/2
    starts=np.c_[r_sep+sg*centers,np.full(count,axis[1])]
    assert shapely.contains(polygon,shapely.points(starts)).all()
    paths=trace_paths(field,starts,step=.001,max_length=8)
    hits=first_hits(paths,wall)
    drift=max(float(np.max(abs((field.values(p)-field.values([p[0]])[0])/span))) for p in paths)
    ranges=[]
    for segment in sorted(set(hits.segment[hits.segment>=0])):
        mask=hits.segment==segment;points=hits.points[mask]
        v=wall[segment+1]-wall[segment];alpha=(points-wall[segment])@v/np.dot(v,v)
        low,high=float(alpha.min()),float(alpha.max())
        row=dict(source=name,segment=int(segment),points=[(wall[segment]+low*v).tolist(),(wall[segment]+high*v).tolist()],
                 sampled_hit_count=int(mask.sum()),source_x_range_mm=[float(centers[mask].min()*1000),float(centers[mask].max()*1000)],
                 sampled_only=True,whole_segment_classification=False)
        ranges.append(row);footprints.append(row)
    display=[]
    for offset in [1,3,5,10,20,40,60,80]:
        i=int(np.argmin(abs(centers-offset/1000)))
        if hits.segment[i]>=0:
            display.append(dict(offset_mm=float(centers[i]*1000),points=clip(paths[i],hits.points[i]).tolist(),segment=int(hits.segment[i])))
    sources[name]=dict(midplane_RZ_m=[float(r_sep),float(axis[1])],source_radial_direction=sg,
                       maximum_offset_mm=limit*1000,source_cells=count,trace_step_m=.001,
                       hit_counts={str(r['segment']):r['sampled_hit_count'] for r in ranges},unresolved_paths=int(np.sum(hits.segment<0)),
                       max_psiN_drift=drift,paths=display,ranges=ranges,
                       heat_distribution_calculated=False,power_split=None)
    print('source',name,json.dumps({k:v for k,v in sources[name].items() if k not in ['paths','ranges']}),flush=True)

# The lower wall arc closes the PFR outline; no heat is inferred in this region.
arcs=[boundary.project(Point(legs[n]['strike_RZ_m'])) for n in ['inner','outer']]
a,b=sorted(arcs)
arc=np.asarray(substring(boundary,a,b).coords)
if np.max(arc[:,1])>xp[1]+.02:
    raise ValueError('selected PFR arc crosses the upper wall')
inner=np.asarray(legs['inner']['points']);outer=np.asarray(legs['outer']['points'])
if np.linalg.norm(arc[0]-inner[-1])>1e-8:arc=arc[::-1]
pfr=np.vstack([inner,arc[1:],outer[-2::-1]])
assert Polygon(pfr).is_valid

result=dict(schema='no-transport-baseline-review.v1',shot_sha256=manifest['shot_sha256'],time_ms=1400,
            baseline_geometry_sha256=manifest['config']['geometry_sha256'],wall=wall.tolist(),lcfs=lcfs.tolist(),axis_RZ_m=axis.tolist(),
            xpoint_RZ_m=xp.tolist(),legs=legs,pfr_polygon=pfr.tolist(),sources=sources,footprints=footprints,
            existing_design_segments=[68,59],requested_design='union of existing floor/shelf and baseline IMP/OMP SOL first-hit wall portions',
            source_interpretation='Actual IMP and OMP launch planes; separate future optimization and normalization, no 50:50 split',
            scope='Baseline topology and sampled first hits; no new heat solve, optimization, or transport',
            limitations=['Sampled footprints are not certified continuous footprint boundaries','All independent-source heat profiles and new optimization remain to be defined',
                         'Private flux region heat is not determined by this direct SOL mapping model'])
from design_region_policy import apply_design_policy
result=apply_design_policy(result)
(OUT/'baseline_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
np.savez_compressed(OUT/'field.npz',R=field.r,Z=field.z,psi=field.psi)
print('WROTE',OUT/'baseline_review.json',flush=True)
