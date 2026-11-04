"""Finite-volume distribution screen; no registration, normalization or global proof."""
import numpy as np
import copy,math
from scipy.spatial import cKDTree

POLICY='regional-distribution-convergence.v3'
LIMITS=dict(linf_relative=.01,l1_relative=.005,a2_relative=.005,power_relative=.005,
            wall_distance_m=50e-6,wall_length_difference_m=50e-6)

def profile(wall,q):
    wall=np.asarray(wall,dtype=float);q=np.asarray(q,dtype=float)
    if wall.ndim!=2 or wall.shape[1]!=2 or q.shape!=(len(wall)-1,) or len(q)<1:
        raise ValueError('invalid cell dimensions')
    ds=np.linalg.norm(np.diff(wall,axis=0),axis=1)
    if not np.all(np.isfinite(wall)) or not np.all(np.isfinite(q)) or np.any(q<0) or np.any(ds<=0) or np.any(wall[:,0]<=0):
        raise ValueError('invalid finite positive geometry / nonnegative q')
    area=np.pi*(wall[:-1,0]+wall[1:,0])*ds
    if np.dot(q,area)<=0:raise ValueError('nonpositive deposited power')
    return dict(wall=wall,q=q,area=area,length=float(ds.sum()),u=np.r_[0,np.cumsum(ds)]/ds.sum())

def from_npz(path,power_MW=1.):
    with np.load(path,allow_pickle=False) as z:
        wall=z['receiving_wall'];power=z['cell_power'];area=np.diff(z['area_edges'])
        if not np.isfinite(power_MW) or power_MW<=0 or not np.all(np.isfinite(power)) or np.any(power<0) or np.any(area<=0):raise ValueError('invalid power/area')
        p=profile(wall,power*power_MW/area)
        if not np.allclose(area,p['area'],rtol=2e-8,atol=1e-13):raise ValueError('area does not match physical wall')
        if abs(float(power.sum())-1)>2e-12:raise ValueError('assigned power not conserved; no renormalization')
        points=z['points'];q=z['q_MW_m2']
        if points.shape!=(len(q),2) or not np.all(np.isfinite(points)) or not np.all(np.isfinite(q)) or np.any(q<0):raise ValueError('invalid continuous samples')
        # Project each saved sample onto neighboring wall segments, retaining its
        # physical origin/orientation. Never align peaks or translate profiles.
        nearest=cKDTree(wall).query(points)[1]
        best=np.full(len(q),np.inf);u=np.zeros(len(q))
        for offset in (-1,0):
            ids=(nearest+offset)%(len(wall)-1);a=wall[ids];d=wall[ids+1]-a
            t=np.clip(np.sum((points-a)*d,axis=1)/np.sum(d*d,axis=1),0,1)
            dist=np.linalg.norm(points-(a+t[:,None]*d),axis=1);mask=dist<best
            best[mask]=dist[mask];u[mask]=(p['u'][ids]+t*(p['u'][ids+1]-p['u'][ids]))[mask]
        if best.max()>1e-7:raise ValueError('continuous sample not on its saved wall')
        order=np.argsort(u);u=u[order];q=q[order]
        unique,first=np.unique(u,return_index=True)
        # Duplicates must have the same function value, not arbitrary tie choice.
        maximum=np.maximum.reduceat(q,first);minimum=np.minimum.reduceat(q,first)
        if np.max(maximum-minimum)>1e-8*max(1.,q.max()):raise ValueError('conflicting continuous samples')
        # The stored wrapped wall kernel is periodic; cover the seam explicitly.
        seam=float(np.interp(0,np.r_[unique[-1]-1,unique,unique[0]+1],np.r_[maximum[-1],maximum,maximum[0]]))
        p['continuous']=(np.r_[0,unique[(unique>0)&(unique<1)],1],np.r_[seam,maximum[(unique>0)&(unique<1)],seam])
        return p

def continuous_compare(a,b,limit=.01):
    ua,qa=a['continuous'];ub,qb=b['continuous'];u=np.unique(np.r_[ua,ub])
    diff=abs(np.interp(u,ua,qa)-np.interp(u,ub,qb));scale=max(qa.max(),qb.max())
    error=float(diff.max()/scale)
    return dict(passed=bool(np.isfinite(error) and error<=limit),linf_relative=error,limit=limit,
        samples=len(u),scope='piecewise linear reconstruction of saved continuous samples; not a continuum bound')

def distribution_compare(a,b,canonical=None):
    """Compare same-support reconstructed q; cell-step differences are diagnostics."""
    c=a if canonical is None else canonical
    step=compare(a,b,c);ua,qa=a['continuous'];ub,qb=b['continuous']
    u=np.unique(np.r_[ua,ub,c['u']]);x=(u[:-1]+u[1:])/2;half=np.diff(u)/2
    # Three-point Gauss integrates q² times R exactly on every linear interval.
    g,w=np.polynomial.legendre.leggauss(3);nodes=x[:,None]+half[:,None]*g
    r=np.interp(nodes.ravel(),c['u'],c['wall'][:,0]).reshape(nodes.shape)
    area=2*np.pi*r*c['length']*half[:,None]*w
    va=np.interp(nodes.ravel(),ua,qa).reshape(nodes.shape);vb=np.interp(nodes.ravel(),ub,qb).reshape(nodes.shape)
    powers=[float(np.sum(v*area)) for v in (va,vb)]
    A=[P*P/float(np.sum(v*v*area)) for P,v in zip(powers,(va,vb))]
    point=continuous_compare(a,b);l1=float(np.sum(abs(va-vb)*area)/(sum(powers)/2))
    metrics=dict(linf_relative=point['linf_relative'],l1_relative=l1,
        a2_relative=abs(A[0]-A[1])/max(A),power_relative=abs(powers[0]-powers[1])/max(powers),
        wall_distance_m=step['wall_distance_m'],wall_length_difference_m=step['wall_length_difference_m'])
    failures=[k for k,val in metrics.items() if not np.isfinite(val) or val>LIMITS[k]]
    # Reject a poor sample reconstruction even when two bad reconstructions agree.
    reconstruction=[abs(P-float(np.dot(p['q'],p['area'])))/float(np.dot(p['q'],p['area'])) for P,p in zip(powers,(a,b))]
    if max(reconstruction)>.001:failures.append('continuous_reconstruction_power')
    if step['a2_relative']>LIMITS['a2_relative']:failures.append('cell_a2_relative')
    return dict(passed=not failures,failures=failures,**metrics,limits=LIMITS,
        reconstructed_A2_m2=A,reconstructed_power_MW=powers,reconstruction_power_relative=reconstruction,
        native_cell_diagnostic=step,sample_comparison=point,
        scope='same physical support, saved-sample piecewise-linear q; not a continuum certificate')

def qualify(validation,distribution):
    from qualification_selection import passes_improvement
    v=copy.deepcopy(validation);old=v['verdict'];v['previous_verdict']=copy.deepcopy(old)
    thresholds=dict(max_peak_ratio_change=.002,max_absolute_peak_relative_change=.005,
        max_region_power_l1_change=.005,max_wall_cdf_change=.002)
    numeric=all(math.isfinite(old.get(k,float('nan'))) and 0<=old[k]<=t for k,t in thresholds.items())
    numeric &= math.isfinite(v['wall_only_refinement']['max_peak_relative_change']) and 0<=v['wall_only_refinement']['max_peak_relative_change']<=.005
    ratios=[l['rows']['selected']['peak_ratio'] for l in v['levels']]
    gain=1-ratios[2]
    improved=all(math.isfinite(r) and 0<r<=1+1e-10 for r in ratios) and gain>2*old['max_peak_ratio_change']+1e-4
    passed=bool(numeric and distribution['passed'])
    v['verdict']={**old,'resolution_screen_passed':passed,'improvement_confirmed':bool(passed and improved),
                  'peak_location_diagnostic_only':True,'distribution_policy':POLICY}
    # Existing strict independent geometry/thermal/feasibility requirements still apply.
    if not passes_improvement(v):
        v['verdict']['improvement_confirmed']=False
        if not improved:v['verdict']['non_improvement']=True
        else:v['verdict']['resolution_screen_passed']=False
    v['distribution']=distribution;v['qualification_policy']=POLICY
    return v

def compare(a,b,canonical=None,limits=None):
    limits=dict(LIMITS if limits is None else limits);c=a if canonical is None else canonical
    u=np.unique(np.r_[a['u'],b['u'],c['u']])
    def coords(p):return np.column_stack([np.interp(u,p['u'],p['wall'][:,j]) for j in range(2)])
    pa,pb,pc=map(coords,(a,b,c))
    distance=max(float(np.linalg.norm(p-pc,axis=1).max()) for p in (pa,pb))
    distance=max(distance,float(np.linalg.norm(pa-pb,axis=1).max()))
    length=max(abs(p['length']-c['length']) for p in (a,b))
    length=max(length,abs(a['length']-b['length']))
    area=np.pi*(pc[:-1,0]+pc[1:,0])*np.linalg.norm(np.diff(pc,axis=0),axis=1)
    def values(p):return p['q'][np.minimum(np.searchsorted(p['u'],u[:-1],side='right')-1,len(p['q'])-1)]
    qa,qb=values(a),values(b);P=[float(np.dot(q,area)) for q in (qa,qb)]
    A=[power**2/float(np.dot(q*q,area)) for power,q in zip(P,(qa,qb))]
    metrics=dict(linf_relative=float(np.max(abs(qa-qb))/max(qa.max(),qb.max())),
        l1_relative=float(np.dot(abs(qa-qb),area)/(sum(P)/2)),
        a2_relative=abs(A[0]-A[1])/max(A),power_relative=abs(P[0]-P[1])/max(P),
        wall_distance_m=distance,wall_length_difference_m=length)
    failures=[k for k,v in metrics.items() if not np.isfinite(v) or v>limits[k]]
    return dict(passed=not failures,failures=failures,**metrics,a2_m2=A,power_MW=P,limits=limits,
                comparison_cells=len(area),scope='saved cell averages on fixed corresponding wall; subcell convergence not proven')
