"""Conservative area-measure boundary-footprint model, not a 2D plasma solver.

Closed connected wall diffusion replaces legacy same-flux-branch clipping.
The finite-width convention has no grazing-incidence projection.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from shapely.geometry import LineString,Polygon
from .broadening import _cdf_interval


def _parameters(center,sigma,period):
    if not all(np.isfinite(x) for x in [center,sigma,period]) or sigma<=0 or period<=0:
        raise ValueError('positive finite sigma and period; finite center required')
    return center%period,sigma/period


def _modes(ratio):
    if ratio>2:return np.array([],float),np.array([],float)
    # First omitted Fourier coefficient <1e-17; successive terms decay faster.
    n=np.arange(1,int(np.ceil(np.sqrt(2*np.log(1e17))/(2*np.pi*ratio)))+1,dtype=float)
    return n,np.exp(-2*(np.pi*n*ratio)**2)


def wrapped_mass(lo,hi,center,sigma,period):
    center,ratio=_parameters(center,sigma,period)
    lo,hi=np.broadcast_arrays(np.asarray(lo,float),np.asarray(hi,float))
    if np.any(~np.isfinite(lo)|~np.isfinite(hi)|(lo<0)|(hi<lo)|(hi>period)):
        raise ValueError('ordered intervals inside one period required')
    if ratio<.12:
        # Integration covers [-P,2P], leaving at least 8.33 sigma at each end.
        value=sum(_cdf_interval(((lo+k*period)-center)/sigma,((hi+k*period)-center)/sigma) for k in [-1,0,1])
    else:
        value=(hi-lo)/period
        n,c=_modes(ratio)
        for mode,coefficient in zip(n,c):
            value=value+coefficient/(np.pi*mode)*(np.sin(2*np.pi*mode*(hi-center)/period)-np.sin(2*np.pi*mode*(lo-center)/period))
    if np.any(value < -2e-14):raise FloatingPointError('negative wrapped Gaussian mass')
    return np.maximum(value,0.)


def wrapped_density(a,center,sigma,period):
    center,ratio=_parameters(center,sigma,period);a=np.asarray(a,float)%period
    if not np.all(np.isfinite(a)):raise ValueError('finite receiving coordinates required')
    if ratio<.12:
        d=a-center
        d=np.where(d>period/2,d-period,np.where(d< -period/2,d+period,d))
        value=sum(np.exp(-.5*((d+k*period)/sigma)**2) for k in [-1,0,1])/(np.sqrt(2*np.pi)*sigma)
    else:
        value=np.ones_like(a);n,c=_modes(ratio)
        for mode,coefficient in zip(n,c):value+=2*coefficient*np.cos(2*np.pi*mode*(a-center)/period)
        value/=period
    if np.any(value < -2e-14/period):raise FloatingPointError('negative wrapped Gaussian density')
    return np.maximum(value,0.)


class AreaBoundary:
    """One closed material boundary; all areas and partial-segment maps exact."""
    def __init__(self,wall,step,solids=None):
        wall=np.asarray(wall,float)
        if (wall.ndim!=2 or wall.shape[1]!=2 or len(wall)<4 or not np.all(np.isfinite(wall))
                or np.any(wall[:,0]<=0) or np.linalg.norm(wall[0]-wall[-1])>1e-10
                or not np.isfinite(step) or step<=0):raise ValueError('one closed positive-R wall and positive step required')
        line=LineString(wall)
        if not line.is_simple or not Polygon(wall).is_valid:raise ValueError('wall must be one simple connected closed component')
        if solids is not None and not solids.is_empty and line.relate_pattern(solids,'T********'):
            raise ValueError('receiving boundary crosses material interior')
        self.original=wall;d=np.diff(wall,axis=0);self.original_lengths=np.linalg.norm(d,axis=1)
        if np.any(self.original_lengths<=0):raise ValueError('zero-length wall segment')
        self.original_areas=np.pi*(wall[:-1,0]+wall[1:,0])*self.original_lengths
        self.original_edges=np.r_[0,np.cumsum(self.original_areas)]
        points=[wall[0]];parents=[]
        for j,(a,b,length) in enumerate(zip(wall[:-1],wall[1:],self.original_lengths)):
            count=max(1,int(np.ceil(length/step)))
            points.extend(np.linspace(a,b,count+1)[1:]);parents.extend([j]*count)
        self.wall=np.asarray(points);self.parent=np.asarray(parents)
        self.delta=np.diff(self.wall,axis=0);self.ds=np.linalg.norm(self.delta,axis=1)
        self.tangent=self.delta/self.ds[:,None]
        self.areas=np.pi*(self.wall[:-1,0]+self.wall[1:,0])*self.ds
        self.edges=np.r_[0,np.cumsum(self.areas)];self.total=float(self.edges[-1])
        self.arc=np.r_[0,np.cumsum(self.ds)]

    def hit_area(self,points,segments):
        points=np.asarray(points,float);segments=np.asarray(segments,int)
        if points.shape!=(len(segments),2) or np.any(segments<0) or np.any(segments>=len(self.original)-1):
            raise ValueError('one valid parent segment per hit')
        a=self.original[segments];d=self.original[segments+1]-a
        fraction=((points-a)*d).sum(axis=1)/(d*d).sum(axis=1)
        if np.any(fraction< -1e-8)|np.any(fraction>1+1e-8)|np.any(np.linalg.norm(a+fraction[:,None]*d-points,axis=1)>1e-7):
            raise ValueError('first hit is not on declared wall segment')
        fraction=np.clip(fraction,0,1)
        value=self.original_edges[segments]+2*np.pi*self.original_lengths[segments]*(a[:,0]*fraction+.5*d[:,0]*fraction**2)
        return value%self.total

    def points_at(self,area):
        area=np.asarray(area,float)%self.total
        cell=np.clip(np.searchsorted(self.edges,area,side='right')-1,0,len(self.ds)-1)
        da=area-self.edges[cell];r=self.wall[cell,0];dr=self.delta[cell,0]
        y=da/(np.pi*self.ds[cell]);t=y/(r+np.sqrt(np.maximum(r*r+dr*y,0)))
        return self.wall[cell]+t[:,None]*self.delta[cell],cell

    def samples(self):
        t=np.array([0.,.25,.5,.75,1.])
        points=(self.wall[:-1,None,:]+t[None,:,None]*self.delta[:,None,:]).reshape(-1,2)
        area=(self.edges[:-1,None]+2*np.pi*self.ds[:,None]*(self.wall[:-1,0,None]*t+.5*self.delta[:,0,None]*t*t)).ravel()
        quadrature=(2*np.pi*points[:,0].reshape(-1,5)*self.ds[:,None]*np.array([7,32,12,32,7])/90).ravel()
        return points,area,np.repeat(np.arange(len(self.ds)),5),quadrature

    def area_at_arc(self,arc):
        arc=np.asarray(arc,float)%self.arc[-1]
        cell=np.clip(np.searchsorted(self.arc,arc,side='right')-1,0,len(self.ds)-1)
        t=(arc-self.arc[cell])/self.ds[cell]
        return self.edges[cell]+2*np.pi*self.ds[cell]*(self.wall[cell,0]*t+.5*self.delta[cell,0]*t*t)


def _refined_peak(area,q,centers,widths,weights,period):
    """Refine sampled local maxima in dimensionless brackets, including the seam.

    Donor-local and regular wall samples bracket modes. Resolution studies remain
    necessary: this numerical search is not a formal global mixture-mode proof.
    """
    narrow=widths/period<.12
    mu=centers[narrow];sigma=widths[narrow];amplitude=weights[narrow]/(np.sqrt(2*np.pi)*sigma)
    broad=list(zip(centers[~narrow],widths[~narrow],weights[~narrow]))
    def density(a):
        a=a%period;d=a-mu
        d=np.where(d>period/2,d-period,np.where(d< -period/2,d+period,d))
        with np.errstate(over='ignore'):
            value=sum(np.dot(amplitude,np.exp(-.5*((d+k*period)/sigma)**2)) for k in [-1,0,1])
        return float(value+sum(w*wrapped_density(a,c,s,period) for c,s,w in broad))
    aa,index=np.unique(area%period,return_index=True);qq=q[index]
    candidates=np.flatnonzero((qq>=np.roll(qq,1))&(qq>=np.roll(qq,-1))&
        ((qq>np.roll(qq,1))|(qq>np.roll(qq,-1)))&(qq>qq.max()*1e-10))
    best=int(np.argmax(q));best_a=float(area[best]%period);best_q=float(q[best])
    for j in candidates:
        lo=aa[j-1] if j else aa[-1]-period
        hi=aa[j+1] if j+1<len(aa) else aa[0]+period
        if hi<=lo:continue
        result=minimize_scalar(lambda t:-density(lo+t*(hi-lo)),bounds=(0.,1.),method='bounded',options={'xatol':1e-11})
        if -result.fun>best_q:best_a=float((lo+result.x*(hi-lo))%period);best_q=float(-result.fun)
    return best_a,best_q,len(candidates)


def deposit_boundary(field,hits,wall,weights,sigma_flux,grad_sep,step,solids=None,region_by_segment=None):
    mesh=AreaBoundary(wall,step,solids)
    w=np.asarray(weights,float);sigma_flux=np.asarray(sigma_flux,float)
    if w.ndim!=1 or np.any(~np.isfinite(w)) or np.any(w<0) or not np.isfinite(w.sum()) or w.sum()<=0:
        raise ValueError('finite nonnegative donor weights with positive total required')
    if sigma_flux.shape!=w.shape or np.shape(hits.segment)!=w.shape:
        raise ValueError('one width and hit per donor required')
    if not np.isfinite(grad_sep) or grad_sep<=0:raise ValueError('positive source gradient required')
    good=(np.asarray(hits.segment)>=0)&(w>0);ids=np.flatnonzero(good)
    centers=np.full(len(w),np.nan);physical=np.full(len(w),np.nan);widths=np.full(len(w),np.nan)
    gradients=np.full(len(w),np.nan)
    if len(ids):
        gradients[ids]=np.linalg.norm(field.gradient(np.asarray(hits.points)[ids]),axis=1)
        if np.any(~np.isfinite(gradients[ids])) or np.any(gradients[ids]<=0):
            raise ValueError('undefined width conversion at zero/nonfinite target gradient')
        physical[ids]=sigma_flux[ids]*grad_sep/gradients[ids]
        widths[ids]=2*np.pi*np.asarray(hits.points)[ids,0]*physical[ids]
        if np.any(~np.isfinite(widths[ids])) or np.any(widths[ids]<=0):
            raise ValueError('positive finite converted footprint width required')
        if np.any(widths[ids]<32*np.finfo(float).eps*mesh.total):
            raise ValueError('footprint width below floating-point coordinate resolution')
        centers[ids]=mesh.hit_area(np.asarray(hits.points)[ids],np.asarray(hits.segment)[ids])
    points,area,sample_cell,quadrature=mesh.samples();regular_count=len(points)
    if len(ids):
        extra=(centers[ids,None]+np.minimum(widths[ids,None],mesh.total)*np.array([-2.,-1.,-.5,0.,.5,1.,2.]))%mesh.total
        ep,ec=mesh.points_at(extra.ravel());points=np.vstack([points,ep]);area=np.r_[area,extra.ravel()];sample_cell=np.r_[sample_cell,ec]
    region=np.zeros(len(mesh.original)-1,int) if region_by_segment is None else np.asarray(region_by_segment,int)
    if region.shape!=(len(mesh.original)-1,) or np.any(region<0):raise ValueError('one nonnegative region per wall segment required')
    cell_region=region[mesh.parent];regions=int(region.max())+1
    q=np.zeros(len(points));power=np.zeros(len(mesh.ds));closure=np.zeros(len(w))
    region_power=np.zeros((len(w),regions));mean_distance=np.full(len(w),np.nan);far_power=0.
    arc_mid=(mesh.arc[:-1]+mesh.arc[1:])/2
    original_arc=np.r_[0,np.cumsum(mesh.original_lengths)]
    for i in ids:
        mass=wrapped_mass(mesh.edges[:-1],mesh.edges[1:],centers[i],widths[i],mesh.total)
        closure[i]=mass.sum()
        power+=w[i]*mass
        q+=w[i]*wrapped_density(area,centers[i],widths[i],mesh.total)
        region_power[i]=w[i]*np.bincount(cell_region,weights=mass,minlength=regions)
        segment=int(hits.segment[i]);hit_arc=original_arc[segment]+np.linalg.norm(np.asarray(hits.points)[i]-mesh.original[segment])
        distance=np.abs(arc_mid-hit_arc);distance=np.minimum(distance,mesh.arc[-1]-distance)
        mean_distance[i]=np.dot(mass,distance)
        if mesh.arc[-1]>.2:
            a,b=mesh.area_at_arc([hit_arc-.1,hit_arc+.1])
            near=wrapped_mass(a,b,centers[i],widths[i],mesh.total) if a<=b else (
                wrapped_mass(a,mesh.total,centers[i],widths[i],mesh.total)+wrapped_mass(0,b,centers[i],widths[i],mesh.total))
            far_power+=w[i]*max(0.,1-float(near))
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(power)):
        raise FloatingPointError('nonfinite boundary heat flux or cell power')
    peak_area,peak_q,peak_brackets=_refined_peak(area,q,centers[ids],widths[ids],w[ids],mesh.total)
    pp,pc=mesh.points_at([peak_area]);points=np.vstack([points,pp]);sample_cell=np.r_[sample_cell,pc];q=np.r_[q,peak_q]
    quadrature_power=float(quadrature@q[:regular_count])
    if not np.isfinite(peak_q) or not np.isfinite(quadrature_power):raise FloatingPointError('nonfinite peak or quadrature power')
    deposited=float(power.sum());unhit=float(w[np.asarray(hits.segment)<0].sum());prescribed=float(w.sum())
    balance=abs(deposited+unhit-prescribed)/prescribed
    maxclosure=float(np.max(abs(closure[ids]-1))) if len(ids) else 0.
    if balance>2e-12 or maxclosure>2e-12:raise FloatingPointError('periodic kernel failed power closure')
    return dict(wall=mesh.wall,cell_parent=mesh.parent,points=points,sample_cell=sample_cell,
        sample_parent=mesh.parent[sample_cell],normals=np.column_stack([-mesh.tangent[sample_cell,1],mesh.tangent[sample_cell,0]]),
        q=q,cell_power=power,areas=mesh.areas,area_edges=mesh.edges,total_area_m2=mesh.total,
        sigma_flux_m=sigma_flux,sigma_physical_m=physical,sigma_area_m2=widths,target_gradient=gradients,
        donor_area_m2=centers,donor_region_power=region_power,donor_mean_wall_distance_m=mean_distance,
        quadrature_integrated_power=quadrature_power,
        ledger=dict(prescribed=prescribed,deposited=deposited,unhit=unhit,unresolved=unhit,balance_error=balance),
        diagnostics=dict(max_donor_closure_error=maxclosure,
            peak_refinement_brackets=peak_brackets,
            wall_distance_estimate_error_bound_m=float(mesh.ds.max()/2),
            broad_kernel_power_fraction=float(w[good&(widths/mesh.total>.1)].sum()/prescribed),
            power_beyond_0p1m_wall_distance=float(far_power/prescribed),
            width_convention='field-normal width converted to area; no grazing-incidence projection',
            receiver='single connected closed boundary; area-measure periodic Gaussian',
            missing_power_renormalized=False))
