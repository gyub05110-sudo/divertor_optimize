"""Additive GPU variant of the frozen deposit/evaluate functions.

CDF powers, geometry, thermal closure and continuous peak refinement are the
preserved CPU algorithms. Only the sampled-density sum is batched on CUDA.
Full CPU/GPU parity is required before use. Original source files are untouched.
"""
import numpy as np
from time import perf_counter
from .boundary_transport import AreaBoundary,wrapped_mass,_refined_peak
from .boundary_evaluator import BoundaryEvaluator
from .gpu_footprint import density_sum

def deposit_accelerated(field,hits,wall,weights,sigma_flux,grad_sep,step,solids=None,region_by_segment=None):
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
        region_power[i]=w[i]*np.bincount(cell_region,weights=mass,minlength=regions)
        segment=int(hits.segment[i]);hit_arc=original_arc[segment]+np.linalg.norm(np.asarray(hits.points)[i]-mesh.original[segment])
        distance=np.abs(arc_mid-hit_arc);distance=np.minimum(distance,mesh.arc[-1]-distance)
        mean_distance[i]=np.dot(mass,distance)
        if mesh.arc[-1]>.2:
            a,b=mesh.area_at_arc([hit_arc-.1,hit_arc+.1])
            near=wrapped_mass(a,b,centers[i],widths[i],mesh.total) if a<=b else (
                wrapped_mass(a,mesh.total,centers[i],widths[i],mesh.total)+wrapped_mass(0,b,centers[i],widths[i],mesh.total))
            far_power+=w[i]*max(0.,1-float(near))
    q=density_sum(area,centers[ids],widths[ids],w[ids],mesh.total)
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


class AcceleratedEvaluator(BoundaryEvaluator):
    def evaluate(self,c):
        start=perf_counter();c=np.asarray(c,float)
        if c.shape!=(self.dimension,) or np.any(~np.isfinite(c)) or np.max(abs(c))>self.bound+1e-12:
            raise ValueError('finite bounded coefficient vector required')
        key=tuple(c)
        if key in self.cache:self.cache.move_to_end(key);return self.cache[key]
        ratios,error=self.precise_lengths(c)
        r=dict(geometry_valid=False,feasible=False,penalty=1000.,errors=[],precise_length_ratios=ratios,length_integration_error=error)
        if max(ratios.values())>1.2+1e-12:r['errors']=['material_length_exceeds_limit']
        else:
            r.update(self.direct.evaluate(c))
            if not r['geometry_valid']:r['errors']=['geometry_constraint_failed']
            else:
                try:
                    hits=self.hits(r);thermal=self.thermal_paths.evaluate(hits)
                    region=np.zeros(len(r['wall'])-1,int)
                    names=['fixed']+[s['name'] for s in self.surfaces]
                    for j,name in enumerate(names[1:],1):region[r['surface_segments'][name]]=j
                    dep=deposit_accelerated(self.field,hits,r['wall'],self.source.weights,thermal['sigma_m'],
                        self.thermal_paths.grad_sep,self.wall_step,self.vessel_solids,region)
                    incidence=self.field.incidence(dep['points'],dep['normals'])
                    powered=dep['cell_power'][dep['sample_cell']]>0
                    angles={};angle_ok=True
                    for s in self.surfaces:
                        mask=np.isin(dep['sample_parent'],r['surface_segments'][s['name']])&powered
                        minimum=float(np.min(incidence[mask])) if np.any(mask) else None
                        passed=s['angle_floor_deg'] is None or minimum is None or (np.isfinite(minimum) and minimum>=s['angle_floor_deg'])
                        angles[s['name']]=dict(minimum_angle_deg=minimum,passed=bool(passed));angle_ok &= passed
                    r.update(deposition=dep,thermal=thermal,hits=hits,region_names=names,
                        receiving_angles=angles,direct_peak_q=float(r['peak_q'][0]),peak_q=float(np.max(dep['q'])))
                    r['feasible']=bool(r['feasible'] and angle_ok and dep['ledger']['unhit']<=.001 and dep['ledger']['balance_error']<=2e-12)
                    if not angle_ok:r['errors'].append('receiving_angle_constraint_failed');r['penalty']+=1000
                    if dep['ledger']['unhit']>.001:r['errors'].append('unhit_power');r['penalty']+=1000
                except (ValueError,FloatingPointError,OverflowError) as exc:
                    r.update(feasible=False,penalty=1e6,errors=[str(exc)])
        r['elapsed_seconds']=perf_counter()-start
        self.cache[key]=r
        while len(self.cache)>6:self.cache.popitem(last=False)
        return r
