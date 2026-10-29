"""Fast fixed-input geometry evaluation with complete boundary footprints."""
from collections import OrderedDict
from dataclasses import asdict
from time import perf_counter
import numpy as np
from scipy.interpolate import BSpline
from shapely.geometry import Point
from .independent_sources import prepare_evaluator,midplane_source
from .reference_transport import ThermalPaths,thermal_summary
from .frozen_mapping import Hits
from .boundary_transport import deposit_boundary


def prepare_boundary(root,source,reference,count=320,step=.0005,nodes=161,wall_step=.001,method='cubic',quadrature=64,bound=.1):
    direct,meta=prepare_evaluator(root,source,count=count,step=step,nodes=nodes,method=method,bound=bound,boundary_mode='continuous')
    # This is a newly owned evaluator. Old source files and saved results are untouched.
    rz=meta['source_RZ_m'];edges=direct.sources[0].edges
    imposed=midplane_source(direct.field,*rz,edges,reference.lambda_q_m,source)
    direct.base.sources=[imposed]
    direct.boundary_guard.h0=np.array([1/(np.linalg.norm(direct.field.magnetic([rz])[0])*imposed.normalization)])
    zero=np.zeros(direct.dimension)
    direct.base.baseline=direct.base.evaluate(zero)
    direct.base.baseline_peaks=np.array([max(direct.base.baseline['maps'][0].q)])
    direct.baseline=direct.evaluate(zero);direct.baseline_peaks=direct.baseline['peak_q']
    meta=dict(meta,widths_mm=[reference.lambda_q_m*1000],source_power_W=reference.power_W)
    return BoundaryEvaluator(direct,meta,reference,wall_step,quadrature)


class BoundaryEvaluator:
    def __init__(self,direct,metadata,reference,wall_step=.001,quadrature=64,thermal_paths=None):
        if len(direct.sources)!=1 or direct.sources[0].width!=reference.lambda_q_m:
            raise ValueError('one source matching supplied reference width required')
        self.direct,self.metadata,self.reference=direct,metadata,reference
        self.source=direct.sources[0];self.wall_step=wall_step;self.cache=OrderedDict()
        self.thermal_paths=thermal_paths or ThermalPaths(direct.field,direct.paths,metadata['source_RZ_m'],self.source,reference,quadrature_order=quadrature)
        if self.thermal_paths.reference!=reference:raise ValueError('thermal cache reference mismatch')
        self.length_rules=[]
        knots=np.r_[np.zeros(4),np.linspace(0,1,6)[1:-1],np.ones(4)]
        for order in [16,32]:
            x,w=np.polynomial.legendre.leggauss(order)
            points=np.concatenate([a+.1*(x+1) for a in np.arange(5)/5])
            weights=np.tile(w*.1,5)
            self.length_rules.append((BSpline(knots,np.eye(8),3).derivative()(points)[:,2:-2],weights))
        self.baseline=self.evaluate(np.zeros(self.dimension))
        self.baseline_peak=self.baseline.get('peak_q',np.nan)

    def __getattr__(self,name):return getattr(self.direct,name)

    def precise_lengths(self,c):
        values={};error=0.
        for j,surface in enumerate(self.surfaces):
            length=self.lengths[surface['name']]
            estimates=[float(w@np.sqrt(length**2+(basis@c[4*j:4*j+4])**2)/length) for basis,w in self.length_rules]
            values[surface['name']]=max(estimates)+abs(estimates[1]-estimates[0])
            error=max(error,abs(estimates[1]-estimates[0]))
        return values,error

    def hits(self,r):
        m=r['maps'][0]
        return Hits(m.segment,m.points,np.array([path.project(Point(point)) if seg>=0 else np.inf for path,point,seg in zip(self.paths,m.points,m.segment)]))

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
                    dep=deposit_boundary(self.field,hits,r['wall'],self.source.weights,thermal['sigma_m'],
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

    def describe(self,r):
        out={k:r[k] for k in ['geometry_valid','feasible','errors','precise_length_ratios','length_integration_error','elapsed_seconds']}
        if 'deposition' not in r:return out
        dep=r['deposition'];k=int(np.argmax(dep['q']))
        good=np.isfinite(dep['sigma_physical_m']);w=self.source.weights
        out.update(peak_q=float(r['peak_q']),peak_MW_m2=float(r['peak_q']*self.reference.power_W/1e6),
            peak_ratio=float(r['peak_q']/self.baseline_peak),peak_RZ_m=dep['points'][k].tolist(),
            ledger=dep['ledger'],diagnostics=dep['diagnostics'],
            area_quadrature_error=float(abs(dep['quadrature_integrated_power']-dep['ledger']['deposited'])),
            receiving_angles=r['receiving_angles'],
            region_power={n:float(dep['donor_region_power'][:,j].sum()) for j,n in enumerate(r['region_names'])},
            sigma_physical_mean_m=float(w[good]@dep['sigma_physical_m'][good]/w[good].sum()),
            mean_wall_transport_distance_m=float(w[good]@dep['donor_mean_wall_distance_m'][good]/w[good].sum()),
            max_wall_displacement_m=float(np.max(np.linalg.norm(r['wall']-self.direct.baseline['wall'],axis=1))),
            thermal=thermal_summary(r['thermal'],w,self.reference))
        return out


def export_boundary(path,e,r):
    if 'deposition' not in r:return
    d=r['deposition']
    np.savez_compressed(path,wall=r['wall'],points=d['points'],q_MW_m2=d['q']*e.reference.power_W/1e6,
        cell_parent=d['cell_parent'],cell_power=d['cell_power'],receiving_wall=d['wall'],area_edges=d['area_edges'],
        sigma_physical_m=d['sigma_physical_m'],sigma_area_m2=d['sigma_area_m2'],target_gradient=d['target_gradient'],
        donor_area_m2=d['donor_area_m2'],donor_region_power=d['donor_region_power'],
        donor_mean_wall_distance_m=d['donor_mean_wall_distance_m'],source_x_m=e.source.x,
        hit_points=r['hits'].points,hit_segment=r['hits'].segment,**{'thermal_'+k:v for k,v in r['thermal'].items()})
