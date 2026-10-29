"""Nominal operating-point adapter; original direct and broadening code preserved."""
from dataclasses import asdict
import numpy as np
from shapely.geometry import Point
from .frozen_mapping import Hits
from .reference_transport import ThermalPaths, thermal_summary
from .transport_deposition import DepositionControls, deposit_sigma
from .broadening_study import summarize as base_summary


class TransportEvaluator:
    def __init__(self,direct,metadata,reference,mode='thermal',wall_step_m=.002,
                 anchor_m=.005,thermal_paths=None,quadrature_order=128):
        if mode not in ('direct','thermal','frozen_transport'):
            raise ValueError('unknown transport control mode')
        index=[i for i,s in enumerate(direct.sources) if abs(s.width-reference.lambda_q_m)<1e-14]
        if len(index)!=1:raise ValueError('nominal lambda must match an existing source width')
        self.direct,self.reference,self.mode=direct,reference,mode
        self.index=index[0];self.sources=[direct.sources[self.index]]
        self.metadata=dict(metadata,source_power_W=reference.power_W,widths_mm=[reference.lambda_q_m*1000])
        self.controls=DepositionControls(wall_step_m,anchor_m)
        self.config=self.controls
        self.thermal_paths=(thermal_paths if thermal_paths is not None else
            ThermalPaths(direct.field,direct.paths,metadata['source_RZ_m'],self.sources[0],reference,
                         quadrature_order=quadrature_order))
        if self.thermal_paths.reference!=reference:
            raise ValueError('cached thermal paths use a different reference point')
        self.path_lengths=self.thermal_paths.lengths
        baseline_hits=self._hits(direct.baseline)
        self.frozen_sigma=self.thermal_paths.evaluate(baseline_hits)['sigma_m'].copy()
        self.baseline=self.evaluate(np.zeros(direct.dimension))
        self.baseline_peaks=self.baseline.get('peak_q')
        self.baseline_admissible=bool(self.baseline['feasible'] and self.baseline_peaks is not None
                                     and np.all(np.isfinite(self.baseline_peaks)) and np.all(self.baseline_peaks>0))

    def __getattr__(self,name):return getattr(self.direct,name)

    def _hits(self,result):
        m=result['maps'][self.index]
        return Hits(m.segment,m.points,np.array([
            p.project(Point(point)) if seg>=0 else np.inf
            for p,point,seg in zip(self.paths,m.points,m.segment)]))

    def evaluate(self,coefficients):
        original=self.direct.evaluate(coefficients);r=dict(original)
        r.update(model_version='reference-two-point-one-way.v1',mode=self.mode,reference=asdict(self.reference))
        if not r['geometry_valid']:return r
        hits=self._hits(original)
        r['maps']=[original['maps'][self.index]]
        r['peak_q']=np.asarray(original['peak_q'])[[self.index]]
        r['surface_power']={k:[v[self.index]] for k,v in original['surface_power'].items()}
        if r.get('boundary_witness'):
            r['boundary_witness']=dict(r['boundary_witness'])
            if 'q' in r['boundary_witness']:
                r['boundary_witness']['q']=[r['boundary_witness']['q'][self.index]]
        # Rebuild single-width accounting; the original aggregates five widths.
        nominal=r['maps'][0]
        r['missing']=float(nominal.unmapped_power)
        r['power_balance_error']=float(nominal.balance_error)
        angle_ok=True;penalty=0.
        for surface in self.surfaces:
            floor=surface['angle_floor_deg']
            if floor is None:continue
            status=nominal.angle_application(floor,r['surface_segments'][surface['name']],nominal.hit_power)
            if status['status']=='target_angle_violation':
                angle_ok=False;penalty+=100*(floor-status['minimum_angle_deg'])
            elif status['status']=='invalid_target_angle':
                angle_ok=False;penalty+=1000
        witness=r.get('boundary_witness')
        boundary_ok=bool(witness and witness['resolved'] and witness['angle_passed'])
        if not boundary_ok:penalty+=1000
        illuminated=any(v[0]>0 for v in r['surface_power'].values())
        r['feasible']=bool(angle_ok and boundary_ok and illuminated and r['missing']<=.001 and r['power_balance_error']<=1e-4)
        r['penalty']=penalty+1000*max(0,r['missing']-.001)+(0 if illuminated else 1000)
        nominal_feasible,nominal_penalty=r['feasible'],r['penalty']
        if self.mode=='direct':
            r['peak_q_W_m2']=self.reference.power_W*r['peak_q']
            return r
        r['direct_maps']=r.pop('maps')
        r['direct_boundary_witness']=r.pop('boundary_witness',None)
        r.pop('peak_q',None)
        try:
            thermal=self.thermal_paths.evaluate(hits)
            sigma=self.frozen_sigma if self.mode=='frozen_transport' else thermal['sigma_m']
            dep=deposit_sigma(self.field,self.paths,self.path_lengths,hits,r['wall'],
                self.sources[0].weights[:,None],self.metadata['source_RZ_m'],
                -1 if self.metadata['source']=='IMP' else 1,self.lcfs,self.vessel_solids,
                self.controls,sigma_m=sigma)
        except (ValueError,FloatingPointError) as error:
            r.update(feasible=False,penalty=1e6,model_errors=[str(error)])
            return r
        r['thermal']=thermal
        r['thermal_diagnostics']=thermal_summary(thermal,self.sources[0].weights,self.reference)
        good=np.isfinite(sigma);w=self.sources[0].weights
        r['deposited_width_statistics']=dict(min_m=float(np.min(sigma[good])),max_m=float(np.max(sigma[good])),
            power_weighted_mean_m=float(np.dot(w[good],sigma[good])/w[good].sum()),
            fixed_at_baseline=self.mode=='frozen_transport')
        r['deposition']=dep;r['peak_q']=np.max(dep['q'],axis=0)
        r['peak_q_W_m2']=self.reference.power_W*r['peak_q']
        r['missing']=float(np.max(dep['ledger']['unresolved']))
        r['power_balance_error']=float(np.max(dep['ledger']['balance_error']))
        incidence=self.field.incidence(dep['points'],dep['normals'])
        powered=np.any(dep['cell_power'][dep['sample_cell']]>0,axis=1)
        angle_ok=True;angle_penalty=0.;angles={}
        for surface in self.surfaces:
            name,floor=surface['name'],surface['angle_floor_deg']
            mask=np.isin(dep['sample_parent'],r['surface_segments'][name])&powered
            minimum=float(np.min(incidence[mask])) if np.any(mask) else None
            passed=floor is None or minimum is None or (np.isfinite(minimum) and minimum>=floor)
            angles[name]=dict(minimum_angle_deg=minimum,floor_deg=floor,
                              receiving_sample_count=int(mask.sum()),passed=bool(passed))
            if not passed:
                angle_ok=False;angle_penalty+=100*max(0.,floor-minimum) if np.isfinite(minimum) else 1000
        r['receiving_angles']=angles
        r['surface_power']={name:np.sum(dep['cell_power'][np.isin(dep['cell_parent'],ids)],axis=0).tolist()
                            for name,ids in r['surface_segments'].items()}
        finite=bool(np.all(np.isfinite(dep['q'])) and np.all(dep['q']>=0))
        illuminated=any(v[0]>0 for v in r['surface_power'].values())
        r['feasible']=bool(nominal_feasible and dep['valid'] and finite and angle_ok and illuminated
                          and r['missing']<=self.controls.missing_tolerance and r['power_balance_error']<=1e-10)
        r['penalty']=nominal_penalty+angle_penalty+1000*max(0,r['missing']-self.controls.missing_tolerance)+(0 if dep['valid'] and finite and illuminated else 1000)
        r['model_errors']=list(dep['reasons'])
        return r


def summarize(evaluator,result):
    row=base_summary(evaluator,result)
    row['reference']=asdict(evaluator.reference)
    row['mode']=evaluator.mode
    if 'peak_q' in row:
        row['peak_q_MW_m2']=(np.asarray(row['peak_q'])*evaluator.reference.power_W/1e6).tolist()
    if 'thermal_diagnostics' in result:row['thermal']=result['thermal_diagnostics']
    if 'deposited_width_statistics' in result:row['deposited_width_statistics']=result['deposited_width_statistics']
    return row
