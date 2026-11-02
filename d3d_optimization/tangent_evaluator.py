"""OMP candidate constraints with a separate immutable comparison baseline."""
from copy import copy
from collections import OrderedDict
import numpy as np
from .boundary_evaluator import BoundaryEvaluator
from .accelerated_transport import AcceleratedEvaluator
from .tangent_direct import TangentMultiSurface,TangentGuard
from .tangent_geometry import Curve

ACTIVE=('omp_entry','shelf','omp_connector','floor')


class TangentEvaluator(BoundaryEvaluator):
    def curves(self,c):
        return {s['name']:self.direct.base.curve(s['name'],c[j*4:j*4+4]) for j,s in enumerate(self.surfaces)}

    def precise_lengths(self,c):
        ratios={}; error=0.
        for name,curve in self.curves(c).items():
            a,b=curve.length_ratio(16),curve.length_ratio(32)
            ratios[name]=max(a,b)+abs(a-b); error=max(error,abs(a-b))
        return ratios,error

    def constraints(self,c):
        c=np.asarray(c,float)
        if c.shape!=(self.dimension,) or not np.all(np.isfinite(c)) or np.max(abs(c))>self.bound+1e-12:
            raise ValueError('finite bounded full coefficients required')
        curves=self.curves(c); k={n:curves[n].max_curvature() for n in ACTIVE}
        floor=curves['floor']; gap,u=floor.minimum_distance(self.lip)
        t=floor.tangent(1); target=self.left_direction
        angle=float(np.degrees(np.arctan2(abs(t[0]*target[1]-t[1]*target[0]),t@target)))
        opposite_zero=all(np.all(c[4*j:4*j+4]==0) for j,s in enumerate(self.surfaces) if s['name'] not in ACTIVE)
        return dict(passed=bool(max(k.values())<=1/self.rmin+1e-7 and gap>=self.gap_min-1e-10 and angle<1e-8 and opposite_zero),
            radius_minimum_m=self.rmin,curvature_max_1_m=k,
            actual_minimum_radius_m={n:1/v if v>0 else None for n,v in k.items()},
            gap_m=gap,gap_minimum_m=self.gap_min,gap_floor_parameter=u,
            gap_floor_point=floor.point(u).tolist(),lip_point=self.lip.tolist(),
            tangent_angle_error_deg=angle,opposite_coefficients_zero=opposite_zero)

    def evaluate(self,c):
        check=self.constraints(c)
        if not check['passed']:
            ratios,error=self.precise_lengths(c)
            return dict(geometry_valid=False,feasible=False,penalty=1000.,errors=['continuous_shape_constraint_failed'],
                precise_length_ratios=ratios,length_integration_error=error,elapsed_seconds=0.,shape_constraints=check)
        method=AcceleratedEvaluator.evaluate if self.backend=='gpu' else BoundaryEvaluator.evaluate
        result=method(self,c); result['shape_constraints']=check
        return result

    def describe(self,r):
        row=super().describe(r)
        if 'shape_constraints' in r: row['shape_constraints']=r['shape_constraints']
        if 'geometry_screen' in r: row['geometry_screen']=r['geometry_screen']
        return row


def constrained_evaluator(original,rmin,backend='cpu',wall_step=None):
    if not np.isfinite(rmin) or rmin<=0 or backend not in ('cpu','gpu'):
        raise ValueError('positive radius and cpu/gpu backend required')
    if original.metadata['source']!='OMP': raise ValueError('this study changes only OMP')
    e=object.__new__(TangentEvaluator); e.__dict__.update(original.__dict__)
    # copy.copy probes __setstate__ before base exists, recursing through the
    # frozen wrapper's __getattr__. Copy its owned attributes explicitly.
    e.direct=object.__new__(type(original.direct)); e.direct.__dict__.update(original.direct.__dict__)
    e.direct.base=copy(original.direct.base)
    e.direct.base.__class__=TangentMultiSurface
    e.direct.boundary_guard=copy(original.direct.boundary_guard)
    e.direct.boundary_guard.__class__=TangentGuard; e.direct.boundary_guard.e=e.direct.base
    e.cache=OrderedDict(); e.rmin=float(rmin); e.backend=backend
    if wall_step is not None: e.wall_step=float(wall_step)
    surface=next(s for s in e.surfaces if s['name']=='floor'); i=surface['segment_index']
    e.left_direction=e.wall[i+2]-e.wall[i+1]
    # Fixed lip is the endpoint before the horizontal overhang above the floor.
    e.lip=np.array([1.37256,-1.3247])
    if np.min(np.linalg.norm(e.wall-e.lip,axis=1))>1e-10:
        raise ValueError('expected fixed lip endpoint missing from supplied geometry')
    a,b=e.targets['floor'][[0,-1]]
    e.gap_min=Curve(a,b,np.zeros(4)).minimum_distance(e.lip)[0]
    # baseline, baseline_peak, and direct.baseline are the ORIGINAL geometry.
    return e
