"""Immutable candidates bind every geometry, thermal and export operation."""
from dataclasses import dataclass
import hashlib,json
import numpy as np
from shapely.geometry import LineString,Polygon
from .imp_tangent import imp_evaluator,expand_imp
from .tangent_evaluator import constrained_evaluator

GEOMETRY_POLICY = 'full-wall-chord-1e-7.v2'
CHORD_TOLERANCE_M = 1e-7


@dataclass(frozen=True)
class Candidate:
    source: str
    radius_mm: int
    vector: tuple

    def __post_init__(self):
        if self.source not in ('IMP','OMP') or self.radius_mm not in (20,50,100):raise ValueError('unsupported source/radius')
        x=np.array(self.vector)
        if x.shape!=((9,) if self.source=='IMP' else (16,)) or not np.all(np.isfinite(x)):
            raise ValueError('finite regional candidate required')
        if max(abs(x[:8] if self.source=='IMP' else x))>.1+1e-12:raise ValueError('coefficient bound exceeded')
        if self.source=='IMP' and not 0<=x[8]<=1:raise ValueError('junction fraction outside [0,1]')

    @property
    def coefficients(self):return expand_imp(self.vector[:8]) if self.source=='IMP' else np.r_[self.vector,np.zeros(8)]
    @property
    def junction_fraction(self):return self.vector[8] if self.source=='IMP' else None
    @property
    def identity(self):return hashlib.sha256(json.dumps(self.record(),sort_keys=True).encode()).hexdigest()
    def record(self):return dict(source=self.source,radius_mm=self.radius_mm,vector=list(self.vector))
    @classmethod
    def from_record(cls,row):return cls(row['source'],row['radius_mm'],tuple(row['vector']))
    def bind(self,e):
        if e.metadata['source']!=self.source or abs(e.rmin*1000-self.radius_mm)>1e-8:
            raise ValueError('candidate/evaluator case mismatch')
        if self.source=='IMP':e.set_junction(self.junction_fraction)
        return self.coefficients


def evaluator(original,mm,backend='gpu'):
    return (imp_evaluator if original.metadata['source']=='IMP' else constrained_evaluator)(original,mm/1000,backend)


def assembled_wall(e,c,chord_tolerance=None):
    curves=e.curves(c); by_index={s['segment_index']:s['name'] for s in e.surfaces}; points=[e.wall[0]]
    for i in range(len(e.wall)-1):
        if i not in by_index:points.append(e.wall[i+1]);continue
        name=by_index[i]; q=curves[name]; count=len(e.targets[name])
        if chord_tolerance is not None:
            # Bound chord error by sup|p''| du²/8. d'' is linear per cubic
            # span, so its supremum occurs at span endpoints (knots).
            second=max(abs(q.d2(np.linspace(0,1,6))))
            # Canonical screen mesh must not depend on the thermal mesh.
            count=max(161,int(np.ceil(np.sqrt(second/(8*chord_tolerance))))+1)
        p=q.point(np.linspace(0,1,count)); p[0]=q.a; p[-1]=q.b
        points.extend(p[1:])
    return np.array(points)


def geometry_screen(e,c,chord_tolerance=CHORD_TOLERANCE_M):
    if chord_tolerance is not None and (not np.isfinite(chord_tolerance) or chord_tolerance<=0):
        raise ValueError('positive finite chord tolerance required')
    wall=assembled_wall(e,c,chord_tolerance)
    line=LineString(wall); errors=[]
    if not line.is_simple:errors.append('wall_self_intersection')
    if not Polygon(wall).is_valid:errors.append('invalid_wall_polygon')
    if line.intersects(Polygon(e.lcfs)):errors.append('lcfs_intersection')
    if e.vessel_solids is not None and line.intersects(e.vessel_solids):errors.append('vessel_intersection')
    return dict(policy=GEOMETRY_POLICY if chord_tolerance==CHORD_TOLERANCE_M else 'diagnostic-only',
                chord_tolerance_m=chord_tolerance,points=len(wall),passed=not errors,errors=errors)


def geometry_check(e,c,chord_tolerance=CHORD_TOLERANCE_M):
    return geometry_screen(e,c,chord_tolerance)['passed']


def evaluate_candidate(e,candidate):
    c=candidate.bind(e); check=e.constraints(c)
    screen=geometry_screen(e,c) if check['passed'] else None
    if screen is not None and not screen['passed']:
        ratios,error=e.precise_lengths(c)
        return dict(geometry_valid=False,feasible=False,penalty=1000.,errors=['complete_wall_invalid'],
                    elapsed_seconds=0.,precise_length_ratios=ratios,length_integration_error=error,
                    shape_constraints=check,geometry_screen=screen)
    result=dict(e.evaluate(c))  # Do not attach policy metadata to the evaluator's cached object.
    if screen is not None:result['geometry_screen']=screen
    return result
