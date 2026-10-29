"""Additive IMP G1 geometry; historical baseline and OMP studies stay frozen."""
from copy import copy
from collections import OrderedDict
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.interpolate import BSpline, PPoly
from scipy.optimize import minimize
from .tangent_geometry import Curve, KNOTS
from .tangent_direct import TangentMultiSurface, TangentGuard
from .tangent_evaluator import TangentEvaluator

ACTIVE=('imp_sloped','imp_vertical')


def junction_bounds(first,second):
    angles=np.unwrap([np.arctan2(first[1],first[0]),np.arctan2(second[1],second[0])])
    # Exclude near-singular normal graphs, while exploring beyond the original
    # corner's angular sector. Both reference dot products stay >= 0.2 |t||v|.
    half=np.arccos(.2)
    return max(angles)-half,min(angles)+half


class EndTangents(Curve):
    def __init__(self,a,b,c,start_tangent=None,end_tangent=None):
        super().__init__(a,b,c,end_tangent)
        if start_tangent is not None:
            t=np.asarray(start_tangent,float)
            if t.shape!=(2,) or not np.all(np.isfinite(t)) or t@self.v<=0:
                raise ValueError('endpoint tangent must point forward along reference')
            coefficients=self.spline.c.copy()
            coefficients[1]=self.length**2*(t@self.n)/(t@self.v)/15
            self.spline=BSpline(KNOTS,coefficients,3)
            self.d1=self.spline.derivative(); self.d2=self.spline.derivative(2)
            pp=PPoly.from_spline(self.spline)
            self.spans=[(float(a),float(b-a),P(pp.c[:,j][::-1]))
                        for j,(a,b) in enumerate(zip(pp.x[:-1],pp.x[1:])) if b>a]


class IMPDirect(TangentMultiSurface):
    def curve(self,name,c):
        a,b=self.targets[name][[0,-1]]
        if name not in ACTIVE: return Curve(a,b,c)
        s=self.targets['imp_sloped']; v=self.targets['imp_vertical']
        first=s[-1]-s[0]; second=v[-1]-v[0]
        lower,upper=junction_bounds(first,second)
        theta=(1-self.junction_fraction)*lower+self.junction_fraction*upper
        common=np.array([np.cos(theta),np.sin(theta)])
        return EndTangents(a,b,c,first if name=='imp_sloped' else common,
                           common if name=='imp_sloped' else second)


class IMPEvaluator(TangentEvaluator):
    def set_junction(self,fraction):
        if not np.isfinite(fraction) or not 0<=fraction<=1:
            raise ValueError('junction fraction must be finite and in [0,1]')
        if self.direct.base.junction_fraction!=float(fraction):
            self.cache.clear()
            self.direct.base.junction_fraction=float(fraction)

    def constraints(self,c):
        c=np.asarray(c,float)
        if c.shape!=(self.dimension,) or not np.all(np.isfinite(c)) or np.max(abs(c))>self.bound+1e-12:
            raise ValueError('finite bounded full coefficients required')
        curves=self.curves(c); k={n:curves[n].max_curvature() for n in ACTIVE}
        s,v=curves['imp_sloped'],curves['imp_vertical']; errors=[]
        for a,b in [(s.tangent(0),s.v),(s.tangent(1),v.tangent(0)),(v.tangent(1),v.v)]:
            errors.append(float(np.degrees(np.arctan2(abs(a[0]*b[1]-a[1]*b[0]),a@b))))
        opposite_zero=bool(np.all(c[:16]==0))
        return dict(passed=bool(max(k.values())<=1/self.rmin+1e-7 and max(errors)<1e-8 and opposite_zero),
                    radius_minimum_m=self.rmin,curvature_max_1_m=k,
                    actual_minimum_radius_m={n:1/x if x>0 else None for n,x in k.items()},
                    tangent_angle_errors_deg=errors,junction_fraction=self.direct.base.junction_fraction,
                    junction_angle_bounds_deg=np.degrees(junction_bounds(s.v,v.v)).tolist(),
                    opposite_coefficients_zero=opposite_zero)


def imp_evaluator(original,rmin,backend='cpu',wall_step=None):
    if not np.isfinite(rmin) or rmin<=0 or backend not in ('cpu','gpu'):
        raise ValueError('positive radius and cpu/gpu backend required')
    if original.metadata['source']!='IMP': raise ValueError('IMP source required')
    e=object.__new__(IMPEvaluator); e.__dict__.update(original.__dict__)
    e.direct=object.__new__(type(original.direct)); e.direct.__dict__.update(original.direct.__dict__)
    e.direct.base=copy(original.direct.base); e.direct.base.__class__=IMPDirect
    e.direct.base.junction_fraction=.5
    e.direct.boundary_guard=copy(original.direct.boundary_guard)
    e.direct.boundary_guard.__class__=TangentGuard; e.direct.boundary_guard.e=e.direct.base
    e.cache=OrderedDict(); e.rmin=float(rmin); e.backend=backend
    if wall_step is not None: e.wall_step=float(wall_step)
    sloped=e.targets['imp_sloped']; vertical=e.targets['imp_vertical']
    if np.linalg.norm(sloped[-1]-vertical[0])>1e-12: raise ValueError('IMP surfaces must share the engineering junction')
    for name,u in [('imp_sloped',0),('imp_vertical',1)]:
        s=next(s for s in e.surfaces if s['name']==name); i=s['segment_index']
        fixed=e.wall[i]-e.wall[i-1] if u==0 else e.wall[i+2]-e.wall[i+1]
        direction=e.targets[name][-1]-e.targets[name][0]
        if fixed@direction<=0 or abs(fixed[0]*direction[1]-fixed[1]*direction[0])>1e-12:
            raise ValueError('IMP outer endpoint must be collinear with adjoining fixed surface')
    return e


def expand_imp(c):
    c=np.asarray(c,float)
    if c.shape!=(8,): raise ValueError('eight IMP coefficients required')
    return np.r_[np.zeros(16),c]


def project_imp(e,target,fraction=.5,maxiter=120,return_best=False):
    """Smooth seed projection, followed by exact continuous acceptance."""
    target=np.asarray(target,float)
    if target.shape!=(8,): raise ValueError('eight target coefficients required')
    u=np.linspace(0,1,501); basis=BSpline(KNOTS,np.eye(8),3)
    b1=basis.derivative()(u)[:,2:6]; b2=basis.derivative(2)(u)[:,2:6]
    def margins(y):
        e.set_junction(y[8]); result=[]
        for j,name in enumerate(ACTIVE):
            q=e.direct.base.curve(name,np.zeros(4)); x=y[j*4:j*4+4]*.02
            first=q.d1(u)+b1@x; second=q.d2(u)+b2@x
            k=q.length*second/(q.length**2+first**2)**1.5
            length=np.trapezoid(np.sqrt(q.length**2+first**2),u)/q.length
            result.extend([.995-e.rmin*k,.995+e.rmin*k,np.array([1.1999-length])])
        return np.concatenate(result)
    y0=np.r_[np.clip(target,-e.bound,e.bound)/.02,fraction]
    result=minimize(lambda y:float(np.sum((y[:8]-target/.02)**2)+.2*(y[8]-fraction)**2),
                    y0,method='SLSQP',bounds=[(-e.bound/.02,e.bound/.02)]*8+[(0,1)],
                    constraints=[dict(type='ineq',fun=margins)],options=dict(maxiter=maxiter,ftol=1e-9))
    x=np.r_[result.x[:8]*.02,result.x[8]]; e.set_junction(x[8]); c=expand_imp(x[:8])
    if e.constraints(c)['passed'] and max(e.precise_lengths(c)[0].values())<=1.2:
        return x
    return x if return_best else None
