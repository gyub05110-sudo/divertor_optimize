"""Continuous IMP design boundary and a zero-power SOL-limit peak witness.

This is a numerical peak guard, not a continuous/global maximum certificate.
The archived first-hit, source, and shape kernels are deliberately unchanged.
"""
from copy import deepcopy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from shapely.geometry import LineString, Point
from .frozen_mapping import first_hits
from .shape import NormalDisplacementBSpline, ReferencePolyline


def continuous_design(review, field, psi_sep, buffer_m=0.):
    result = deepcopy(review)
    wall = np.asarray(review['wall'])
    a, b = wall[69:71]
    f = lambda s: float(field.values([a+s*(b-a)])[0]-psi_sep)
    grid = np.linspace(0, 1, 257)
    roots = [brentq(f, l, h, xtol=1e-14) for l,h in zip(grid[:-1],grid[1:]) if f(l)*f(h)<0]
    if len(roots)!=1:
        raise ValueError('expected one continuous IMP strike on original wall 69')
    strike = a+roots[0]*(b-a)
    if np.linalg.norm(strike-np.asarray(review['legs']['inner']['strike_RZ_m']))>1e-4:
        raise ValueError('continuous root is inconsistent with the inner leg')
    start = roots[0]-buffer_m/np.linalg.norm(b-a)
    if start < 0:
        raise ValueError('IMP buffer extends beyond original segment 69')
    va, vb = wall[70:72]
    cutoff = review['optimization_policy']['imp_maximum_Z_m']
    end = va+(cutoff-va[1])/(vb[1]-va[1])*(vb-va)
    for row in result['optimization_footprints']:
        if row['segment']==69: row['points']=[(a+start*(b-a)).tolist(), b.tolist()]
        if row['segment']==70: row['points']=[va.tolist(), end.tolist()]
    return result, strike


def extended_leg(field, xp, psi_sep, source, spacing):
    hessian=np.array([[field.spline.ev(xp[1],xp[0],dy=2),field.spline.ev(xp[1],xp[0],dx=1,dy=1)],
                      [field.spline.ev(xp[1],xp[0],dx=1,dy=1),field.spline.ev(xp[1],xp[0],dx=2)]])
    eigen,basis=np.linalg.eigh(hessian)
    if not eigen[0]<0<eigen[1]: raise ValueError('X-point must be a saddle')
    rays=[]
    for sign in [-1,1]:
        v=np.sqrt(-eigen[0])*basis[:,1]+sign*np.sqrt(eigen[1])*basis[:,0]
        v/=np.linalg.norm(v)
        if v[1]>0:v=-v
        rays.append(v)
    direction=sorted(rays,key=lambda v:v[0])[0 if source=='IMP' else 1]
    start=np.asarray(xp)+.0001*direction
    for _ in range(8):
        g=field.gradient([start])[0]
        start-=(field.values([start])[0]-psi_sep)*g/np.dot(g,g)
    g=field.gradient([start])[0]
    sign=np.sign(np.dot([-g[1],g[0]],direction))
    def rhs(s,p):
        g=field.gradient([p])[0]; t=np.array([-g[1],g[0]])
        return sign*t/np.linalg.norm(t)
    def margin(s,p):
        return min(p[0]-field.r[0],field.r[-1]-p[0],p[1]-field.z[0],field.z[-1]-p[1])-.002
    margin.terminal=True; margin.direction=-1
    sol=solve_ivp(rhs,(0,3),start,method='DOP853',rtol=2e-10,atol=1e-12,
                  max_step=spacing,events=margin,dense_output=True)
    if not sol.success:raise ValueError('separatrix leg integration failed')
    s=np.r_[np.arange(0,sol.t[-1],spacing),sol.t[-1]]
    path=sol.sol(s).T
    if np.max(abs(field.values(path)-psi_sep))>1e-9:
        raise ValueError('separatrix leg flux drift too large')
    return LineString(path)


class BoundaryGuard:
    def __init__(self, evaluator, xp, psi_sep, source, source_rz, spacing):
        self.e=evaluator; self.psi_sep=psi_sep
        self.path=extended_leg(evaluator.field,xp,psi_sep,source,spacing)
        self.h0=np.array([1/(np.linalg.norm(evaluator.field.magnetic([source_rz])[0])*s.normalization)
                          for s in evaluator.sources])

    def witness(self, result, coefficients):
        hit=first_hits([self.path],result['wall'])
        k=int(hit.segment[0])
        if k<0: return dict(resolved=False,reason='no_separatrix_first_hit')
        original=int(result['original_segment'][k])
        surface=next(((j,s) for j,s in enumerate(self.e.surfaces) if s['segment_index']==original),None)
        if surface:
            j,surface=surface
            target=self.e.targets[surface['name']]
            ref=ReferencePolyline([target[0],target[-1]])
            spline=NormalDisplacementBSpline(coefficients[4*j:4*j+4])
            point=lambda s: spline.deform(ref,[s])[0]
            tangent=lambda s: spline.deformed_tangent(ref,[s])[0]
            ids=result['surface_segments'][surface['name']]
            index=k-int(ids[0]); n=len(ids)
            lo=max(0,(index-2)/n); hi=min(1,(index+3)/n)
        else:
            a,b=result['wall'][k:k+2]
            point=lambda s: a+s*(b-a)
            tangent=lambda s: b-a
            lo,hi=0.,1.
        def residual(s): return float(self.e.field.values([point(s)])[0]-self.psi_sep)
        grid=np.linspace(lo,hi,13)
        vals=[residual(s) for s in grid]
        roots=[float(s) for s,v in zip(grid,vals) if abs(v)<1e-14]
        roots += [brentq(residual,l,h,xtol=1e-14) for l,h,u,v in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]) if u*v<0]
        if not roots:return dict(resolved=False,reason='first_hit_flux_root_not_bracketed')
        s=min(roots,key=lambda s:np.linalg.norm(point(s)-hit.points[0]))
        p=point(s); t=tangent(s); normal=np.array([-t[1],t[0]])/np.linalg.norm(t)
        delta=float(np.linalg.norm(p-hit.points[0]))
        # Local root refinement must remain near the observed first intersection.
        if delta>2e-3 or abs(residual(s))>1e-11:
            return dict(resolved=False,reason='first_hit_root_refinement_moved_too_far')
        magnetic=self.e.field.magnetic([p])[0]
        q=self.h0*abs(np.dot(magnetic[:2],normal))
        angle=float(np.degrees(np.arcsin(np.clip(abs(np.dot(magnetic[:2],normal))/np.linalg.norm(magnetic),0,1))))
        floor=surface.get('angle_floor_deg') if surface else None
        return dict(resolved=True,point=p.tolist(),q=q.tolist(),angle_deg=angle,
                    angle_passed=bool(floor is None or angle>=floor),surface=surface['name'] if surface else 'fixed',
                    split_segment=original,psi_residual=residual(s),refinement_distance_m=delta,
                    normal_RZ=normal.tolist(),zero_power_limit=True)


class BoundaryAwareEvaluator:
    def __init__(self, evaluator, xp, psi_sep, source, source_rz, spacing):
        self.base=evaluator
        self.boundary_guard=BoundaryGuard(evaluator,xp,psi_sep,source,source_rz,spacing)
        self.baseline=self.evaluate(np.zeros(evaluator.dimension))
        if not self.baseline['feasible']:
            raise ValueError('continuous-boundary baseline is infeasible')
        self.baseline_peaks=self.baseline['peak_q']

    def __getattr__(self, name):return getattr(self.base,name)

    def evaluate(self, coefficients):
        r=self.base.evaluate(coefficients)
        if not r['geometry_valid']:return r
        w=self.boundary_guard.witness(r,coefficients)
        r['boundary_witness']=w
        peak=np.array([max(m.q) for m in r['maps']])
        if w['resolved']:
            peak=np.maximum(peak,w['q'])
            if not w['angle_passed']:
                r['feasible']=False; r['penalty']+=100*(1.5-w['angle_deg'])
        else:
            r['feasible']=False; r['penalty']+=1000
        r['peak_q']=peak
        return r
