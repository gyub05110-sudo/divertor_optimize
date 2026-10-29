"""Continuous cubic geometry with an optional oriented endpoint tangent.

Stationary polynomial roots cover entire spans, including both one-sided knot
limits; the wall mesh is not used to enforce curvature or lip clearance.
"""
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.interpolate import BSpline, PPoly

KNOTS=np.r_[np.zeros(4),np.arange(1,5)/5,np.ones(4)]


def real_roots(poly, width):
    roots=poly.trim(tol=1e-15).roots()
    return [float(x.real) for x in roots if abs(x.imag)<1e-8 and -1e-12<=x.real<=width+1e-12]


class Curve:
    def __init__(self,a,b,c,end_tangent=None):
        self.a=np.asarray(a,float); self.b=np.asarray(b,float); c=np.asarray(c,float)
        if self.a.shape!=(2,) or self.b.shape!=(2,) or c.shape!=(4,) or not np.all(np.isfinite(np.r_[self.a,self.b,c])):
            raise ValueError('finite endpoints and four coefficients required')
        self.v=self.b-self.a; self.length=float(np.linalg.norm(self.v))
        if self.length<=0: raise ValueError('nonzero segment required')
        self.n=np.array([-self.v[1],self.v[0]])/self.length
        slope=0.
        if end_tangent is not None:
            t=np.asarray(end_tangent,float)
            if t.shape!=(2,) or not np.all(np.isfinite(t)) or t@self.v<=0:
                raise ValueError('endpoint tangent must point forward along reference')
            slope=self.length**2*(t@self.n)/(t@self.v)
        self.spline=BSpline(KNOTS,np.r_[0.,0.,c,-slope/15,0.],3)
        self.d1=self.spline.derivative(); self.d2=self.spline.derivative(2)
        pp=PPoly.from_spline(self.spline)
        self.spans=[(float(a),float(b-a),P(pp.c[:,j][::-1])) for j,(a,b) in enumerate(zip(pp.x[:-1],pp.x[1:])) if b>a]

    def point(self,u):
        u=np.asarray(u); return self.a+u[...,None]*self.v+self.spline(u)[...,None]*self.n

    def tangent(self,u):
        return self.v+np.asarray(self.d1(u))[...,None]*self.n

    def curvature(self,u):
        return self.length*self.d2(u)/(self.length**2+self.d1(u)**2)**1.5

    def max_curvature(self):
        maximum=0.
        for a,width,d in self.spans:
            p=d.deriv(); q=p.deriv(); r=q.deriv()
            station=r*(self.length**2+p*p)-3*p*q*q
            x=np.clip([0.,width]+real_roots(station,width),0,width)
            maximum=max(maximum,float(np.max(abs(self.length*q(x)/(self.length**2+p(x)**2)**1.5))))
        return maximum

    def minimum_distance(self,point):
        point=np.asarray(point,float); best=(np.inf,0.)
        for a,width,d in self.spans:
            xy=[P([self.a[j]-point[j]+a*self.v[j],self.v[j]])+self.n[j]*d for j in range(2)]
            station=xy[0]*xy[0].deriv()+xy[1]*xy[1].deriv()
            x=np.clip([0.,width]+real_roots(station,width),0,width)
            d2=xy[0](x)**2+xy[1](x)**2; k=int(np.argmin(d2))
            if d2[k]<best[0]: best=(float(d2[k]),float(a+x[k]))
        return float(np.sqrt(best[0])),best[1]

    def length_ratio(self,order=32):
        x,w=np.polynomial.legendre.leggauss(order)
        u=np.concatenate([a+.1*(x+1) for a in np.arange(5)/5])
        return float(np.tile(w*.1,5)@np.sqrt(self.length**2+self.d1(u)**2)/self.length)
