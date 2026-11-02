"""Positive prefix integrals and monotone thermal closure enclosures."""
from fractions import Fraction as F
from .enclosure_math import (I,D,DOWN,UP,interval,hull,intersect,norm,power,positive_root)
from .reference_transport import E_CHARGE,M_DEUTERIUM


class PathThermal:
    def __init__(self,model,last_edge):
        p=model['path']
        if not 0<=last_edge<len(p)-1:raise ValueError('invalid last path edge')
        self.h=[norm([interval(b)-interval(a) for a,b in zip(p[j],p[j+1])]) for j in range(last_edge+1)]
        if any(x.lo<=0 for x in self.h):raise ValueError('nonpositive path edge')
        self.metric=[interval(x) for x in model['metric'][:last_edge+2]]
        self.integrand=[interval(x) for x in model['integrand'][:last_edge+2]]
        self.mapping=[interval(x) for x in model['mapping'][:last_edge+2]]
        if any(len(v)!=last_edge+2 for v in (self.metric,self.integrand,self.mapping)):
            raise ValueError('missing prefix thermal node')
        if any(x.lo<1 for x in self.metric) or any(x.lo<=0 for x in self.integrand+self.mapping):
            raise ValueError('nonpositive or invalid prefix metric')
        self.L=[I(0)];self.J=[I(0)]
        for j,h in enumerate(self.h):
            self.L.append(self.L[-1]+h*(self.metric[j]+self.metric[j+1])/2)
            self.J.append(self.J[-1]+h*(self.integrand[j]+self.integrand[j+1])/2)

    def at(self,j,t):
        t=interval(t)
        if not 0<=j<len(self.h) or t.lo<0 or t.hi>1:raise ValueError('invalid path fraction')
        def integral(values,prefix):
            def endpoint(x):
                x=I(x)
                return prefix[j]+self.h[j]*(x*(1-x/2)*values[j]+x.square()*values[j+1]/2)
            return I(endpoint(t.lo).lo,endpoint(t.hi).hi)
        # All inverse thermal coordinates lie before the first hit. Including
        # the complete last segment is a safe relaxation, not a sampled bound.
        return integral(self.metric,self.L),integral(self.integrand,self.J),hull(self.mapping[:j+2])

    def average(self,J,y,E,panels=32):
        """Validated integral of mapping(F^-1(p(u))) over the entire unit interval.

        Enumerate every compatible inverse segment. Positivity bounds the
        fractional inverse without a floating root or sampled-extrema claim.
        """
        if not isinstance(panels,int) or not 1<=panels<=512:raise ValueError('thermal panels 1..512 required')
        S=(y.square()+E).sqrt()+y;ratio=E/S.square()
        cuts=[D(0)]+[DOWN.divide(D(j),D(panels)) for j in range(1,panels)]+[D(1)]
        total=I(0)
        for a,b in zip(cuts,cuts[1:]):
            u=I(a,b)
            partial=intersect(J*(1-u)*(1+u*ratio),I(0,J.hi))
            values=[]
            for j,h in enumerate(self.h):
                if partial.hi<self.J[j].lo or partial.lo>self.J[j+1].hi:continue
                delta=intersect(partial-self.J[j],I(0,(self.J[j+1]-self.J[j]).hi))
                if delta is None:continue
                t=intersect(delta/(h*hull(self.integrand[j:j+2])),I(0,1))
                if t is None:continue
                q=self.mapping[j]+t*(self.mapping[j+1]-self.mapping[j])
                values.append(intersect(q,hull(self.mapping[j:j+2])))
            if not values:raise ArithmeticError('positive inverse integral lost all branches')
            total+=hull(values)*(I(b)-a)
        return intersect(total,hull(self.mapping))


def sheath_root(E,rhs):
    """Enclose every positive root y³+E*y=rhs in a positive input box."""
    E=interval(E);rhs=interval(rhs)
    if E.lo<=0 or rhs.lo<=0:raise ValueError('positive thermal root parameters required')
    def root(e,r):
        lo=D(0);hi=min(positive_root(I(r),3).hi,UP.divide(r,e))
        for _ in range(160):
            mid=DOWN.divide(DOWN.add(lo,hi),D(2))
            if not lo<mid<hi:break
            residual=power(I(mid),3)+I(e)*mid-r
            if residual.hi<0:lo=mid
            elif residual.lo>0:hi=mid
            elif F(mid)**3+F(e)*F(mid)<=F(r):lo=mid
            else:hi=mid
        return I(lo,hi)
    return I(root(E.hi,rhs.lo).lo,root(E.lo,rhs.hi).hi)


def thermal_width(qu,field_ratio,J,average,reference):
    qu,field_ratio,J,average=map(interval,(qu,field_ratio,J,average))
    n,chi,gamma,kappa=[interval(reference[x]) for x in ('density_m3','chi_perp_m2_s','gamma','kappa0')]
    if any(x.lo<=0 for x in (qu,field_ratio,J,average,n,chi,gamma,kappa)):
        raise ValueError('strictly positive finite thermal inputs required')
    e=interval(E_CHARGE);mass=interval(M_DEUTERIUM)
    C=gamma*e*n/2*(2*e/mass).sqrt()
    H=qu*field_ratio/C;E=I('3.5')*qu*J/kappa
    y=sheath_root(E,power(H,3)*H.sqrt())
    upstream=y.square()+E
    factor=positive_root(upstream,7)
    sigma_factor=(3*e*n*chi/kappa).sqrt()*factor*J/(upstream.sqrt()+y)
    sigma=sigma_factor*average
    if sigma.lo<=0:raise ValueError('positive spreading lower bound unresolved')
    return {'sigma_m':sigma,'sigma_factor':sigma_factor,'Tt_eV':power(positive_root(y,7),4),
            'Tu_eV':factor.square(),'target_y':y,'energy':E}
