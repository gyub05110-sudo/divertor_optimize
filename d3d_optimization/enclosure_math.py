"""Directed arithmetic for explicitly frozen mathematical models.

No binary floating-point function is used as an interval proof oracle.
"""
from decimal import Decimal as D,InvalidOperation
from fractions import Fraction
from math import factorial
from .transport_peak_bounds import I, DOWN, UP, rational_interval


def interval(value):
    try:
        result=value if isinstance(value,I) else I(*value) if isinstance(value,(list,tuple)) else I(value)
    except (TypeError,InvalidOperation) as exc:
        raise ValueError('finite scalar or interval required') from exc
    if not result.lo.is_finite() or not result.hi.is_finite():
        raise ValueError('finite interval required')
    return result


def record(x):
    return [str(x.lo),str(x.hi)]


def hull(values):
    values=list(values)
    if not values:raise ValueError('empty enclosure')
    return I(min(v.lo for v in values),max(v.hi for v in values))


def intersect(a,b):
    lo=max(a.lo,b.lo);hi=min(a.hi,b.hi)
    return I(lo,hi) if lo<=hi else None


def power(a,n):
    if n<0 or not isinstance(n,int):raise ValueError('nonnegative integer power required')
    out=I(1)
    for _ in range(n):out=out*a
    return out


def positive_root(a,n):
    """Root bound by directed sign tests; uncertain signs retain a bracket."""
    a=interval(a)
    if a.lo<0 or not isinstance(n,int) or n<1:raise ValueError('nonnegative root argument required')
    def endpoint(x,upper):
        if not x:return D(0)
        lo=D(0);hi=D(10)**((x.adjusted()//n)+1)
        # The decimal power above is an exact power of ten.
        for _ in range(180):
            mid=DOWN.divide(DOWN.add(lo,hi),D(2))
            if not lo<mid<hi:break
            y=power(I(mid),n)
            if y.hi<x:lo=mid
            elif y.lo>x:hi=mid
            else:
                # Exact comparison settles the ambiguous rounded residual.
                if Fraction(mid)**n <= Fraction(x):lo=mid
                else:hi=mid
        return hi if upper else lo
    return I(endpoint(a.lo,False),endpoint(a.hi,True))


def sin_cos(x):
    x=interval(x)
    if x.lo < -8 or x.hi>8:raise ValueError('angle enclosure requires [-8,8] radians')
    mid=DOWN.divide(DOWN.add(x.lo,x.hi),D(2))
    m=Fraction(mid)
    # Taylor degree 81/80, with Lagrange remainder using |f^(n)|<=1.
    s=sum(((-1)**k*m**(2*k+1)/factorial(2*k+1) for k in range(41)),Fraction())
    c=sum(((-1)**k*m**(2*k)/factorial(2*k) for k in range(41)),Fraction())
    rs=abs(m)**82/factorial(82);rc=abs(m)**81/factorial(81)
    radius=max(UP.subtract(mid,x.lo),UP.subtract(x.hi,mid))
    def widen(a,r):
        q=rational_interval(a-r,a+r)+I(radius.copy_negate(),radius)
        return intersect(q,I(-1,1))
    return widen(s,rs),widen(c,rc)


def derivative(knots,controls,degree):
    t=[interval(x) for x in knots]
    if any(x.lo!=x.hi for x in t):raise ValueError('exact spline knots required')
    if degree<1:raise ValueError('positive spline degree required')
    result=[]
    for i in range(len(controls)-1):
        den=t[i+degree+1]-t[i+1]
        if den.lo<=0:raise ValueError('unsupported repeated interior derivative knot')
        result.append(degree*(controls[i+1]-controls[i])/den)
    return knots[1:-1],result,degree-1


def spline_range(knots,controls,degree,u):
    bounds=[interval(x) for x in knots]
    if any(x.lo!=x.hi for x in bounds):raise ValueError('exact spline knots required')
    t=[x.lo for x in bounds];u=interval(u)
    c=[interval(x) for x in controls]
    if (len(t)!=len(c)+degree+1 or degree<0 or any(a>b for a,b in zip(t,t[1:]))
            or u.lo<t[degree] or u.hi>t[len(c)]):
        raise ValueError('invalid spline or out-of-domain parameter')
    values=[]
    for span in range(degree,len(c)):
        if t[span]>=t[span+1] or u.hi<t[span] or u.lo>t[span+1]:continue
        v=I(max(u.lo,t[span]),min(u.hi,t[span+1]))
        work=c[span-degree:span+1]
        for r in range(1,degree+1):
            for j in range(degree,r-1,-1):
                i=span-degree+j
                alpha=(v-t[i])/(I(t[i+degree-r+1])-t[i])
                work[j]=(1-alpha)*work[j-1]+alpha*work[j]
        # B-spline positivity provides an independent convex-hull enclosure.
        values.append(intersect(work[degree],hull(c[span-degree:span+1])))
    return hull(values)


def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),I(0))
def norm(a):return sum((v.square() for v in a),I(0)).sqrt()


class SplineField:
    """Exact exported tensor B-spline coefficients, Z rows / R columns."""
    def __init__(self,model):
        self.r=list(model['r_knots']);self.z=list(model['z_knots'])
        if any(interval(x).lo!=interval(x).hi for x in self.r+self.z):
            raise ValueError('exact spline knots required')
        self.kr=int(model['degree_r']);self.kz=int(model['degree_z'])
        self.c=[[interval(x) for x in row] for row in model['coefficients']]
        if (len(self.c)!=len(self.z)-self.kz-1 or
                any(len(row)!=len(self.r)-self.kr-1 for row in self.c)):
            raise ValueError('tensor coefficient shape mismatch')
        self.K=interval(model['toroidal_K'])
        self.dr=[derivative(self.r,row,self.kr)[1] for row in self.c]
        self.dz=[]
        for i in range(len(self.c)-1):
            den=interval(self.z[i+self.kz+1])-interval(self.z[i+1])
            if den.lo<=0:raise ValueError('invalid field derivative knots')
            self.dz.append([self.kz*(b-a)/den for a,b in zip(self.c[i],self.c[i+1])])

    def _tensor(self,c,r,z,kr,kz,p):
        # Evaluate only Z basis-support rows needed by the requested range.
        values=[];tz=[interval(x).lo for x in z];zz=p[1]
        for span in range(kz,len(c)):
            if tz[span]>=tz[span+1] or zz.hi<tz[span] or zz.lo>tz[span+1]:continue
            local=[spline_range(r,c[j],kr,p[0]) for j in range(span-kz,span+1)]
            work=local[:];v=I(max(zz.lo,tz[span]),min(zz.hi,tz[span+1]))
            for q in range(1,kz+1):
                for j in range(kz,q-1,-1):
                    i=span-kz+j
                    alpha=(v-tz[i])/(I(tz[i+kz-q+1])-tz[i])
                    work[j]=(1-alpha)*work[j-1]+alpha*work[j]
            values.append(intersect(work[kz],hull(local)))
        return hull(values)

    def gradient(self,p):
        p=list(map(interval,p))
        for x,t,k in zip(p,[self.r,self.z],[self.kr,self.kz]):
            if x.lo<interval(t[k]).lo or x.hi>interval(t[-k-1]).hi:
                raise ValueError('field enclosure outside grid')
        return [self._tensor(self.dr,self.r[1:-1],self.z,self.kr-1,self.kz,p),
                self._tensor(self.dz,self.r,self.z[1:-1],self.kr,self.kz-1,p)]

    def quantities(self,p):
        p=list(map(interval,p))
        if p[0].lo<=0:raise ValueError('nonpositive R')
        g=norm(self.gradient(p))
        if g.lo<=0:raise ValueError('gradient enclosure includes zero')
        magnetic=(g.square()+self.K.square()).sqrt()/p[0]
        return g,magnetic
