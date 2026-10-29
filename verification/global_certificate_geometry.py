"""Small auditable whole-box enclosures and exact rational geometry.

Only basic finite binary64 operations are enclosed. No ODE, spline-field,
GEOS predicate, transcendental or physical-model certification is implied.
"""
from fractions import Fraction as F
import math
import numpy as np

ROUND_PAD=F(1,2**44)

def outward(values):
    if not all(math.isfinite(v) for v in values):raise ValueError('finite interval arithmetic required')
    return math.nextafter(min(values),-math.inf),math.nextafter(max(values),math.inf)

def iadd(a,b):return outward([a[0]+b[0],a[1]+b[1]])

def imul(a,b):return outward([x*y for x in a for y in b])

def displacement_interval(x, coefficients):
    """Enclose the actual de Boor recurrence over all coefficient intervals.

    alpha and 1-alpha are coefficient-independent rounded constants, just as
    in the production code. Every multiply/add involving a coefficient is
    outward enclosed. The active count is the production four-control case.
    """
    if len(coefficients)!=4 or not 0<=x<=1:raise ValueError('four-control spline on [0,1] required')
    if any(not math.isfinite(v) for pair in coefficients for v in pair) or any(a>b for a,b in coefficients):
        raise ValueError('finite ordered coefficient intervals required')
    if x in (0.,1.):return (0.,0.)
    knots=np.r_[np.zeros(4),np.linspace(0.,1.,6)[1:-1],np.ones(4)]
    full=[(0.,0.)]*2+list(coefficients)+[(0.,0.)]*2
    span=min(max(int(np.searchsorted(knots,x,side='right')-1),3),7)
    work=list(full[span-3:span+1])
    for recursion in range(1,4):
        for local in range(3,recursion-1,-1):
            index=span-3+local
            denominator=knots[index+4-recursion]-knots[index]
            alpha=float(0. if denominator==0 else (x-knots[index])/denominator)
            beta=1.-alpha
            work[local]=iadd(imul((beta,beta),work[local-1]),imul((alpha,alpha),work[local]))
    return work[3]

def point(p):return tuple(v if isinstance(v,F) else F(float(v)) for v in p)
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]

def hull(points):
    ps=sorted(set(point(p) for p in points))
    if len(ps)<=1:return ps
    def half(seq):
        result=[]
        for p in seq:
            while len(result)>1 and cross(sub(result[-1],result[-2]),sub(p,result[-1]))<=0:result.pop()
            result.append(p)
        return result
    return half(ps)[:-1]+half(reversed(ps))[:-1]

def contains(poly,p):
    p=point(p)
    if len(poly)<3:
        return p in poly or (len(poly)==2 and cross(sub(poly[1],poly[0]),sub(p,poly[0]))==0
                            and all(min(a,b)<=v<=max(a,b) for a,b,v in zip(poly[0],poly[1],p)))
    return all(cross(sub(b,a),sub(p,a))>=0 for a,b in zip(poly,poly[1:]+poly[:1]))

def gap(a,b,axis):
    axis=point(axis); pa=[dot(p,axis) for p in a];pb=[dot(p,axis) for p in b]
    return max(min(pb)-max(pa),min(pa)-max(pb))

def check_separator(a,b,axis):return bool(a and b and gap(a,b,axis)>0)

def separator(a,b):
    # Floating axes are only proposals. Their projections and final comparisons
    # are exact rational operations; an inaccurate axis cannot create a proof.
    axes=[(1.,0.),(0.,1.)]
    # Try the small edge envelope first, before the potentially long path hull.
    for poly in (b,a):
        for p,q in zip(poly,poly[1:]+poly[:1]):
            d=sub(q,p);axes.append((float(d[1]),-float(d[0])))
    af=np.asarray(a,dtype=float);bf=np.asarray(b,dtype=float)
    for axis in axes:
        if axis==(0.,0.):continue
        pa=af@axis;pb=bf@axis
        # This filter can only miss a certificate. It cannot authorize one:
        # every accepted witness still passes the exact rational predicate.
        if max(min(pb)-max(pa),min(pa)-max(pb))>0 and check_separator(a,b,axis):return axis
    return None

def node_envelope(position,normal,interval):
    """Enclose fl(position + fl(d*normal)) for all enclosed d.

    With |position|<=4, |normal|<=2, |d|<=1/4, the multiply/add roundoff is
    below 2^-49 in either coordinate under IEEE nearest rounding, including
    subnormal cases. The larger explicit 2^-44 pad is used here. Positions,
    normals and d are the frozen outputs of coefficient-independent operations.
    """
    if (max(abs(float(v)) for v in position)>4 or max(abs(float(v)) for v in normal)>2
            or max(abs(float(v)) for v in interval)>.25):
        raise ValueError('roundoff-envelope magnitude assumptions violated')
    p,n=point(position),point(normal)
    vertices=[]
    for d in map(F,interval):
        q=[p[j]+d*n[j] for j in range(2)]
        vertices += [(q[0]+sx*ROUND_PAD,q[1]+sy*ROUND_PAD) for sx in (-1,1) for sy in (-1,1)]
    return hull(vertices)

def segment_intersection(p,q,a,b):
    p,q,a,b=map(point,(p,q,a,b))
    if any(max(min(p[j],q[j]),min(a[j],b[j]))>min(max(p[j],q[j]),max(a[j],b[j])) for j in range(2)):
        return None
    r,s=sub(q,p),sub(b,a);den=cross(r,s);ap=sub(a,p)
    if den==0:
        if cross(ap,r)==0:raise ValueError('collinear overlap is not certified')
        return None
    t,u=cross(ap,s)/den,cross(ap,r)/den
    if 0<=t<=1 and 0<=u<=1:return t,(p[0]+t*r[0],p[1]+t*r[1])
    return None

def first_fixed_intersection(path,fixed):
    path=[point(p) for p in path]
    fixed=[(int(i),[point(p) for p in ends]) for i,ends in fixed]
    for j,(p,q) in enumerate(zip(path[:-1],path[1:])):
        hits=[]
        for i,(a,b) in fixed:
            hit=segment_intersection(p,q,a,b)
            if hit is not None:hits.append((hit[0],i,hit[1]))
        if hits:
            t,i,location=min(hits)
            return dict(fixed_index=i,path_edge=j,path_fraction=t,point=location,prefix=path[:j+1]+[location])
    raise ValueError('stored path has no fixed-wall intersection')
