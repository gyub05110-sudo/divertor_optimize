"""Replayable P/area lower bounds for an exact frozen spline-wall model.

This is a relaxation, not a certificate for the floating-point transport code.
All inequalities are rational except square roots, whose rational upper bounds
are checked by squaring. No sampled peak is accepted as a certified upper bound.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import copy
import hashlib
import json
import math
import sys

SCHEMA = 'restricted-area-lower-bound.v1'


def rational(value):
    try:
        return F(value)
    except (ValueError, OverflowError, TypeError) as exc:
        raise ValueError('finite rational value required') from exc


@dataclass(frozen=True, init=False)
class Interval:
    lo: F
    hi: F

    def __init__(self, lo, hi=None):
        lo, hi = rational(lo), rational(lo if hi is None else hi)
        if lo > hi:
            raise ValueError('reversed interval')
        object.__setattr__(self,'lo',lo)
        object.__setattr__(self,'hi',hi)

    def __add__(self, other):
        other = other if isinstance(other,Interval) else Interval(other)
        return Interval(self.lo+other.lo,self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi,-self.lo)

    def __sub__(self, other):
        return self + -(other if isinstance(other,Interval) else Interval(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        other = other if isinstance(other,Interval) else Interval(other)
        products = [a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(products),max(products))

    __rmul__ = __mul__

    def abs_upper(self):
        return max(abs(self.lo),abs(self.hi))


@lru_cache(None)
def pi_bounds():
    def atan_bounds(x, terms):
        total=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
        next_term=(-1)**terms*x**(2*terms+1)/F(2*terms+1)
        return min(total,total+next_term),max(total,total+next_term)
    a,b=atan_bounds(F(1,5),24)
    c,d=atan_bounds(F(1,239),8)
    # Machin identity, with alternating-series remainder bounds.
    return 16*a-4*d,16*b-4*c


def sqrt_upper(value):
    value=rational(value)
    if value < 0:raise ValueError('negative square root')
    if not value:return F(0)
    try:approximate=float(value)
    except OverflowError:approximate=math.inf
    if not math.isfinite(approximate) or approximate<sys.float_info.min:
        # Avoid underflow-driven nextafter loops or overflow. Scale to an
        # approximately 54-bit integer square root, using only exact integers.
        shift=(value.numerator.bit_length()-value.denominator.bit_length())//2-54
        scale=F(2)**shift
        scaled=value/(scale*scale)
        root=math.isqrt(scaled.numerator//scaled.denominator)
        if root*root*scaled.denominator<scaled.numerator:root+=1
        return root*scale
    result=math.sqrt(approximate)
    # Also handles a float conversion that rounds down. The final inequality is exact.
    while F(result)**2 < value:
        result=math.nextafter(result,math.inf)
        if not math.isfinite(result):raise ValueError('square root overflow')
    return F(result)


def spline_range(knots, controls, u, derivative=0):
    """Exact interval de Boor on one positive knot span of a clamped cubic."""
    t=list(map(rational,knots)); c=list(controls); degree=3
    if derivative not in (0,1):raise ValueError('only displacement and first derivative supported')
    if derivative:
        c=[(c[i+1]-c[i])*(F(degree)/(t[i+degree+1]-t[i+1]))
           for i in range(len(c)-1)]
        t=t[1:-1];degree-=1
    middle=(u.lo+u.hi)/2
    span=next((i for i in range(degree,len(c)) if t[i]<=middle<t[i+1]),None)
    if span is None or u.lo<t[span] or u.hi>t[span+1]:
        raise ValueError('interval must lie within one positive knot span')
    work=c[span-degree:span+1]
    for r in range(1,degree+1):
        for j in range(degree,r-1,-1):
            i=span-degree+j
            alpha=(u-t[i])*(1/(t[i+degree-r+1]-t[i]))
            work[j]=(1-alpha)*work[j-1]+alpha*work[j]
    return work[degree]


def validate_model(model):
    wall=model['wall']
    if len(wall)<2 or any(len(p)!=2 for p in wall):raise ValueError('2D wall required')
    for p in wall:
        if rational(p[0])<=0:raise ValueError('positive R anchors required')
        rational(p[1])
    if rational(model['power_MW'])<=0:raise ValueError('positive power required')
    seen=set()
    for q in model['curves']:
        i=q['segment'];t=list(map(rational,q['knots']));c=q['controls']
        if not isinstance(i,int) or i in seen or not 0<=i<len(wall)-1:
            raise ValueError('unique valid curve segment required')
        seen.add(i)
        if len(c)!=8 or len(t)!=12 or t[:4]!=[0]*4 or t[-4:]!=[1]*4 or any(a>=b for a,b in zip(t[3:8],t[4:9])):
            raise ValueError('clamped cubic with four simple interior knots required')
        if len(q['normal'])!=2:raise ValueError('2D frozen normal required')
        for x in c+q['normal']:rational(x)
        if rational(c[0])!=0 or rational(c[-1])!=0:raise ValueError('fixed endpoints required')


def _curve_area(wall,q,control,domain,subcells):
    a=list(map(rational,wall[q['segment']]))
    b=list(map(rational,wall[q['segment']+1]))
    v=[y-x for x,y in zip(a,b)];n=list(map(rational,q['normal']))
    controls=[Interval(x) for x in q['controls']]
    if control is not None:controls[control]=domain
    knots=list(map(rational,q['knots']));area=F(0)
    for left,right in zip(knots[3:8],knots[4:9]):
        for j in range(subcells):
            u=Interval(left+(right-left)*F(j,subcells),left+(right-left)*F(j+1,subcells))
            d=spline_range(knots,controls,u)
            dp=spline_range(knots,controls,u,1)
            radial=a[0]+u*v[0]+n[0]*d
            if radial.hi<=0:raise ValueError('no positive R enclosed')
            speed=sqrt_upper(sum(((v[k]+n[k]*dp).abs_upper()**2 for k in range(2)),F(0)))
            area+=2*pi_bounds()[1]*radial.hi*speed*(u.hi-u.lo)
    return area


def area_upper(model, variable, domain, subcells=4):
    validate_model(model)
    if not isinstance(subcells,int) or not 1<=subcells<=128:raise ValueError('subcells must be 1..128')
    curves={q['segment']:q for q in model['curves']}
    if variable is not None:
        if len(variable)!=2 or variable[0] not in curves or variable[1] not in range(2,6):
            raise ValueError('variable must identify an interior spline control')
        if not isinstance(domain,Interval):raise ValueError('interval domain required')
    area=F(0)
    for i,(a,b) in enumerate(zip(model['wall'],model['wall'][1:])):
        if i in curves:
            control=variable[1] if variable is not None and variable[0]==i else None
            area+=_curve_area(model['wall'],curves[i],control,domain,subcells)
        else:
            a=list(map(rational,a));b=list(map(rational,b))
            length=sqrt_upper(sum(((y-x)**2 for x,y in zip(a,b)),F(0)))
            area+=pi_bounds()[1]*(a[0]+b[0])*length
    if area<=0:raise ValueError('positive receiving area required')
    return area


def build_bounds(model, variable, domain, splits=16, subcells=4):
    """Best-bound bisection; retain ALL children, no numerical pruning."""
    domain=Interval(*domain)
    if domain.lo==domain.hi:raise ValueError('nonzero domain width required')
    if not isinstance(splits,int) or not 0<=splits<=256:raise ValueError('splits must be 0..256')
    def node(identity, interval, inherited=F(0)):
        area=area_upper(model,variable,interval,subcells)
        bound=max(inherited,rational(model['power_MW'])/area)
        return {'id':identity,'domain':[str(interval.lo),str(interval.hi)],
                'area_upper_m2':str(area),'lower_bound_MW_m2':str(bound)}
    leaves=[node(1,domain)];history=[]
    for iteration in range(splits+1):
        weakest=min(leaves,key=lambda n:(F(n['lower_bound_MW_m2']),n['id']))
        history.append({'split_count':iteration,'lower_bound_MW_m2':weakest['lower_bound_MW_m2']})
        if iteration==splits:break
        interval=Interval(*weakest['domain']);mid=(interval.lo+interval.hi)/2
        inherited=F(weakest['lower_bound_MW_m2']);identity=weakest['id']
        leaves.remove(weakest)
        leaves.extend([node(2*identity,Interval(interval.lo,mid),inherited),
                       node(2*identity+1,Interval(mid,interval.hi),inherited)])
    frozen=copy.deepcopy(model)
    return {'schema':SCHEMA,'model':frozen,'model_sha256':hashlib.sha256(json.dumps(frozen,sort_keys=True,allow_nan=False).encode()).hexdigest(),
            'variable':list(variable),'domain':[str(domain.lo),str(domain.hi)],
            'splits':splits,'subcells':subcells,'leaves':sorted(leaves,key=lambda n:n['id']),'history':history,
            'scope':'exact frozen continuous spline wall; nonnegative heat flux; exact full power deposition; positive R',
            'global_optimum_proven':False}


def replay(record):
    expected=build_bounds(record['model'],record['variable'],record['domain'],record['splits'],record['subcells'])
    if record!=expected:raise ValueError('bound replay mismatch')
    return True


def compare_numerical(record, peak):
    replay(record)
    if not math.isfinite(peak) or peak<=0:raise ValueError('positive finite numerical peak required')
    lower=F(record['history'][-1]['lower_bound_MW_m2'])
    return {'numerical_candidate_peak_MW_m2':peak,'conditional_model_lower_bound_MW_m2':float(lower),
            'uncertified_difference_MW_m2':float(F(peak)-lower),
            'certified_upper_bound_MW_m2':None,'certified_gap_MW_m2':None,'global_optimum_proven':False,
            'missing_obligations':['rigorous feasibility for the same exact model',
                                   'field tracing, interpolation and heat solver error enclosures',
                                   'continuous heat-flux peak upper enclosure',
                                   'production wall versus exact frozen spline model enclosure']}
