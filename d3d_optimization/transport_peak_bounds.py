"""Directed enclosures of infinite periodic Gaussian mixtures.

The inputs are frozen constants or explicitly assumed parameter intervals.
This module does NOT certify field tracing, thermal widths, shape feasibility,
or coverage of a geometry parameter box. No numerical sampled extrema become
certified transport-parameter intervals.
"""
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING
import copy
import hashlib
import json
from .optimality_bounds import pi_bounds

DOWN=Context(prec=42,rounding=ROUND_FLOOR)
UP=Context(prec=42,rounding=ROUND_CEILING)
SCHEMA='periodic-kernel-peak-enclosure.v1'
FAR_EXP_UPPER=UP.next_plus(UP.exp(D(-100)))


def decimal(x):
    value=D.from_float(x) if isinstance(x,float) else D(x)
    if not value.is_finite():raise ValueError('finite decimal input required')
    return value


class I:
    __slots__=('lo','hi')

    def __init__(self,lo,hi=None):
        self.lo=decimal(lo);self.hi=decimal(lo if hi is None else hi)
        if self.lo>self.hi:raise ValueError('reversed interval')

    def __add__(self,other):
        other=as_interval(other)
        return I(DOWN.add(self.lo,other.lo),UP.add(self.hi,other.hi))

    __radd__=__add__

    def __neg__(self):return I(self.hi.copy_negate(),self.lo.copy_negate())
    def __sub__(self,other):return self+-as_interval(other)
    def __rsub__(self,other):return as_interval(other)+-self

    def __mul__(self,other):
        other=as_interval(other)
        pairs=[(a,b) for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return I(min(DOWN.multiply(a,b) for a,b in pairs),max(UP.multiply(a,b) for a,b in pairs))

    __rmul__=__mul__

    def __truediv__(self,other):
        other=as_interval(other)
        if other.lo<=0<=other.hi:raise ValueError('division interval includes zero')
        return self*I(DOWN.divide(D(1),other.hi),UP.divide(D(1),other.lo))

    def square(self):
        low=D(0) if self.lo<=0<=self.hi else min(self.lo.copy_abs(),self.hi.copy_abs())
        high=max(self.lo.copy_abs(),self.hi.copy_abs())
        return I(DOWN.multiply(low,low),UP.multiply(high,high))

    def sqrt(self):
        if self.lo<0:raise ValueError('negative square root')
        return I(max(D(0),DOWN.next_minus(DOWN.sqrt(self.lo))),UP.next_plus(UP.sqrt(self.hi)))

    def exp(self):
        # Decimal exp and sqrt are correctly rounded to nearest; the adjacent
        # representable numbers therefore enclose them regardless of rounding mode.
        return I(max(D(0),DOWN.next_minus(DOWN.exp(self.lo))),UP.next_plus(UP.exp(self.hi)))


def as_interval(value):return value if isinstance(value,I) else I(value)


def rational_interval(lo,hi):
    return I(DOWN.divide(D(lo.numerator),D(lo.denominator)),
             UP.divide(D(hi.numerator),D(hi.denominator)))


class KernelBounds:
    def __init__(self,model):
        self.period=I(*model['period_m2'])
        if self.period.lo<=0:raise ValueError('positive period required')
        if not model['donors']:raise ValueError('at least one donor required')
        self.donors=[];self.tail_upper=D(0);self.fixed=self.period.lo==self.period.hi
        pi=rational_interval(*pi_bounds());normalization=(2*pi).sqrt()
        for row in model['donors']:
            if ('center_fraction' in row)==('center_m2' in row):
                raise ValueError('exactly one center representation required')
            weight=I(*row['weight_MW']);sigma=I(*row['sigma_m2'])
            center=I(*row.get('center_fraction',row.get('center_m2')))
            mu=center if 'center_fraction' in row else center/self.period
            if 'center_m2' in row:
                if center.lo<0 or center.hi>self.period.lo:
                    raise ValueError('all area centers must lie within every allowed period')
                # Validation proves the true quotient is in [0,1]. Intersect
                # outward division with that exact domain, including c=A.
                mu=I(max(D(0),mu.lo),min(D(1),mu.hi))
            if weight.lo<0 or sigma.lo<=0 or mu.lo<0 or mu.hi>1:
                raise ValueError('nonnegative weight, positive sigma, center fraction in [0,1] required')
            self.fixed &= all(v.lo==v.hi for v in (weight,center,sigma))
            if weight.hi==0:continue
            amplitude=weight/(normalization*sigma)
            wraps=int(UP.divide(UP.multiply(D(10),sigma.hi),self.period.lo).to_integral_value(rounding=ROUND_CEILING))+1
            if wraps>128:raise ValueError('image-series budget exceeded; broad Fourier verifier required')
            # For t,mu in [0,1], omitted images |k|>=K+1 have
            # distance >= K*P. Their successive exponential ratios are <= r.
            ratio=I(self.period.lo)/I(sigma.hi)
            first=(-I('.5')*(wraps*ratio).square()).exp()
            decay=(-I('.5')*(2*wraps+1)*ratio.square()).exp()
            tail=2*I(amplitude.hi)*I(first.hi)/(1-I(decay.hi))
            self.tail_upper=UP.add(self.tail_upper,tail.hi)
            uniform=None
            if DOWN.divide(sigma.lo,self.period.hi)>=D('.12'):
                # Fourier series: sum_{n>=1} exp(-a*n*n) <= exp(-a)/(1-exp(-3a)).
                # |cos|<=1 gives a spatially uniform enclosure, especially tight
                # for very broad footprints where image-wise intervals are loose.
                exponent=-2*pi.square()*(sigma/self.period).square()
                mass=I(exponent.exp().hi)/(1-I((3*exponent).exp().hi))
                uniform=weight/self.period*(1+I(-1,1)*2*mass)
            self.donors.append((mu,sigma,amplitude,wraps,uniform,tail.hi))
        if not self.donors:raise ValueError('positive total upper power required')

    def density(self,t):
        if t.lo<0 or t.hi>1:raise ValueError('receiving coordinate must be in [0,1]')
        total=I(0)
        for mu,sigma,amplitude,wraps,uniform,tail in self.donors:
            donor=I(0)
            for k in range(-wraps,wraps+1):
                distance=self.period*(t-mu+k)/sigma
                exponent=-I('.5')*distance.square()
                # Rigorous cheap upper for far tails, never discard power.
                if exponent.hi < -100:
                    value=I(0,FAR_EXP_UPPER)
                else:value=exponent.exp()
                donor+=amplitude*value
            donor+=I(0,tail)
            if uniform is not None:donor=I(max(donor.lo,uniform.lo),min(donor.hi,uniform.hi))
            total+=donor
        return total


def certify(model,tolerance='0.00001',max_splits=1000):
    tolerance=decimal(tolerance)
    if tolerance<=0 or not isinstance(max_splits,int) or not 0<=max_splits<=10000:
        raise ValueError('positive tolerance and split budget 0..10000 required')
    evaluator=KernelBounds(model)
    def midpoint(a,b):
        # Endpoints are dyadic with <=10000 divisions. Bound brackets remain
        # exhaustive even if midpoint rounds: the SAME decimal joins both sides.
        return DOWN.divide(DOWN.add(a,b),D(2))
    def leaf(identity,a,b,parent_upper=None):
        value=evaluator.density(I(a,b))
        return {'id':identity,'interval':[str(a),str(b)],
                'upper':str(min(value.hi,parent_upper) if parent_upper is not None else value.hi)}
    leaves=[leaf(1,D(0),D(1))]
    strongest=max(evaluator.donors,key=lambda row:row[2].hi)[0]
    probes=[D(0),D('.5'),midpoint(strongest.lo,strongest.hi)]
    lower=max(evaluator.density(I(t)).lo for t in probes)
    history=[];stop_reason='split_budget'
    for iteration in range(max_splits+1):
        active=max(leaves,key=lambda row:(D(row['upper']),-row['id']))
        upper=D(active['upper']);gap=UP.subtract(upper,lower)
        history.append({'splits':iteration,'lower':str(lower),'upper':str(upper)})
        if gap<=tolerance:stop_reason='tolerance';break
        if iteration==max_splits:break
        a,b=map(D,active['interval']);mid=midpoint(a,b)
        if not a<mid<b:stop_reason='coordinate_precision';break
        lower=max(lower,evaluator.density(I(mid)).lo)
        leaves.remove(active)
        leaves.extend([leaf(2*active['id'],a,mid,upper),leaf(2*active['id']+1,mid,b,upper)])
    frozen=copy.deepcopy(model)
    return {'schema':SCHEMA,'model':frozen,
            'model_sha256':hashlib.sha256(json.dumps(frozen,sort_keys=True,allow_nan=False).encode()).hexdigest(),
            'tolerance':str(tolerance),'max_splits':max_splits,'stop_reason':stop_reason,
            'scope':'fixed infinite periodic kernel' if evaluator.fixed else 'conditional interval-kernel family',
            'peak_bracket_MW_m2':[str(lower),str(upper)],'gap_upper_MW_m2':str(gap),
            'tolerance_met':gap<=tolerance,'tail_upper_MW_m2':str(evaluator.tail_upper),
            'leaves':sorted(leaves,key=lambda row:row['id']),'history':history,
            'geometry_global_optimum_proven':False,'transport_parameter_enclosure_proven':False}


def replay(record):
    if record!=certify(record['model'],record['tolerance'],record['max_splits']):
        raise ValueError('kernel certificate replay mismatch')
    return True
