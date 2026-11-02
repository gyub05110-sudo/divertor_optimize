from decimal import Decimal as D,localcontext
from fractions import Fraction as F
import copy
import math
import numpy as np
import pytest
from d3d_optimization import transport_peak_bounds as pb
from d3d_optimization.boundary_transport import wrapped_density


def kernel(center=.37,sigma=.04):
    return {'period_m2':[1.,1.],'donors':[{'weight_MW':[1.,1.],
            'center_fraction':[center,center],'sigma_m2':[sigma,sigma]}]}


def test_decimal_interval_contains_exact_operations():
    x=pb.I('0.1','0.3');y=pb.I('-2','3')
    assert (x*y).lo<=D('-.6') and (x*y).hi>=D('.9')
    assert F((x/pb.I('3')).lo)<=F(1,30)
    with pytest.raises(ValueError):x/pb.I('-1','1')
    assert pb.I(0).exp().lo<=1<=pb.I(0).exp().hi
    assert F(pb.I(.1).lo)==F(.1)==F(pb.I(.1).hi)


def test_independent_high_precision_image_sum_is_enclosed():
    model=kernel(.37,.04);result=pb.KernelBounds(model).density(pb.I(.39))
    with localcontext() as context:
        context.prec=80
        pi=D('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628')
        center=D.from_float(.37);sigma=D.from_float(.04);x=D.from_float(.39)
        reference=sum(((-((x-center+k)/sigma)**2/2).exp() for k in range(-4,5)),D(0))/(sigma*(2*pi).sqrt())
    assert result.lo<=reference<=result.hi


@pytest.mark.parametrize('center,sigma',[(0.,.04),(.37,.04),(.99,.3),(.5,1.2)])
def test_fixed_peak_bracket_encloses_numerical_global_peak(center,sigma):
    model=kernel(center,sigma)
    result=pb.certify(model,tolerance='0.00001',max_splits=250)
    peak=float(wrapped_density(center,center,sigma,1.))
    lo,hi=map(float,result['peak_bracket_MW_m2'])
    assert lo-1e-12<=peak<=hi+1e-12
    assert result['tolerance_met']
    assert hi-lo<=1.000001e-5
    assert result['geometry_global_optimum_proven'] is False
    assert pb.replay(result)
    intervals=sorted(result['leaves'],key=lambda n:D(n['interval'][0]))
    assert D(intervals[0]['interval'][0])==0 and D(intervals[-1]['interval'][1])==1
    assert all(a['interval'][1]==b['interval'][0] for a,b in zip(intervals,intervals[1:]))
    bad=copy.deepcopy(result);bad['peak_bracket_MW_m2'][1]='0'
    with pytest.raises(ValueError):pb.replay(bad)


def test_interval_donors_bound_every_member_not_just_midpoint():
    model=kernel();model['period_m2']=[.99,1.01]
    model['donors'][0]['center_fraction']=[.36,.38]
    model['donors'][0]['sigma_m2']=[.035,.045]
    result=pb.certify(model,tolerance='.001',max_splits=10)
    lo,hi=map(float,result['peak_bracket_MW_m2'])
    for period in [.99,1.01]:
        for center in [.36,.38]:
            for sigma in [.035,.045]:
                peak=float(wrapped_density(center*period,center*period,sigma,period))
                assert lo<=peak<=hi
    assert not result['tolerance_met']
    assert result['scope']=='conditional interval-kernel family'


def test_tail_and_two_donor_intervals_enclose_density():
    model=kernel(.99,.15);model['donors'].append(kernel(.25,.03)['donors'][0])
    evaluator=pb.KernelBounds(model)
    for a in np.linspace(0,.99,20):
        b=min(1.,a+.01);interval=evaluator.density(pb.I(float(a),float(b)))
        x=np.linspace(a,b,11)
        values=wrapped_density(x,.99,.15,1)+wrapped_density(x,.25,.03,1)
        assert float(interval.lo)<=min(values)+1e-12
        assert float(interval.hi)>=max(values)-1e-12
    assert evaluator.tail_upper>=0


@pytest.mark.parametrize('key,value',[('sigma_m2',[0,1]),('center_fraction',[-.1,.1]),('weight_MW',[-1,1])])
def test_invalid_input_rejected(key,value):
    model=kernel();model['donors'][0][key]=value
    with pytest.raises(ValueError):pb.certify(model)


def test_exact_area_centers_are_preserved_without_float_division():
    model={'period_m2':[73.3,73.3],'donors':[{'weight_MW':[1,1],
            'center_m2':[71.1,71.1],'sigma_m2':[.1,.1]}]}
    result=pb.certify(model,tolerance='.0001',max_splits=200)
    assert result['scope']=='fixed infinite periodic kernel'
    assert result['tolerance_met']
    assert float(result['peak_bracket_MW_m2'][0])<=1/math.sqrt(2*math.pi)/.1+1e-12


def test_exact_area_center_at_period_is_a_valid_seam():
    result=pb.certify({'period_m2':[73.3,73.3],'donors':[
        {'weight_MW':[1,1],'center_m2':[73.3,73.3],'sigma_m2':[.1,.1]}]},
        tolerance='.0001',max_splits=20)
    assert result['scope']=='fixed infinite periodic kernel'
    assert result['tolerance_met']


def test_frozen_adapter_and_three_coefficient_grid():
    from types import SimpleNamespace
    from verification.run_transport_peak_pilot import freeze_kernel,candidates
    from d3d_optimization.campaign_geometry import Candidate
    e=SimpleNamespace(reference=SimpleNamespace(power_W=1e6),source=SimpleNamespace(weights=np.array([.2,.8])))
    r={'deposition':{'donor_area_m2':np.array([7.,8.]),'sigma_area_m2':np.array([.1,.2]),'total_area_m2':10.}}
    model=freeze_kernel(e,r)
    assert model['donors'][0]['center_m2']==[7.,7.]
    assert len(model['donors'])==2
    seed=Candidate('IMP',20,tuple([.01]*8+[.4]))
    rows=list(candidates(seed))
    assert len(rows)==125 and len({x.identity for x in rows})==125
    assert any(x==seed for x in rows)
    assert all(x.vector[3:]==seed.vector[3:] for x in rows)
