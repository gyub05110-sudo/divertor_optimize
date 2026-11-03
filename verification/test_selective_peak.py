import copy
import importlib.util
import math
import pytest


def kernel(centers=(.25,.75),sigma=.02):
    return {'period_m2':[1,1],'donors':[
        {'weight_MW':[.5,.5],'center_fraction':[x,x],'sigma_m2':[sigma,sigma]} for x in centers]}


def test_module_exists():
    assert importlib.util.find_spec('d3d_optimization.selective_peak') is not None


def test_separated_peaks_rank_below_coincident_with_replay():
    from d3d_optimization.selective_peak import certify,replay,compare
    a=certify(kernel(),tolerance='.01',max_splits=100,seconds=10,problem_id='same')
    b=certify(kernel((.5,.5)),tolerance='.01',max_splits=100,seconds=10,problem_id='same')
    assert a['tolerance_met'] and b['tolerance_met']
    lo,hi=map(float,a['peak_bracket_MW_m2'])
    # Binary64 analytic check is not a directed-rounding proof oracle.
    assert lo-1e-12<=.5/(math.sqrt(2*math.pi)*.02)<=hi+1e-12
    assert compare(a,b,problem_id_a='same',problem_id_b='same')['status']=='a_lower'
    assert replay(a)
    bad=copy.deepcopy(a);bad['peak_bracket_MW_m2']=['0','0']
    with pytest.raises(ValueError,match='replay'):replay(bad)


def test_broad_sigma_and_zero_budget_keep_valid_bound():
    from d3d_optimization.selective_peak import certify
    m=kernel(sigma=1000000)
    r=certify(m,seconds=2,max_splits=10,tolerance='.001')
    lo,hi=map(float,r['peak_bracket_MW_m2'])
    assert lo<=1<=hi and hi-lo<.001
    m['donors'][0]['sigma_m2']=[.02,1000000]
    r=certify(m,seconds=0,max_splits=100)
    assert r['stop_reason']=='time_budget' and r['completed_splits']==0
    assert r['peak_bracket_MW_m2'][0]!='NaN'


def test_near_ties_and_model_identity_are_not_false_winners():
    from d3d_optimization.selective_peak import certify,compare
    a=certify(kernel(),seconds=0,problem_id='same')
    assert compare(a,a,problem_id_a='same',problem_id_b='same')['status']=='unresolved'
    with pytest.raises(ValueError,match='same problem'):
        compare(a,a,problem_id_a='16 donors',problem_id_b='320 donors')
    with pytest.raises(ValueError):certify(kernel(),seconds=float('nan'))


def test_comparison_rejects_relabelled_or_forged_certificate():
    from d3d_optimization.selective_peak import certify,compare
    a=certify(kernel(),seconds=0,problem_id='A')
    b=copy.deepcopy(a);b['problem_id']='B'
    with pytest.raises(ValueError):compare(a,b,problem_id_a='A',problem_id_b='A')
    b=copy.deepcopy(a);b['peak_bracket_MW_m2']=['0','0']
    with pytest.raises(ValueError):compare(a,b,problem_id_a='A',problem_id_b='A')


def test_broad_narrow_union_contains_both_endpoint_models():
    from d3d_optimization.selective_peak import SpatialBounds
    from d3d_optimization.transport_peak_bounds import I
    m=kernel();m['donors'][0]['sigma_m2']=['.01','1000000']
    family=SpatialBounds(m).density(I('.25'))
    for width in ('.01','.2','1000000'):
        point=copy.deepcopy(m);point['donors'][0]['sigma_m2']=[width,width]
        value=SpatialBounds(point).density(I('.25'))
        assert family.lo<=value.lo and family.hi>=value.hi


def test_coordinate_precision_returns_retained_evidence():
    from d3d_optimization.selective_peak import certify,replay
    m={'period_m2':['1','1'],'donors':[{'weight_MW':['1','1'],
       'center_fraction':['.3','.3'],'sigma_m2':['1e-60','1e-60']}]}
    r=certify(m,tolerance='1e-100',max_splits=1000,seconds=10)
    assert r['stop_reason']=='coordinate_precision' and not r['tolerance_met']
    assert replay(r)
