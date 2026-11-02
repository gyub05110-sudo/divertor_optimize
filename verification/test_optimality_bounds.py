"""Bounds must enclose functions, cover the domain, and avoid false certificates."""
from fractions import Fraction as F
import copy
import math
import subprocess
import sys
import numpy as np
import pytest
from scipy.interpolate import BSpline
from d3d_optimization import optimality_bounds as ob


def cylinder():
    # Open meridian is sufficient for testing the area integral.
    return {'wall': [[2., 0.], [2., 3.]], 'curves': [], 'power_MW': 1.}


def spline_model():
    return {'wall': [[2., 0.], [2., 3.]], 'power_MW': 1., 'curves': [
        {'segment': 0, 'normal': [1., 0.], 'knots': [0.]*4+[.2,.4,.6,.8]+[1.]*4,
         'controls': [0.,0.,.01,-.02,.03,0.,0.,0.]}]}


def test_pi_and_sqrt_are_rational_enclosures():
    lo, hi = ob.pi_bounds()
    assert lo < F('3.14159265358979323846264338327950288419716939937510') < hi
    for value in [F(0), F(2), F(1234567, 13), F(1,10**40)]:
        upper = ob.sqrt_upper(value)
        assert upper*upper >= value
        assert float(upper) == pytest.approx(math.sqrt(float(value)))


def test_sqrt_does_not_depend_on_float_range():
    for value in [F(10)**1000,F(10)**-1000]:
        upper=ob.sqrt_upper(value)
        assert upper*upper>=value
        assert upper*upper<value*(1+F(1,10**14))


def test_subnormal_sqrt_terminates_without_nextafter_walk():
    code="from fractions import Fraction as F; from d3d_optimization.optimality_bounds import sqrt_upper; x=F('6e-324'); y=sqrt_upper(x); assert y*y>=x"
    subprocess.run([sys.executable,'-c',code],check=True,timeout=3)


def test_interval_spline_encloses_actual_values_and_derivatives():
    q = spline_model()['curves'][0]
    controls = [ob.Interval(x) for x in q['controls']]
    controls[2] = ob.Interval(-.04,.05)
    for value in [-.04,0.,.05]:
        c = q['controls'].copy(); c[2] = value
        spline = BSpline(q['knots'],c,3)
        for a,b in zip(q['knots'][3:8],q['knots'][4:9]):
            for derivative in [0,1]:
                result = ob.spline_range(q['knots'],controls,ob.Interval(a,b),derivative)
                values = spline(np.linspace(a,b,53),nu=derivative)
                assert min(values) >= float(result.lo)-1e-12
                assert max(values) <= float(result.hi)+1e-12


def test_cylinder_area_upper_and_lower_bound():
    area = ob.area_upper(cylinder(),None,None,subcells=2)
    assert float(area) >= 12*math.pi
    assert float(area) == pytest.approx(12*math.pi,rel=1e-14)


def test_partition_replay_and_no_false_optimality_certificate():
    result = ob.build_bounds(spline_model(),[0,2],[-.04,.05],splits=5,subcells=2)
    assert ob.replay(result)
    assert len(result['leaves']) == 6
    leaves=sorted(result['leaves'],key=lambda n:F(n['domain'][0]))
    assert leaves[0]['domain'][0]==result['domain'][0]
    assert leaves[-1]['domain'][1]==result['domain'][1]
    assert all(a['domain'][1]==b['domain'][0] for a,b in zip(leaves,leaves[1:]))
    history = [F(h['lower_bound_MW_m2']) for h in result['history']]
    assert history == sorted(history)
    for coefficients in np.linspace(-.04,.05,11):
        q = spline_model()['curves'][0]; c=q['controls'].copy();c[2]=coefficients
        u=np.linspace(0,1,10001); d=BSpline(q['knots'],c,3)
        area=np.trapezoid(2*np.pi*(2+d(u))*np.sqrt(9+d(u,nu=1)**2),u)
        assert float(history[-1]) <= 1/area
    report=ob.compare_numerical(result,2.4)
    assert report['certified_gap_MW_m2'] is None
    assert report['global_optimum_proven'] is False
    bad=copy.deepcopy(result);bad['history'][-1]['lower_bound_MW_m2']='3'
    with pytest.raises(ValueError):ob.replay(bad)
    bad=copy.deepcopy(result);bad['leaves'].pop()
    with pytest.raises(ValueError):ob.replay(bad)


@pytest.mark.parametrize('domain',[[1,0],[float('nan'),1]])
def test_invalid_domain_rejected(domain):
    with pytest.raises(ValueError):ob.build_bounds(spline_model(),[0,2],domain)


def test_freeze_actual_imp_family_preserves_candidate_and_control_mapping():
    from verification.run_optimality_bounds import freeze_family
    from verification.test_campaign_consistency import geometry_example
    e,candidate=geometry_example('IMP/R020')
    frozen,variable=freeze_family(e,candidate,2,1.)
    assert frozen['power_MW']==1.
    curve=next(q for q in frozen['curves'] if q['segment']==variable[0])
    assert variable[1]==4
    assert curve['controls'][variable[1]]==candidate.vector[2]
    area=ob.area_upper(frozen,variable,ob.Interval(candidate.vector[2]),subcells=4)
    wall=observation_wall(e,candidate)
    numerical_area=np.sum(np.pi*(wall[:-1,0]+wall[1:,0])*np.linalg.norm(np.diff(wall,axis=0),axis=1))
    assert float(area)>=numerical_area


def observation_wall(e,candidate):
    from d3d_optimization.campaign_geometry import assembled_wall
    return assembled_wall(e,candidate.bind(e),1e-7)


def test_source_snapshot_is_self_checked_and_never_overwritten(tmp_path):
    from verification.run_optimality_bounds import snapshot_code,digest
    root=tmp_path/'repo';root.mkdir()
    for directory in ['d3d_optimization','verification']:
        (root/directory).mkdir();(root/directory/'example.py').write_text('VALUE = 1\n')
    destination=tmp_path/'snapshot'
    hashes=snapshot_code(root,destination)
    assert len(hashes)==2
    assert all(digest(destination/p)==h for p,h in hashes.items())
    (root/'verification/example.py').write_text('VALUE = 2\n')
    with pytest.raises(FileExistsError):snapshot_code(root,destination)
    assert all(digest(destination/p)==h for p,h in hashes.items())
