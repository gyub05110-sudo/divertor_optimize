import numpy as np
import pytest
from scipy.special import ndtr
from d3d_optimization.flux_footprint import branch_power

def test_cell_probabilities_and_unresolved_power_use_absolute_input_weights():
    d=branch_power(np.array([-1.,0.,2.]), [.5], [1.], [2.])
    np.testing.assert_allclose(d['cell_power'],2*np.diff(ndtr([-1.5,-.5,1.5])))
    assert d['below']>0 and d['above']>0
    assert d['cell_power'].sum()+d['below']+d['above']==pytest.approx(2.)
    rev=branch_power([2.,0.,-1.],[.5],[1.],[2.])
    np.testing.assert_allclose(rev['cell_power'][::-1],d['cell_power'])
    assert rev['below']==pytest.approx(d['below'])

def test_ambiguous_geometry_is_rejected_and_zero_width_has_explicit_domain():
    with pytest.raises(ValueError):branch_power([0.,1.,.5],[.5],[.1],[1.])
    with pytest.raises(ValueError):branch_power([0.,1.],[.5],[0.],[1.])
    with pytest.raises(ValueError):branch_power([0.,1.],[.5],[.1],[-1.])

def test_resolved_small_tail_is_positive_after_nine_standard_deviations():
    d=branch_power([9.,10.,11.],[0.],[1.],[1.])
    assert d['cell_power'][0]==pytest.approx(ndtr(-9.)-ndtr(-10.),rel=2e-12,abs=0)

def test_horizontal_limit_and_tilted_cell_mass():
    from verification.nf_incidence_limit import compare_target
    control=compare_target(0.,0.,.002)
    tilted=compare_target(2.,0.,.002)
    assert control['total_variation']<.005
    assert tilted['total_variation']>.25
    assert abs(tilted['reference_peak']/control['reference_peak']-1/np.sqrt(5))<.004

