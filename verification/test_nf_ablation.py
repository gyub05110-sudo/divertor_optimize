from types import SimpleNamespace
import numpy as np
import pytest
from verification.nf_width_ablation import BaselineWidths

def test_ablation_substitutes_only_width_without_mutating_thermal_result():
    stored={'sigma_m':np.array([.2,.3]),'Tu_eV':np.array([10.,12.])}
    actual=SimpleNamespace(evaluate=lambda hits,**kwargs:stored,grad_sep=1.)
    frozen=BaselineWidths(actual,[.1,.4]);r=frozen.evaluate(None)
    np.testing.assert_array_equal(r['sigma_m'],[.1,.4])
    np.testing.assert_array_equal(stored['sigma_m'],[.2,.3])
    assert r['Tu_eV'] is stored['Tu_eV'] and frozen.grad_sep==1.
    with pytest.raises(ValueError):BaselineWidths(actual,[.1,np.nan])

