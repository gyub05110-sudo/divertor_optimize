import numpy as np
import pytest
from cpd_optimization.evaluator import midplane_separatrix, deadline_remaining


class CircularFlux:
    r = np.linspace(.5,4.,50)
    def values(self,p):
        p=np.asarray(p)
        return (p[:,0]-2.)**2+p[:,1]**2


def test_midplane_roots_use_axis_and_selected_flux_not_d3d_guess():
    field=CircularFlux()
    assert midplane_separatrix(field,[2.,.2],1.,'IMP') == pytest.approx(2.-np.sqrt(.96))
    assert midplane_separatrix(field,[2.,.2],1.,'OMP') == pytest.approx(2.+np.sqrt(.96))
    with pytest.raises(ValueError): midplane_separatrix(field,[2.,.2],1.,'both')


def test_deadline_never_restarts_elapsed_budget():
    assert deadline_remaining(100., 10., now=104.) == 6.
    assert deadline_remaining(100., 10., now=111.) == 0.
