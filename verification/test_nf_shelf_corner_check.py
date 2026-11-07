import numpy as np
import pytest

from verification.nf_shelf_corner_check import anchored_shelf, first_passing_distance


def test_native_shelf_basis_preserves_both_positions_and_tangents():
    result = anchored_shelf([1.6452,-1.2448], [1.40117,-1.24477])
    assert result['basis_count']==4
    assert result['max_position_difference_m'] < 1e-14
    assert result['max_tangent_difference'] < 1e-14


def test_first_passing_distance_distinguishes_below_and_unresolved():
    assert first_passing_distance(lambda d: 1+d, 1.5, 1) == pytest.approx(.5)
    assert first_passing_distance(lambda d: np.zeros_like(d)+1.6, 1.5, 1) == 0
    assert first_passing_distance(lambda d: np.zeros_like(d)+1., 1.5, 1) is None
