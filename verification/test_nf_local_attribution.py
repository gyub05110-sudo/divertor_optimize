"""Interaction and physical-coordinate regressions for frozen-profile attribution."""
from itertools import product

import numpy as np

from verification.nf_local_attribution import shapley, wall_location


def test_donor_effects_split_interactions_over_all_orders():
    # The 6ab term is shared between a and b; the 12abc term equally among all three.
    # The two donors have opposite signed responses, preventing accidental pooling.
    scale = np.array([1., -2.])
    cube = {(a,b,c): scale*(9. + a - 2*b + 3*c + 6*a*b + 12*a*b*c)
            for a,b,c in product(range(2), repeat=3)}
    np.testing.assert_allclose(shapley(cube), np.array([8., 5., 7.])[:,None]*scale)


def test_partial_conical_area_at_physical_arclength():
    wall = np.array([[1., 0.], [2., 0.], [2., 1.]])
    profile = {'wall': wall, 'u': np.array([0., .5, 1.]), 'length': 2.}
    edges = np.pi*np.array([0., 3., 7.])
    mid = wall_location(profile, edges, .25)
    assert mid['cell'] == 0 and mid['s_m'] == .5
    np.testing.assert_allclose(mid['RZ_m'], [1.5, 0.])
    np.testing.assert_allclose(mid['area_m2'], 1.25*np.pi)
    assert wall_location(profile, edges, .5)['cell'] == 1
