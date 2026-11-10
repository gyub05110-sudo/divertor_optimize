"""Input-only tests for the native Stage 5 handoff."""
from importlib import import_module, util

import numpy as np
import pytest

from d3d_optimization.frozen_field import FrozenField, exponential_source


@pytest.fixture
def bridge():
    assert util.find_spec('verification.stage5_inputs'), 'Stage 5 native input bridge is required'
    return import_module('verification.stage5_inputs')


@pytest.mark.parametrize('current,polarity', [(-93., 1.), (122., -1.)])
def test_native_eqdsk_serializes_signed_field_flux_grid_axis_and_lcfs(bridge, current, polarity):
    r, z = np.linspace(1, 3, 7), np.linspace(-1, 1, 5)
    rr, zz = np.meshgrid(r, z)
    field = FrozenField(r, z, polarity*(rr**2 + 2*zz + .3), current)
    axis = np.array([1.6, .17])
    lcfs = np.array([[1.25, -.3], [1.8, -.5], [2.2, .4], [1.4, .5]])
    wall = np.array([[1.1, -.8], [2.8, -.8], [2.8, .8], [1.1, .8], [1.1, -.8]])
    level = polarity * 5.17
    lines = bridge.geqdsk(field, axis, level, lcfs, wall, 'test').decode().splitlines()
    assert lines[0].split()[-2:] == ['7', '5']
    counts = lines.index('4 5')
    values = np.array([float(line[i:i+16]) for line in lines[1:counts] for i in range(0,len(line),16)])
    np.testing.assert_allclose(values[5:10], [*axis, polarity*(axis[0]**2+2*axis[1]+.3),level, .0288*current/axis[0]])
    np.testing.assert_allclose(values[20:27], .0288*current)
    np.testing.assert_allclose(values[20+4*len(r):20+4*len(r)+rr.size].reshape(zz.shape),field.psi)
    vertices = [float(line[i:i+16]) for line in lines[counts+1:] for i in range(0,len(line),16)]
    np.testing.assert_allclose(np.reshape(vertices,(-1,2)), np.vstack((lcfs,wall)))


def test_donors_are_offset_from_case_separatrix_with_original_watts_and_downward_sign(bridge):
    r, z = np.linspace(1, 3, 7), np.linspace(-1, 1, 5)
    rr, zz = np.meshgrid(r, z)
    field = FrozenField(r, z, -rr**2 + zz, -123.)
    crossing, zaxis, wall_r = 2.13, -.027, 2.19
    wall = np.array([[1.1,-.8],[wall_r,-.8],[wall_r,.8],[1.1,.8],[1.1,-.8]])
    edges = (wall_r-crossing)*np.linspace(0,1,9)**2
    src = exponential_source(field,crossing,zaxis,edges,.003)
    watts = 3e6*src.weights
    table = bridge.donor_table(field, crossing, zaxis, wall, .003, 3e6, src.x, watts)
    np.testing.assert_array_equal(table[:,4], crossing+src.x)
    np.testing.assert_array_equal(table[:,5], zaxis)
    np.testing.assert_array_equal(table[:,6], 1)
    np.testing.assert_array_equal(table[:,8], watts)
    with pytest.raises(ValueError,match='assigned power'):
        bridge.donor_table(field,crossing,zaxis,wall,.003,3e6,src.x,watts[::-1])
    with pytest.raises(ValueError,match='offset'):
        bridge.donor_table(field,crossing,zaxis,wall,.003,3e6,src.x+.01,watts)


def test_panel_keeps_all_labels_even_if_two_roles_have_same_winner(bridge):
    roles = bridge.requested_roles('70481becc0')
    assert roles[:7] == bridge.requested_roles('3264fdd71f')
    assert set(roles[:7]) == {(30713,r) for r in ('baseline','hotspot','full','floor-only','floor-entry','floor-connector','floor-shelf')}
    assert set(roles[7:]) == {(30757, 'hotspot'),(30757,'full')}


def test_source_cdf_preserves_monotonicity_when_float_sum_creeps_above_one(bridge):
    weights = np.r_[np.full(100, .01), 1e-16]
    original = weights.copy()
    assert np.cumsum(weights)[-2] > 1.
    cdf = bridge.native_cdf(weights)
    assert np.all(np.diff(cdf) >= 0) and cdf[-1] == 1.
    np.testing.assert_array_equal(weights,original)
    np.testing.assert_allclose(cdf,np.cumsum(weights),rtol=0,atol=2e-15)
