"""An upper LCFS saddle can explain a lower-divertor input rejection without accepting the frame."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from d3d_optimization.frozen_field import FrozenField
from d3d_optimization.user_input import load_user_input


def topology():
    assert importlib.util.find_spec("d3d_optimization.user_topology"), "upper LCFS evidence is missing"
    from d3d_optimization import user_topology
    return user_topology


def analytic_upper_null():
    # psi = (R-2)^2 + Z^2/2 - Z^3/3 gives an elliptic axis at (2,0)
    # and an upper saddle at (2,1), with both inward branches given analytically.
    r, z = np.linspace(1., 3., 21), np.linspace(-1., 1.5, 26)
    rr, zz = np.meshgrid(r, z)
    field = FrozenField(r, z, (rr-2)**2 + zz**2/2 - zz**3/3, 1.)
    levels = np.linspace(-.5, 1., 63)
    radius = np.sqrt(np.maximum(0., 1/6 - levels**2/2 + levels**3/3))
    ring = np.r_[np.c_[2+radius[1:-1], levels[1:-1]], [[2, 1]],
                 np.c_[2-radius[-2:0:-1], levels[-2:0:-1]], [[2, -.5]]]
    return field, ring, np.array([2., 0.])


def test_upper_saddle_requires_two_lcfs_branches_and_elliptic_axis():
    diagnose = topology().upper_saddle_evidence
    field, lcfs, axis = analytic_upper_null()
    result = diagnose(field, lcfs, axis)
    assert result["classification"] == "upper_lcfs_saddle"
    assert len(set(result["neighbour_branch_indices"])) == 2
    assert max(result["neighbour_branch_angles_deg"]) < 2
    np.testing.assert_allclose(result["saddle_RZ_m"], [2, 1], atol=1e-8)
    assert result["axis_hessian_determinant"] > 0 and result["saddle_hessian_determinant"] < 0
    # Both chords now point along the same branch.
    invalid = lcfs.copy()
    i = np.argmax(lcfs[:, 1])
    invalid[i+1] = .5*(lcfs[i]+lcfs[i-1])
    assert diagnose(field, invalid, axis) is None
    assert diagnose(field, lcfs, np.array([2.1, .2])) is None


def test_smooth_elliptic_contour_has_no_upper_saddle():
    field, _, axis = analytic_upper_null()
    angle = np.linspace(0, 2*np.pi, 100)
    ring = axis + .2*np.c_[np.cos(angle), np.sin(angle)]
    assert topology().upper_saddle_evidence(field, ring, axis) is None


def test_nearby_small_gradient_is_not_a_coincident_upper_saddle():
    field, lcfs, axis = analytic_upper_null()
    # With small flux units an absolute residual can pass even when the exact
    # saddle is outside the stated 10-micrometre matching region.
    scaled = FrozenField(field.r, field.z, .001*field.psi, 1.)
    lcfs[np.argmax(lcfs[:, 1]), 1] += 2e-5
    assert topology().upper_saddle_evidence(scaled, lcfs, axis) is None


@pytest.mark.parametrize("ident", ["8d1575d8c4", "a7e51f17f7", "8bb901c2a1"])
def test_previous_lower_x_failure_remains_rejected_with_upper_lcfs_evidence(ident):
    shot = Path(f"D:/DivertorOptimization/studies/nf-stage3-20261108/inputs/raw/d3d_shot_{ident}.parquet")
    if not shot.is_file():
        pytest.skip("requires archived Stage 3 development discharge")
    with pytest.raises(ValueError, match="upper LCFS saddle") as caught:
        load_user_input(shot, 1400.)
    assert str(caught.value.__cause__) in str(caught.value)
    evidence = json.loads(json.dumps(caught.value.evidence, allow_nan=False))
    assert evidence["saddle_distance_to_lcfs_top_m"] < 1e-6
    assert 0 < max(evidence["neighbour_branch_angles_deg"]) < 8


def development_shot():
    path = Path("D:/DivertorOptimization/studies/nf-stage3-20261108/inputs/raw/d3d_shot_8d1575d8c4.parquet")
    if not path.is_file():
        pytest.skip("requires archived Stage 3 development discharge")
    return path


def test_cli_writes_topology_evidence_and_original_reason(tmp_path):
    from optimize_divertor import main
    destination = tmp_path/"upper-null"
    assert main(["--equilibrium", str(development_shot()), "--time-ms", "1400",
                 "--output", str(destination), "--budget-seconds", "60"]) == 2
    result = json.loads((destination/"result.json").read_bytes())
    assert result["status"] == "input_error" and "sampling outside EFIT grid" in result["reason"]
    assert result["equilibrium_topology_diagnostic"]["classification"] == "upper_lcfs_saddle"


def test_secondary_topology_diagnostic_cannot_hide_the_original_error(monkeypatch):
    def singular(*_):
        raise np.linalg.LinAlgError("ill-conditioned local fit")
    monkeypatch.setattr(topology(), "upper_saddle_evidence", singular)
    with pytest.raises(ValueError, match="lower-divertor topology unresolved.*sampling outside EFIT grid") as caught:
        load_user_input(development_shot(), 1400.)
    assert caught.value.evidence is None
