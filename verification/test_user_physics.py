from functools import lru_cache
from importlib.util import find_spec
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from d3d_optimization.campaign_geometry import Candidate, evaluator
from d3d_optimization.case import load_case
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.tangent_geometry import Curve


ROOT = Path(__file__).resolve().parents[1]
SHOT = ROOT / "input/D3D/d3d/d3d_shot_1a3155370c.parquet"
GEOM = ROOT / "input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json"
LEVEL = ("cubic", 32, .002, 81, .004)


class PhysicalInput:
    time_ms = 1400
    equilibrium_sha256 = "test"
    provenance = {"geometry_sha256": "test"}

    def __init__(self, toroidal_scale=1.):
        self.review = json.loads((ROOT / "baseline_review/baseline_review.json").read_text("utf8"))
        self.geometry = json.loads(GEOM.read_text("utf8"))
        self.toroidal_scale = toroidal_scale

    @lru_cache(maxsize=2)
    def field(self, method):
        field, lcfs, axis = load_case(SHOT, self.time_ms, method)
        field.tf_current_ka *= self.toroidal_scale
        return field, lcfs, axis


def physics_api():
    assert find_spec("d3d_optimization.user_physics"), "user field must feed production physics directly"
    from d3d_optimization import user_physics
    return user_physics


def test_fixed_shelf_limit_uses_anchored_normal_and_both_endpoints():
    api = physics_api()
    inputs = PhysicalInput(toroidal_scale=5.)
    limit = api.shelf_limits(inputs)
    assert limit["admissible"] is False
    assert len(limit["endpoints"]) == 2
    for end in limit["endpoints"]:
        field = inputs.field("cubic")[0]
        b = field.magnetic([end["RZ_m"]])[0]
        normal = np.array(end["normal_RZ"])
        expected = np.degrees(np.arcsin(abs(b[:2] @ normal) / np.linalg.norm(b)))
        maximum = np.degrees(np.arcsin(np.linalg.norm(b[:2]) / np.linalg.norm(b)))
        assert end["actual_deg"] == pytest.approx(expected)
        assert end["max_deg"] == pytest.approx(maximum)
    assert limit["actual_deg"] == min(end["actual_deg"] for end in limit["endpoints"])


@pytest.mark.parametrize("toroidal_scale, expected", [(1., True), (5., False)])
def test_user_width_measurable_original_and_distinct_replacement(toroidal_scale, expected, monkeypatch):
    api = physics_api()
    # The historical constructor initially evaluates five unrelated widths and requires
    # their union to pass. They have no role in a single requested operating point.
    def historical_width_forbidden(*args, **kwargs):
        pytest.fail("historical width constructor was called")
    monkeypatch.setattr("d3d_optimization.multisurface.exponential_source", historical_width_forbidden)
    inputs = PhysicalInput(toroidal_scale)
    ref = ReferencePoint(lambda_q_m=.0041)
    original = api.prepare_user(inputs, ref, LEVEL)
    baseline = original.baseline
    assert baseline["physics_measurable"] and baseline["reference_only"] and not baseline["feasible"]
    assert original.baseline_peak > 0 and len(original.sources) == 1
    assert original.sources[0].width == ref.lambda_q_m
    assert baseline["deposition"]["ledger"]["balance_error"] < 2e-12
    assert baseline["shelf_exposure"]["receiving"]["passed"] is expected
    immutable_wall = baseline["wall"].copy()
    e = evaluator(original, 20, "cpu")
    result = api.evaluate_user_candidate(e, Candidate("OMP", 20, tuple(np.zeros(16))))
    assert result["feasible"] is expected
    assert result["shape_constraints"]["passed"] and result["geometry_screen"]["passed"]
    assert not np.array_equal(result["wall"], immutable_wall)
    np.testing.assert_array_equal(baseline["wall"], immutable_wall)
    assert result["shelf_exposure"]["receiving"]["passed"] is expected
    assert result["shelf_exposure"]["receiving"]["checked_zero_numeric_power_cells"] >= 0
    assert result["shelf_exposure"]["direct"]["subfloor_power_fraction"] >= 0


def test_direct_uses_exact_curved_normal_and_dark_cells_still_constrain():
    api = physics_api()
    curve = Curve([1., 0.], [2., 0.], [0., .04, -.06, 0.])
    u = np.array([.3, .5, .7])
    pts = curve.point(u)

    class UniformField:
        def magnetic(self, points):
            return np.tile([1., 0., 1.], (len(points), 1))

    field = UniformField()
    heat = SimpleNamespace(segment=np.array([1, 1, 3]), points=pts,
                           hit_power=np.array([.2, .3, .5]))
    verts = curve.point([0., .25, .50, .75, 1.])
    result = dict(maps=[heat], surface_segments={"shelf": np.array([0, 1, 2, 3])},
                  deposition=dict(wall=verts, cell_parent=np.arange(4),
                                  # No numeric deposition at either failed endpoint.
                                  cell_power=np.array([0., .5, .5, 0.]),
                                  sample_parent=np.array([1]), points=curve.point([.5])))
    exposure, angles = api._shelf_exposure(field, curve, 1.5, result)
    tangent = curve.tangent(u)
    expected = np.degrees(np.arcsin(abs(tangent[:,1]) / np.linalg.norm(tangent, axis=1) / np.sqrt(2)))
    np.testing.assert_allclose(angles, expected)
    assert exposure["direct"]["minimum_angle_deg"] == pytest.approx(min(expected))
    assert exposure["direct"]["subfloor_power_fraction"] == pytest.approx(heat.hit_power[expected<1.5].sum())
    assert exposure["receiving"]["passed"] is False
    assert exposure["receiving"]["checked_zero_numeric_power_cells"] == 2
    assert exposure["receiving"]["endpoint_angles_deg"] == [0., 0.]


@pytest.mark.parametrize("center", [.002, .998])
def test_near_endpoint_subfloor_minimum_between_passing_samples(center):
    from d3d_optimization.user_physics import _shelf_exposure
    curve = Curve([1., 0.], [2., 0.], np.zeros(4))

    class ShallowTrough:
        def magnetic(self, points):
            u = points[:, 0] - 1.
            poloidal = np.tan(np.radians(1.4999)) + .9 * (u - center)**2
            return np.column_stack([np.zeros_like(u), poloidal, np.ones_like(u)])

    grid = np.linspace(0., 1., 161)
    wall = curve.point(grid)
    donor = SimpleNamespace(segment=np.array([0]), points=wall[[0]], hit_power=np.array([1.]))
    native = dict(maps=[donor], surface_segments={"shelf": np.arange(160)},
                  deposition=dict(wall=wall, cell_parent=np.arange(160),
                                  cell_power=np.full(160, 1/160), points=wall,
                                  sample_parent=np.minimum(np.arange(161),159)))
    exposure, _ = _shelf_exposure(ShallowTrough(), curve, 1.5, native)
    row = exposure["receiving"]
    assert min(row["endpoint_angles_deg"]) > 1.5
    assert not row["passed"]
    assert row["minimum_angle_deg"] == pytest.approx(1.4999, abs=1e-10)
    assert row["minimum_u"] == pytest.approx(center, abs=1e-7)


def test_two_flux_roots_on_fixed_imp_segment_use_resolved_strike_without_excluding_omp():
    from d3d_optimization import user_physics as api
    from d3d_optimization.continuous_boundary import continuous_design
    from d3d_optimization.independent_sources import split_design_wall

    review = PhysicalInput().review
    a, b = np.asarray(review["wall"][69:71])
    v = b-a
    strike = a+.241*v
    review["legs"]["inner"] = dict(segment=69, strike_RZ_m=strike.tolist(), refinement_m=1e-8)

    class TwoCrossings:
        def values(self, points):
            u = (np.asarray(points)-a) @ v / (v @ v)
            return (u-.241)*(u-.783)

    with pytest.raises(ValueError, match="expected one continuous IMP strike"):
        continuous_design(review, TwoCrossings(), 0.)
    assert hasattr(api, "_engineering_regions"), "fixed IMP label cannot veto lower OMP geometry"
    design, identity = api._engineering_regions(review, TwoCrossings(), 0.)
    assert identity == "resolved_inner_first_hit_on_69"
    original = split_design_wall(review)[1]
    updated = split_design_wall(design)[1]
    assert [(s["name"],s["endpoints_RZ_m"]) for s in updated[:4]] == [
        (s["name"],s["endpoints_RZ_m"]) for s in original[:4]]
    assert len(updated) == 6
    np.testing.assert_allclose(updated[4]["endpoints_RZ_m"][0], strike, atol=1e-14)
