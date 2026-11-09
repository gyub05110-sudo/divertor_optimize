"""Fixed-material exclusions use the same closed polygon and wall as candidate screening."""
from dataclasses import replace
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from shapely.geometry import LineString, Polygon

from d3d_optimization.campaign_geometry import CHORD_TOLERANCE_M, Candidate, assembled_wall
from d3d_optimization.frozen_field import FrozenField
from d3d_optimization.independent_sources import split_design_wall
from d3d_optimization.tangent_direct import TangentMultiSurface
from d3d_optimization.user_input import REFERENCE, UserInput, load_user_input


def check(user):
    assert importlib.util.find_spec("d3d_optimization.user_precheck"), "fixed-wall precheck is missing"
    from d3d_optimization.user_precheck import fixed_wall_precheck
    record = fixed_wall_precheck(user)
    return json.loads(json.dumps(record, allow_nan=False))


@pytest.fixture
def material():
    reference = json.loads(REFERENCE.read_bytes())
    # Select the supplied fixed inner reporting portions.
    review = dict(reference, legs={"inner": {"segment": 70}})
    r, z = np.linspace(.8, 2.6, 9), np.linspace(-1.5, 1.5, 9)
    field = FrozenField(r, z, np.zeros((9, 9)), 1.)
    return UserInput(0., "", review, {}, {}, field, np.empty((0, 2)), np.array([1.5, 0.]))


def with_plasma(user, polygon):
    return replace(user, _lcfs=np.array(polygon.exterior.coords))


def square(center, radius):
    return Polygon(np.array(center) + radius * np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]]))


def candidate_wall(user, vector):
    wall, surfaces, _ = split_design_wall(user.review)
    targets = {s["name"]: np.linspace(*s["endpoints_RZ_m"], 9) for s in surfaces}
    e = SimpleNamespace(wall=wall, surfaces=surfaces, targets=targets)
    e.curves = lambda coeff: {s["name"]: TangentMultiSurface.curve(e, s["name"], coeff[4*j:4*j+4])
                              for j, s in enumerate(surfaces)}
    c = Candidate("OMP", 50, tuple(vector)).coefficients
    return LineString(assembled_wall(e, c, CHORD_TOLERANCE_M))


def test_movable_shelf_interior_does_not_prove_exclusion(material):
    _, surfaces, _ = split_design_wall(material.review)
    endpoints = next(s["endpoints_RZ_m"] for s in surfaces if s["name"] == "shelf")
    plasma = square(np.mean(endpoints, axis=0), .0001)
    user = with_plasma(material, plasma)
    assert candidate_wall(user, np.zeros(16)).intersects(plasma)
    assert not candidate_wall(user, np.r_[np.zeros(4), .03*np.ones(4), np.zeros(8)]).intersects(plasma)
    record = check(user)
    assert not record["excluded"]
    assert record["fixed_segment_intersections"] == record["fixed_adjustable_endpoint_intersections"] == []
    assert record["adjustable_endpoints_tested"] == 8


def test_closed_polygon_contact_at_movable_endpoint_is_exclusion(material):
    _, surfaces, _ = split_design_wall(material.review)
    coord = next(s["endpoints_RZ_m"][1] for s in surfaces if s["name"] == "shelf")
    plasma = Polygon(np.array(coord) + [[0., 0.], [-.0002, -.0002], [.0002, -.0002]])
    user = with_plasma(material, plasma)
    record = check(user)
    assert record["excluded"] and candidate_wall(user, np.full(16, .02)).intersects(plasma)
    endpoint, = record["fixed_adjustable_endpoint_intersections"]
    assert (endpoint["region"], endpoint["endpoint"], endpoint["original_material_segment"]) == ("shelf", 1, 59)
    assert endpoint["witness_RZ_m"] == coord and not endpoint["witness_strictly_interior"]
    assert endpoint["overlap_length_m"] == 0.


@pytest.mark.parametrize("name,material_id", [("imp_sloped", 69), ("imp_vertical", 70)])
def test_inactive_imp_splines_are_fixed_material(material, name, material_id):
    _, surfaces, _ = split_design_wall(material.review)
    face = next(s for s in surfaces if s["name"] == name)
    plasma = square(np.mean(face["endpoints_RZ_m"], axis=0), .0001)
    user = with_plasma(material, plasma)
    record = check(user)
    hit, = record["fixed_segment_intersections"]
    assert record["excluded"] and candidate_wall(user, np.full(16, .02)).intersects(plasma)
    assert (hit["region"], hit["original_material_segment"], hit["split_segment"]) == (name, material_id, face["segment_index"])
    assert hit["witness_strictly_interior"] and hit["overlap_length_m"] > .0001


def test_prior_b5e07dde39_has_only_fixed_original_segments_four_and_five():
    shot = Path("D:/DivertorOptimization/studies/nf-stage3-20261108/inputs/raw/d3d_shot_b5e07dde39.parquet")
    if not shot.exists():
        pytest.skip("requires archived Stage 3 development discharge")
    assert sha256(shot.read_bytes()).hexdigest() == "a783f1ae24d01a159ad1123e24dd51f99ba6f60dd3f7577851b28d1229b0b749"
    user = load_user_input(shot, 1400.)
    record = check(user)
    assert record["excluded"] and record["fixed_adjustable_endpoint_intersections"] == []
    assert record["fixed_segments_tested"] == 77
    assert [hit["original_material_segment"] for hit in record["fixed_segment_intersections"]] == [4, 5]
    np.testing.assert_allclose([hit["overlap_length_m"] for hit in record["fixed_segment_intersections"]],
                               [.014904089065212898, .0018863851784908271], rtol=0, atol=1e-12)
