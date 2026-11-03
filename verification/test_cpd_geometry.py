"""Geometric contracts for independent CPD material components (metres)."""
from pathlib import Path
import numpy as np
import pytest
from shapely.geometry import Point

import d3d_optimization  # enables the repository's isolated scientific dependencies
from cpd_optimization.geometry import build_geometry, load_geometry, sample_three_point_arc

SOURCE = Path(__file__).resolve().parents[1] / 'input/CPD_5/CPD_5_nodes.yaml'


@pytest.mark.parametrize('start,mid,end', [
    ((1., 0.), (0., 1.), (-1., 0.)),
    ((1., 0.), (0., -1.), (-1., 0.)),
    ((1., 0.), (-1., 0.), (0., 1.)),  # clockwise major arc
])
def test_arc_passes_exactly_through_three_points_on_correct_circle(start, mid, end):
    points = sample_three_point_arc(start, mid, end, max_chord_m=0.07)
    np.testing.assert_array_equal(points[0], start)
    np.testing.assert_array_equal(points[-1], end)
    assert any(np.array_equal(p, mid) for p in points)
    np.testing.assert_allclose(np.linalg.norm(points, axis=1), 1., atol=2e-14)
    assert np.max(np.linalg.norm(np.diff(points, axis=0), axis=1)) <= 0.07
    if mid == (0., 1.):
        assert points[:, 1].min() >= -1e-14
    if mid == (0., -1.):
        assert points[:, 1].max() <= 1e-14
    if mid == (-1., 0.):
        assert points[:, 1].min() < -.99


def test_degenerate_arc_is_explicit_error():
    with pytest.raises(ValueError, match='collinear|degenerate'):
        sample_three_point_arc((0, 0), (1, 0), (2, 0))


def square(x0, y0, x1, y1):
    return {'nodes': [dict(node_id=i+1, x=x, z=z,
                           operation='moveTo' if i == 0 else 'lineTo')
                      for i, (x, z) in enumerate([(x0,y0),(x1,y0),(x1,y1),(x0,y1)])]}


def detached_example():
    return {'metadata': {'units': 'm'}, 'components': {
        'Blanket': {'contours': {'outer': square(0, 2, 1, 3)}},
        'Divertor': {'contours': {'outer': square(0, 0, 1, 1)}},
        'Shield': {'contours': {'outer': square(-2, 0, -1, 3)}},
        'VV': {'contours': {'outer': square(-4, -4, 4, 4),
                            'inner': square(-3, -3, 3, 3)}}}}


def test_gap_ray_hits_actual_vessel_not_artificial_bridge():
    geometry = build_geometry(detached_example())
    hit = geometry.first_hit((2, 1.5), (-3.5, 1.5))
    assert hit.component == 'Shield'
    np.testing.assert_allclose(hit.point_m, (-1, 1.5))
    assert hit.distance_m == pytest.approx(3)
    assert hit.face.start_node_id == 2
    assert hit.face.end_node_id == 3
    assert geometry.first_hit((2, 1.5), (-.5, 1.5)) is None
    vessel_hit = geometry.ray_first_hit((2, 1.5), (1, 0), max_distance_m=10)
    assert vessel_hit.component == 'VV'
    assert vessel_hit.face.contour == 'inner'
    np.testing.assert_allclose(vessel_hit.point_m, (3, 1.5))


def test_contacts_tangent_collinear_and_start_inside():
    geometry = build_geometry(detached_example())
    hit = geometry.first_hit((2, 1), (-.5, 1))
    np.testing.assert_allclose(hit.point_m, (1, 1))
    assert hit.component == 'Divertor'
    hit = geometry.first_hit((.5, .5), (2, .5))
    assert hit.distance_m == 0
    assert hit.face is None  # the origin is in material, not on a face
    assert geometry.first_hit((2, 1.5), (2, 1.5)) is None


def test_invalid_polygon_is_rejected_without_silent_repair():
    data = detached_example()
    nodes = data['components']['Divertor']['contours']['outer']['nodes']
    nodes[1]['x'], nodes[1]['z'] = 1, 1
    nodes[2]['x'], nodes[2]['z'] = 1, 0
    with pytest.raises(ValueError, match='Divertor.*invalid'):
        build_geometry(data)


def test_real_cpd_units_independent_solids_and_provenance():
    geometry = load_geometry(SOURCE)
    assert set(geometry.solids) == {'Blanket', 'Divertor', 'Shield', 'VV'}
    assert all(s.is_valid and s.area > 0 for s in geometry.solids.values())
    assert len(geometry.solids['VV'].interiors) == 1
    assert not geometry.solids['VV'].covers(Point(4, 0))
    assert geometry.solids['VV'].covers(Point(1.95, 0))
    contour = geometry.contours[('Divertor', 'outer')]
    np.testing.assert_allclose(contour.nodes[0].point_m, (2.055, -1.533))
    np.testing.assert_array_equal(contour.points_m[0], contour.points_m[-1])
    assert len(contour.segment_faces) == len(contour.points_m)-1
    assert len(contour.faces) == 18
    assert contour.faces[-1].start_node_id == 18
    assert contour.faces[-1].end_node_id == 1
    assert contour.faces[-1].is_closing
    assert geometry.solids['Blanket'].distance(geometry.solids['Divertor']) > .015
    # The 20 mm inboard slot really remains open.
    assert geometry.first_hit((2.8, -1.523), (2.4, -1.523)) is None
    for face in contour.faces:
        if face.operation == 'threePointArc':
            assert any(np.array_equal(p, face.midpoint_m) for p in face.points_m)


def test_movable_divertor_rebuild_preserves_original_and_arc_midpoints():
    geometry = load_geometry(SOURCE)
    original = geometry.contours[('Divertor', 'outer')].nodes[9].point_m
    moved = geometry.with_divertor_nodes({10: (original[0]+.01, original[1])})
    np.testing.assert_array_equal(geometry.contours[('Divertor', 'outer')].nodes[9].point_m, original)
    assert moved.solids['Divertor'].is_valid
    assert moved.solids['Divertor'].symmetric_difference(geometry.solids['Divertor']).area > 0
    assert moved.solids['Blanket'].equals_exact(geometry.solids['Blanket'], 0)
    np.testing.assert_array_equal(moved.contours[('Divertor', 'outer')].nodes[8].midpoint_m,
                                  geometry.contours[('Divertor', 'outer')].nodes[8].midpoint_m)


def test_unknown_units_fail_explicitly():
    data = detached_example()
    data['metadata']['units'] = 'inches'
    with pytest.raises(ValueError, match='units'):
        build_geometry(data)
