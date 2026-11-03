"""Contracts for fixed-face CPD normal displacements and geometric rejection."""
from pathlib import Path
import importlib

import numpy as np
import pytest
from shapely.geometry import Polygon

import d3d_optimization
from cpd_optimization.geometry import load_geometry, SourceFace

SOURCE = Path(__file__).resolve().parents[1] / 'input/CPD_5/CPD_5_nodes.yaml'


def module():
    return importlib.import_module('cpd_optimization.deformation')


def test_zero_preserves_source_geometry_and_fixed_faces():
    g = load_geometry(SOURCE, max_chord_m=.003)
    family = module().DivertorDeformation(g)
    before = g.contours['Divertor', 'outer'].points_m.copy()
    result = family.deform(np.zeros(family.dimension))
    assert family.dimension == 40
    assert 'Divertor.outer:17->18' in family.slices
    assert result['geometry_valid'], result['errors']
    np.testing.assert_allclose(result['wall'], before, rtol=0, atol=2e-14)
    assert len(result['segment_faces']) == len(result['wall'])-1
    assert all(abs(r-1) < 1e-13 for r in result['length_ratios'].values())
    changed = family.deform(np.full(family.dimension, .0001))
    for face in g.contours['Divertor', 'outer'].faces:
        if face.face_id not in family.slices:
            np.testing.assert_array_equal(changed['face_points'][face.face_id], face.points_m)
    np.testing.assert_array_equal(g.contours['Divertor', 'outer'].points_m, before)


def test_circle_derivatives_endpoints_and_length_include_normal_curvature():
    face = SourceFace('Divertor', 'outer', 8, 9, 'threePointArc',
                      np.array([[3.,0.],[2.,1.]]), (2+2**-.5,2**-.5))
    curve = module().FaceCurve(face)
    c = np.array([.04, .06, .02, -.01])
    p, first, second = curve.evaluate(np.array([0.,1.]), c)
    p0, first0, _ = curve.evaluate(np.array([0.,1.]), np.zeros(4))
    np.testing.assert_allclose(p, p0, atol=1e-14)
    np.testing.assert_allclose(first, first0, atol=1e-14)
    t = np.linspace(.01,.99,91); h=1e-5
    pp, dd, d2 = curve.evaluate(t,c)
    plus = curve.evaluate(t+h,c)[0]; minus=curve.evaluate(t-h,c)[0]
    np.testing.assert_allclose(dd,(plus-minus)/(2*h),atol=1e-7)
    np.testing.assert_allclose(d2,(plus-2*pp+minus)/h**2,atol=2e-5)
    metrics=curve.metrics(c)
    dense=curve.evaluate(np.linspace(0,1,100001),c)[0]
    exact=np.linalg.norm(np.diff(dense,axis=0),axis=1).sum()
    assert metrics['length_m'] == pytest.approx(exact,rel=2e-9)
    assert curve.reference_length_m == pytest.approx(np.pi/2)
    assert metrics['length_ratio'] == pytest.approx(exact/(np.pi/2),rel=2e-9)


def test_displacement_bound_and_candidate_chords():
    g=load_geometry(SOURCE,max_chord_m=.003)
    family=module().DivertorDeformation(g)
    c=np.linspace(-.02,.02,family.dimension)
    result=family.deform(c)
    for face_id, sl in family.slices.items():
        curve=family.curves[face_id]
        t=np.linspace(0,1,1001)
        displacement=np.linalg.norm(curve.evaluate(t,c[sl])[0]-curve.evaluate(t,np.zeros(4))[0],axis=1)
        assert displacement.max() <= .1+1e-12
        assert np.linalg.norm(np.diff(result['face_points'][face_id],axis=0),axis=1).max() <= .003*(1+1e-10)
    with pytest.raises(ValueError,match='bounded'):
        family.deform(np.full(family.dimension,.10001))


def test_collision_and_lcfs_intrusion_are_explicit_rejections():
    g=load_geometry(SOURCE,max_chord_m=.003)
    family=module().DivertorDeformation(g)
    base=family.deform(np.zeros(family.dimension))
    # A supplied LCFS occupying material must never be silently clipped.
    rejected=family.deform(np.zeros(family.dimension),lcfs=base['solid'])
    assert not rejected['geometry_valid']
    assert 'lcfs_intrusion' in rejected['errors']
    # Independent obstacle geometry exercises collision detection without relying
    # on the particular source's narrow-gap deformation direction.
    g.solids['Blanket']=Polygon(base['solid'].exterior)
    rejected=module().DivertorDeformation(g).deform(np.zeros(family.dimension))
    assert not rejected['geometry_valid']
    assert 'collision_Blanket' in rejected['errors']


def test_evaluator_api_and_finer_active_arc_sampling():
    g=load_geometry(SOURCE,max_chord_m=.003)
    family=module().DivertorDeformation(g)
    result=family.evaluate(np.zeros(family.dimension),max_chord_m=.001)
    assert result['feasible']
    assert len(result['segment_face_ids']) == len(result['wall'])-1
    assert set(result['segment_face_ids']) == set(range(1,19))
    for name,curve in family.curves.items():
        if curve.arc:
            assert np.linalg.norm(np.diff(result['face_points'][name],axis=0),axis=1).max() <= .001*(1+1e-10)
    assert result['min_curvature_radius_m'] >= .020


def test_large_short_face_deformation_rejected_by_curvature_or_length():
    g=load_geometry(SOURCE,max_chord_m=.003)
    family=module().DivertorDeformation(g)
    shortest=min(family.curves,key=lambda key:family.curves[key].reference_length_m)
    c=np.zeros(family.dimension)
    c[family.slices[shortest]]=[.1,-.1,.1,-.1]
    result=family.deform(c)
    assert not result['geometry_valid']
    assert any(e.startswith(('curvature_radius_below_limit:','material_length_exceeds_limit:')) for e in result['errors'])
