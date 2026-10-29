import importlib.util
import numpy as np
import pytest


def api():
    assert importlib.util.find_spec('d3d_optimization.tangent_geometry') is not None, 'new constrained geometry module missing'
    from d3d_optimization.tangent_geometry import Curve
    return Curve


def test_oriented_endpoint_tangent_and_unchanged_endpoints():
    Curve=api(); a=np.array([1.4277,-1.35782]); b=np.array([1.1645,-1.35894])
    out=np.array([1.0882576947643479,-1.2853220897608946])-b
    curve=Curve(a,b,[.01,-.02,.03,.04],out)
    np.testing.assert_allclose(curve.point([0,1]),[a,b],atol=1e-15)
    t=curve.tangent(1); assert np.dot(t,out)>0
    assert abs(np.linalg.det(np.array([t,out])))<1e-14


def test_old_spline_parity_without_imposed_tangent():
    Curve=api()
    from d3d_optimization.shape import NormalDisplacementBSpline,ReferencePolyline
    c=[.01,-.02,.03,-.04]; a=[1.,-1.]; b=[1.3,-1.02]; u=np.linspace(0,1,101)
    old=NormalDisplacementBSpline(c); ref=ReferencePolyline([a,b]); new=Curve(a,b,c)
    np.testing.assert_allclose(new.point(u),old.deform(ref,u),atol=1e-14)
    np.testing.assert_allclose(new.tangent(u),old.deformed_tangent(ref,u),atol=1e-13)


def test_continuous_curvature_and_gap_against_dense_samples():
    Curve=api(); rng=np.random.default_rng(229)
    u=np.linspace(0,1,20001); tip=np.array([1.12,-.97])
    for _ in range(6):
        curve=Curve([1,-1],[1.3,-1],rng.normal(0,.02,4))
        peak=curve.max_curvature(); sampled=np.max(abs(curve.curvature(u)))
        assert sampled <= peak*(1+1e-10)
        assert abs(peak-sampled)<1e-4
        gap=curve.minimum_distance(tip)[0]; sampled=np.min(np.linalg.norm(curve.point(u)-tip,axis=1))
        assert gap <= sampled+1e-13
        assert sampled-gap<2e-7


def test_curvature_sign_symmetry_and_invalid_inputs():
    Curve=api(); c=np.array([.01,-.02,.03,-.01])
    a=Curve([1,0],[1.3,0],c); b=Curve([1,0],[1.3,0],-c)
    assert a.max_curvature()==pytest.approx(b.max_curvature())
    with pytest.raises(ValueError): Curve([1,0],[1,0],c)
    with pytest.raises(ValueError): Curve([1,0],[1.3,0],c,[-1,0])
    with pytest.raises(ValueError): Curve([1,0],[1.3,0],[np.nan]*4)
