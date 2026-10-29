import importlib.util
from pathlib import Path
import numpy as np
import pytest


def api():
    assert importlib.util.find_spec('d3d_optimization.imp_tangent'), 'IMP tangent adapter missing'
    from d3d_optimization import imp_tangent
    return imp_tangent


def test_both_endpoint_tangents_and_exact_curvature():
    m=api(); a=[1.,-.2]; b=[.94,-.1]; start=[-.1,.1]; end=[-.01,.1]
    q=m.EndTangents(a,b,[.001,-.001,.002,-.002],start,end)
    np.testing.assert_allclose(q.point([0,1]),[a,b],atol=1e-15)
    for u,t in [(0,start),(1,end)]:
        t=np.array(t); got=q.tangent(u)
        assert got@t>0
        assert abs(got[0]*t[1]-got[1]*t[0])<1e-14
    dense=np.max(abs(q.curvature(np.linspace(0,1,50001))))
    assert q.max_curvature()>=dense-1e-8
    assert q.max_curvature()-dense<1e-4
    with pytest.raises(ValueError): m.EndTangents(a,b,np.zeros(4),[1,-1],end)


@pytest.fixture(scope='module')
def original():
    from d3d_optimization.boundary_evaluator import prepare_boundary
    from d3d_optimization.reference_transport import ReferencePoint
    return prepare_boundary(Path(__file__).resolve().parents[1],'IMP',ReferencePoint())


def test_three_imp_connections_and_omp_preservation(original):
    m=api(); e=m.imp_evaluator(original,.02,'cpu')
    old=original.baseline['wall'].copy()
    for fraction in [0,.3,.5,.8,1]:
        e.set_junction(fraction); curves=e.curves(np.zeros(24))
        s,v=curves['imp_sloped'],curves['imp_vertical']
        for first,second in [(s.tangent(0),s.v),(s.tangent(1),v.tangent(0)),(v.tangent(1),v.v)]:
            assert first@second>0
            assert abs(first[0]*second[1]-first[1]*second[0])<1e-14
        for name in ['floor','shelf','omp_entry','omp_connector']:
            expected=original.targets[name]
            np.testing.assert_allclose(curves[name].point(np.linspace(0,1,len(expected))),expected,atol=5e-15)
    np.testing.assert_array_equal(original.baseline['wall'],old)
    assert e.baseline is original.baseline
    with pytest.raises(ValueError): e.set_junction(-.01)
    with pytest.raises(ValueError): e.set_junction(np.nan)


def test_junction_change_invalidates_cache_and_does_not_mutate_original(original):
    e=api().imp_evaluator(original,.02,'cpu')
    e.cache['sentinel']={}; e.set_junction(.7)
    assert not e.cache
    assert not hasattr(original.direct.base,'junction_fraction')
    c=np.zeros(24); c[0]=.001
    assert not e.constraints(c)['passed']


def test_geometry_projected_seed_and_gpu_parity(original):
    m=api(); e=m.imp_evaluator(original,.05,'gpu')
    x=m.project_imp(e,np.zeros(8),.5)
    assert x is not None, '50 mm geometry feasible seed required'
    c=m.expand_imp(x[:8]); e.set_junction(x[8])
    check=e.constraints(c)
    assert check['passed']
    assert max(e.precise_lengths(c)[0].values())<=1.2
    gpu=e.evaluate(c)
    cpu=m.imp_evaluator(original,.05,'cpu'); cpu.set_junction(x[8]); r=cpu.evaluate(c)
    assert gpu['feasible'] and r['feasible']
    assert gpu['peak_q']==pytest.approx(r['peak_q'],rel=1e-10)
    assert gpu['deposition']['ledger']['balance_error']<2e-12
    assert max(check['tangent_angle_errors_deg'])<1e-8
    np.testing.assert_allclose(gpu['deposition']['q'],r['deposition']['q'],rtol=2e-12,atol=1e-13)
    first=gpu['wall'].copy(); first_q=gpu['deposition']['q'].copy()
    e.set_junction(.65); other=e.curves(c)['imp_sloped'].point(np.linspace(0,1,161))
    assert np.max(abs(other-gpu['targets']['imp_sloped']))>1e-4
    e.evaluate(c); e.set_junction(x[8]); again=e.evaluate(c)
    np.testing.assert_array_equal(again['wall'],first)
    np.testing.assert_array_equal(again['deposition']['q'],first_q)
