import importlib.util
from pathlib import Path
import numpy as np


def test_affine_shape_propagates_without_changing_comparison_baseline():
    assert importlib.util.find_spec('d3d_optimization.tangent_evaluator') is not None, 'tangent-aware evaluator missing'
    from d3d_optimization.boundary_evaluator import prepare_boundary
    from d3d_optimization.reference_transport import ReferencePoint
    from d3d_optimization.tangent_evaluator import constrained_evaluator
    root=Path(__file__).resolve().parents[1]
    old=prepare_boundary(root,'OMP',ReferencePoint(),count=32,step=.002,nodes=161,wall_step=.003)
    wall=old.baseline['wall'].copy(); peak=old.baseline_peak
    e=constrained_evaluator(old,.02,backend='cpu')
    assert e.baseline_peak==peak
    np.testing.assert_array_equal(e.baseline['wall'],wall)
    np.testing.assert_array_equal(old.baseline['wall'],wall)
    c=np.zeros(24); curve=e.direct.base.curve('floor',c[12:16])
    r=e.direct.base.evaluate(c)
    np.testing.assert_allclose(r['targets']['floor'],curve.point(np.linspace(0,1,161)),atol=1e-14)
    assert not np.array_equal(r['targets']['floor'],old.baseline['targets']['floor'])
    lengths,_=e.precise_lengths(c)
    assert lengths['floor']>1
    report=e.constraints(c)
    assert report['tangent_angle_error_deg']<1e-8
    full=e.evaluate(c)
    assert full['feasible']
    witness=full['boundary_witness']; assert witness['resolved'] and witness['surface']=='floor'
    from scipy.optimize import minimize_scalar
    point=np.asarray(witness['point'])
    u=minimize_scalar(lambda u: np.sum((curve.point(u)-point)**2),bounds=(0,1),method='bounded',
        options={'xatol':1e-14}).x
    assert np.linalg.norm(curve.point(u)-point)<1e-8
    t=curve.tangent(u); normal=np.array([-t[1],t[0]])/np.linalg.norm(t)
    np.testing.assert_allclose(normal,witness['normal_RZ'],atol=1e-7)
    np.testing.assert_array_equal(old.baseline['wall'],wall)
    assert old.baseline_peak==peak and old.direct.boundary_guard.e is old.direct.base
    from verification.validate_tangent_curvature import geometry_audit
    assert geometry_audit(e,c,full)['independent_geometry_passed']
    narrow=c.copy();narrow[12]=-.02
    assert e.constraints(narrow)['gap_m']<e.gap_min
    assert not e.evaluate(narrow)['feasible']
    moved_imp=c.copy();moved_imp[16]=1e-4
    assert not e.constraints(moved_imp)['opposite_coefficients_zero']
    assert not e.evaluate(moved_imp)['feasible']
    e.rmin=.1
    assert not e.constraints(c)['passed']  # stricter radius must not bypass the affine curve
    assert e.direct.base.sources is old.direct.base.sources
