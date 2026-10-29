import inspect
import json
from pathlib import Path
import numpy as np
import pytest
from shapely.geometry import LineString
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.case import separatrix
from d3d_optimization.frozen_mapping import first_hits

ROOT = Path(__file__).resolve().parents[1]

def continuous(source='IMP', count=48, nodes=41):
    assert 'boundary_mode' in inspect.signature(prepare_evaluator).parameters, 'continuous boundary policy missing'
    return prepare_evaluator(ROOT, source, count=count, nodes=nodes, boundary_mode='continuous')

def test_continuous_imp_boundary_covers_strike_and_vertex_without_sample_gaps():
    e, meta = continuous()
    surfaces = {s['original_segment']:np.asarray(s['endpoints_RZ_m']) for s in meta['surfaces']}
    review = json.loads((ROOT/'baseline_review/baseline_review.json').read_text())
    _, _, psi = separatrix(e.field, e.lcfs, review['axis_RZ_m'])
    assert abs(e.field.values([surfaces[69][0]])[0]-psi)<1e-12
    np.testing.assert_allclose(surfaces[69][1], surfaces[70][0], atol=1e-14)
    assert surfaces[70][-1,1] == pytest.approx(-1.15)
    assert {53,54,55}.isdisjoint(surfaces)
    assert LineString(e.wall).hausdorff_distance(LineString(review['wall']))<1e-12

def test_clamped_imp_strike_cannot_claim_sampled_peak_reduction():
    e, meta = continuous()
    c = np.zeros(e.dimension)
    j = next(j for j,s in enumerate(e.surfaces) if s['original_segment']==69)
    c[4*j:4*j+4] = .001
    r = e.evaluate(c)
    assert r['feasible']
    assert r['boundary_witness']['resolved']
    np.testing.assert_allclose(r['boundary_witness']['q'], e.baseline['boundary_witness']['q'], rtol=1e-9)
    assert max(r['peak_q']/e.baseline_peaks)>=1-1e-9
    assert all(abs(s.weights.sum()-1)<1e-14 for s in e.sources)

@pytest.mark.parametrize('source',['IMP','OMP'])
def test_boundary_candidate_keeps_full_first_hit_accounting(source):
    e, meta = continuous(source)
    r = e.evaluate(np.full(e.dimension,.0002))
    assert r['geometry_valid'] and r['boundary_witness']['resolved']
    hit = first_hits(e.paths, r['wall'])
    np.testing.assert_array_equal(hit.segment, r['maps'][0].segment)
    np.testing.assert_allclose(hit.points, r['maps'][0].points, atol=1e-12)
    np.testing.assert_allclose(e.baseline['peak_q']/e.baseline_peaks, 1., atol=1e-14)
    assert abs(r['boundary_witness']['psi_residual'])<1e-11

def test_unresolved_witness_fails_closed():
    e, _ = continuous()
    e.boundary_guard.path = LineString([[2.5,0.],[2.5,.1]])
    r=e.evaluate(np.zeros(e.dimension))
    assert not r['feasible']
    assert not r['boundary_witness']['resolved']

def test_buffer_requires_continuous_mode_and_finite_nonnegative_distance():
    assert 'imp_buffer_m' in inspect.signature(prepare_evaluator).parameters
    for value in [-.001, float('nan'), .5]:
        with pytest.raises(ValueError):
            prepare_evaluator(ROOT,'IMP',count=8,nodes=9,boundary_mode='continuous',imp_buffer_m=value)
    with pytest.raises(ValueError):
        prepare_evaluator(ROOT,'IMP',count=8,nodes=9,imp_buffer_m=.005)

def test_search_summary_uses_boundary_peak_not_only_sampled_rays():
    from types import SimpleNamespace
    from optimize_independent_sources import summary
    e=SimpleNamespace(baseline_peaks=np.ones(5))
    r=dict(geometry_valid=True,feasible=True,maps=[SimpleNamespace(q=[.8])]*5,
           peak_q=np.ones(5),angles={},surface_power={},missing=0.,power_balance_error=0.)
    assert summary(e,r)['worst_all_wall_peak_ratio']==1.

def test_numerical_acceptance_rejects_grid_dependent_gain():
    import importlib.util
    spec=importlib.util.find_spec('verification.validate_continuous_run')
    assert spec is not None, 'multi-resolution acceptance missing'
    from verification.validate_continuous_run import acceptance
    levels=[dict(feasible=True,audit_passed=True,peak_ratios=[x]*5) for x in [.98,.99,1.,1.,1.]]
    verdict=acceptance(levels)
    assert verdict['numerically_verified_feasible']
    assert not verdict['improvement_confirmed']
    levels=[dict(feasible=True,audit_passed=True,peak_ratios=[x]*5) for x in [.65,.645,.644,.6435,.6438]]
    assert acceptance(levels)['improvement_confirmed']
    levels[3]['audit_passed']=False
    assert not acceptance(levels)['numerically_verified_feasible']
