from pathlib import Path
import numpy as np
import pytest
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.broadening import BroadeningConfig
from d3d_optimization.broadened_evaluator import BroadenedEvaluator

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module', params=['IMP', 'OMP'])
def direct(request):
    return prepare_evaluator(ROOT, request.param, count=48, nodes=41, boundary_mode='continuous')


def test_zero_spread_is_exact_dispatch(direct):
    e, meta = direct
    b = BroadenedEvaluator(e, meta, BroadeningConfig(0.))
    c = np.full(e.dimension, .0001)
    actual, expected = b.evaluate(c), e.evaluate(c)
    np.testing.assert_array_equal(actual['peak_q'], expected['peak_q'])
    assert actual['feasible'] == expected['feasible']
    for a, x in zip(actual['maps'], expected['maps']):
        np.testing.assert_array_equal(a.q, x.q)
        np.testing.assert_array_equal(a.segment_power, x.segment_power)


def test_positive_spread_has_independent_conservative_ledger(direct):
    e, meta = direct
    before = e.baseline['peak_q'].copy()
    weights = np.column_stack([s.weights.copy() for s in e.sources])
    b = BroadenedEvaluator(e, meta, BroadeningConfig(.0001, wall_step_m=.005))
    r = b.baseline
    assert r['geometry_valid']
    assert r['model_version'] == 'connection-deposition-surrogate.v1'
    np.testing.assert_array_equal(e.baseline['peak_q'], before)
    np.testing.assert_array_equal(np.column_stack([s.weights for s in e.sources]), weights)
    ledger = r['deposition']['ledger']
    np.testing.assert_allclose(ledger['prescribed'], np.ones(5), atol=1e-13)
    np.testing.assert_allclose(ledger['deposited'] + ledger['unresolved'], np.ones(5), atol=1e-12)
    assert np.all(r['deposition']['q'] >= 0)
    np.testing.assert_array_equal(r['peak_q'], r['deposition']['q'].max(axis=0))
    # No old unsmeared SOL-limit maximum is inserted in positive-spread peaks.
    assert 'direct_boundary_witness' in r and 'boundary_witness' not in r


def test_invalid_candidate_fails_without_deposition(direct):
    e, meta = direct
    b = BroadenedEvaluator(e, meta, BroadeningConfig(.0001, wall_step_m=.005))
    r = b.evaluate(np.full(e.dimension, e.bound))
    if not r['geometry_valid']:
        assert not r['feasible'] and 'deposition' not in r


def test_blocked_mass_cannot_be_accepted_by_adapter(direct, monkeypatch):
    e, meta = direct
    b = BroadenedEvaluator(e, meta, BroadeningConfig(.0001, wall_step_m=.005))
    import d3d_optimization.broadened_evaluator as module
    original = module.build_deposition
    def broken(*args, **kwargs):
        result = original(*args, **kwargs)
        result['ledger']['unresolved'][:] = .01
        result['valid'] = False
        result['reasons'] = ['missing_power']
        return result
    monkeypatch.setattr(module, 'build_deposition', broken)
    r = b.evaluate(np.zeros(e.dimension))
    assert not r['feasible'] and r['missing'] >= .01


def test_zero_q_grazing_sample_on_powered_cell_cannot_evade_angle_floor(direct, monkeypatch):
    e, meta = direct
    if meta['source'] != 'OMP':
        pytest.skip('OMP illuminates the angle-constrained shelf')
    b = BroadenedEvaluator(e, meta, BroadeningConfig(.0001, wall_step_m=.005))
    assert b.baseline['feasible']
    import d3d_optimization.broadened_evaluator as module
    original, incidence = module.build_deposition, e.field.incidence
    chosen = {}
    def grazing(*args, **kwargs):
        dep = original(*args, **kwargs)
        shelf = b.baseline['surface_segments']['shelf']
        powered = np.any(dep['cell_power'][dep['sample_cell']] > 0, axis=1)
        ids = np.flatnonzero(np.isin(dep['sample_parent'], shelf) & powered)
        assert len(ids)
        chosen['index'], chosen['points'] = ids[0], dep['points']
        dep['q'][ids[0]] = 0.
        return dep
    def incidence_with_zero(points, normals):
        angles = incidence(points,normals)
        if points is chosen.get('points'):
            angles[chosen['index']] = 0.
        return angles
    monkeypatch.setattr(module, 'build_deposition', grazing)
    monkeypatch.setattr(e.field, 'incidence', incidence_with_zero)
    r = b.evaluate(np.zeros(e.dimension))
    assert not r['feasible']
    assert not r['receiving_angles']['shelf']['passed']
