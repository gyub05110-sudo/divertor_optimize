import copy
from decimal import Decimal

import pytest

from d3d_optimization.selective_peak import replay


def kernel(centers=(.25, .75), sigma=.02, power=.5):
    return {'period_m2': [1, 1], 'donors': [
        {'weight_MW': [power, power], 'center_fraction': [x, x],
         'sigma_m2': [sigma, sigma]} for x in centers]}


def test_builder_one_operation_monotonic_replay_and_owned_source():
    from d3d_optimization.adaptive_peak import PeakBuilder
    source = kernel()
    builder = PeakBuilder(source, 'same')
    original = copy.deepcopy(source)
    source['donors'][0]['weight_MW'][0] = 900
    last = builder.certificate()
    assert replay(last) and last['model'] == original
    for count in range(1, 13):
        assert builder.advance()
        assert replay(last)  # Later advances cannot mutate an exported snapshot.
        current = builder.certificate()
        assert len(current['operations']) == count
        assert replay(current)
        lo, hi = map(Decimal, current['peak_bracket_MW_m2'])
        old_lo, old_hi = map(Decimal, last['peak_bracket_MW_m2'])
        assert old_lo <= lo <= hi <= old_hi
        last = current
    last['model']['donors'].clear()
    last['operations'].clear()
    assert replay(builder.certificate())
    assert builder.certificate()['model'] == original


def test_builder_initialization_is_lazy_and_zero_cap_retains_certificate(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    def forbidden(*args):
        pytest.fail('zero operations must not construct the evaluator')
    monkeypatch.setattr(ap.selective_peak, 'SpatialBounds', forbidden)
    builder = ap.PeakBuilder(kernel(), 'same', max_operations=0)
    assert not builder.advance()
    assert replay(builder.certificate())


def test_separated_peak_rank_resolves_early_and_stops():
    from d3d_optimization.adaptive_peak import rank
    result = rank({'spread': kernel(), 'coincident': kernel((.5, .5))},
                  {'spread': 'same', 'coincident': 'same'},
                  [('spread', 'coincident')], seconds=10, max_operations=128)
    pair = result['comparisons'][0]
    assert pair['status'] == 'a_lower'
    assert result['stop_reason'] == 'all_comparisons_resolved'
    assert 0 < result['completed_operations'] < 128
    assert pair['provisional_resolved_after_operation'] == result['completed_operations']
    assert result['replay_count'] == 2
    assert all(replay(cert) for cert in result['certificates'].values())
    assert result['end_to_end_seconds'] >= result['generation_seconds'] + result['replay_seconds']


def test_shared_baseline_replayed_once_and_scheduler_is_deterministic(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    models = {'baseline': kernel(), 'b': kernel(), 'c': kernel()}
    calls = []
    real_replay = ap.selective_peak.replay
    def counted(record):
        calls.append(record)
        return real_replay(record)
    monkeypatch.setattr(ap.selective_peak, 'replay', counted)
    result = ap.rank(models, dict.fromkeys(models, 'same'),
                     [('baseline', 'b'), ('baseline', 'c')], seconds=10, max_operations=8)
    assert len(calls) == result['replay_count'] == 3
    assert result['completed_operations'] == 8
    assert result['stop_reason'] == 'operation_budget'
    assert result['schedule'][0]['label'] == 'baseline'
    assert all(p['status'] == 'unresolved' for p in result['comparisons'])
    again = ap.rank(models, dict.fromkeys(models, 'same'),
                    [('baseline', 'b'), ('baseline', 'c')], seconds=10, max_operations=8)
    assert [(s['label'], s['operation']) for s in result['schedule']] == [
        (s['label'], s['operation']) for s in again['schedule']]


@pytest.mark.parametrize('seconds,operations', [(0, 50), (10, 0)])
def test_zero_budget_valid_unresolved(seconds, operations):
    from d3d_optimization.adaptive_peak import rank
    models = {'a': kernel(), 'b': kernel()}
    result = rank(models, dict.fromkeys(models, 'same'), [('a', 'b')],
                  seconds=seconds, max_operations=operations)
    assert result['completed_operations'] == 0
    assert result['comparisons'][0]['status'] == 'unresolved'
    assert result['comparisons'][0]['provisional_resolved_after_operation'] is None
    assert all(replay(cert) for cert in result['certificates'].values())


def test_zero_budget_can_prove_universally_separated_peaks():
    from d3d_optimization.adaptive_peak import rank
    result = rank({'a': kernel(power=0), 'b': kernel()}, {'a': 'same', 'b': 'same'},
                  [('a', 'b')], seconds=0)
    assert result['comparisons'][0]['status'] == 'a_lower'
    assert result['comparisons'][0]['provisional_resolved_after_operation'] == 0


@pytest.mark.parametrize('pairs,identities', [
    ([('a', 'a')], {'a': 'same', 'b': 'same'}),
    ([('a', 'missing')], {'a': 'same', 'b': 'same'}),
    ([('a', 'b'), ('b', 'a')], {'a': 'same', 'b': 'same'}),
    ([('a', 'b'), ('a', 'b')], {'a': 'same', 'b': 'same'}),
    ([('a',)], {'a': 'same', 'b': 'same'}),
    ([('a', 'b')], {'a': 'same'}),
    ([('a', 'b')], {'a': '', 'b': ''}),
    ([('a', 'b')], {'a': '320', 'b': '1280'}),
])
def test_invalid_pairs_and_identities_fail_before_numerics(pairs, identities, monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    def forbidden(*args, **kwargs):
        pytest.fail('validate pair identities before numerical work')
    monkeypatch.setattr(ap, 'PeakBuilder', forbidden)
    with pytest.raises(ValueError):
        ap.rank({'a': kernel(), 'b': kernel()}, identities, pairs)


@pytest.mark.parametrize('kwargs', [
    {'seconds': -1}, {'seconds': float('nan')}, {'seconds': float('inf')},
    {'max_operations': -1}, {'max_operations': 1.5}, {'max_operations': True},
    {'max_operations': 10001},
])
def test_invalid_budget(kwargs):
    from d3d_optimization.adaptive_peak import rank
    with pytest.raises(ValueError):
        rank({'a': kernel(), 'b': kernel()}, {'a': 'same', 'b': 'same'}, [('a', 'b')], **kwargs)


def test_forged_final_certificate_is_rejected(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    real_certificate = ap.PeakBuilder.certificate
    def forged(self, *args, **kwargs):
        cert = real_certificate(self, *args, **kwargs)
        cert['peak_bracket_MW_m2'] = ['0', '0']
        return cert
    monkeypatch.setattr(ap.PeakBuilder, 'certificate', forged)
    with pytest.raises(ValueError, match='replay'):
        ap.rank({'a': kernel(), 'b': kernel()}, {'a': 'same', 'b': 'same'}, [('a', 'b')], seconds=0)


def test_final_replay_failure_propagates(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    def reject(record):
        raise ValueError('independent replay rejected')
    monkeypatch.setattr(ap.selective_peak, 'replay', reject)
    with pytest.raises(ValueError, match='independent replay'):
        ap.rank({'a': kernel(), 'b': kernel()}, {'a': 'same', 'b': 'same'}, [('a', 'b')], seconds=0)


def test_near_tie_keeps_overlap_and_no_provisional_winner():
    from d3d_optimization.adaptive_peak import rank
    result = rank({'a': kernel(), 'b': kernel(power=.500000000001)},
                  {'a': 'same', 'b': 'same'}, [('a', 'b')], seconds=10, max_operations=8)
    assert result['comparisons'][0]['status'] == 'unresolved'
    assert result['comparisons'][0]['provisional_resolved_seconds'] is None
    assert result['completed_operations'] == 8


def test_resolved_pairs_leave_scheduler_and_distinct_identity_groups_work():
    from d3d_optimization.adaptive_peak import rank
    models = {'zero': kernel(power=0), 'a': kernel(), 'b': kernel(), 'c': kernel()}
    result = rank(models, {'zero': 'IMP', 'a': 'IMP', 'b': 'OMP', 'c': 'OMP'},
                  [('zero', 'a'), ('b', 'c')], seconds=10, max_operations=5)
    assert result['comparisons'][0]['status'] == 'a_lower'
    assert result['comparisons'][0]['provisional_resolved_after_operation'] == 0
    assert {op['label'] for op in result['schedule']} <= {'b', 'c'}
    assert result['replay_count'] == 4


def test_rank_owns_models_identities_and_pairs_during_generation(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    models = {'a': kernel(), 'b': kernel()}
    identities = {'a': 'same', 'b': 'same'}
    pairs = [('a', 'b')]
    originals = copy.deepcopy(models)
    real_advance = ap.PeakBuilder.advance
    def mutate_sources(self):
        models['a']['donors'].clear()
        identities['a'] = 'other'
        pairs.clear()
        return real_advance(self)
    monkeypatch.setattr(ap.PeakBuilder, 'advance', mutate_sources)
    result = ap.rank(models, identities, pairs, seconds=10, max_operations=4)
    assert len(result['comparisons']) == 1
    assert result['comparisons'][0]['problem_id'] == 'same'
    for label, record in result['certificates'].items():
        assert record['model'] == originals[label]
        assert replay(record)


def test_valid_certificate_for_wrong_problem_or_model_cannot_replace_owned_evidence(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    from d3d_optimization.selective_peak import certify
    replacement = certify(kernel(power=0), seconds=0, problem_id='same')
    assert replay(replacement)
    monkeypatch.setattr(ap.PeakBuilder, 'certificate', lambda self, **kw: copy.deepcopy(replacement))
    with pytest.raises(ValueError, match='identity'):
        ap.rank({'a': kernel(), 'b': kernel()}, {'a': 'same', 'b': 'same'}, [('a', 'b')], seconds=0)


def test_coordinate_exhaustion_retains_replayable_partition():
    from d3d_optimization.adaptive_peak import PeakBuilder
    builder = PeakBuilder(kernel(centers=(.3,), sigma='1e-60'), 'same', max_operations=1000)
    for _ in range(1000):
        if builder.advance() is None:
            break
    record = builder.certificate()
    assert record['stop_reason'] == 'coordinate_precision'
    assert 0 < len(record['operations']) < 1000
    assert replay(record)


@pytest.mark.parametrize('center', [[-1, 0], [0, 2], ['NaN', 'NaN']])
def test_zero_budget_still_validates_centers(center):
    from d3d_optimization.adaptive_peak import PeakBuilder
    model = kernel()
    model['donors'][0]['center_fraction'] = center
    with pytest.raises(ValueError):
        PeakBuilder(model, 'same', max_operations=0)


def test_atomic_budget_overrun_still_replays_and_reports_all_time(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    clock = [0.]
    monkeypatch.setattr(ap, 'perf_counter', lambda: clock[0])
    operation = ap.selective_peak._operation
    def slow_operation(*args):
        operation(*args)
        clock[0] += 2.
    monkeypatch.setattr(ap.selective_peak, '_operation', slow_operation)
    result = ap.rank({'a': kernel(), 'b': kernel()}, {'a': 'same', 'b': 'same'},
                     [('a', 'b')], seconds=1, max_operations=100)
    assert result['completed_operations'] == 1
    assert result['stop_reason'] == 'time_budget'
    assert result['generation_seconds'] == 2
    assert result['budget_overrun_seconds'] == 1
    assert result['largest_atomic_operation_seconds'] == 2
    assert result['replay_count'] == 2
    assert result['replay_seconds'] == 2
    assert result['end_to_end_seconds'] == 4


def test_false_replay_result_is_not_trusted(monkeypatch):
    from d3d_optimization import adaptive_peak as ap
    monkeypatch.setattr(ap.selective_peak, 'replay', lambda record: False)
    with pytest.raises(ValueError, match='replay failed'):
        ap.rank({'a': kernel(), 'b': kernel()}, {'a': 'same', 'b': 'same'}, [('a', 'b')], seconds=0)
