"""Deadline-bounded, independent IMP/OMP CPD geometry searches.

Preparation and baseline acceptance precede the campaign clock. Search and
320/640-donor qualification share that clock. Individual numerical calls cannot
be preempted safely: measured conservative cost guards avoid starting calls near
the deadline, and any unexpected overrun is recorded without truncating writes.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
from time import perf_counter

import numpy as np


SOURCES = ('IMP', 'OMP')
SEARCH_COUNT = 160
QUALIFICATION_COUNTS = (320, 640)


def numerical_settings(count):
    """Joint donor, trajectory, receiver and thermal-quadrature refinement."""
    if count == 160:
        return dict(path_step=.002, wall_step=.003, arc_step=.003, quadrature=64)
    if count == 320:
        return dict(path_step=.001, wall_step=.001, arc_step=.001, quadrature=128)
    if count == 640:
        return dict(path_step=.0005, wall_step=.0005, arc_step=.0005, quadrature=256)
    raise ValueError('Supported donor counts are 160, 320 and 640')


def _jsonable(value):
    if isinstance(value, dict): return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray): return _jsonable(value.tolist())
    if isinstance(value, np.generic): return _jsonable(value.item())
    if isinstance(value, float) and not math.isfinite(value): return None
    if isinstance(value, Path): return str(value)
    return value


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    with temporary.open('w', encoding='utf-8') as stream:
        json.dump(_jsonable(value), stream, indent=2, allow_nan=False)
        stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, path)


def _finite_number(value):
    return isinstance(value, (int, float, np.number)) and np.isfinite(value)


def _common_reasons(description):
    reasons = []
    if not description.get('feasible', False): reasons.append('evaluation_infeasible')
    for key in ('peak_MW_m2', 'area_90_m2'):
        if not _finite_number(description.get(key)) or description[key] <= 0:
            reasons.append('invalid_'+key)
    ledger = description.get('ledger', {})
    required = ('prescribed', 'deposited', 'unhit', 'unresolved', 'balance_error')
    if any(not _finite_number(ledger.get(key)) or ledger[key] < 0 for key in required) or ledger.get('prescribed', 0) <= 0:
        reasons.append('invalid_power_ledger')
    else:
        prescribed = ledger['prescribed']
        balance = abs(ledger['deposited']+ledger['unresolved']-prescribed)/prescribed
        if max(balance, ledger['balance_error']) > 2e-12:
            reasons.append('power_balance_failed')
        # This receiver model has no unresolved route other than missing hits.
        if abs(ledger['unresolved']-ledger['unhit'])/prescribed > 2e-12:
            reasons.append('other_unresolved_power')
    return reasons


def baseline_reasons(description):
    reasons = _common_reasons(description)
    ledger = description.get('ledger', {})
    if 'invalid_power_ledger' not in reasons and ledger['unhit']/ledger['prescribed'] > 1e-8:
        reasons.append('baseline_unhit_power')
    hit = description.get('divertor_first_hit_power_fraction')
    if not _finite_number(hit) or hit < .99:
        reasons.append('baseline_divertor_first_hit_below_0p99')
    return reasons


def candidate_reasons(description, baseline):
    """Fixed-power global-peak eligibility, independent of optimizer decisions."""
    reasons = _common_reasons(description)
    ledger = description.get('ledger', {})
    base = baseline.get('ledger', {})
    if 'invalid_power_ledger' not in reasons:
        if abs(ledger['prescribed']-base['prescribed'])/base['prescribed'] > 2e-12:
            reasons.append('assigned_power_changed')
        if ledger['unhit']/ledger['prescribed'] > base['unhit']/base['prescribed']+1e-10:
            reasons.append('unhit_power_increased')
    if _finite_number(description.get('area_90_m2')) and description['area_90_m2'] < .99*baseline['area_90_m2']:
        reasons.append('footprint_area_decreased')
    return reasons


def qualification_verdict(rows):
    """Require both paired resolutions and gain above a conservative error estimate."""
    rows = {int(key): value for key, value in rows.items()}
    complete = all(count in rows and 'baseline' in rows[count] and 'candidate' in rows[count]
                   for count in QUALIFICATION_COUNTS)
    result = dict(completed=complete, qualified_improvement=False, reasons=[])
    if not complete:
        result['reasons'] = ['qualification_incomplete']
        return result
    for count in QUALIFICATION_COUNTS:
        row = rows[count]
        result['reasons'].extend(f'{count}:{reason}' for reason in baseline_reasons(row['baseline']))
        result['reasons'].extend(f'{count}:{reason}' for reason in candidate_reasons(row['candidate'], row['baseline']))
    if result['reasons']: return result
    b0, b1 = (rows[count]['baseline']['peak_MW_m2'] for count in QUALIFICATION_COUNTS)
    c0, c1 = (rows[count]['candidate']['peak_MW_m2'] for count in QUALIFICATION_COUNTS)
    gain = b1-c1
    uncertainty = abs(b1-b0)+abs(c1-c0)
    result.update(qualified_count=640, peak_gain_MW_m2=gain, peak_gain_fraction=gain/b1,
        estimated_numerical_uncertainty_MW_m2=uncertainty,
        estimated_numerical_uncertainty_fraction=uncertainty/b1,
        uncertainty_method='sum of absolute baseline and candidate combined finite-resolution peak differences (320/1mm/128 to 640/0.5mm/256); not a rigorous bound',
        area_90_ratio=rows[640]['candidate']['area_90_m2']/rows[640]['baseline']['area_90_m2'])
    if gain <= uncertainty+1e-12*max(b1, c1):
        result['reasons'].append('gain_not_above_numerical_uncertainty')
    else:
        result['qualified_improvement'] = True
    return result


class PatternProposals:
    """Seeded coordinate/pattern search; every DOF is tried before focusing."""
    def __init__(self, evaluator, baseline_result, seed):
        family = evaluator.deformation
        self.dimension, self.bound = family.dimension, family.bound
        self.slices = list(family.slices.values())
        self.face_names = list(family.slices)
        self.lengths = np.array([family.curves[name].reference_length_m for name in self.face_names])
        self.face_steps = np.clip(.01*self.lengths, .001, .005)
        self.rng = np.random.default_rng(seed)
        self.iteration = 0; self.stagnation = 0; self.scale = 1.
        self.face_probabilities = self._face_probabilities(evaluator, baseline_result)

    def _face_probabilities(self, evaluator, result):
        # Keep a nonzero exploration probability on every original active face.
        power = np.zeros(len(self.slices))
        try:
            face_ids = np.asarray(result['divertor_face_ids'])
            dep = next(d for d in result['deposition']['depositions'] if d['loop_index'] == evaluator.divertor_loop)
            cell_faces = face_ids[dep['cell_parent']]
            for j, name in enumerate(self.face_names):
                start_node = int(name.rsplit(':', 1)[1].split('->')[0])
                power[j] = np.asarray(dep['cell_power'])[cell_faces == start_node].sum()
        except (KeyError, AttributeError, StopIteration, ValueError, IndexError):
            pass
        return .25/len(power)+(.75*power/power.sum() if power.sum() > 0 else np.full(len(power), .75/len(power)))

    def next(self, best):
        c = np.array(best, float, copy=True)
        i = self.iteration; self.iteration += 1
        if i < 2*self.dimension:
            index = (i//2) % self.dimension
            face = next(j for j, sl in enumerate(self.slices) if sl.start <= index < sl.stop)
            c[index] += (1 if i % 2 == 0 else -1)*self.face_steps[face]
        else:
            # Every fifth proposal uses all faces uniformly; otherwise use
            # baseline deposited-power evidence to prioritize useful faces.
            probabilities = None if i % 5 == 0 else self.face_probabilities
            face = int(self.rng.choice(len(self.slices), p=probabilities))
            sl = self.slices[face]
            step = self.face_steps[face]*self.scale
            long_faces = np.flatnonzero(self.lengths >= .5)
            if i % 20 == 19 and len(long_faces):
                face = int(self.rng.choice(long_faces)); sl = self.slices[face]
                c[sl] = self.rng.uniform(-self.bound, self.bound, sl.stop-sl.start)
            elif i % 5 == 4:
                c[sl] += self.rng.normal(0, step, sl.stop-sl.start)
            else:
                index = int(self.rng.integers(sl.start, sl.stop))
                c[index] += self.rng.choice([-1., 1.])*step*self.rng.choice([.25, .5, 1., 2., 4.])
        return np.clip(c, -self.bound, self.bound)

    def feedback(self, accepted):
        self.stagnation = 0 if accepted else self.stagnation+1
        if self.stagnation >= 2*self.dimension:
            self.scale = max(self.scale*.7, .01)
            self.stagnation = 0


class CampaignRunner:
    """Factory(source, donor_count) and monotonic clock are injectable for tests."""
    def __init__(self, factory, output, *, seconds=10800., reserve=900., seed=20261103,
                 clock=perf_counter, progress=False):
        if not np.isfinite(seconds) or seconds <= 0 or not np.isfinite(reserve) or not 0 < reserve < seconds:
            raise ValueError('positive duration and 0 < reserve < duration required')
        self.factory = factory; self.output = Path(output)
        self.seconds = float(seconds); self.reserve = float(reserve); self.seed = int(seed)
        self.clock = clock; self.progress = progress; self.started = None
        self.evaluators = {}; self.proposals = {}; self.costs = {}; self.io_cost = 0.
        self.last_progress = -math.inf
        self.report = dict(status='preparing_baselines', campaign_started=False,
            created_utc=datetime.now(timezone.utc).isoformat(), requested_seconds=self.seconds,
            qualification_reserve_seconds=self.reserve, seed=self.seed,
            search_count=SEARCH_COUNT, qualification_counts=list(QUALIFICATION_COUNTS),
            numerical_settings={n:numerical_settings(n) for n in (SEARCH_COUNT,)+QUALIFICATION_COUNTS},
            independent_sources=True, sources={}, search_evaluations=0,
            algorithm='seeded coordinate and face-pattern search; broad long-face exploration',
            objective='minimum global receiving-surface peak at fixed assigned power',
            footprint_guard='area_90_m2 >= 0.99 of same-resolution immutable baseline',
            deadline_policy='conservative measured-cost guards; no unsafe mid-evaluation termination',
            deadline_events=[])

    def _elapsed(self):
        return 0. if self.started is None else max(0., self.clock()-self.started)

    def _checkpoint(self):
        start = self.clock()
        self.report['campaign_elapsed_seconds'] = self._elapsed()
        self.report['campaign_remaining_seconds'] = max(0., self.seconds-self._elapsed())
        self.report['updated_utc'] = datetime.now(timezone.utc).isoformat()
        atomic_json(self.output/'study.json', self.report)
        status = {key: self.report[key] for key in ('status','campaign_started','campaign_elapsed_seconds',
            'campaign_remaining_seconds','search_evaluations','updated_utc')}
        status['sources'] = {name: dict(evaluations=data.get('evaluations', 0),
            best_peak_MW_m2=data.get('best_search', {}).get('peak_MW_m2'),
            qualified_improvement=data.get('qualification', {}).get('qualified_improvement', False))
            for name, data in self.report['sources'].items()}
        atomic_json(self.output/'status.json', status)
        self.io_cost = max(self.io_cost, self.clock()-start)
        if self.progress and (self.clock()-self.last_progress >= 60 or self.report['status'] != 'searching'):
            print(json.dumps(status), flush=True); self.last_progress = self.clock()

    def _timed(self, source, count, kind, operation):
        start = self.clock()
        try: return operation()
        finally:
            elapsed = max(1e-6, self.clock()-start)
            key = (source, count, kind)
            self.costs[key] = max(self.costs.get(key, 0.), elapsed)

    def _estimate(self, source, count, kind):
        choices = [(abs(count-n), n, cost) for (s,n,k),cost in self.costs.items() if s == source and k == kind]
        if not choices: return .01
        _, previous, cost = min(choices)
        ratio = max(1., count/previous)
        # Cost scales with donor count and receiving/trajectory refinement;
        # subsequent measurements at a finer resolution replace this forecast.
        key = 'wall_step' if kind == 'evaluate' else 'path_step'
        refinement = max(1., numerical_settings(previous)[key]/numerical_settings(count)[key])
        return 1.5*cost*ratio*refinement

    def _can_start(self, source, count, kind, cutoff):
        margin = self._estimate(source, count, kind)+max(2*self.io_cost, min(.25, .005*self.seconds))
        if self.clock()+margin <= cutoff: return True
        self.report['deadline_events'].append(dict(stage=self.report['status'], source=source,
            count=count, operation=kind, remaining_seconds=max(0., cutoff-self.clock()),
            estimated_required_seconds=margin))
        return False

    def _export(self, evaluator, filename, result):
        started = self.clock()
        path = self.output/filename
        temporary = path.with_name(path.stem+'.tmp.npz')
        try:
            evaluator.export(temporary, result)
            os.replace(temporary, path)
        finally:
            self.io_cost = max(self.io_cost, self.clock()-started)
        return str(path.resolve())

    def _record_candidate(self, row):
        row.update(equilibrium_sha256=self.report.get('equilibrium_sha256'),
                   campaign_elapsed_seconds=self._elapsed())
        with (self.output/'candidates.jsonl').open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(_jsonable(row), allow_nan=False)+'\n')
            stream.flush(); os.fsync(stream.fileno())

    def _check_metadata(self, evaluator):
        for key in ('equilibrium_sha256', 'geometry_sha256'):
            value = evaluator.metadata.get(key)
            if not value: raise ValueError('Missing '+key)
            if key in self.report and value != self.report[key]:
                raise ValueError('Fixed '+key+' changed between evaluators')
            self.report[key] = value
        fixed = self.report.setdefault('fixed_inputs', {})
        source = evaluator.metadata.get('source')
        source_fixed = fixed.setdefault('sources', {}).setdefault(source, {})
        for key in ('reference', 'source_RZ_m', 'source_support_m'):
            container = fixed if key == 'reference' else source_fixed
            if key not in evaluator.metadata:
                if key in container: raise ValueError('Fixed '+key+' missing from later evaluator')
                continue
            value = _jsonable(evaluator.metadata[key])
            if key in container and value != container[key]:
                raise ValueError('Fixed '+key+' changed between evaluators')
            container[key] = deepcopy(value)

    def _prepare(self):
        preparation_start = self.clock()
        for index, source in enumerate(SOURCES):
            evaluator = self._timed(source, SEARCH_COUNT, 'prepare', lambda: self.factory(source, SEARCH_COUNT))
            self._check_metadata(evaluator)
            baseline = self._timed(source, SEARCH_COUNT, 'evaluate', lambda: evaluator.evaluate(None))
            desc = evaluator.describe(baseline)
            reasons = baseline_reasons(desc)
            data = dict(metadata=deepcopy(evaluator.metadata), baseline=deepcopy(desc),
                baseline_acceptance_reasons=reasons, evaluations=0, accepted_updates=0,
                qualification=qualification_verdict({}), qualification_rows={})
            self.report['sources'][source] = data
            if 'deposition' in baseline:
                data['baseline_artifact'] = self._export(evaluator, f'baseline-{source}-160.npz', baseline)
            if reasons:
                self.report['status'] = 'baseline_rejected'
                self.report['preparation_seconds'] = self.clock()-preparation_start
                self._checkpoint(); return False
            self.evaluators[source] = evaluator
            data['best_coefficients'] = np.zeros(evaluator.deformation.dimension).tolist()
            data['deformation_basis'] = deepcopy(getattr(evaluator.deformation, 'basis_description', None))
            data['best_search'] = deepcopy(desc)
            data['best_search_artifact'] = self._export(evaluator, f'best-{source}-search-000000.npz', baseline)
            self.proposals[source] = PatternProposals(evaluator, baseline, self.seed+index)
            data['proposal_face_probabilities'] = self.proposals[source].face_probabilities.tolist()
            self._checkpoint()
        self.report['preparation_seconds'] = self.clock()-preparation_start
        return True

    def _search(self):
        cutoff = self.started+self.seconds-self.reserve
        self.report['status'] = 'searching'; self._checkpoint()
        turn = 0
        while True:
            source = SOURCES[turn % len(SOURCES)]; turn += 1
            if not self._can_start(source, SEARCH_COUNT, 'evaluate', cutoff): break
            evaluator = self.evaluators[source]; data = self.report['sources'][source]
            coefficients = self.proposals[source].next(data['best_coefficients'])
            accepted = False; artifact = None
            try:
                result = self._timed(source, SEARCH_COUNT, 'evaluate', lambda: evaluator.evaluate(coefficients))
                desc = evaluator.describe(result)
                reasons = candidate_reasons(desc, data['baseline'])
                accepted = not reasons and desc['peak_MW_m2'] < data['best_search']['peak_MW_m2']
                if accepted:
                    artifact = self._export(evaluator, f"best-{source}-search-{data['evaluations']+1:06d}.npz", result)
                    data.update(best_coefficients=coefficients.tolist(), best_search=deepcopy(desc),
                                best_search_artifact=artifact, best_found_elapsed_seconds=self._elapsed())
                    data['accepted_updates'] += 1
            except Exception as exc:
                accepted = False
                desc = dict(feasible=False, errors=[f'{type(exc).__name__}: {exc}'])
                reasons = ['candidate_evaluation_failed']
            data['evaluations'] += 1; self.report['search_evaluations'] += 1
            self.proposals[source].feedback(accepted)
            self._record_candidate(dict(source=source, count=SEARCH_COUNT,
                evaluation=data['evaluations'], coefficients=coefficients.tolist(),
                accepted=accepted, acceptance_reasons=reasons, description=desc,
                artifact=artifact if accepted else None))
            self._checkpoint()
        self.report['search_elapsed_seconds'] = self._elapsed()

    def _qualify(self):
        self.report['status'] = 'qualifying'; self._checkpoint()
        cutoff = self.started+self.seconds
        for count in QUALIFICATION_COUNTS:
            for source in SOURCES:
                data = self.report['sources'][source]
                if not self._can_start(source, count, 'prepare', cutoff): return
                row = {}; data['qualification_rows'][str(count)] = row
                try:
                    evaluator = self._timed(source, count, 'prepare', lambda: self.factory(source, count))
                    self._check_metadata(evaluator)
                    row['metadata'] = deepcopy(evaluator.metadata)
                    if not self._can_start(source, count, 'evaluate', cutoff): return
                    baseline = self._timed(source, count, 'evaluate', lambda: evaluator.evaluate(None))
                    row['baseline'] = evaluator.describe(baseline)
                    if 'deposition' in baseline:
                        row['baseline_artifact'] = self._export(evaluator, f'baseline-{source}-{count}.npz', baseline)
                    reasons = baseline_reasons(row['baseline'])
                    if reasons:
                        row['errors'] = reasons; self._checkpoint(); continue
                    if not self._can_start(source, count, 'evaluate', cutoff): return
                    coefficients = np.asarray(data['best_coefficients'])
                    candidate = self._timed(source, count, 'evaluate', lambda: evaluator.evaluate(coefficients))
                    row['candidate'] = evaluator.describe(candidate)
                    row['coefficients'] = coefficients.tolist()
                    if 'deposition' in candidate:
                        row['candidate_artifact'] = self._export(evaluator, f'candidate-{source}-{count}.npz', candidate)
                    self._record_candidate(dict(source=source, count=count, stage='qualification',
                        coefficients=coefficients.tolist(), description=row['candidate'],
                        acceptance_reasons=candidate_reasons(row['candidate'], row['baseline'])))
                except Exception as exc:
                    row['errors'] = [f'{type(exc).__name__}: {exc}']
                data['qualification'] = qualification_verdict(data['qualification_rows'])
                self._checkpoint()

    def run(self):
        self.output.mkdir(parents=True, exist_ok=True)
        if (self.output/'study.json').exists():
            raise FileExistsError('Choose a new output directory; saved studies are immutable')
        self._checkpoint()
        try:
            if not self._prepare(): return self.report
            self.started = self.clock()
            self.report.update(campaign_started=True, campaign_started_utc=datetime.now(timezone.utc).isoformat())
            self._search(); self._qualify()
            for data in self.report['sources'].values():
                data['qualification'] = qualification_verdict(data['qualification_rows'])
            complete = all(data['qualification']['completed'] for data in self.report['sources'].values())
            self.report['status'] = 'completed' if complete else 'time_expired_without_qualification'
            if not complete and any(row.get('errors') for data in self.report['sources'].values() for row in data['qualification_rows'].values()):
                self.report['status'] = 'qualification_failed'
            self.report['deadline_overrun_seconds'] = max(0., self._elapsed()-self.seconds)
            if self.report['deadline_overrun_seconds'] > 0:
                self.report['status'] = 'completed_after_deadline' if complete else 'time_expired_without_qualification'
        except KeyboardInterrupt:
            self.report['status'] = 'interrupted'
        except Exception as exc:
            self.report['status'] = 'preparation_failed' if self.started is None else 'campaign_failed'
            self.report['failure'] = f'{type(exc).__name__}: {exc}'
        finally:
            self._checkpoint()
        return self.report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes', default='input/CPD_5/CPD_5_nodes.yaml')
    parser.add_argument('--equilibrium', required=True, help='Validated frozen equilibrium NPZ')
    parser.add_argument('--output', required=True, help='New study directory, normally on D:')
    parser.add_argument('--cache-dir')
    parser.add_argument('--seconds', type=float, default=10800.)
    parser.add_argument('--reserve', type=float, default=900.)
    parser.add_argument('--seed', type=int, default=20261103)
    parser.add_argument('--smoke-seconds', type=float, help='Explicit short smoke budget; overrides --seconds')
    args = parser.parse_args()
    seconds = args.seconds if args.smoke_seconds is None else args.smoke_seconds
    if args.smoke_seconds is not None and args.reserve >= seconds:
        args.reserve = .3*seconds
    from cpd_optimization.evaluator import CPDEvaluator
    cache = Path(args.cache_dir) if args.cache_dir else Path(args.output).parent/'path_cache'
    def factory(source, count):
        return CPDEvaluator(args.nodes, args.equilibrium, source, count=count,
            **numerical_settings(count), max_path=30.,
            cache_dir=cache)
    report = CampaignRunner(factory, args.output, seconds=seconds, reserve=args.reserve,
                            seed=args.seed, progress=True).run()
    print(json.dumps(_jsonable(dict(status=report['status'], output=str(Path(args.output).resolve()),
        campaign_elapsed_seconds=report['campaign_elapsed_seconds'],
        qualification={source:data.get('qualification') for source,data in report['sources'].items()})), indent=2), flush=True)
    return 0 if report['status'] == 'completed' else 2


if __name__ == '__main__':
    raise SystemExit(main())
