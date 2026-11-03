"""Incremental, replayed peak rankings for declared finite kernel families.

Scheduling is a runtime heuristic only. Final decisions use independently
replayed certificates, never provisional decisions or caller verification flags.
Feasibility, continuous transport, and global optimality remain separate.
"""
import copy
import math
from time import perf_counter

from . import selective_peak
from .enclosure_math import interval
from .geometry_transport_certificate import digest
from .transport_peak_bounds import D, DOWN, UP, I


MAX_OPERATIONS = 10000


def _operation_budget(value):
    if type(value) is not int or not 0 <= value <= MAX_OPERATIONS:
        raise ValueError('integer operation budget in 0..10000 required')


def _time_budget(value):
    try:
        valid = math.isfinite(value) and value >= 0
    except (TypeError, OverflowError):
        valid = False
    if not valid:
        raise ValueError('finite nonnegative time budget required')


def _identity(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('nonempty problem identity required')


class PeakBuilder:
    """Own a kernel copy and append one existing proof operation per advance.

    ``advance`` returns the completed operation (a detached dict), or None
    after the operation cap or coordinate precision is exhausted. The expensive
    spatial evaluator is initialized lazily within the first atomic advance.
    """

    def __init__(self, model, problem_id, *, max_operations=MAX_OPERATIONS):
        started = perf_counter()
        _identity(problem_id)
        _operation_budget(max_operations)
        self._model = copy.deepcopy(model)
        self._problem_id = problem_id
        self._max_operations = max_operations
        self._state = selective_peak._initial(self._model)
        # The universal bound does not inspect centers. Validate them even when
        # the generation budget is zero, without constructing SpatialBounds.
        area = interval(self._model['period_m2'])
        for row in self._model['donors']:
            fraction, absolute = row.get('center_fraction'), row.get('center_m2')
            if (fraction is None) == (absolute is None):
                raise ValueError('one center representation required')
            center = interval(fraction if fraction is not None else absolute)
            if center.lo < 0 or center.hi > (D(1) if fraction is not None else area.lo):
                raise ValueError('center outside period')
        self._model_sha256 = digest(self._model)
        self._identity_sha256 = digest({'model': self._model, 'problem_id': problem_id})
        self._evaluator = None
        self._points = None
        self._operations = []
        self._splits = 0
        self._stop_reason = 'incremental'
        self._elapsed = perf_counter() - started
        self._largest_atomic = self._elapsed

    def bracket(self):
        """Immutable current endpoints; these are provisional until replay."""
        return self._state['lower'], max(x['upper'] for x in self._state['leaves'])

    def advance(self):
        if len(self._operations) >= self._max_operations:
            self._stop_reason = 'operation_budget'
            return None
        if self._stop_reason == 'coordinate_precision':
            return None
        started = perf_counter()
        try:
            if self._evaluator is None:
                self._evaluator = selective_peak.SpatialBounds(self._model)
                ordered = sorted(self._model['donors'], key=lambda row: UP.divide(
                    I(*row['weight_MW']).hi, I(*row['sigma_m2']).lo), reverse=True)
                area = I(*self._model['period_m2'])
                points = [D(0), D('.5')]
                for row in ordered[:8]:
                    center = (I(*row['center_fraction']) if 'center_fraction' in row
                              else I(*row['center_m2']) / area)
                    points.append(max(D(0), min(D(1), DOWN.divide(
                        DOWN.add(center.lo, center.hi), D(2)))))
                self._points = list(dict.fromkeys(points))
            index = len(self._operations)
            if index < len(self._points):
                op = {'kind': 'probe', 'point': str(self._points[index])}
            else:
                active = max(self._state['leaves'], key=lambda x: (x['upper'], -x['id']))
                midpoint = DOWN.divide(DOWN.add(active['a'], active['b']), D(2))
                if not active['a'] < midpoint < active['b']:
                    self._stop_reason = 'coordinate_precision'
                    return None
                op = {'kind': 'split', 'id': active['id']}
            selective_peak._operation(self._state, self._evaluator, op)
            self._operations.append(op)
            self._splits += op['kind'] == 'split'
            return dict(op)
        finally:
            elapsed = perf_counter() - started
            self._elapsed += elapsed
            self._largest_atomic = max(self._largest_atomic, elapsed)

    def certificate(self, *, seconds_budget=0, stop_reason=None):
        """Export a detached existing-v1 certificate; this does not imply replay."""
        _time_budget(seconds_budget)
        seconds_budget = float(seconds_budget)
        evidence = selective_peak._evidence(self._state)
        return dict(
            schema='selective-spatial-peak.v1', model=copy.deepcopy(self._model),
            model_sha256=self._model_sha256, problem_id=self._problem_id,
            identity_sha256=self._identity_sha256, tolerance='.01',
            max_splits=self._max_operations, max_operations=self._max_operations,
            seconds_budget=seconds_budget, stop_reason=stop_reason or self._stop_reason,
            completed_splits=self._splits, operations=copy.deepcopy(self._operations),
            **evidence, tolerance_met=D(evidence['gap_upper_MW_m2']) <= D('.01'),
            elapsed_seconds=self._elapsed,
            budget_overrun_seconds=max(0., self._elapsed - seconds_budget),
            largest_atomic_operation_seconds=self._largest_atomic,
            scope='declared finite periodic-kernel family only',
            geometry_global_optimum_proven=False)


def _status(a, b):
    la, ua = a
    lb, ub = b
    return 'a_lower' if ua < lb else 'b_lower' if ub < la else 'unresolved'


def _validate_inputs(models, problem_ids, pairs):
    if not isinstance(models, dict) or not models:
        raise ValueError('nonempty labeled model dictionary required')
    if not isinstance(problem_ids, dict):
        raise ValueError('problem identity dictionary required')
    for label in models:
        if not isinstance(label, str) or not label.strip():
            raise ValueError('nonempty model labels required')
        if label not in problem_ids:
            raise ValueError('missing problem identity')
        _identity(problem_ids[label])
    if not isinstance(pairs, (list, tuple)):
        raise ValueError('list of comparison pairs required')
    seen = set()
    checked = []
    for pair in pairs:
        if (not isinstance(pair, (list, tuple)) or len(pair) != 2
                or any(not isinstance(x, str) or x not in models for x in pair)):
            raise ValueError('two known labels required per pair')
        a, b = pair
        key = frozenset(pair)
        if a == b or key in seen:
            raise ValueError('self and duplicate comparisons are not allowed')
        if problem_ids[a] != problem_ids[b]:
            raise ValueError('same problem identity required within each pair')
        seen.add(key)
        checked.append((a, b))
    return checked


def rank(models, problem_ids, pairs, seconds=136, max_operations=4096):
    """Share one builder/replay per label and refine unresolved pairs only.

    The shared generation deadline includes validation, initialization,
    scheduling, and provisional decisions. Initialization is also reported as
    a subcomponent. Final certificate output and mandatory replay are outside
    that deadline and included in end-to-end time. Deadline checks occur between
    atomic advances; one slow initialization/advance can overrun the budget.
    """
    started = perf_counter()
    _time_budget(seconds)
    seconds = float(seconds)
    _operation_budget(max_operations)
    pairs = _validate_inputs(models, problem_ids, pairs)
    identities = {label: problem_ids[label] for label in models}
    init_start = perf_counter()
    builders = {label: PeakBuilder(model, identities[label], max_operations=max_operations)
                for label, model in sorted(models.items())}
    initialization_seconds = perf_counter() - init_start
    largest_atomic = max(builder._largest_atomic for builder in builders.values())
    unresolved = set(range(len(pairs)))
    provisional = {}
    exhausted = set()
    schedule = []

    def update_provisional():
        for index in sorted(unresolved):
            a, b = pairs[index]
            if _status(builders[a].bracket(), builders[b].bracket()) != 'unresolved':
                provisional[index] = (len(schedule), perf_counter() - started)
                unresolved.remove(index)

    update_provisional()
    while unresolved:
        if perf_counter() - started >= seconds:
            stop_reason = 'time_budget'
            break
        if len(schedule) >= max_operations:
            stop_reason = 'operation_budget'
            break
        degrees = {}
        for index in sorted(unresolved):
            for label in pairs[index]:
                if label not in exhausted:
                    degrees[label] = degrees.get(label, 0) + 1
        if not degrees:
            stop_reason = 'coordinate_precision'
            break
        # Directed arithmetic is unnecessary for this heuristic, but keeps it
        # deterministic and independent of the caller's decimal context.
        def score(label):
            lo, hi = builders[label].bracket()
            return UP.multiply(UP.subtract(hi, lo), D(degrees[label]))
        selected = max(sorted(degrees), key=score)
        # Scheduling can itself consume the remaining budget.
        if perf_counter() - started >= seconds:
            stop_reason = 'time_budget'
            break
        operation = builders[selected].advance()
        largest_atomic = max(largest_atomic, builders[selected]._largest_atomic)
        if operation is None:
            exhausted.add(selected)
        else:
            schedule.append(dict(index=len(schedule) + 1, label=selected,
                                 operation=operation, elapsed_seconds=perf_counter() - started))
            update_provisional()
    else:
        stop_reason = 'all_comparisons_resolved'
    generation_seconds = perf_counter() - started

    output_started = perf_counter()
    certificates = {label: builder.certificate(seconds_budget=seconds, stop_reason=stop_reason)
                    for label, builder in builders.items()}
    output_seconds = perf_counter() - output_started
    replay_started = perf_counter()
    for label, record in certificates.items():
        if selective_peak.replay(record) is not True:
            raise ValueError('independent peak replay failed')
        if (record['identity_sha256'] != builders[label]._identity_sha256
                or record.get('problem_id') != identities[label]):
            raise ValueError('final certificate identity differs from owned kernel')
    replay_seconds = perf_counter() - replay_started
    # Only these independently replayed records determine returned decisions.
    output_started = perf_counter()
    comparisons = []
    for index, (a, b) in enumerate(pairs):
        resolved_op, resolved_seconds = provisional.get(index, (None, None))
        comparisons.append(dict(
            a=a, b=b, problem_id=identities[a],
            status=_status(tuple(map(D, certificates[a]['peak_bracket_MW_m2'])),
                           tuple(map(D, certificates[b]['peak_bracket_MW_m2']))),
            provisional_resolved_after_operation=resolved_op,
            provisional_resolved_seconds=resolved_seconds))
    result = dict(
        schema='adaptive-peak-ranking.v1', certificates=certificates,
        comparisons=comparisons, schedule=schedule, completed_operations=len(schedule),
        replay_count=len(certificates), stop_reason=stop_reason,
        seconds_budget=seconds, max_operations=max_operations,
        initialization_seconds=initialization_seconds, generation_seconds=generation_seconds,
        replay_seconds=replay_seconds, budget_overrun_seconds=max(0., generation_seconds - seconds),
        largest_atomic_operation_seconds=largest_atomic,
        scope='declared finite periodic-kernel peak ranking only; feasibility is separate',
        geometry_global_optimum_proven=False, continuous_transport_proven=False)
    result['output_seconds'] = output_seconds + perf_counter() - output_started
    result['end_to_end_seconds'] = perf_counter() - started
    return result
