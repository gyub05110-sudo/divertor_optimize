"""Live, cold-start anytime search; historical and production modules stay intact.

The supervisor supplies the process launch clock and writes the returned record.
Only completed, eligible fine events can replace the original-wall fallback.
SearchCase owns adaptive proposals, seeding, objective feedback and restarts.
"""
from collections import Counter
from dataclasses import asdict
import json
import os
from pathlib import Path
from time import perf_counter_ns, time

import numpy as np

from d3d_optimization.broadening_study import save
from d3d_optimization.campaign_geometry import Candidate, evaluator
from d3d_optimization.campaign_search import SearchCase
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_base_budget import measurement, prep, score
from verification.nf_base_live import draw
from verification.validate_boundary_transport import LEVELS

POLICIES = ("coarse_adaptive", "fine_adaptive", "random_fine")
CHECKPOINTS = (200, 400, 600)
SOURCE = "OMP"
RADIUS = 20
BATCH = 4
POPULATION = 64
PROJECTED_SEEDS = 12
POPULATION_TOLERANCE = 1e-12
NS_PER_S = 1_000_000_000


def _measurement(value):
    """Every score, including a geometric rejection, has the same audit keys."""
    return dict(eligible=bool(value.get("eligible", value.get("feasible", False))),
                errors=list(value.get("errors", [])), ledger=value.get("ledger", {}),
                peak=value.get("peak", value.get("peak_MW_m2")), A2=value.get("A2"))


def snapshot(record, t):
    """Pure reconstruction at t, without trusting a saved optimizer/incumbent.

    Counts include completed on-time events only. ``fine`` and ``coarse`` count
    fresh scores; cached requests have their own count. Late/failed totals are
    restricted to events that have ended by t, even if t exceeds the admission
    budget. A late value can therefore be audited but never win.
    """
    t = float(t)
    cutoff = min(t, float(record.get("budget_s", 600.)))
    counts = dict(fine=0, coarse=0, eligible_fine=0, eligible_coarse=0,
                  cached=0, proposals=0, rejected=0, late=0, failed=0)
    baseline = record.get("baseline")
    ready = record.get("baseline_ready_s")
    if baseline is None or ready is None:
        baseline_events = [e for e in record.get("events", [])
                           if e["stage"] == "baseline" and e.get("status") == "completed"
                           and e.get("measurement") is not None]
        if baseline_events:
            baseline_event = min(baseline_events, key=lambda e: e["end_s"])
            baseline, ready = baseline_event["measurement"], baseline_event["end_s"]
    available = bool(baseline and baseline.get("eligible") and ready is not None and ready <= cutoff)
    incumbent = (dict(kind="original_wall", id=None, candidate=None,
                      peak=baseline["peak"], A2=baseline["A2"]) if available else None)
    first = None
    rows = sorted(enumerate(record.get("events", [])),
                  key=lambda item: (item[1].get("end_s", float("inf")), item[0]))
    for _, event in rows:
        end = event.get("end_s")
        if end is None or end > t:
            continue
        status = event.get("status")
        if status in ("late", "failed"):
            counts[status] += 1
        if status != "completed" or end > cutoff:
            continue
        stage = event["stage"]
        if stage in ("proposal", "draw"):
            counts["proposals"] += 1
        if event.get("rejection"):
            counts["rejected"] += 1
        if stage not in ("fine", "coarse"):
            continue
        if event.get("cached"):
            counts["cached"] += 1
            continue
        counts[stage] += 1
        measured = event.get("measurement", {})
        if not measured.get("eligible"):
            continue
        counts["eligible_" + stage] += 1
        if stage != "fine" or not available:
            continue
        candidate = Candidate.from_record(event["candidate"])
        peak = measured["peak"]
        if first is None and peak <= .95 * baseline["peak"]:
            first = end
        if peak < incumbent["peak"]:
            incumbent = dict(kind="new_design", id=candidate.identity,
                             candidate=candidate.record(), peak=peak, A2=measured["A2"])
    return dict(t_s=t, available=available, incumbent=incumbent,
                first_5pct_s=first, counts=counts)


class Deadline(Exception):
    """An in-flight operation is retained, but no further action may start."""


class _DeadlineConfig(dict):
    """Let native should_stop check this arm's monotonic clock between seeds.

    Its existing epoch comparison receives an always-future or already-past
    deadline. Other config values and every seed/proposal rule are unchanged;
    wall-clock adjustments cannot extend the process-start admission budget.
    """
    def __init__(self, arm, **values):
        super().__init__(values)
        self.arm = arm

    def __getitem__(self, key):
        if key in ("search_deadline_epoch", "deadline_epoch"):
            return 0. if self.arm.clock() >= self.arm.deadline else float("inf")
        return super().__getitem__(key)


class _MeasuredEvaluator:
    """Keep SearchCase's real score/parity logic and add the full measurement."""
    def __init__(self, evaluator, guard):
        self.wrapped = evaluator
        self.guard = guard
        self.last_measurement = None

    def __getattr__(self, name):
        return getattr(self.wrapped, name)

    def constraints(self, coefficients):
        self.guard()
        return self.wrapped.constraints(coefficients)

    def evaluate(self, coefficients):
        self.guard()
        return self.wrapped.evaluate(coefficients)

    def describe(self, result):
        self.last_measurement = _measurement(measurement(result, self.reference))
        return dict(self.wrapped.describe(result), measurement=self.last_measurement)


def _population_evidence(population, startup, archive_count, inherited_count):
    """Record affine diversity without changing native seed or restart rules."""
    vectors = np.asarray(population, dtype=float)
    if vectors.ndim != 2 or not len(vectors) or not np.isfinite(vectors).all():
        raise ValueError("initial population must be a finite nonempty matrix")
    axis_min, axis_max = vectors.min(axis=0), vectors.max(axis=0)
    span = axis_max-axis_min
    center = vectors.mean(axis=0)
    # OMP's native coefficient bounds are fixed at [-0.1, 0.1] m.
    # Use their widths, never this sample's observed span, as the rank scale.
    scale = np.full(vectors.shape[1], .2)
    singular_values = np.linalg.svd((vectors-center)/scale, compute_uv=False)
    rank = int(np.count_nonzero(singular_values > POPULATION_TOLERANCE))
    full_space = vectors.shape[1] == 16 and bool(np.all(span > POPULATION_TOLERANCE)) and rank == 16
    return dict(vectors=vectors.tolist(), count=len(vectors), dimension=vectors.shape[1],
                axis_min=axis_min.tolist(), axis_max=axis_max.tolist(), axis_span=span.tolist(),
                center=center.tolist(), scale=scale.tolist(), rank=rank,
                singular_values=singular_values.tolist(),
                rank_method="centered_and_fixed_bound_width_scaled", rank_tolerance=POPULATION_TOLERANCE,
                axis_span_tolerance=POPULATION_TOLERANCE, full_space_passed=full_space,
                full_space_required=startup, startup=startup,
                archive_count=archive_count, inherited_count=inherited_count,
                archive_origin="current_run" if archive_count or inherited_count else "none")


class _TimedSearchCase(SearchCase):
    """Instrument the existing optimizer, including its run_slice restarts."""
    def __init__(self, *args, arm, level, **kwargs):
        self.arm, self.level = arm, level
        self.proposal_event = None
        super().__init__(*args, **kwargs)
        self.e = _MeasuredEvaluator(self.e, arm.guard)
        self.cpu = _MeasuredEvaluator(self.cpu, arm.guard)

    def initialize(self, inherited=()):
        event = self.arm.begin("initialize", level=self.level, restart=self.state["restart"])
        previous_optimizer = self.state.get("optimizer")
        startup = self.state["restart"] == 0
        try:
            inherited = tuple(inherited)
            archive_count, inherited_count = len(self.state["archive"]), len(inherited)
            event.update(archive_count=archive_count, inherited_count=inherited_count,
                         full_space_required=startup)
            if startup and (archive_count or inherited_count or previous_optimizer is not None):
                raise ValueError("cold initialization requires no prior optimizer, archive or inherited shapes")
            try:
                super().initialize(inherited)
            finally:
                optimizer = self.state.get("optimizer")
                if optimizer is not None and optimizer is not previous_optimizer:
                    event["initial_population"] = _population_evidence(
                        optimizer["population"], startup, archive_count, inherited_count)
            evidence = event["initial_population"]
            if startup and (not evidence["full_space_passed"] or evidence["count"] != POPULATION):
                raise ValueError("initial population requires 64 members, all 16 axis spans > 1e-12, "
                                 "and centered fixed-bound-scaled affine rank 16")
            # Later restarts can inherit current-run archive members. Keep the
            # native policy, record its diversity, and enforce the gate only at
            # cold startup. This does not assert later evaluated coverage.
        except BaseException as exc:
            self.arm.fail(event, exc)
            raise
        self.arm.finish(event)

    def checkpoint(self, force=True):
        self.arm.guard()
        return super().checkpoint(force=force)

    def score(self, x, stage):
        candidate = Candidate(self.source, self.mm, tuple(map(float, x)))
        row = self.arm.register(candidate)
        optimizer = self.state.get("optimizer")
        context = dict(search_stage=stage, generation=None, index=None)
        if stage == "differential_evolution":
            index = optimizer["index"]
            context.update(generation=optimizer["generation"], index=index,
                           proposal_kind="population_initialization" if optimizer["energies"][index] is None else "de_trial")
        elif stage == "coordinate_polish":
            context.update(generation=optimizer["generation"], index=optimizer["index"],
                           proposal_kind="coordinate_polish")
        else:
            context["proposal_kind"] = "external_validation" if stage == "coarse_promotion" else stage
        if self.proposal_event is not None:
            self.proposal_event.update(id=row["id"], candidate=row["candidate"], **context)
            self.arm.finish(self.proposal_event)
            self.proposal_event = None
        self.e.last_measurement = None
        event = self.arm.begin(self.level, id=row["id"], candidate=row["candidate"], **context)
        try:
            scored = super().score(x, stage)
            measured = scored.get("measurement") or self.e.last_measurement or _measurement(scored)
            event.update(measurement=_measurement(measured), cached=bool(scored.get("cached")),
                         physics_evaluated=bool(scored.get("physics_evaluated")),
                         objective=scored["objective"])
        except BaseException as exc:
            if self.e.last_measurement is not None:
                event["measurement"] = self.e.last_measurement
            else:
                event["measurement"] = _measurement(dict(errors=[f"{type(exc).__name__}: {exc}"]))
            self.arm.fail(event, exc)
            raise
        self.arm.finish(event)
        # finish() raises before SearchCase.step can accept a late objective.
        if self.level == "coarse" and event["measurement"]["eligible"] and not event["cached"]:
            if row["id"] not in self.arm.coarse_seen:
                self.arm.coarse_seen.add(row["id"])
                self.arm.group.append(event)
        return scored

    def step(self):
        self.proposal_event = self.arm.begin("proposal", level=self.level)
        try:
            super().step()
        except BaseException as exc:
            if self.proposal_event is not None and "status" not in self.proposal_event:
                self.arm.fail(self.proposal_event, exc)
            raise
        finally:
            self.proposal_event = None
        self.arm.guard()
        if self.level == "coarse" and len(self.arm.group) == BATCH:
            group = self.arm.group
            self.arm.group = []
            winner = min(group, key=lambda e: (e["measurement"]["peak"], e["start_s"]))
            self.arm.record["promotions"].append(dict(
                ids=[e["id"] for e in group], selected=winner["id"],
                at_s=self.arm.elapsed(self.arm.clock())))
            self.arm.fine.score(winner["candidate"]["vector"], "coarse_promotion")
        self.arm.progress()


class _Arm:
    def __init__(self, record, out, started_ns, clock_ns):
        self.record, self.out, self.started, self.clock = record, Path(out), started_ns, clock_ns
        self.deadline = started_ns + round(record["budget_s"] * NS_PER_S)
        self.by_id, self.coarse_seen, self.group = {}, set(), []
        self.last_save = self.started
        self.fine = None
        self.cases = {}
        self.stream_started = False

    def elapsed(self, ns):
        return (ns - self.started) / NS_PER_S

    def guard(self):
        if self.clock() >= self.deadline:
            raise Deadline

    def begin(self, stage, **labels):
        start = self.clock()
        if start >= self.deadline:
            raise Deadline
        event = dict(stage=stage, start_s=self.elapsed(start), **labels)
        self.record["events"].append(event)
        return event

    def finish(self, event):
        end = self.clock()
        event.update(end_s=self.elapsed(end), status="late" if end > self.deadline else "completed")
        self.append_event(event)
        if end > self.deadline:
            raise Deadline

    def fail(self, event, exc):
        event.update(end_s=self.elapsed(self.clock()),
                     status="late" if isinstance(exc, Deadline) else "failed")
        if not isinstance(exc, Deadline):
            event["error"] = f"{type(exc).__name__}: {exc}"
        self.append_event(event)

    def append_event(self, event):
        """Append once after completion; subsequent work pays this I/O cost."""
        payload = json.dumps(event, allow_nan=False) + "\n"
        self.out.mkdir(parents=True, exist_ok=True)
        mode = "a" if self.stream_started else "x"
        with (self.out/"events.jsonl").open(mode, encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        self.stream_started = True

    def timed(self, stage, operation, observe=None, **labels):
        event = self.begin(stage, **labels)
        try:
            result = operation()
            if observe is not None:
                observe(result, event)
        except BaseException as exc:
            self.fail(event, exc)
            raise
        self.finish(event)
        return result, event

    def register(self, candidate):
        # Candidate.identity is recomputed from its immutable record every time.
        candidate = Candidate.from_record(candidate.record())
        key = candidate.identity
        if key not in self.by_id:
            row = dict(i=len(self.by_id), id=key, candidate=candidate.record(),
                       first_seen_s=self.elapsed(self.clock()))
            self.by_id[key] = row
            self.record["candidates"].append(row)
        return self.by_id[key]

    def progress(self):
        self.guard()
        if self.clock() - self.last_save < 2 * NS_PER_S:
            return
        self.record["checkpoints"] = {str(t): snapshot(self.record, t) for t in CHECKPOINTS
                                      if t <= self.elapsed(self.clock())}
        self.timed("progress_io", lambda: save(self.out/"anytime_progress.json", self.record))
        self.last_save = self.clock()


def _cold_inputs(root, out):
    """Reject the historical sources read by SearchCase.warm_vectors."""
    forbidden = [root/"campaign_protected_candidates.json", root/"campaign_warm_candidates.json"]
    forbidden += [root/"results/omp-tangent-curvature-20261029"/f"R{mm:03d}"/"manifest.json"
                  for mm in (100, 50, 20)]
    present = [str(path) for path in forbidden if path.exists()]
    if present:
        raise ValueError("cold input root contains historical shape initialization: " + ", ".join(present))
    for level in ("fine", "coarse"):
        folder = out/("search_" + level)
        if folder.exists() and any(folder.iterdir()):
            raise ValueError("anytime arm requires an empty search directory: " + str(folder))
    if (out/"events.jsonl").exists():
        raise ValueError("anytime arm event stream already exists: " + str(out/"events.jsonl"))
    out.mkdir(parents=True, exist_ok=True)


def run_arm(policy, input_root, seed, out, started_ns, budget_s=600, clock_ns=perf_counter_ns):
    """Run one arm until the process-start admission deadline, with no size cap."""
    if policy not in POLICIES:
        raise ValueError(policy)
    if not np.isfinite(budget_s) or budget_s <= 0:
        raise ValueError("positive finite budget required")
    root, out = Path(input_root), Path(out)
    record = dict(schema="nf-anytime-search.v1", policy=policy, seed=int(seed),
                  budget_s=float(budget_s), startup_s=(clock_ns()-started_ns)/NS_PER_S,
                  events=[], candidates=[], rejected=Counter(), completed_draws=0,
                  promotions=[], baseline=None, baseline_ready_s=None, checkpoints={},
                  error=None, source=SOURCE, radius_mm=RADIUS,
                  adaptive_config=dict(population=POPULATION, projected_seeds=PROJECTED_SEEDS,
                                       initial_population_gate="startup_only_full_affine_dimension"),
                  levels=dict(coarse=LEVELS[0], fine=LEVELS[2]))
    arm = _Arm(record, out, started_ns, clock_ns)
    try:
        arm.timed("cold_inputs", lambda: _cold_inputs(root, out))
        reference = ReferencePoint()
        record["reference"] = asdict(reference)
        (original, _), _ = arm.timed("prepare_fine", lambda: prep(root, SOURCE, reference, 2))
        baseline, baseline_event = arm.timed(
            "baseline", lambda: _measurement(measurement(original.baseline, reference)),
            observe=lambda value, event: event.update(measurement=value))
        record.update(baseline=baseline, baseline_ready_s=baseline_event["end_s"])
        remaining = max(0., (arm.deadline-clock_ns())/NS_PER_S)
        config = _DeadlineConfig(arm, seed=int(seed), population=POPULATION, projected_seeds=PROJECTED_SEEDS, backend="gpu",
                                 fingerprint=f"nf-anytime-search.v1:{root.resolve()}:{policy}:{seed}",
                                 search_deadline_epoch=time()+remaining, deadline_epoch=time()+remaining)

        def build(original, level):
            case = _TimedSearchCase(original, RADIUS, out/("search_"+level), root,
                                   config, arm=arm, level=level)
            arm.cases[level] = case
            return case

        arm.fine, _ = arm.timed("build_fine", lambda: build(original, "fine"))
        if policy == "coarse_adaptive":
            (low, _), _ = arm.timed("prepare_coarse", lambda: prep(root, SOURCE, reference, 0))
            case, _ = arm.timed("build_coarse", lambda: build(low, "coarse"))
        else:
            case = arm.fine
        if policy != "random_fine":
            arm.guard()
            # Preserve the production run_slice restart thresholds and inherited
            # current-run archive; staging has no historical initialization.
            case.run_slice(max(0., (arm.deadline-clock_ns())/NS_PER_S))
            record["stop"] = "deadline" if clock_ns() >= arm.deadline else "search_stopped"
        else:
            rng = np.random.default_rng(seed)
            while True:
                def observe_draw(value, event):
                    candidate, reason = value
                    record["completed_draws"] += 1
                    if reason:
                        record["rejected"][reason] += 1
                        event["rejection"] = reason
                    else:
                        row = arm.register(candidate)
                        event.update(id=row["id"], candidate=row["candidate"])

                (candidate, reason), _ = arm.timed(
                    "draw", lambda: draw(case.e, rng, record["candidates"]), observe=observe_draw)
                if reason is None:
                    case.score(candidate.vector, "projected_random")
                arm.progress()
    except Deadline:
        record["stop"] = "deadline"
    except Exception as exc:
        record.update(stop="failed", error=f"{type(exc).__name__}: {exc}")
    if not record["events"]:
        # Import/startup alone may exhaust the budget. Preserve an empty raw
        # stream for the supervisor without starting any search operation.
        try:
            out.mkdir(parents=True, exist_ok=True)
            with (out/"events.jsonl").open("x", encoding="utf-8"):
                pass
        except Exception as exc:
            record.update(stop="failed", error=f"{type(exc).__name__}: {exc}")
    record["rejected"] = dict(record["rejected"])
    record["parity"] = {level: arm.cases[level].state.get("parity") if level in arm.cases else None
                        for level in ("fine", "coarse")}
    record["checkpoints"] = {str(t): snapshot(record, t) for t in CHECKPOINTS}
    record["final"] = snapshot(record, budget_s)
    record["elapsed_s"] = arm.elapsed(clock_ns())
    return record
