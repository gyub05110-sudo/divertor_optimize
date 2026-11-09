"""Serial, budgeted OMP search. Only completed fine scores select the returned wall."""
from collections import Counter, deque
from heapq import heappop, heappush
from math import isfinite
from time import perf_counter

import numpy as np

from .campaign_geometry import Candidate
from .campaign_search import scales
from .tangent_evaluator import ACTIVE
from .user_surface_policy import REGION_FAMILIES, SEARCH_FAMILIES
from run_tangent_curvature_study import project_floor


COARSE_BATCH = 4
ELITE_SIZE = 8


def _measurement(raw):
    value = raw.get("peak_q")
    peak = float(value) if value is not None and np.ndim(value) == 0 else None
    peak = peak if peak is not None and isfinite(peak) and peak > 0. else None
    eligible = bool(raw.get("feasible") and not raw.get("reference_only") and peak is not None)
    errors = list(raw.get("errors", []))
    if raw.get("feasible") and peak is None:
        errors.append("invalid_whole_wall_peak")
    return peak, eligible, errors


class _Proposals:
    """Projected cold restarts and signed mutations around this run's elite walls."""

    def __init__(self, geometry, radius_mm, seed, family):
        self.geometry = geometry
        self.radius_mm = radius_mm
        self.rng = np.random.default_rng(seed)
        self.frozen = np.repeat([family != "full" and name != "floor" and REGION_FAMILIES[name] != family
                                 for name in ACTIVE], 4)
        self.scale = scales(geometry)

    def draw(self, index, coarse_elite, fine_elite):
        rng, bound = self.rng, self.geometry.bound
        parent = None
        if index < 4 or index % 5 == 0 or not coarse_elite:
            # All families consume the same random draws and share their origin.
            # The first projected tangent floor is distinct from the original wall.
            target = np.zeros(4) if index == 0 else rng.uniform(-.85*bound, .85*bound, 4)
            x = np.zeros(16)
            perturbation = rng.uniform(-1., 1., 12) * self.scale[:12] * rng.uniform(.05, .5)
            if index % 2:
                x[:12] = perturbation
            kind = "projected_origin" if index == 0 else "projected_restart"
        else:
            pool = fine_elite if fine_elite and index % 3 else coarse_elite
            # Usually refine the leader. Revisit other low-peak basins as well.
            rank = rng.integers(min(4, len(pool))) if index % 4 == 0 else 0
            parent = pool[rank][2]
            x = np.array(parent.vector)
            step = rng.uniform(-1., 1., 16) * self.scale * rng.choice([.12, .4, 1.2])
            if index % 3 == 0:
                # Carry forward successful directions in the projected family.
                x += rng.uniform(.3, .8) * (x - np.array(coarse_elite[-1][2].vector))
            if index % 3 == 1:
                # Deliberately give all families some pure floor moves.
                step[:12] = 0.
            x += step
            target = x[12:]
            kind = "local_projected"
        x[self.frozen] = 0.
        floor = project_floor(self.geometry, np.clip(target, -bound, bound))
        if floor is None:
            return None, kind, parent
        x[12:] = floor
        candidate = Candidate("OMP", self.radius_mm, tuple(map(float, np.clip(x, -bound, bound))))
        return candidate, kind, parent


def _retain(elite, entry):
    elite.append(entry)
    # The proposal index retains the first occurrence of an exact peak tie.
    elite.sort(key=lambda row: row[:2])
    del elite[ELITE_SIZE:]


def search_peak(coarse_evaluator, fine_evaluator, *, radius_mm, seed, deadline, started,
                baseline, observe, family="full", evaluate=None):
    """Find the lowest eligible fine peak measured before an absolute perf_counter deadline.

    Pass the native fine reference result with within_budget=False if completed late.
    A reference-only wall supplies no eligible fallback. Returned rich arrays are
    retained only for the selected design. Counts for coarse/fine include attempts,
    while eligible_coarse/eligible_fine count only on-time eligible scores.
    """
    if family not in SEARCH_FAMILIES:
        raise ValueError(f"expected resolved screened search family, got {family}")
    if not isfinite(started) or not isfinite(deadline) or deadline < started:
        raise ValueError("expected finite process start and deadline in perf_counter seconds")
    counts = dict.fromkeys(("proposals", "projection_rejected", "duplicates", "coarse", "fine",
                            "eligible_coarse", "eligible_fine", "failed", "late"), 0)
    errors = Counter()
    selected = None
    if baseline is not None and baseline.get("within_budget", True) and _measurement(baseline)[1]:
        selected = baseline
    winner = None
    first_feasible = None
    overrun = 0.
    seen = set()
    pending = []
    coarse_elite, fine_elite = [], []
    coarse_durations, fine_durations = deque(maxlen=4), deque(maxlen=4)
    screens_since_fine = 0
    proposals = None
    termination = "deadline"
    # No physics import or proposal initialization after a setup-exhausted budget.
    if perf_counter() < deadline:
        if evaluate is None:
            from .user_physics import evaluate_user_candidate as evaluate
        proposals = _Proposals(coarse_evaluator, radius_mm, seed, family)

    def emit(stage, begin, end, status, **fields):
        nonlocal overrun
        on_time = end <= deadline
        if not on_time:
            counts["late"] += 1
            overrun = max(overrun, end-deadline)
        if status == "failed":
            counts["failed"] += 1
        errors.update(fields.get("errors", []))
        observe(dict(stage=stage, start_s=begin-started, end_s=end-started, duration_s=end-begin,
                     status=status if on_time or status == "failed" else "late",
                     within_budget=on_time, **fields))

    while (begin := perf_counter()) < deadline:
        # Confirm the first eligible shape immediately. Later keep a rolling
        # shortlist of all unconfirmed coarse scores, including previous runners-up.
        # Near the deadline spend the remaining time on the best pending fine scores.
        reserve = (max(fine_durations, default=0.) + 2*max(coarse_durations, default=0.))
        promote = pending and (not counts["fine"] or screens_since_fine >= COARSE_BATCH
                               or deadline-begin <= reserve)
        if promote:
            coarse_peak, index, candidate, kind = heappop(pending)
            level, evaluator = "fine", fine_evaluator
            screens_since_fine = 0
        else:
            index = counts["proposals"]
            counts["proposals"] += 1
            try:
                candidate, kind, parent = proposals.draw(index, coarse_elite, fine_elite)
            except Exception as exc:
                emit("proposal", begin, perf_counter(), "failed", index=index,
                     errors=[f"{type(exc).__name__}: {exc}"])
                if isinstance(exc, (ValueError, ArithmeticError)):
                    continue
                termination = "failed"
                break
            end = perf_counter()
            identity = candidate.identity if candidate is not None else None
            if candidate is None:
                counts["projection_rejected"] += 1
                reason = "floor_projection_failed"
            elif identity in seen:
                counts["duplicates"] += 1
                reason = "duplicate_candidate"
            else:
                reason = None
            emit("proposal", begin, end, "rejected" if reason else "completed", index=index,
                 id=identity, candidate=candidate.record() if candidate is not None else None,
                 proposal_kind=kind, parent_id=parent.identity if parent is not None else None,
                 errors=[reason] if reason else [])
            begin = perf_counter()
            if begin >= deadline or reason:
                continue
            seen.add(identity)
            coarse_peak = None
            level, evaluator = "coarse", coarse_evaluator
            screens_since_fine += 1

        fields = dict(id=candidate.identity, index=index, candidate=candidate.record(),
                      proposal_kind=kind, incumbent_updated=False)
        if level == "fine":
            fields["coarse_peak_q"] = coarse_peak
        begin = perf_counter()
        if begin >= deadline:
            break
        counts[level] += 1
        try:
            raw = evaluate(evaluator, candidate)
            peak, eligible, reasons = _measurement(raw)
        except Exception as exc:
            emit(level, begin, perf_counter(), "failed", **fields, eligible=False, peak_q=None,
                 errors=[f"{type(exc).__name__}: {exc}"])
            if isinstance(exc, (ValueError, ArithmeticError)):
                continue
            termination = "failed"
            break
        end = perf_counter()
        fields.update(peak_q=peak, eligible=eligible, errors=reasons)
        if "geometry_valid" in raw:
            fields["geometry_valid"] = bool(raw["geometry_valid"])
        if "shape_constraints" in raw:
            fields["shape_passed"] = bool(raw["shape_constraints"]["passed"])
        if "geometry_screen" in raw:
            fields["whole_wall_passed"] = bool(raw["geometry_screen"]["passed"])
        if "shelf_exposure" in raw:
            fields["shelf_exposure"] = raw["shelf_exposure"]
        (fine_durations if level == "fine" else coarse_durations).append(end-begin)
        if eligible and end <= deadline:
            counts["eligible_" + level] += 1
            if level == "coarse":
                entry = (peak, index, candidate, kind)
                heappush(pending, entry)
                _retain(coarse_elite, entry)
            else:
                _retain(fine_elite, (peak, index, candidate, coarse_peak, end-started))
                if first_feasible is None:
                    first_feasible = end-started
                if selected is None or peak < selected["peak_q"]:
                    winner, selected = candidate, raw
                    fields["incumbent_updated"] = True
        emit(level, begin, end, "completed", **fields)

    finalists = [dict(id=candidate.identity, candidate=candidate.record(), index=index,
                      peak_q=peak, coarse_peak_q=coarse, end_s=end)
                 for peak, index, candidate, coarse, end in fine_elite]
    status = "qualified" if winner is not None else ("fallback" if selected is not None else "no_design")
    return dict(status=status, termination=termination, winner=winner, selected=selected,
                counts=counts, errors=dict(errors),
                first_feasible_elapsed_s=first_feasible, inflight_overrun_s=overrun,
                elapsed_s=perf_counter()-started, fine_finalists=finalists)
