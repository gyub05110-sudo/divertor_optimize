"""Cold 32-member OMP differential evolution at the supplied input's fine level."""
from collections import Counter
from math import isfinite
from time import perf_counter

import numpy as np

from .campaign_geometry import Candidate
from .campaign_search import bounds, retain, scales
from .campaign_state import accept, new_population, propose
from .user_search import _measurement, _retain
from run_tangent_curvature_study import project_floor


POPULATION = 32
PROJECTED_SEEDS = 12


def _initialize(fine, seed, archive, inherited, deadline):
    """Cold OMP branch of campaign_search.seed_population, with no historical file access."""
    rng = np.random.default_rng(seed)
    centers = [np.asarray(row["candidate"]["vector"]) for row in archive[:12]]
    centers.extend(inherited)
    for j in range(PROJECTED_SEEDS):
        if perf_counter() >= deadline:
            return None
        target = np.zeros(4) if j == 0 else rng.uniform(-.08, .08, 4)
        floor = project_floor(fine, target)
        if floor is not None:
            centers.append(np.r_[np.zeros(12), floor])
    if not centers:
        centers = [np.zeros(16)]
    initial = centers[:POPULATION]
    scale, box = scales(fine), bounds("OMP")
    while len(initial) < POPULATION:
        i = len(initial)
        amplitude = (.05, .2, .7, 2.)[i % 4]
        x = centers[i % len(centers)] + rng.uniform(-1., 1., 16) * scale * amplitude
        if i % 13 == 0:
            x = rng.uniform(box[:, 0], box[:, 1])
        initial.append(np.clip(x, box[:, 0], box[:, 1]))
    return new_population(initial, seed)


def _objective(fine, raw, peak, eligible):
    """Native SearchCase geometry residuals with the supplied-input eligibility predicate."""
    if eligible:
        return peak / fine.baseline_peak
    shape = raw.get("shape_constraints", {})
    if shape and not shape["passed"]:
        violation = max(0., fine.rmin * max(shape["curvature_max_1_m"].values()) - 1.)
        violation += max(0., (fine.gap_min - shape["gap_m"]) / .03)
        return 10. + violation + max(0., max(raw["precise_length_ratios"].values()) - 1.2)
    if not raw.get("geometry_valid", False):
        excess = max(0., max(raw["precise_length_ratios"].values()) - 1.2)
        collision = "material_length_exceeds_limit" not in raw["errors"]
        return 10. + excess + collision
    return 30. + min(raw["penalty"], 1000.) / 1000.


def search_direct(fine, *, radius_mm, seed, started, deadline, baseline, observe, evaluate=None):
    """Direct fine comparator; the clock, admission, and rich selected result match search_peak."""
    if not isfinite(started) or not isfinite(deadline) or deadline < started:
        raise ValueError("expected finite process start and deadline in perf_counter seconds")
    if evaluate is None:
        from .user_physics import evaluate_user_candidate as evaluate
    counts = dict.fromkeys(("proposals", "projection_rejected", "duplicates", "coarse", "fine",
                            "eligible_coarse", "eligible_fine", "failed", "late"), 0)
    selected = winner = best = first_feasible = None
    if baseline is not None and baseline.get("within_budget", True) and _measurement(baseline)[1]:
        selected = baseline
    errors = Counter()
    archive, finalists, inherited = [], [], []
    objectives = {}
    restart = last_improvement = 0
    best_infeasible = optimizer = None
    overrun = 0.
    termination = "deadline"
    scale, box = scales(fine), bounds("OMP")

    def emit(stage, begin, end, status, **fields):
        nonlocal overrun
        on_time = end <= deadline
        if not on_time:
            overrun = max(overrun, end-deadline)
            counts["late"] += 1
        if status == "failed":
            counts["failed"] += 1
        errors.update(fields.get("errors", []))
        observe(dict(stage=stage, start_s=begin-started, end_s=end-started, duration_s=end-begin,
                     within_budget=on_time, status=status if on_time or status == "failed" else "late",
                     **fields))

    while (begin := perf_counter()) < deadline:
        if optimizer is None:
            initial_seed = seed + radius_mm + 10000 + 100003 * restart
            fields = dict(restart=restart, population=POPULATION, projected_seeds=PROJECTED_SEEDS,
                          seed=initial_seed, archive_count=len(archive), inherited_count=len(inherited))
            try:
                optimizer = _initialize(fine, initial_seed, archive, inherited, deadline)
                if optimizer is not None:
                    vectors = np.asarray(optimizer["population"])
                    rank = int(np.linalg.matrix_rank((vectors-vectors.mean(axis=0))/.2, tol=1e-12))
                    fields.update(rank=rank, vectors=vectors.tolist(), axis_span=np.ptp(vectors, axis=0).tolist())
                    if restart == 0 and (rank != 16 or min(fields["axis_span"]) <= 1e-12):
                        raise ValueError("cold 32-member OMP population lacks full affine rank or axis span")
            except Exception as exc:
                emit("initialize", begin, perf_counter(), "failed", **fields, errors=[f"{type(exc).__name__}: {exc}"])
                termination = "failed"
                break
            emit("initialize", begin, perf_counter(), "completed", **fields, errors=[])
            last_improvement = counts["proposals"]
            continue

        index, target = counts["proposals"], optimizer["index"]
        generation = optimizer["generation"]
        polish = best is not None and index % 19 == 0
        if polish:
            x = np.array(best.vector)
            j = (index // 19) % (2 * len(x))
            factor = (.25, .08, .02, .005)[(index // (38 * len(x))) % 4]
            x[j//2] += (-1 if j % 2 else 1) * factor * scale[j//2]
            x = np.clip(x, box[:, 0], box[:, 1])
            kind = "coordinate_polish"
        else:
            x = propose(optimizer, box)
            kind = "population_initialization" if generation == 0 else "de_trial"
        candidate = Candidate("OMP", radius_mm, tuple(map(float, x)))
        identity = candidate.identity
        fields = dict(id=identity, candidate=candidate.record(), index=index, target_index=target,
                      generation=generation, restart=restart, proposal_kind=kind)
        duplicate = identity in objectives
        counts["proposals"] += 1
        counts["duplicates"] += int(duplicate)
        emit("proposal", begin, perf_counter(), "rejected" if duplicate else "completed",
             **fields, parent_id=best.identity if polish else None,
             errors=["duplicate_candidate"] if duplicate else [])
        begin = perf_counter()
        if begin >= deadline:
            break
        if duplicate:
            score = objectives[identity]
        else:
            counts["fine"] += 1
            fields.update(incumbent_updated=False, coarse_peak_q=None, cached=False)
            try:
                raw = evaluate(fine, candidate)
                peak, eligible, reasons = _measurement(raw)
                score = _objective(fine, raw, peak, eligible)
            except Exception as exc:
                emit("fine", begin, perf_counter(), "failed", **fields, eligible=False, peak_q=None,
                     physics_evaluated=False, errors=[f"{type(exc).__name__}: {exc}"])
                if not isinstance(exc, (ValueError, ArithmeticError)):
                    termination = "failed"
                    break
                # Keep a failed target from blocking every following DE member.
                score = 1e6
            else:
                end = perf_counter()
                fields.update(peak_q=peak, eligible=eligible, errors=reasons,
                              physics_evaluated="deposition" in raw, objective=score)
                if "geometry_valid" in raw:
                    fields["geometry_valid"] = bool(raw["geometry_valid"])
                if "shape_constraints" in raw:
                    fields["shape_passed"] = bool(raw["shape_constraints"]["passed"])
                if "geometry_screen" in raw:
                    fields["whole_wall_passed"] = bool(raw["geometry_screen"]["passed"])
                if "shelf_exposure" in raw:
                    fields["shelf_exposure"] = raw["shelf_exposure"]
                if end <= deadline:
                    if eligible:
                        counts["eligible_fine"] += 1
                        if first_feasible is None:
                            first_feasible = end-started
                        if best is None or score < objectives[best.identity]:
                            best = candidate
                            last_improvement = counts["proposals"]
                        if selected is None or peak < selected["peak_q"]:
                            winner, selected = candidate, raw
                            fields["incumbent_updated"] = True
                        _retain(finalists, (peak, index, candidate, end-started))
                        archive = retain(archive, dict(identity=identity, candidate=candidate.record(),
                                                       objective=score), scale)
                    elif best_infeasible is None or score < best_infeasible[0]:
                        best_infeasible = (score, x.copy())
                emit("fine", begin, end, "completed", **fields)
                if end > deadline:
                    break
            objectives[identity] = score
        if not polish:
            accept(optimizer, x, score)
        stalled = counts["proposals"] - last_improvement
        if optimizer["generation"] >= 100 or (optimizer["generation"] >= 15 and stalled > 25*POPULATION):
            restart += 1
            if not archive and best_infeasible is not None:
                inherited.append(best_infeasible[1])
            optimizer = None

    rows = [dict(id=c.identity, candidate=c.record(), index=i, peak_q=peak, coarse_peak_q=None, end_s=end)
            for peak, i, c, end in finalists]
    status = "qualified" if winner is not None else ("fallback" if selected is not None else "no_design")
    return dict(status=status, termination=termination, winner=winner, selected=selected,
                counts=counts, errors=dict(errors), first_feasible_elapsed_s=first_feasible,
                inflight_overrun_s=overrun, elapsed_s=perf_counter()-started, fine_finalists=rows,
                population=POPULATION, projected_seeds=PROJECTED_SEEDS)
