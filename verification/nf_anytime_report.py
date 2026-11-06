"""Read-only evidence aggregation and Korean report for the full-space study.

Run only after the serial live timing campaign has ended. This module never
calls a physical evaluator, edits a frozen input, or substitutes a missing
outcome. Matplotlib/NumPy are loaded only for post-trial analysis and rendering.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
import math
import os
from pathlib import Path
from statistics import median


CODE_ROOT = Path(__file__).resolve().parents[1]
SEEDS = (20611, 20627, 20639)
POLICIES = ("coarse_adaptive", "fine_adaptive", "random_fine")
CHECKPOINTS = (200, 400, 600)
SETTINGS = ("cubic640", "cubic1280", "quintic640", "quintic1280")
BASELINE_ID = "original_wall"
NATIVE_RTOL = 1e-10
RANK_TOLERANCE = 1e-12
GEOMETRY_ERRORS = frozenset(("continuous_shape_constraint_failed", "complete_wall_invalid",
                             "material_length_exceeds_limit", "geometry_constraint_failed"))
LABELS = dict(coarse_adaptive="저비용 적응 탐색", fine_adaptive="정밀 적응 탐색",
              random_fine="투영 무작위 정밀 탐색")
COLORS = dict(coarse_adaptive="#27649B", fine_adaptive="#C46B22", random_fine="#8D7417")
STYLES = dict(coarse_adaptive="-", fine_adaptive="--", random_fine=":")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _identity(candidate):
    return hashlib.sha256(json.dumps(candidate, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _order(first, second):
    """Backend-replay tolerance diagnostic, never the primary win rule."""
    if abs(first-second) <= NATIVE_RTOL*max(abs(first), abs(second)):
        return 0
    return -1 if first < second else 1


def _exact_order(first, second):
    return -1 if first < second else 1 if first > second else 0


def _completed(event, cutoff):
    return event.get("status") == "completed" and event.get("end_s", math.inf) <= cutoff


def _baseline(arm, cutoff):
    rows = [e for e in arm["events"] if e.get("stage") == "baseline"
            and _completed(e, cutoff) and e.get("measurement", {}).get("eligible")]
    require(len(rows) <= 1, "multiple completed original-wall baselines")
    return rows[0] if rows else None


def reconstruct(arm, t):
    """Independent admission reconstruction; never read a saved incumbent."""
    cutoff = min(float(t), float(arm["budget_s"]))
    result = dict(available=False, kind=None, id=None, candidate=None, peak=None, A2=None,
                  baseline=None, baseline_ready_s=None, first_5pct_s=None)
    if not arm.get("native_pool_eligible", True):
        return result
    base = _baseline(arm, cutoff)
    if base is None:
        return result
    state = base["measurement"]
    result.update(available=True, kind="original_wall", peak=state["peak"], A2=state["A2"],
                  baseline=state, baseline_ready_s=base["end_s"])
    for event in sorted(arm["events"], key=lambda e: e.get("end_s", math.inf)):
        if (event.get("stage") != "fine" or not _completed(event, cutoff)
                or event.get("cached") or not event.get("measurement", {}).get("eligible")
                or event["end_s"] < base["end_s"]):
            continue
        measured = event["measurement"]
        identity = _identity(event["candidate"])
        require(identity == event["id"], "candidate identity differs from event coefficients")
        if result["first_5pct_s"] is None and measured["peak"] <= .95*state["peak"]:
            result["first_5pct_s"] = event["end_s"]
        if measured["peak"] < result["peak"]:
            result.update(kind="new_design", id=identity, candidate=event["candidate"],
                          peak=measured["peak"], A2=measured["A2"])
    return result


def anytime_curve(arm):
    """A step series begins at actual baseline readiness, never at time zero."""
    final = reconstruct(arm, arm["budget_s"])
    if not final["available"]:
        return []
    base = final["baseline"]
    points = [dict(t_s=final["baseline_ready_s"], peak=base["peak"], A2=base["A2"],
                   peak_ratio=1., id=None, kind="original_wall")]
    for event in sorted(arm["events"], key=lambda e: e.get("end_s", math.inf)):
        if (event.get("stage") != "fine" or not _completed(event, arm["budget_s"])
                or event.get("cached") or not event.get("measurement", {}).get("eligible")
                or event["end_s"] < points[0]["t_s"]):
            continue
        state = event["measurement"]
        if state["peak"] < points[-1]["peak"]:
            points.append(dict(t_s=event["end_s"], peak=state["peak"], A2=state["A2"],
                               peak_ratio=state["peak"]/base["peak"], id=event["id"], kind="new_design"))
    horizon = arm["budget_s"]
    if arm.get("stop") != "deadline" and arm.get("elapsed_s") is not None:
        horizon = min(horizon, arm["elapsed_s"])
    if points[-1]["t_s"] < horizon:
        points.append(dict(points[-1], t_s=horizon))
    return points


def first_improvement(arm):
    result = reconstruct(arm, arm["budget_s"])
    observed = result["first_5pct_s"] is not None
    certified = arm.get("native_pool_eligible", True)
    elapsed = arm.get("elapsed_s", arm.get("process", {}).get("elapsed_s"))
    censor = min(arm["budget_s"], elapsed) if elapsed is not None else None
    status, stop = arm.get("audit_status"), arm.get("stop")
    if observed:
        reason = "observed"
    elif status == "process_failed":
        reason = "process_failed_uncertified"
    elif status == "arm_failed" or stop == "failed":
        reason = "arm_failed"
    elif stop == "deadline" and censor is not None and censor >= arm["budget_s"]:
        reason = "deadline_without_5pct_improvement"
    elif stop == "search_stopped":
        reason = "search_stopped_before_budget" if censor is not None and censor < arm["budget_s"] else "search_stopped"
    else:
        reason = (stop or status or "unknown_stop")+"_without_5pct_improvement"
    return dict(event_observed=observed, time_s=result["first_5pct_s"],
                censor_time_s=None if observed else censor, reason=reason,
                certified=certified, baseline_available=result["available"])


def late_result_diagnostic(arm):
    """Describe discarded outcomes without altering admission or first-hit times."""
    incumbent = reconstruct(arm, arm["budget_s"])
    late = [event for event in arm["events"] if event.get("status") == "late"]
    rows = []
    for event in late:
        if event.get("stage") != "fine":
            continue
        measured = event.get("measurement")
        identity, candidate = event.get("id"), event.get("candidate")
        if candidate is not None:
            require(identity == _identity(candidate), "late candidate identity differs from coefficients")
        genuine = not event.get("cached") and event.get("physics_evaluated") is True
        eligible = bool(measured and measured.get("eligible"))
        observables = bool(measured and all(isinstance(measured.get(k), (int, float)) and
                          math.isfinite(measured[k]) and measured[k] > 0 for k in ("peak", "A2")))
        valid = genuine and eligible and observables
        reason = ("cached_result" if event.get("cached") else "measurement_missing" if not measured else
                  "deposition_result_not_confirmed" if not genuine else "ineligible" if not eligible else
                  "observables_missing" if not observables else "incumbent_unavailable" if not incumbent["available"] else "compared")
        comparable = valid and incumbent["available"]
        diagnostic = _order(measured["peak"], incumbent["peak"]) if comparable else None
        rows.append(dict(id=identity, candidate=candidate, start_s=event.get("start_s"), end_s=event.get("end_s"),
                         measurement=measured, cached=bool(event.get("cached")), physics_evaluated=event.get("physics_evaluated"),
                         genuine_eligible_fine_observation=valid, comparison_available=comparable,
                         comparison_reason=reason,
                         strict_improvement=measured["peak"] < incumbent["peak"] if comparable else None,
                         relative_improvement=1-measured["peak"]/incumbent["peak"] if comparable else None,
                         replay_tolerance_ordering={-1: "late_lower", 0: "within_replay_tolerance", 1: "incumbent_lower"}.get(diagnostic)))
    return dict(cutoff_s=arm["budget_s"], incumbent_at_cutoff=incumbent,
                late_coarse_count=sum(e.get("stage") == "coarse" for e in late), late_fine_count=len(rows),
                genuine_eligible_late_fine_count=sum(r["genuine_eligible_fine_observation"] for r in rows),
                compared_fine_count=sum(r["comparison_available"] for r in rows),
                strict_improvement_count=sum(r["strict_improvement"] is True for r in rows),
                improvement_beyond_replay_tolerance_count=sum(r["replay_tolerance_ordering"] == "late_lower" for r in rows),
                replay_diagnostic_rtol=NATIVE_RTOL, fine_results=rows,
                changes_admitted_incumbent=False, changes_first_improvement_time=False,
                scope="discarded genuine eligible fine observations versus the admitted fine incumbent; coarse outcomes are counted only")


def vector_coverage(vectors):
    """Affine coverage of observed coefficient vectors, with fixed 0.2 m scale."""
    import numpy as np
    unique = list(dict.fromkeys(tuple(float(value) for value in row) for row in vectors))
    if not unique:
        return dict(request_count=len(vectors), unique_count=0, rank=0, axis_min=None, axis_max=None,
                    axis_span=None, varied_axes=0, scale=[.2]*16, singular_values=[],
                    rank_tolerance=RANK_TOLERANCE)
    x = np.asarray(unique, float)
    require(x.ndim == 2 and x.shape[1] == 16 and np.isfinite(x).all(), "invalid observed OMP vectors")
    span = np.ptp(x, axis=0)
    singular = np.linalg.svd((x-x.mean(axis=0))/.2, compute_uv=False)
    return dict(request_count=len(vectors), unique_count=len(unique),
                rank=int(np.count_nonzero(singular > RANK_TOLERANCE)),
                axis_min=x.min(axis=0).tolist(), axis_max=x.max(axis=0).tolist(), axis_span=span.tolist(),
                varied_axes=int(np.count_nonzero(span > RANK_TOLERANCE)), scale=[.2]*16,
                singular_values=singular.tolist(), rank_tolerance=RANK_TOLERANCE)


def _score_counts(events, level, cutoff):
    requested = [e for e in events if e.get("stage") == level and e.get("start_s", math.inf) <= cutoff]
    ended = [e for e in requested if e.get("end_s", math.inf) <= cutoff]
    completed = [e for e in ended if e.get("status") == "completed"]
    fresh = [e for e in completed if not e.get("cached")]
    geometry = [e for e in fresh if GEOMETRY_ERRORS.intersection(e.get("measurement", {}).get("errors", []))]
    deposition = [e for e in fresh if e.get("physics_evaluated") is True]
    unknown = [e for e in fresh if not e.get("physics_evaluated") and e not in geometry]
    return dict(requests_started=len(requested), requests_completed=len(completed),
                fresh_requests_completed=len(fresh), in_flight=len(requested)-len(ended),
                late_requests=sum(e.get("status") == "late" for e in ended),
                failed_requests=sum(e.get("status") == "failed" for e in ended),
                deposition_result_evaluations=len(deposition),
                cached_requests=len(completed)-len(fresh), geometry_rejections=len(geometry),
                other_non_deposition_results=len(unknown),
                jointly_eligible=sum(bool(e.get("measurement", {}).get("eligible")) for e in fresh),
                proposal_kinds=dict(Counter(e.get("proposal_kind", "unrecorded") for e in completed)))


def checkpoint_metrics(arm, t):
    cutoff = min(float(t), arm["budget_s"])
    active = "coarse" if arm["policy"] == "coarse_adaptive" else "fine"
    completed = [e for e in arm["events"] if _completed(e, cutoff)]
    proposed = [e["candidate"]["vector"] for e in completed if "candidate" in e and
                ((e.get("stage") == "proposal" and e.get("level") == active)
                 or (e.get("stage") == "draw" and active == "fine"))]
    scored = [e for e in completed if e.get("stage") == active and not e.get("cached")]
    deposition = [e["candidate"]["vector"] for e in scored if e.get("physics_evaluated") is True]
    eligible = [e["candidate"]["vector"] for e in scored if e.get("measurement", {}).get("eligible")]
    choice = reconstruct(arm, cutoff)
    choice.update(cutoff_s=t, case_id=arm["case_id"], seed=arm["seed"], policy=arm["policy"],
                  active_level=active, counts={level: _score_counts(arm["events"], level, cutoff)
                                               for level in ("coarse", "fine")},
                  draw_rejections=dict(Counter(e["rejection"] for e in completed if e.get("rejection"))),
                  coverage={name: vector_coverage(values) for name, values in
                            (("proposed", proposed), ("deposition", deposition), ("eligible", eligible))},
                  evidence_certified=arm.get("native_pool_eligible", True),
                  peak_reduction_pct=None, A2_change_pct=None, q_down=False, q_down_A2_up=False,
                  q_down_beyond_replay_tolerance=False, q_down_A2_up_beyond_replay_tolerance=False)
    if choice["available"]:
        base = choice["baseline"]
        down = choice["peak"] < base["peak"]
        diagnostic_down = _order(choice["peak"], base["peak"]) < 0
        choice.update(peak_reduction_pct=100*(1-choice["peak"]/base["peak"]),
                      A2_change_pct=100*(choice["A2"]/base["A2"]-1), q_down=down,
                      q_down_A2_up=down and choice["A2"] > base["A2"],
                      q_down_beyond_replay_tolerance=diagnostic_down,
                      q_down_A2_up_beyond_replay_tolerance=diagnostic_down and _order(choice["A2"], base["A2"]) > 0)
    return choice


def paired_comparisons(arms, checkpoints=CHECKPOINTS):
    groups, rows = {}, []
    for arm in arms:
        key = (arm["case_id"], arm["seed"])
        require(arm["policy"] not in groups.setdefault(key, {}), "duplicate paired arm")
        groups[key][arm["policy"]] = arm
    for (case, seed), group in sorted(groups.items()):
        for cutoff in checkpoints:
            for first, second in combinations(POLICIES, 2):
                a = reconstruct(group[first], cutoff) if first in group else dict(available=False)
                b = reconstruct(group[second], cutoff) if second in group else dict(available=False)
                row = dict(case_id=case, seed=seed, cutoff_s=cutoff, first_policy=first, second_policy=second,
                           both_available=a["available"] and b["available"], ordering="unavailable",
                           first_peak=a.get("peak"), second_peak=b.get("peak"), difference=None,
                           relative_difference=None, replay_tolerance_ordering=None,
                           replay_diagnostic_rtol=NATIVE_RTOL)
                if row["both_available"]:
                    order, diagnostic = _exact_order(a["peak"], b["peak"]), _order(a["peak"], b["peak"])
                    row.update(ordering={-1: "first_lower", 0: "exact_tie", 1: "second_lower"}[order],
                               replay_tolerance_ordering={-1: "first_lower", 0: "within_replay_tolerance", 1: "second_lower"}[diagnostic],
                               difference=a["peak"]-b["peak"], relative_difference=a["peak"]/b["peak"]-1)
                rows.append(row)
    return rows


def _range(values):
    return [min(values), max(values)] if values else None


def stage_costs(arm):
    stages = Counter()
    for event in arm["events"]:
        if event.get("end_s") is not None:
            stages[event["stage"]] += event["end_s"]-event["start_s"]
    startup, elapsed = arm.get("startup_s"), arm.get("elapsed_s")
    gaps = elapsed-startup-sum(stages.values()) if startup is not None and elapsed is not None else None
    if gaps is not None:
        require(gaps >= -1e-6, "recorded stage durations exceed the worker clock")
        gaps = max(0., gaps)
    parity = [row for row in arm.get("parity", {}).values() if row is not None]
    cpu = [row["cpu_seconds"] for row in parity if row.get("cpu_seconds") is not None]
    gpu = [row["gpu_seconds"] for row in parity if row.get("gpu_seconds") is not None]
    scores = [e for e in arm["events"] if e.get("stage") in ("fine", "coarse")]
    process_elapsed = arm.get("process", {}).get("elapsed_s")
    base = _baseline(arm, arm["budget_s"])
    return dict(startup_s=startup, baseline_ready_s=base["end_s"] if base else None,
                elapsed_s=elapsed, process_elapsed_s=process_elapsed, admission_budget_s=arm["budget_s"],
                stage_seconds=dict(stages), clock_gaps_seconds=gaps,
                clock_gap_definition="charged event append/logging, optimizer bookkeeping and other uninstrumented process time; not removed from the budget",
                preparation_and_build_seconds=sum(stages[k] for k in ("cold_inputs", "prepare_fine", "prepare_coarse", "baseline", "build_fine", "build_coarse")),
                initialization_seconds=stages["initialize"],
                admitted_scoring_event_seconds=sum(e["end_s"]-e["start_s"] for e in scores if _completed(e, arm["budget_s"])),
                recorded_scoring_event_seconds=sum(e["end_s"]-e["start_s"] for e in scores if e.get("end_s") is not None),
                overrun_seconds=max(0., elapsed-arm["budget_s"]) if elapsed is not None else None,
                process_overrun_seconds=max(0., process_elapsed-arm["budget_s"]) if process_elapsed is not None else None,
                process_finalization_and_observation_seconds=max(0., process_elapsed-elapsed)
                if process_elapsed is not None and elapsed is not None else None,
                reported_cpu_parity_seconds=sum(cpu) if cpu else None,
                reported_gpu_parity_evaluation_seconds=sum(gpu) if gpu else None,
                parity_payload_count=len(parity), parity_is_already_inside_score_events=True,
                exact_total_solver_invocations_known=False)


def summarize_live(arms):
    metrics = [checkpoint_metrics(arm, t) for arm in arms for t in CHECKPOINTS]
    case_ids = sorted({arm["case_id"] for arm in arms})
    groups = []
    for case in case_ids:
        for policy in POLICIES:
            for cutoff in CHECKPOINTS:
                selected = [r for r in metrics if (r["case_id"], r["policy"], r["cutoff_s"]) == (case, policy, cutoff)]
                valid = [r for r in selected if r["available"]]
                gains = [r["peak_reduction_pct"] for r in valid]
                active = "coarse" if policy == "coarse_adaptive" else "fine"
                counts = Counter()
                kinds = Counter()
                draw_rejections = Counter()
                for row in selected:
                    current = row["counts"][active]
                    counts.update({k: v for k, v in current.items() if k != "proposal_kinds"})
                    kinds.update(current["proposal_kinds"])
                    draw_rejections.update(row["draw_rejections"])
                groups.append(dict(case_id=case, policy=policy, cutoff_s=cutoff,
                                   seed_count=len(selected), available_count=len(valid),
                                   uncertified_arm_count=sum(not r["evidence_certified"] for r in selected),
                                   peak_reduction_pct_median=median(gains) if gains else None,
                                   peak_reduction_pct_seed_range=_range(gains),
                                   q_down_count=sum(r["q_down"] for r in selected),
                                   q_down_A2_up_count=sum(r["q_down_A2_up"] for r in selected),
                                   q_down_beyond_replay_tolerance_count=sum(r["q_down_beyond_replay_tolerance"] for r in selected),
                                   active_counts_sum=dict(counts), proposal_kind_counts_sum=dict(kinds),
                                   draw_rejections_sum=dict(draw_rejections),
                                   coverage_rank_seed_ranges={pool: _range([r["coverage"][pool]["rank"] for r in selected])
                                                              for pool in ("proposed", "deposition", "eligible")}))
    arm_rows = []
    for arm in arms:
        recorded_horizon = max([arm.get("elapsed_s") or 0.] +
                               [e.get("end_s", 0.) for e in arm["events"]])
        arm_rows.append(dict(case_id=arm["case_id"], seed=arm["seed"], policy=arm["policy"],
                             audit_status=arm["audit_status"], stop=arm.get("stop"),
                             elapsed_s=arm.get("elapsed_s"), process=arm.get("process"),
                             first_5pct=first_improvement(arm), curve=anytime_curve(arm),
                             late_result_diagnostic=late_result_diagnostic(arm),
                             costs=stage_costs(arm),
                             parity=arm.get("parity", {}),
                             all_recorded_counts={level: _score_counts(arm["events"], level, recorded_horizon)
                                                  for level in ("coarse", "fine")},
                             late_scores=[e for e in arm["events"] if e.get("stage") in ("fine", "coarse")
                                          and e.get("status") == "late"]))
    fine_phase = [r for r in metrics if r["policy"] == "fine_adaptive" and r["cutoff_s"] == 600]
    phase_kinds = Counter()
    for row in fine_phase:
        phase_kinds.update(row["counts"]["fine"]["proposal_kinds"])
    phase = dict(arm_count=len(fine_phase), certified_arm_count=sum(r["evidence_certified"] for r in fine_phase),
                 zero_completed_de_trial_arm_count=sum(r["counts"]["fine"]["proposal_kinds"].get("de_trial", 0) == 0
                                                       for r in fine_phase),
                 completed_proposal_kinds=dict(phase_kinds),
                 only_initialization_and_polish=bool(phase_kinds) and set(phase_kinds) <= {"population_initialization", "coordinate_polish"},
                 arms=[dict(case_id=r["case_id"], seed=r["seed"], evidence_certified=r["evidence_certified"],
                            completed_proposal_kinds=r["counts"]["fine"]["proposal_kinds"]) for r in fine_phase])
    return dict(equilibrium_count=len(case_ids), seed_repeats_per_equilibrium=len(SEEDS),
                arm_count=len(arms), independent_discharge_count_claim=False,
                arms=arm_rows, checkpoint_metrics=metrics, groups=groups,
                fine_adaptive_phase_at_600=phase,
                paired_comparisons=paired_comparisons(arms),
                count_definitions=dict(
                    requests_started="score events begun by checkpoint; includes unfinished requests",
                    requests_completed="score events finished on time by checkpoint; includes cache hits",
                    deposition_result_evaluations="fresh completed events with physics_evaluated=true, meaning a deposition result was returned; not all solver entries or CPU/GPU invocations",
                    geometry_rejections="fresh completed scores containing an explicit registered geometry error",
                    eligible="joint geometry and physical eligibility; not geometry-only feasibility",
                    coverage="actual completed active-level proposals with recorded vectors, returned deposition results, or jointly eligible scores; rejected random draws without recorded vectors cannot enter rank; unique vectors; centered affine rank scaled by fixed 0.2 m bounds",
                    seed_ranges="descriptive within-equilibrium seed ranges, not confidence intervals"))


def exact_equal_k_expectation(values, k, baseline):
    require(isinstance(k, int) and 0 <= k <= len(values), "invalid equal-k allocation")
    if k == 0:
        return float(baseline)
    ordered = sorted(baseline if value is None else min(value, baseline) for value in values)
    denominator = math.comb(len(ordered), k)
    return math.fsum(value*math.comb(len(ordered)-i-1, k-1)/denominator
                     for i, value in enumerate(ordered[:len(ordered)-k+1]))


def panel_metrics(rows, baseline):
    universe = [r for r in rows if r["coarse"]["eligible"]]
    winners = [min(universe[i:i+4], key=lambda r: (r["coarse"]["peak"], r["i"]))
               for i in range(0, len(universe)-3, 4)]
    best = lambda pool: min([baseline["peak"]] + [r["fine"]["peak"] for r in pool if r["fine"]["eligible"]])
    selected, reference = best(winners), best(rows)
    expected = exact_equal_k_expectation([r["fine"]["peak"] if r["fine"]["eligible"] else None
                                         for r in universe], len(winners), baseline["peak"])
    return dict(panel_size=len(rows), coarse_eligible_count=len(universe),
                fine_requests_k=len(winners), uniform_population_size=len(universe),
                uniform_subset_count=math.comb(len(universe), len(winners)),
                promoted_ids=[r["id"] for r in winners], selected_peak=selected, reference_peak=reference,
                regret_relative=selected/reference-1, uniform_expected_peak=expected,
                uniform_expected_regret_relative=expected/reference-1,
                equal_wall_time_comparison=False,
                scope="locked finite panel, original-wall fallback, common coarse-eligible universe, uniform equal-k sampling without replacement")


def panel_overlap(overlap, panel_ids):
    if overlap is None:
        return dict(available=False, held_out_validation_claim=False,
                    reason="live-candidate overlap provenance is not recorded; overlap is unknown, not zero")
    for live_key, intersection_key in (("live_proposed", "proposed_intersection"),
                                       ("live_on_time_fine_observed", "on_time_fine_intersection")):
        live_ids, intersection = overlap[live_key+"_ids"], overlap[intersection_key+"_ids"]
        require(len(live_ids) == len(set(live_ids)) == overlap[live_key+"_count"], "invalid live overlap population")
        expected = set(panel_ids).intersection(live_ids)
        require(set(intersection) == expected and len(intersection) == len(expected) == overlap[intersection_key+"_count"],
                "panel/live overlap counts or identities differ")
    require(not overlap["held_out_validation_claim"] and not overlap["live_outcomes_used_for_generation"],
            "panel overlap provenance contradicts preregistered selection")
    return dict(available=True, **{key: value for key, value in overlap.items()
                                  if key not in ("live_proposed_ids", "live_on_time_fine_observed_ids")})


def _consistent_panel_observations(first, second):
    """Repeated native records must agree; this does not establish physical accuracy."""
    consistent = first["eligible"] == second["eligible"]
    for key in ("peak", "A2"):
        a, b = first.get(key), second.get(key)
        a_present = isinstance(a, (int, float)) and math.isfinite(a) and a > 0
        b_present = isinstance(b, (int, float)) and math.isfinite(b) and b > 0
        consistent = consistent and a_present == b_present
        if a_present and b_present:
            consistent = consistent and math.isclose(a, b, rel_tol=NATIVE_RTOL, abs_tol=0.)
    require(consistent, "duplicate panel observation differs in eligibility or native observables")


def panel_diagnostics(panels):
    """Post-hoc shape reuse from existing sealed panel outcomes only."""
    require(bool(panels), "no finite panels for shape diagnostics")
    require(len({(p["case_id"], p["seed"]) for p in panels}) == len(panels), "duplicate case/seed panel")
    cases, shapes = {}, {}
    for panel in sorted(panels, key=lambda p: (p["case_id"], p["seed"])):
        case = cases.setdefault(panel["case_id"], dict(panels=0, rows=0, candidates={}, baselines=[]))
        case["panels"] += 1
        provenance = dict(seed=panel["seed"], panel_sha256=panel["panel_sha256"],
                          result_sha256=panel["result_sha256"])
        base = panel["baseline"]
        require(base["eligible"] and all(isinstance(base.get(k), (int, float)) and
                math.isfinite(base[k]) and base[k] > 0 for k in ("peak", "A2")), "invalid panel original-wall baseline")
        if case["baselines"]:
            _consistent_panel_observations(case["baselines"][0]["measurement"], base)
        case["baselines"].append(dict(**provenance, row_i=None, measurement=base))
        require(len({r["id"] for r in panel["rows"]}) == len(panel["rows"]), "duplicate shape inside a fixed panel")
        for row in panel["rows"]:
            identity, candidate, measured = row["id"], row["candidate"], row["fine"]
            require(identity == _identity(candidate), "panel shape identity differs from coefficients")
            require(identity not in shapes or shapes[identity] == candidate, "same panel ID has different shapes")
            shapes[identity] = candidate
            if measured["eligible"]:
                require(all(isinstance(measured.get(k), (int, float)) and math.isfinite(measured[k]) and
                            measured[k] > 0 for k in ("peak", "A2")), "eligible panel observables are missing")
            observations = case["candidates"].setdefault(identity, [])
            if observations:
                _consistent_panel_observations(observations[0]["measurement"], measured)
            observations.append(dict(**provenance, row_i=row["i"], measurement=measured))
            case["rows"] += 1

    def representative(observations):
        # Exact minimum preserves the requested best-observed finite-union rule.
        return min(observations, key=lambda o: (o["measurement"]["peak"], o["seed"], o["row_i"])) \
            if observations[0]["measurement"]["eligible"] else observations[0]

    choices, transfer_rows, per_case = [], [], []
    for case_id, case in cases.items():
        selected_id, selected = BASELINE_ID, case["baselines"][0]
        observations = case["baselines"]
        for identity, pool in sorted(case["candidates"].items()):
            observation = representative(pool)
            if observation["measurement"]["eligible"] and observation["measurement"]["peak"] < selected["measurement"]["peak"]:
                selected_id, selected, observations = identity, observation, pool
        measured = selected["measurement"]
        choices.append(dict(case_id=case_id, id=selected_id, candidate=shapes.get(selected_id),
                            kind="original_wall" if selected_id == BASELINE_ID else "finite_panel_design",
                            peak=measured["peak"], A2=measured["A2"], measurement=measured,
                            selected_observation=selected, observations=observations))
        per_case.append(dict(case_id=case_id, panel_count=case["panels"], fine_measurement_rows=case["rows"],
                             unique_shape_count=len(case["candidates"]), shape_ids=sorted(case["candidates"]),
                             duplicate_measurement_rows=case["rows"]-len(case["candidates"])))
    for source in choices:
        for destination_id, case in cases.items():
            baseline_observation = case["baselines"][0]
            base = baseline_observation["measurement"]
            pool = case["baselines"] if source["id"] == BASELINE_ID else case["candidates"].get(source["id"], [])
            selected = (pool[0] if source["id"] == BASELINE_ID else representative(pool)) if pool else None
            measured = selected["measurement"] if selected else None
            valid = bool(measured and measured["eligible"])
            status = ("baseline" if source["id"] == BASELINE_ID else "eligible" if valid else
                      "ineligible" if measured is not None else "unavailable")
            transfer_rows.append(dict(source_case_id=source["case_id"], destination_case_id=destination_id,
                                      id=source["id"], status=status, measurement=measured,
                                      peak=measured.get("peak") if measured else None, A2=measured.get("A2") if measured else None,
                                      baseline=base, baseline_observation=baseline_observation,
                                      selected_observation=selected, observations=pool,
                                      peak_reduction_pct=100*(1-measured["peak"]/base["peak"]) if valid else None,
                                      A2_change_pct=100*(measured["A2"]/base["A2"]-1) if valid else None,
                                      exact_order_vs_destination_wall=_exact_order(measured["peak"], base["peak"]) if valid else None))
    occurrences = Counter(identity for case in cases.values() for identity in case["candidates"])
    common = sorted(set.intersection(*(set(case["candidates"]) for case in cases.values())))
    counts = dict(panel_count=len(panels), equilibrium_count=len(cases),
                  fine_measurement_rows=sum(case["rows"] for case in cases.values()),
                  candidate_case_count=sum(len(case["candidates"]) for case in cases.values()),
                  global_unique_shape_count=len(shapes), global_shape_ids=sorted(shapes),
                  common_to_all_shape_count=len(common), common_to_all_shape_ids=common,
                  shared_across_equilibria_shape_count=sum(count >= 2 for count in occurrences.values()),
                  original_wall_in_shape_counts=False, per_case=per_case)
    return dict(shape_counts=counts, posthoc_transfer=dict(posthoc=True, new_physics_evaluations=0,
        held_out_validation_claim=False, global_optimum_claim=False, live_checkpoint_selection=False,
        source_choices=choices, rows=transfer_rows, duplicate_consistency_rtol=NATIVE_RTOL,
        selection_rule="strict best eligible observed fine peak in each finite panel union plus original wall; wall retained on exact ties",
        duplicate_rule="check eligibility and available native peak/A2; preserve all observations, use exact minimum eligible peak",
        baseline_rule="first fixed seed in each equilibrium after checking repeated original-wall observations",
        scope="post-hoc reuse of the same shape ID in already scored destination panels; missing and ineligible outcomes are not imputed"))


def observable_native_replay_matched(replay):
    """A vacuous or eligibility-only match cannot verify an observable pair."""
    return bool(replay["all_available_observations_matched"] and replay["previously_scored"]
                and replay["compared_count"] > 0 and replay["observable_compared_count"] > 0)


def native_replay_coverage(candidates):
    replay = [row["native_reproducibility"] for row in candidates]
    eligibility_only = [row for row in replay if row["previously_scored"]
                        and row["compared_count"] > 0 and row["observable_compared_count"] == 0
                        and row["eligibility_only_count"] == row["compared_count"]]
    return dict(locked_candidate_count=len(replay),
                previously_scored_candidate_count=sum(row["previously_scored"] for row in replay),
                unobserved_candidate_count=sum(not row["previously_scored"] for row in replay),
                compared_observation_count=sum(row["compared_count"] for row in replay),
                compared_observable_count=sum(row["observable_compared_count"] for row in replay),
                eligibility_only_observation_count=sum(row["eligibility_only_count"] for row in replay),
                candidates_with_observable_comparisons=sum(row["observable_compared_count"] > 0 for row in replay),
                eligibility_only_candidate_count=len(eligibility_only),
                other_nonobservable_candidate_count=sum(row["previously_scored"] and
                    row["observable_compared_count"] == 0 and row not in eligibility_only for row in replay),
                observable_replay_matched_candidate_count=sum(observable_native_replay_matched(row) for row in replay))


def setting_pair(first, second, *, same_design=False):
    rows, signs, first_values, second_values, diagnostic_signs = {}, [], [], [], []
    for name in SETTINGS:
        a, b = first.get(name), second.get(name)
        valid = bool(a and b and a["eligible"] and b["eligible"])
        sign = _exact_order(a["peak"], b["peak"]) if valid else None
        rows[name] = dict(both_eligible=valid, first=a, second=b, ordering=sign,
                          relative_difference=a["peak"]/b["peak"]-1 if valid else None)
        if valid:
            signs.append(sign)
            diagnostic_signs.append(_order(a["peak"], b["peak"]))
            first_values.append(a["peak"])
            second_values.append(b["peak"])
    complete = len(signs) == len(SETTINGS)
    a_range, b_range = _range(first_values), _range(second_values)
    relation, ordering, diagnostic_relation = ("incomplete_or_ineligible",)*3
    if complete:
        relation = ("first_below_second" if a_range[1] < b_range[0] else
                    "second_below_first" if b_range[1] < a_range[0] else "overlap_or_touch")
        diagnostic_relation = ("first_below_second" if _order(a_range[1], b_range[0]) < 0 else
                               "second_below_first" if _order(b_range[1], a_range[0]) < 0 else "overlap_or_touch")
        if same_design:
            ordering = "same_design"
        elif all(sign == 0 for sign in signs):
            ordering = "exact_tie_all_settings"
        elif all(sign < 0 for sign in signs):
            ordering = "first_lower_all_settings"
        elif all(sign > 0 for sign in signs):
            ordering = "second_lower_all_settings"
        elif min(signs) < 0 < max(signs):
            ordering = "setting_order_reversal"
        else:
            ordering = "first_not_higher_with_ties" if min(signs) < 0 else "second_not_higher_with_ties"
    return dict(settings=rows, all_settings_eligible=complete, paired_ordering=ordering,
                first_peak_range=a_range, second_peak_range=b_range, range_relation=relation,
                replay_tolerance_range_relation=diagnostic_relation,
                replay_tolerance_setting_orderings=dict(zip(SETTINGS, diagnostic_signs)) if complete else None,
                fine_ordering_resolved_by_ranges=relation in ("first_below_second", "second_below_first"),
                same_design=same_design, replay_diagnostic_rtol=NATIVE_RTOL,
                primary_range_rule="exact finite-setting envelopes; strict separation, closed overlap",
                confidence_interval=False, continuum_error_bound=False)


def require_terminal_audit(frozen, audit):
    expected = {(case["case_id"], seed, policy) for case in frozen["cases"]
                for seed in frozen["seeds"] for policy in frozen["policies"]}
    keys = [(row["case_id"], row["seed"], row["policy"]) for row in audit["rows"]]
    require(len(expected) == 27 and len(keys) == 27 and set(keys) == expected,
            "report requires exactly 27 registered audited arms")
    require(all(row["status"] in ("completed", "arm_failed", "process_failed") for row in audit["rows"]),
            "all 27 arms must be terminal and independently audited")


def load_live(study):
    from verification.nf_anytime_reference import load_audited_arms, validate_frozen
    from verification.nf_anytime_campaign import validate_initialization
    frozen = validate_frozen(study)
    audit = read_json(study/"audit.json")
    require_terminal_audit(frozen, audit)
    arms, hashes = load_audited_arms(study, frozen)
    for arm in arms:
        prefix = Path("live")/arm["case_id"]/str(arm["seed"])/arm["policy"]
        arm["process"] = read_json(study/prefix/"process.json")
        if arm["native_pool_eligible"]:
            record = read_json(study/prefix/"result.json")
            validate_initialization(record)
            for key in ("elapsed_s", "startup_s", "stop", "error", "parity", "baseline_ready_s"):
                arm[key] = record.get(key)
            for cutoff in CHECKPOINTS:
                actual, certified = reconstruct(arm, cutoff), arm["checkpoints"][str(cutoff)]
                require(actual["available"] == certified["available"], "report reconstruction availability differs")
                if actual["available"]:
                    require(all(actual[k] == certified[k] for k in ("kind", "id", "candidate", "peak", "A2", "first_5pct_s")),
                            "independent report reconstruction differs from audited checkpoints")
    return frozen, arms, hashes


def _sealed(study, kind, name, hashes):
    relative = Path(name)
    require(not relative.is_absolute() and ".." not in relative.parts, "unsafe post-trial evidence path")
    path = study/kind/relative
    expected = path.with_name(path.name+".sha256").read_text(encoding="ascii").strip()
    actual = digest(path)
    require(actual == expected, "post-trial byte seal differs: "+str(path))
    hashes[(Path(kind)/relative).as_posix()] = actual
    return json.loads(path.read_text(encoding="utf-8"))


def _post_context(study, kind, hashes):
    if not (study/kind/"summary.json").exists():
        return None
    protocol = _sealed(study, kind, "protocol.json", hashes)
    require(protocol["campaign_sha256"] == digest(study/"frozen.json"), "post-trial campaign hash differs")
    for name, expected in protocol["source_hashes"].items():
        snapshot = study/kind/"code-snapshot"/name
        require(digest(snapshot) == expected and digest(CODE_ROOT/name) == expected,
                "post-trial code provenance differs: "+name)
        require(snapshot.with_name(snapshot.name+".sha256").read_text(encoding="ascii").strip() == expected,
                "post-trial source snapshot seal differs")
    live = _sealed(study, kind, "live-evidence.json", hashes)
    for name, expected in live.items():
        require(digest(study/name) == expected, "post-trial live evidence differs: "+name)
    summary = _sealed(study, kind, "summary.json", hashes)
    require(summary["protocol_sha256"] == hashes[kind+"/protocol.json"], "post-trial summary protocol differs")
    return summary


def load_reference(study, frozen, hashes):
    from verification.nf_anytime_reference import validate_candidate
    summary = _post_context(study, "reference", hashes)
    if summary is None:
        return dict(available=False, reason="고정 패널 사후 계산이 완료되지 않았거나 실행되지 않았습니다.")
    expected = {(case["case_id"], seed) for case in frozen["cases"] for seed in SEEDS}
    require({(row["case_id"], row["seed"]) for row in summary["rows"]} == expected and len(summary["rows"]) == 9,
            "reference summary does not contain exactly nine registered panels")
    rows, diagnostic_panels = [], []
    for row in summary["rows"]:
        prefix = row["case_id"]+"/"+str(row["seed"])
        result = _sealed(study, "reference", prefix+"/result.json", hashes)
        require(row["result_sha256"] == hashes["reference/"+prefix+"/result.json"], "panel result hash differs")
        panel = _sealed(study, "reference", prefix+"/panel.json", hashes)
        require(result["panel_sha256"] == hashes["reference/"+prefix+"/panel.json"], "locked panel hash differs")
        require([(r["id"], r["candidate"]) for r in panel["rows"]] ==
                [(r["id"], r["candidate"]) for r in result["rows"]], "panel identities changed before scoring")
        for item in result["rows"]:
            validate_candidate(item["id"], item["candidate"])
        for name, expected_hash in result["score_hashes"].items():
            scored = _sealed(study, "reference", name, hashes)
            require(hashes["reference/"+name] == expected_hash, "panel score hash differs")
            original = next(r for r in result["rows"] if r["id"] == scored["id"])
            require(scored["measurement"] == original[scored["level"]], "panel score and result differ")
        calculated = panel_metrics(result["rows"], result["baseline"])
        require(math.isclose(calculated["uniform_expected_peak"], row["uniform_expected_peak"], rel_tol=1e-14),
                "independent exact equal-k expectation differs")
        require(calculated["promoted_ids"] == row["promoted_ids"] and
                calculated["selected_peak"] == row["selected"]["peak"] and
                calculated["reference_peak"] == row["reference_best"]["peak"], "independent panel selection differs")
        overlap = result.get("live_overlap")
        require(overlap == panel.get("live_overlap") == row.get("live_overlap"), "panel overlap provenance differs")
        if overlap is not None:
            require(overlap["live_evidence_sha256"] == hashes["reference/live-evidence.json"], "panel overlap live evidence hash differs")
        rows.append(dict(case_id=row["case_id"], seed=row["seed"], complete=result["complete"],
                         live_overlap=panel_overlap(overlap, [r["id"] for r in panel["rows"]]), **calculated))
        diagnostic_panels.append(dict(case_id=row["case_id"], seed=row["seed"], rows=result["rows"],
                                      baseline=result["baseline"], panel_sha256=result["panel_sha256"],
                                      result_sha256=row["result_sha256"]))
    return dict(available=True, rows=rows, costs=summary["costs"], equal_wall_time_comparison=False,
                **panel_diagnostics(diagnostic_panels))


def load_uncertainty(study, frozen, hashes):
    from verification.nf_anytime_uncertainty import assess_case
    summary = _post_context(study, "uncertainty", hashes)
    if summary is None:
        return dict(available=False, reason="고정 후보의 수치 설정 사후 계산이 완료되지 않았거나 실행되지 않았습니다.")
    selection = _sealed(study, "uncertainty", "selection.json", hashes)
    require(summary["selection_sha256"] == hashes["uncertainty/selection.json"], "uncertainty selection hash differs")
    require({r["case_id"] for r in selection["cases"]} == {r["case_id"] for r in frozen["cases"]},
            "uncertainty selection cases differ")
    cases = []
    for selected in selection["cases"]:
        case = selected["case_id"]
        quality = _sealed(study, "uncertainty", case+"/quality.json", hashes)
        claimed = next(row for row in summary["cases"] if row["case_id"] == case)
        require(claimed["quality_sha256"] == hashes["uncertainty/"+case+"/quality.json"], "quality hash differs")
        measured = {r["id"]: {} for r in selected["candidates"]}
        locked = {r["id"]: r["candidate"] for r in selected["candidates"]}
        for name, expected in quality["score_hashes"].items():
            score = _sealed(study, "uncertainty", name, hashes)
            require(hashes["uncertainty/"+name] == expected, "uncertainty score hash differs")
            require(score["selection_sha256"] == hashes["uncertainty/selection.json"] and score["case_id"] == case,
                    "uncertainty score selection or equilibrium differs")
            require(score["id"] in locked and score["candidate"] == locked[score["id"]]
                    and score["setting"] in SETTINGS, "uncertainty score candidate or setting differs")
            measured[score["id"]][score["setting"]] = score["measurement"]
        require(all(set(values) == set(SETTINGS) for values in measured.values()),
                "completed uncertainty summary lacks a registered setting score")
        calculated = assess_case(selected, measured)
        require(all(quality[key] == value for key, value in calculated.items()), "recomputed uncertainty assessment differs")
        pairs, nearest = [], []
        coverage = native_replay_coverage(calculated["candidates"])
        require(all(coverage[key] == value for key, value in calculated["native_replay_coverage"].items()),
                "native replay coverage differs from recomputed candidate observations")
        native_passed = {row["id"]: observable_native_replay_matched(row["native_reproducibility"])
                         for row in calculated["candidates"]}
        choice_groups = {}
        neighbors = {r["selected_id"]: r["neighbor_id"] for r in selected["nearest_observed"]}
        for choice in selected["checkpoint_choices"]:
            choice_groups.setdefault((choice["seed"], choice["cutoff_s"]), []).append(choice)
            neighbor = neighbors.get(choice["id"]) if choice["available"] else None
            row = dict(case_id=case, seed=choice["seed"], cutoff_s=choice["cutoff_s"], policy=choice["policy"],
                       selected_id=choice["id"], neighbor_id=neighbor, available=neighbor is not None)
            if row["available"]:
                row.update(setting_pair(measured[choice["id"]], measured[neighbor]))
                row["native_replay_verified_for_pair"] = native_passed[choice["id"]] and native_passed[neighbor]
                row["resolved_for_live_choices"] = row["fine_ordering_resolved_by_ranges"] and row["native_replay_verified_for_pair"]
            nearest.append(row)
        for (seed, cutoff), choices in sorted(choice_groups.items()):
            by_policy = {choice["policy"]: choice for choice in choices}
            for first, second in combinations(POLICIES, 2):
                a, b = by_policy[first], by_policy[second]
                row = dict(case_id=case, seed=seed, cutoff_s=cutoff, first_policy=first, second_policy=second,
                           first_id=a["id"], second_id=b["id"], available=a["available"] and b["available"])
                if row["available"]:
                    row.update(setting_pair(measured[a["id"]], measured[b["id"]], same_design=a["id"] == b["id"]))
                    row["native_replay_verified_for_pair"] = native_passed[a["id"]] and native_passed[b["id"]]
                    row["resolved_for_live_choices"] = row["fine_ordering_resolved_by_ranges"] and row["native_replay_verified_for_pair"]
                pairs.append(row)
        cases.append(dict(case_id=case, candidates=calculated["candidates"],
                          native_reproducibility_passed=calculated["native_reproducibility_passed"],
                          all_available_native_observations_matched=calculated["all_available_native_observations_matched"],
                          native_replay_coverage=coverage, native_replay_scope=calculated["native_replay_scope"],
                          checkpoint_choices=selected["checkpoint_choices"],
                          method_checkpoint_pairs=pairs, nearest_checkpoint_pairs=nearest,
                          observed_native_pool_size=selected["observed_native_pool_size"]))
    coverage = {key: sum(case["native_replay_coverage"][key] for case in cases)
                for key in cases[0]["native_replay_coverage"]}
    matched = all(case["all_available_native_observations_matched"] for case in cases)
    require(all(coverage[key] == value for key, value in summary["native_replay_coverage"].items()) and
            matched == summary["all_available_native_observations_matched"], "global native replay coverage or match status differs")
    return dict(available=True, cases=cases, costs=summary["costs"],
                native_reproducibility_passed=matched, all_available_native_observations_matched=matched,
                native_replay_coverage=coverage, native_replay_scope=summary["native_replay_scope"],
                range_scope="finite numerical settings within the same transport model; not confidence intervals or continuum error bounds")


def rejection_case_metrics(selected, scores):
    """Independently classify native outcomes; a failed score is unknown, not negative."""
    from verification.nf_anytime_reference import finite_positive, validate_measurement
    identities = [row["id"] for row in selected["candidates"]]
    require(len(set(identities)) == len(identities) and set(scores) == set(identities) and BASELINE_ID in scores,
            "rejection native score identities differ from selection")
    for scored in scores.values():
        validate_measurement(scored["measurement"])
        require(scored["status"] in ("completed", "failed", "not_evaluated_preparation_failed"),
                "rejection native score status is unknown")

    def compare(scored, reference, available):
        state = scored["measurement"]
        peak_known = bool(available and finite_positive(state.get("peak")) and finite_positive(reference.get("peak")))
        area_known = bool(available and finite_positive(state.get("A2")) and finite_positive(reference.get("A2")))
        eligible = bool(available and reference.get("eligible", True) and scored["status"] == "completed" and state["eligible"])
        return dict(reference=reference, reference_available=available,
                    peak_ratio=state["peak"]/reference["peak"] if peak_known else None,
                    peak_reduction_pct=100*(1-state["peak"]/reference["peak"]) if peak_known else None,
                    A2_change_pct=100*(state["A2"]/reference["A2"]-1) if area_known else None,
                    strictly_lower_observed_peak=state["peak"] < reference["peak"] if peak_known else None,
                    eligible_peak_improvement=bool(eligible and peak_known and state["peak"] < reference["peak"]),
                    diagnostic_only=True)

    baseline, rows = scores[BASELINE_ID], []
    for candidate in selected["candidates"]:
        if candidate["id"] == BASELINE_ID:
            continue
        scored = scores[candidate["id"]]
        eligibility = scored["measurement"]["eligible"] if scored["status"] == "completed" else None
        comparisons, seen = [], set()
        for observation in candidate["qualifying_observations"]:
            key = observation["seed"], observation["policy"]
            if key not in seen:
                incumbent = observation["incumbent_600"]
                comparisons.append(dict(seed=key[0], policy=key[1], **compare(scored, incumbent, incumbent["available"])))
                seen.add(key)
        rows.append(dict(candidate, fine=scored, fine_eligible=eligibility,
                         coarse_eligibility_false_negative=eligibility,
                         diagnostic_status="fine_evaluation_failed" if eligibility is None else
                                           "fine_eligible_false_negative" if eligibility else "fine_ineligible",
                         vs_baseline=compare(scored, baseline["measurement"],
                                             baseline["status"] == "completed" and baseline["measurement"]["eligible"]),
                         vs_live_coarse_incumbents=comparisons))
    return dict(case_id=selected["case_id"], baseline=baseline, rows=rows,
                false_negative_count=sum(row["coarse_eligibility_false_negative"] is True for row in rows),
                fine_ineligible_count=sum(row["coarse_eligibility_false_negative"] is False for row in rows),
                fine_evaluation_failed_count=sum(row["coarse_eligibility_false_negative"] is None for row in rows),
                post_hoc=True, primary_winners_changed=False)


def load_rejections(study, frozen, arms, hashes):
    kind = "screen-rejections"
    summary = _post_context(study, kind, hashes)
    if summary is None:
        return dict(available=False, reason="관측된 coarse 열결과 탈락의 별도 사후 native 진단이 완료되지 않았거나 실행되지 않았습니다.")
    from verification.nf_anytime_rejections import protocol_details, select_rejections, source_paths
    protocol = _sealed(study, kind, "protocol.json", hashes)
    selection = _sealed(study, kind, "selection.json", hashes)
    live_hashes = _sealed(study, kind, "live-evidence.json", hashes)
    require(protocol["details"] == protocol_details() and
            set(protocol["source_hashes"]) == {p.relative_to(CODE_ROOT).as_posix() for p in source_paths()} and
            protocol["reference"] == frozen["reference"] and protocol["source"] == "OMP" and
            protocol["radius_mm"] == 20 and protocol["backend"] == "gpu", "rejection native protocol differs")
    metadata = dict(protocol_sha256=hashes[kind+"/protocol.json"], selection_sha256=hashes[kind+"/selection.json"])
    require(all(summary[key] == value for key, value in metadata.items()) and
            selection["protocol_sha256"] == metadata["protocol_sha256"], "rejection selection or protocol hash differs")
    require(selection["live_evidence_hashes"] == live_hashes and
            all(hashes.get(name) == expected for name, expected in live_hashes.items()), "rejection raw evidence differs from audited main records")
    canonical = select_rejections(arms, live_hashes)
    require(all(selection[key] == value for key, value in canonical.items()), "rejection locked selection differs from audited raw events")
    require(summary["post_hoc"] is True and summary["primary_winners_changed"] is False,
            "rejection summary changed the primary comparison scope")
    selected_ids = {case["case_id"] for case in canonical["cases"]}
    require(selected_ids <= {case["case_id"] for case in frozen["cases"]} and
            len(summary["cases"]) == len(selected_ids) and {case["case_id"] for case in summary["cases"]} == selected_ids,
            "rejection result cases differ from selection")
    cases = []
    count_keys = ("false_negative_count", "fine_ineligible_count", "fine_evaluation_failed_count")
    for selected in canonical["cases"]:
        case = selected["case_id"]
        result = _sealed(study, kind, case+"/result.json", hashes)
        claimed = next(row for row in summary["cases"] if row["case_id"] == case)
        require(claimed["result_sha256"] == hashes[kind+"/"+case+"/result.json"] and
                all(result[key] == value for key, value in metadata.items()), "rejection result hash or selection provenance differs")
        expected = {case+"/scores/"+row["id"]+".cubic1280.json": row for row in selected["candidates"]}
        require(set(result["score_hashes"]) == set(expected), "rejection native score set differs")
        scores = {}
        for name, candidate in expected.items():
            scored = _sealed(study, kind, name, hashes)
            require(result["score_hashes"][name] == hashes[kind+"/"+name], "rejection native score hash differs")
            identity = dict(id=candidate["id"], candidate=candidate["candidate"], case_id=case,
                            setting="cubic1280", setting_tuple=["cubic", 1280, .000125, 641, .00025], **metadata)
            require(all(scored[key] == value for key, value in identity.items()), "rejection native score provenance differs")
            scores[candidate["id"]] = scored
        calculated = rejection_case_metrics(selected, scores)
        require(all(result[key] == value for key, value in calculated.items()) and
                all(claimed[key] == calculated[key] for key in count_keys), "rejection independently reconstructed outcomes differ")
        cases.append(dict(calculated, result_sha256=claimed["result_sha256"], score_hashes=result["score_hashes"]))
    totals = {key: sum(case[key] for case in cases) for key in count_keys}
    require(all(summary[key] == value for key, value in totals.items()) and
            all(summary[key] == canonical[key] for key in ("selected_candidate_count", "qualifying_observation_count", "complete_live_coverage")),
            "rejection aggregate counts differ from observed native statuses")
    return dict(available=True, cases=cases, costs=summary["costs"], post_hoc=True, primary_winners_changed=False,
                selected_candidate_count=canonical["selected_candidate_count"],
                qualifying_observation_count=canonical["qualifying_observation_count"],
                complete_live_coverage=canonical["complete_live_coverage"], uncertified_arms=canonical["uncertified_arms"],
                examined_case_ids=canonical["examined_case_ids"], examined_arm_count=canonical["examined_arm_count"],
                native_setting="cubic1280", scope=canonical["scope"], **totals)


def make_summary(study):
    frozen, arms, hashes = load_live(study)
    amendment = read_json(study/"amendment.json")
    require(frozen["population"] == 64 and frozen["projected_seeds"] == 12 and
            amendment["old_population"] == 12 and amendment["new_population"] == 64 and
            amendment["projected_seeds"] == 12 and amendment["comparison_outcomes_used"] is False,
            "full-space protocol amendment differs from the report scope")
    hashes["amendment.json"] = digest(study/"amendment.json")
    prior = Path(amendment["prior_study"])
    require(digest(prior/"frozen.json") == amendment["prior_protocol_sha256"] and
            digest(prior/"interruption.json") == amendment["interruption_sha256"], "aborted-run amendment evidence changed")
    reference = load_reference(study, frozen, hashes)
    uncertainty = load_uncertainty(study, frozen, hashes)
    rejections = load_rejections(study, frozen, arms, hashes)
    return dict(schema="nf-anytime-report.v1", generated_utc=datetime.now(timezone.utc).isoformat(),
                study=str(study), source_sha256=digest(__file__), input_hashes=hashes,
                frozen_case_ids=[row["case_id"] for row in frozen["cases"]],
                seeds=list(SEEDS), policies=list(POLICIES), checkpoints_s=list(CHECKPOINTS),
                reference_inputs=frozen["reference"], amendment=amendment,
                live=summarize_live(arms), finite_panel=reference, numerical_settings=uncertainty,
                screen_rejections=rejections,
                interpretation=dict(original_wall="supplied DIII-D original material wall", source="OMP",
                                    already_studied_equilibria=True, holdout_claim=False,
                                    global_optimum_claim=False, independent_seed_discharge_claim=False,
                                    numerical_ranges_are_confidence_intervals=False,
                                    continuum_error_bound=False,
                                    remaining_scope=["FIREFLY cross-check", "continuous first-hit transport",
                                                     "1800 s extension", "IMP and other operating points"]))


def _fmt(value, digits=3):
    return "미확인" if value is None else f"{value:.{digits}f}"


def _span(values, digits=2):
    return "미확인" if values is None else f"{values[0]:.{digits}f}–{values[1]:.{digits}f}"


def _table(headers, rows):
    clean = lambda value: str(value).replace("|", "/").replace("\n", " ")
    return "\n".join(["| "+" | ".join(headers)+" |", "| "+" | ".join("---" for _ in headers)+" |"] +
                     ["| "+" | ".join(clean(v) for v in row)+" |" for row in rows])


def panel_diagnostics_text(panel):
    counts, transfer = panel["shape_counts"], panel["posthoc_transfer"]
    statuses = dict(baseline="기준벽", eligible="적격", ineligible="부적격", unavailable="미관측")
    return (f"패널 fine 평가 기록 {counts['fine_measurement_rows']}행은 "
            f"서로 다른 형상 {counts['global_unique_shape_count']}개에 해당한다. "
            f"평형별 고유 형상 수의 합은 {counts['candidate_case_count']}이며 "
            f"모든 평형에 공통인 형상은 {counts['common_to_all_shape_count']}개다. "
            "이 수에는 원래 기준벽을 넣지 않았다. 평가 기록 수를 독립 형상 수나 총 열해석 호출 수로 해석하지 않는다.\n\n"+
            _table(["평형", "패널 수", "fine 기록 행", "고유 형상", "평형 내 반복 행"],
                   [(r["case_id"], r["panel_count"], r["fine_measurement_rows"], r["unique_shape_count"],
                     r["duplicate_measurement_rows"]) for r in counts["per_case"]])+
            "\n\n### 사후 입력 간 재사용\n\n"
            "각 출발 평형의 세 고정 패널 합집합(최대 72개 후보)과 원래 기준벽에서 관측 fine 피크가 가장 낮은 형상을 고른 뒤, "
            "같은 ID가 도착 평형의 기존 패널에 있을 때만 그 결과를 읽었다. 출발 선택이 기준벽이면 도착 평형의 기준벽을 쓴다. "
            "동률이면 기준벽을 유지하며, 반복 ID의 적격성과 관측값을 native 재현 허용치 1e-10으로 대조한 후 원값 최저를 보존한다. "
            "기준벽의 반복 관측도 대조하고 고정 시드 순서의 첫 기록을 기준으로 쓴다. "
            "이 허용치는 물리적 불확실성 경계가 아니다. 추가 열계산이나 미관측·부적격 값의 대입은 없다. "
            "음의 피크 감소율은 도착 기준벽보다 나빠졌다는 뜻이다. "
            "이는 사후 유한 집합 진단이며 holdout 전이 검증·전역 최적성·실시간 600초 선택 결과가 아니다. "
            "선택 ID 전체와 원본 패널 해시·시드·행 번호·반복 관측은 summary.json에 있다.\n\n"+
            _table(["출발 평형", "도착 평형", "선택 ID 앞 12자", "상태", "피크 MW/m²", "A2 m²", "피크 감소 %", "A2 변화 %"],
                   [(r["source_case_id"], r["destination_case_id"],
                     "원래 기준벽" if r["id"] == BASELINE_ID else r["id"][:12], statuses[r["status"]],
                     _fmt(r["peak"]), _fmt(r["A2"]), _fmt(r["peak_reduction_pct"]), _fmt(r["A2_change_pct"]))
                    for r in transfer["rows"]]))


def native_replay_intro(numerical):
    rows = []
    for case in numerical["cases"]:
        coverage = case["native_replay_coverage"]
        rows.append((case["case_id"], coverage["locked_candidate_count"],
                     coverage["previously_scored_candidate_count"], coverage["candidates_with_observable_comparisons"],
                     coverage["eligibility_only_candidate_count"], coverage["unobserved_candidate_count"],
                     coverage["other_nonobservable_candidate_count"], coverage["compared_observation_count"],
                     coverage["compared_observable_count"], coverage["eligibility_only_observation_count"]))
    return (f"대조 가능한 이전 native 관측이 모두 일치했는지: {numerical['all_available_native_observations_matched']}. "
            "고정 후보 수에는 각 평형의 기준벽도 포함된다. 기존 native 관측이 없는 후보는 새로 계산했더라도 재현 일치 사례가 아니다. "
            "예를 들어 coarse-only 초과 후보는 이전 정밀 관측이 없을 수 있다. 적격성만 비교한 경우도 피크/A2 관측값 일치와 구분한다. "
            "부적격 결과라도 양쪽에 유한 양수 관측값이 있으면 대조하며 관측값 존재 여부가 달라지면 실패로 남긴다.\n\n"+
            _table(["평형", "고정 후보", "이전 관측 후보", "관측값 대조 후보", "적격성만 후보", "이전 관측 없음",
                    "기타 관측값 없음", "대조 관측 건", "피크/A2 대조 값", "적격성만 관측 건"], rows)+
            "\n\n관측 건수는 동일 후보의 여러 저장 결과를 포함하며 후보 수와 다르다. "
            "표의 '관측값 대조 후보'는 실제 대조가 있었다는 수이며 모두 일치했다는 개별 성공 수로 대신 쓰지 않는다.")


def fine_adaptive_phase_text(phase):
    kinds = phase["completed_proposal_kinds"]
    count = kinds.get("de_trial", 0)
    observed = (f"정밀 적응 탐색 {phase['zero_completed_de_trial_arm_count']}/{phase['arm_count']}개 실행에서 "
                "600초 안의 DE trial 완료가 없었다. ")
    if count == 0:
        observed += "전체 DE trial 완료 0건이다. "
    else:
        observed += f"전체 DE trial 완료는 {count}건이다. "
    observed += (f"완료 요청은 초기 집단 평가 {kinds.get('population_initialization', 0)}건, "
                 f"좌표 polish {kinds.get('coordinate_polish', 0)}건이며 캐시 요청을 포함한다. ")
    if phase["only_initialization_and_polish"]:
        observed += "기록된 완료 단계는 초기 집단 평가와 좌표 polish뿐이다. "
    if phase["certified_arm_count"] != phase["arm_count"]:
        observed += f"이 중 원시 event가 인증된 실행은 {phase['certified_arm_count']}개다. "
    return observed+("따라서 이번 비교는 native 기본 population=64/projected_seeds=12와 600초 예산의 설정 비교다. "
                     "방법별 최적 튜닝이나 충분한 시간에서의 점근적 최적화 우열을 검증한 결과로 확장하지 않는다.")


def late_result_text(live):
    rows, improved = [], []
    for policy in POLICIES:
        diagnostics = [a["late_result_diagnostic"] for a in live["arms"] if a["policy"] == policy]
        rows.append((LABELS[policy], *[sum(d[key] for d in diagnostics) for key in
                     ("late_coarse_count", "late_fine_count", "genuine_eligible_late_fine_count", "compared_fine_count",
                      "strict_improvement_count", "improvement_beyond_replay_tolerance_count")]))
    for arm in live["arms"]:
        diagnostic = arm["late_result_diagnostic"]
        for row in diagnostic["fine_results"]:
            if row["strict_improvement"]:
                improved.append((arm["case_id"], arm["seed"], LABELS[arm["policy"]], row["id"][:12],
                                 _fmt(row["end_s"]), format(diagnostic["incumbent_at_cutoff"]["peak"], ".12g"),
                                 format(row["measurement"]["peak"], ".12g"), format(100*row["relative_improvement"], ".8g"),
                                 row["replay_tolerance_ordering"]))
    text = ("시한을 넘긴 coarse/fine 점수는 보존만 하며 승자·곡선·최초 5% 시각에 넣지 않는다. "
            "아래 비교는 캐시가 아니고 침적 결과가 확인된 전체 적격 late fine만 600초 채택 fine incumbent와 대조한다. "
            "부분 계산·관측값 누락·부적격 결과와 채택 기준이 없는 실행은 비교 불가로 남긴다. "
            "coarse 값은 개수만 세며 fine 피크와 수치 비교하지 않는다. 원값의 엄격한 개선과 상대 1e-10 재현 허용차 진단을 구분한다. "
            "후자는 물리적 불확실성이나 늦은 결과의 채택 근거가 아니다.\n\n"+
            _table(["방법", "late coarse", "late fine", "새 적격 fine", "대조 가능", "엄격 개선", "재현 허용차 밖 개선"], rows))
    if improved:
        text += ("\n\n600초 incumbent보다 원값 피크가 낮았던 late fine만 아래에 보인다. 작은 차이도 표시하며 개선율은 600초 incumbent 대비다.\n\n"+
                 _table(["평형", "시드", "방법", "ID 앞 12자", "완료 초", "600초 피크", "late 피크", "추가 감소 %", "재현 허용차 진단"], improved))
    return text


def numerical_lead(numerical):
    if not numerical["available"]:
        return ""
    descriptions = []
    for other in ("fine_adaptive", "random_fine"):
        pairs = [row for case in numerical["cases"] for row in case["method_checkpoint_pairs"]
                 if row["cutoff_s"] == 600 and row["first_policy"] == "coarse_adaptive" and row["second_policy"] == other]
        relations = Counter(row.get("range_relation", "incomplete_or_ineligible") for row in pairs)
        paired_lower = sum(row.get("paired_ordering") == "first_lower_all_settings" for row in pairs)
        verified_lower = sum(row.get("range_relation") == "first_below_second" and
                             row.get("native_replay_verified_for_pair", False) for row in pairs)
        descriptions.append(f"저비용 / {LABELS[other]}: 저비용 쪽으로 범위 분리 {relations['first_below_second']}/{len(pairs)}, "
                            f"중첩·접촉 {relations['overlap_or_touch']}, 반대 분리 {relations['second_below_first']}, "
                            f"미확인·부적격 {relations['incomplete_or_ineligible']}; "
                            f"같은 설정마다 저비용 피크가 낮음 {paired_lower}/{len(pairs)}. "
                            f"저비용 쪽 분리 중 이전 native 관측 재현까지 확인된 쌍은 {verified_lower}개다.")
    return ("600초 선택의 네 수치 설정 범위를 교차 비교하면 "+" ".join(descriptions)+
            " 같은 설정별 방향 일치와 설정 전체 범위의 분리는 별도 기준이다. "
            "원값의 작은 승패만으로 미세 우열을 확정하지 않으며, 이 범위는 신뢰구간·연속체 오차 상한이 아니다.")


def rejections_text(rejections):
    if not rejections["available"]:
        return rejections["reason"]
    cases, details = [], []
    labels = {True: "fine 적격/false negative", False: "fine 부적격", None: "알 수 없음"}
    for case in rejections["cases"]:
        cases.append((case["case_id"], len(case["rows"]), case["false_negative_count"],
                      case["fine_ineligible_count"], case["fine_evaluation_failed_count"]))
        for row in case["rows"]:
            observations = row["qualifying_observations"]
            coarse_values = [o["event"]["measurement"] for o in observations]
            coarse_errors = sorted({error for state in coarse_values for error in state.get("errors", [])})
            fine = row["fine"]["measurement"]
            fine_errors = list(fine.get("errors", []))
            if row["fine"].get("error"):
                fine_errors.append(row["fine"]["error"])
            comparisons = "; ".join(str(r["seed"])+": "+_fmt(r["peak_reduction_pct"])+"%"
                                    for r in row["vs_live_coarse_incumbents"])
            details.append((case["case_id"], row["id"][:12], labels[row["coarse_eligibility_false_negative"]],
                            "; ".join(_fmt(s.get("peak"))+" / "+_fmt(s.get("A2")) for s in coarse_values),
                            _fmt(fine.get("peak"))+" / "+_fmt(fine.get("A2")),
                            _fmt(row["vs_baseline"]["peak_reduction_pct"]), comparisons,
                            "; ".join(coarse_errors) or "없음", "; ".join(fine_errors) or "없음"))
    return (f"인증된 main event 중 시간 내 완료·비캐시·침적 결과 반환·전체 부적격인 coarse 관측 "
            f"{rejections['qualifying_observation_count']}건에서 평형별 중복을 제거한 "
            f"{rejections['selected_candidate_count']}개 후보를 추가로 확인했다. "
            "열유속 값·오류 종류·이전 fine 결과로 대상을 좁히지 않았고, 목록은 추가 계산 전에 고정했다. "
            "native cubic1280만의 조건부 사후 진단이며 전역 false-negative 비율, 사전 등록 endpoint, 동일 시간 성능 비교가 아니다. "
            "fine 적격이면 기준벽 개선 여부와 관계없이 coarse 적격성 false negative로 세며, 계산 실패는 알 수 없음으로 남긴다. "
            "주요 600초 승자·곡선·개선 시각은 바뀌지 않는다.\n\n"+
            _table(["평형", "선택 후보", "확인된 false negative", "fine 부적격", "알 수 없음"], cases)+
            "\n\n피크/A2의 단위는 각각 MW/m²와 m²다. 비율은 관측값의 진단용 비교이며 부적격·실패 점수를 채택하지 않는다. "
            "600초 비교는 해당 coarse 탈락을 관측한 시드의 고정 coarse-adaptive 승자를 기준으로 한다.\n\n"+
            _table(["평형", "ID 앞 12자", "상태", "coarse 피크/A2", "native fine 피크/A2", "기준벽 대비 피크 감소 %",
                    "시드별 600초 대비", "coarse 오류", "fine 오류"], details)+
            "\n\n전체 coarse 관측·원시 해시·native score 상태·오류와 A2 변화, 기준벽 및 시드별 600초 원값 비교는 "
            "summary.json의 screen_rejections에 보존한다. 추가 native 계산 비용도 이 항목의 costs에 따로 기록하며 600초 예산에 섞지 않는다.")


def korean_report(summary):
    live = summary["live"]
    final = [r for r in live["checkpoint_metrics"] if r["cutoff_s"] == 600]
    observed = sum(a["first_5pct"]["event_observed"] for a in live["arms"])
    lead_pairs = []
    for other in ("fine_adaptive", "random_fine"):
        pairs = [r for r in live["paired_comparisons"] if r["cutoff_s"] == 600 and
                 r["first_policy"] == "coarse_adaptive" and r["second_policy"] == other]
        counts = Counter(r["ordering"] for r in pairs)
        lead_pairs.append(f"저비용 적응 / {LABELS[other]}의 명목 피크 비교는 "
                          f"{counts['first_lower']}/{counts['exact_tie']}/{counts['second_lower']}"
                          f"(첫째가 낮음/정확히 같음/둘째가 낮음), 미확인 {counts['unavailable']}쌍이다.")
    lines = ["# DIII-D 고정 입력에서 시간별 탐색 성능과 수치 설정 민감도", "",
             "600초의 같은 평형·시드별 비교에서 "+" ".join(lead_pairs)+
             " 각 비교는 세 평형 안의 세 시드, 총 아홉 쌍이며 독립 방전 아홉 개가 아니다. "
             "작은 native 차이도 사전 고정한 엄격한 < 규칙의 명목 승패에 남기고, 수치 설정에서의 분리 여부는 별도로 판정한다.", "",
             numerical_lead(summary["numerical_settings"]), "",
             _table(["방법", "600초 사용 가능", "기준 대비 q 감소", "q↓·A2↑"],
                    [(LABELS[policy], sum(r["available"] for r in final if r["policy"] == policy),
                      sum(r["q_down"] for r in final if r["policy"] == policy),
                      sum(r["q_down_A2_up"] for r in final if r["policy"] == policy)) for policy in POLICIES]), "",
             f"세 평형의 시드 반복을 포함한 {live['arm_count']}개 실행을 종료 상태와 원시 event로 대조했다. "
             f"600초에 기준벽보다 피크가 낮은 실행은 {sum(r['q_down'] for r in final)}/{len(final)}, "
             f"피크 감소와 A2 증가가 함께 확인된 실행은 {sum(r['q_down_A2_up'] for r in final)}/{len(final)}이다. "
             f"최초 5% 감소는 {observed}/{len(live['arms'])}개 실행에서 시간 내 관측됐다. "
             "분모는 방법별 실행 수이며 독립 방전 수는 아니다.", "",
             "같은 모델의 관측 후보 중 시간 내 확인한 최저 피크를 비교한다. 전역 최적이나 유일한 최적 형상을 뜻하지 않는다. "
             "세 평형은 이전 연구에서 이미 사용한 입력이며 미사용 holdout으로 부르지 않는다. "
             "OMP 1 MW와 제공된 DIII-D 원래 물질벽을 고정 기준으로 유지했다.", "",
             fine_adaptive_phase_text(live["fine_adaptive_phase_at_600"]), "",
             "## 실제 완료 시각의 best-so-far", "",
             "![평형·시드별 시간에 따른 정밀 피크 비율](anytime_curves.png)", "",
             "각 선은 그 방법의 기준벽 준비가 끝난 실제 시각에서 시작한다. 세로축은 정밀 피크/해당 기준벽 피크다. "
             "적격 정밀 평가가 시간 안에 끝났을 때만 갱신하며 저비용 점수와 초과된 결과는 채택하지 않는다. "
             "준비·임포트·CPU 대조·기록 비용을 포함한 프로세스 시작 기준 시각이다. 실패 실행의 선은 실제 종료 시각에서 멈추며, "
             "그전에 인증된 마지막 결과는 이후 checkpoint 표에서도 유지할 수 있다.", "",
             "## 200/400/600초 비교", "",
             "각 행은 같은 평형·방법의 세 시드 반복이다. 중앙값과 시드 범위는 기술 통계이며 신뢰구간이 아니다. "
             "q 감소/A2 증가와 명목 승패는 원값의 엄격한 < 및 > 규칙을 그대로 쓴다. "
             "상대 1e-10 backend 재현 허용차 진단은 별도 필드이며 물리적 수치 불확실성으로 해석하지 않는다.", "",
             _table(["평형", "방법", "초", "사용 가능/3", "피크 감소 중앙값 %", "시드 범위 %", "q 감소", "q↓·A2↑"],
                    [(g["case_id"], LABELS[g["policy"]], g["cutoff_s"], g["available_count"],
                      _fmt(g["peak_reduction_pct_median"]), _span(g["peak_reduction_pct_seed_range"]),
                      g["q_down_count"], g["q_down_A2_up_count"]) for g in live["groups"]]), "",
             "같은 평형·시드·시간의 쌍만 비교했다. 아래는 세 시드의 첫 방법 피크가 낮음/정확히 같음/둘째 방법 피크가 낮음/미확인 개수다. "
             "작은 native 점수 차이가 추가 수치 설정에서도 분리되는지는 뒤의 판정과 구분한다.", ""]
    pair_rows = []
    for case in summary["frozen_case_ids"]:
        for cutoff in CHECKPOINTS:
            for first, second in combinations(POLICIES, 2):
                rows = [r for r in live["paired_comparisons"] if (r["case_id"], r["cutoff_s"], r["first_policy"], r["second_policy"]) == (case, cutoff, first, second)]
                counts = Counter(r["ordering"] for r in rows)
                pair_rows.append((case, cutoff, LABELS[first]+" / "+LABELS[second],
                                  "/".join(str(counts[k]) for k in ("first_lower", "exact_tie", "second_lower", "unavailable"))))
    lines += [_table(["평형", "초", "방법 쌍", "첫째/정확동률/둘째/미확인"], pair_rows), "",
              "## 최초 5% 개선과 우측 검열", "",
              "미관측 실행의 개선 시각을 600초 성공으로 대입하지 않았다. 실패 실행은 실제 관측 종료 시각에서 검열하며, "
              "감독자 중단으로 인증된 checkpoint가 없는 실행은 성능 곡선에서 제외하고 실패로 남긴다.", "",
              _table(["평형", "시드", "방법", "최초 5% 초", "검열 시각 초", "상태"],
                     [(a["case_id"], a["seed"], LABELS[a["policy"]], _fmt(a["first_5pct"]["time_s"]),
                       _fmt(a["first_5pct"]["censor_time_s"]), a["first_5pct"]["reason"]) for a in live["arms"]]), "",
              "## 시한 초과 결과 진단", "", late_result_text(live), "",
              "## 실제 비용과 준비·대조 비용", "",
              "아래는 세 평형·세 시드에서 관측한 방법별 최소–최대 초다. 600초는 채택 시한이며 실제 종료 시간이 아니다. "
              "종료 중인 평가의 초과분과 최종 직렬화/부모 관측 차이를 보존했다. CPU parity 시간은 scoring event에 이미 포함되므로 "
              "더해서 이중 계산하지 않는다. 준비 시간을 빼서 warm-start 성능으로 바꾸지 않는다. "
              "감독자 timeout에서도 부모가 관측한 프로세스 종료 시각과 초과분은 포함한다. 작업자 최종 시각이 없으면 미확인으로 남기며 부모 시각으로 대체하지 않는다. "
              "실행별 stage 합계와 로그/optimizer bookkeeping 등 채택 시간에 포함된 시각 공백은 summary.json에 있다.", ""]
    cost_rows = []
    for policy in POLICIES:
        costs = [a["costs"] for a in live["arms"] if a["policy"] == policy]
        def cost_range(key, mark_unknown=False):
            values = [c[key] for c in costs if c.get(key) is not None]
            text = _span(_range(values))
            missing = len(costs)-len(values)
            return text+f" (미확인 {missing}/{len(costs)})" if mark_unknown and missing else text
        cost_rows.append((LABELS[policy], cost_range("startup_s"), cost_range("baseline_ready_s"),
                          cost_range("preparation_and_build_seconds"), cost_range("reported_cpu_parity_seconds"),
                          cost_range("admitted_scoring_event_seconds"), cost_range("elapsed_s", True),
                          cost_range("overrun_seconds", True), cost_range("process_elapsed_s", True),
                          cost_range("process_overrun_seconds", True)))
    lines += [_table(["방법", "시작 준비", "기준 준비 완료", "준비/구축 합", "보고된 CPU parity", "시간 내 scoring event 합",
                      "작업자 종료", "작업자 초과", "부모 관측 종료", "부모 관측 초과"], cost_rows), "",
              "## 요청·침적 결과·캐시·기하 탈락", "",
              "아래 수는 각 방법의 활성 탐색 수준에서 동일 평형의 세 시드를 합한 값이다. 저비용 적응 탐색의 활성 수준은 coarse이며 "
              "정밀 승격 요청은 summary.json의 fine 항목에 별도로 남긴다. '침적 결과'는 physics_evaluated=true인 새 후보 평가다. "
              "이 값은 침적 결과가 반환되었다는 뜻이며 총 열계산 진입 횟수나 CPU/GPU 호출 수가 아니다. "
              "첫 parity 대조의 CPU/GPU 시간은 별도 payload로 보존했다. '기하 탈락'은 명시적 기하 오류만 센다. "
              "단계별 초기/DE/좌표 수는 캐시를 포함한 완료 요청 수이며 새 침적 결과 수와 다르다. 생성 탈락에는 투영·벽 검사·근접중복 등의 이유별 합계가 포함된다. "
              "감독자 중단 실행의 남아 있는 원시 개수는 미인증 진단값으로 표시하며 성능 checkpoint에는 사용하지 않는다.", ""]
    count_rows = []
    coverage_rows = []
    for g in live["groups"]:
        c, kinds = g["active_counts_sum"], g["proposal_kind_counts_sum"]
        count_rows.append((g["case_id"], LABELS[g["policy"]], g["cutoff_s"], c.get("requests_started", 0),
                           c.get("requests_completed", 0), c.get("deposition_result_evaluations", 0),
                           c.get("cached_requests", 0), c.get("geometry_rejections", 0),
                           c.get("other_non_deposition_results", 0),
                           sum(g["draw_rejections_sum"].values()),
                           "/".join(str(kinds.get(k, 0)) for k in ("population_initialization", "de_trial", "coordinate_polish"))))
        coverage_rows.append((g["case_id"], LABELS[g["policy"]], g["cutoff_s"],
                              *[_span(g["coverage_rank_seed_ranges"][k], 0) for k in ("proposed", "deposition", "eligible")]))
    lines += [_table(["평형", "방법", "초", "시작", "완료", "침적 결과", "캐시", "기하 탈락", "기타 비침적", "생성 탈락", "초기/DE/좌표"], count_rows), "",
              "## 시간 내 실제 다룬 벡터의 차원", "",
              "후보 벡터를 중심화하고 고정 계수 경계 폭 0.2 m로 나눈 행렬의 rank(절대 특이값 기준 1e-12)를 다시 계산했다. "
              "제안/침적 결과/전체 적격 후보를 분리한다. 아래는 세 시드의 rank 범위이며 각 축 최소·최대·범위와 특이값은 summary.json에 있다. "
              "제안 rank는 벡터가 저장된 완료 제안에 한한다. 무작위 생성에서 벡터 없이 탈락 이유만 남은 시도는 rank에 넣을 수 없다. "
              "초기 64개 집단의 rank=16만으로 시간 내 64개 평가 완료나 DE 단계 진입, 전 형상 공간 탐색을 주장하지 않는다.", "",
              _table(["평형", "방법", "초", "제안 rank", "침적 결과 rank", "전체 적격 rank"], coverage_rows), "",
              "## 고정 유한 패널의 선택 손실", ""]
    panel = summary["finite_panel"]
    if not panel["available"]:
        lines.append(panel["reason"])
    else:
        lines += ["추가 패널 열계산 전에 기하 규칙으로 목록을 고정한 패널의 전수 정밀 최저와 비교한 손실이다. "
                  "동일 시드·생성 규칙 때문에 실시간 random_fine 후보와 겹칠 수 있으며, 새 holdout 후보나 독립 실험으로 세지 않는다. "
                  "선택 규칙은 사전 고정했고 후보 열결과를 목록 선택에 사용하지 않았다. 무작위 기대값은 공통 coarse 적격 모집단에서 "
                  "동일 정밀 요청 수 k개를 중복 없이 고르는 정확한 조합 기대값이며, fine 부적격도 요청을 소비한다. "
                  "이는 동일 벽시계 시간 비교나 실시간 random_fine 성능을 대신하지 않는다.", "",
                  _table(["평형", "시드", "패널 n", "coarse 적격", "k", "선택 손실 %", "무작위 기대 손실 %", "live 제안/정밀 중복", "24개 완성"],
                         [(r["case_id"], r["seed"], r["panel_size"], r["coarse_eligible_count"], r["fine_requests_k"],
                           _fmt(100*r["regret_relative"]), _fmt(100*r["uniform_expected_regret_relative"]),
                           (str(r["live_overlap"]["proposed_intersection_count"])+"/"+str(r["live_overlap"]["on_time_fine_intersection_count"]))
                           if r["live_overlap"]["available"] else "미확인", r["complete"]) for r in panel["rows"]]), "",
                  "중복은 동일 평형의 모든 live 시드·방법과 대조했다. 시간 내 정밀 관측 중복은 적격/부적격을 모두 포함한다. "
                  "과거 전체 형상 archive와의 중복은 검사하지 않았으므로 과거 후보와도 겹칠 수 있다.", "",
                  panel_diagnostics_text(panel)]
    lines += ["", "## 수치 설정과 미세 우열", ""]
    numerical = summary["numerical_settings"]
    if not numerical["available"]:
        lines.append(numerical["reason"])
    else:
        lines += [native_replay_intro(numerical), "",
                  "cubic640/cubic1280/quintic640/quintic1280을 같은 기준벽과 짝지어 비교했다. "
                  "설정마다 같은 방향이어도 전체 범위가 겹치면 미세 우열은 미확정으로 둔다. "
                  "이전 관측값의 native 재현 대조를 확인하지 못한 후보 관련 실시간 우열은 확정하지 않는다. "
                  "범위가 분리된 경우에도 이번 유한 수치 설정에서의 분리이며 확률적 신뢰구간·연속체 오차 상한이 아니다.", ""]
        numerical_rows = []
        selected_rows = []
        for case in numerical["cases"]:
            candidates = {r["id"]: r for r in case["candidates"]}
            for choice in case["checkpoint_choices"]:
                if choice["cutoff_s"] != 600:
                    continue
                detail = candidates.get(choice["id"]) if choice["available"] else None
                reductions = ([r["peak_reduction_pct"] for r in detail["settings"].values() if r["paired_eligible"]]
                              if detail else [])
                area_changes = ([r["A2_change_pct"] for r in detail["settings"].values() if r["paired_eligible"]]
                                if detail else [])
                selected_rows.append((case["case_id"], choice["seed"], LABELS[choice["policy"]],
                                      _span(_range(reductions)), _span(_range(area_changes)),
                                      detail["all_settings_eligible"] if detail else "미확인",
                                      detail["A2_decreased_in_any_setting"] if detail else "미확인",
                                      detail["candidate_max_below_baseline_min"] if detail else "미확인"))
            for cutoff in CHECKPOINTS:
                method = [r for r in case["method_checkpoint_pairs"] if r["cutoff_s"] == cutoff]
                neighbor = [r for r in case["nearest_checkpoint_pairs"] if r["cutoff_s"] == cutoff]
                numerical_rows.append((case["case_id"], cutoff,
                                       sum(r.get("resolved_for_live_choices", False) for r in method),
                                       sum(r.get("range_relation") == "overlap_or_touch" for r in method),
                                       sum(r.get("paired_ordering") == "setting_order_reversal" for r in method),
                                       sum(r.get("same_design", False) for r in method),
                                       sum(r.get("resolved_for_live_choices", False) for r in neighbor),
                                       sum(r["available"] for r in neighbor)))
        lines += [_table(["평형", "초", "방법쌍 범위분리", "중첩/접촉", "설정별 순서역전", "동일형상", "근접피크쌍 분리", "근접쌍 있음"], numerical_rows), "",
                  "방법 쌍은 각 평형·시드·checkpoint에서 모두 대조했다. 근접 후보는 같은 평형에서 시간 내 관측한 적격 정밀 후보 중 "
                  "native 피크가 가장 가까운 별개 후보이며, 형상 공간의 모든 이웃을 뜻하지 않는다. 설정별 paired ordering, "
                  "A2 변화·감소와 부적격 상태는 summary.json의 후보 및 쌍별 원값에 보존했다.", "",
                  "600초 선택 형상의 같은 설정 기준벽 대비 변화는 다음과 같다. 범위는 실제 확인한 설정의 최솟값–최댓값이며 "
                  "전체 설정 적격 여부를 별도로 표시한다.", "",
                  _table(["평형", "시드", "방법", "피크 감소 % 범위", "A2 변화 % 범위", "전 설정 적격", "A2 감소 있음", "기준 범위와 분리"], selected_rows)]
    lines += ["", "## 관측된 coarse 적격성 탈락의 별도 사후 진단", "",
              rejections_text(summary.get("screen_rejections", dict(available=False, reason="추가 사후 진단 자료 없음"))), "",
              "## 설정 정정과 남은 범위", "",
              "처음 population=12/projected_seeds=12 부분 실행에서는 성공한 바닥 투영 12개가 초기 집단을 모두 채워 "
              "DE 집단의 0–11번 계수가 고정됐다. 첫 방법 간 비교가 완료되기 전에 중단했으며 원래 생산 기본값 64/12로 정정했다. "
              "열유속 비교 결과로 조정하지 않았고 이전 부분 실행·스냅샷은 별도 보존하며 이번 27개 결과에 섞지 않았다. "
              "정정 기록과 이전 프로토콜/중단 파일의 해시를 대조했다.", "",
              "FIREFLY 별도 교차검사, 연속 first-hit 수송, 1800초 확장, IMP 및 다른 운전점은 이 단계의 완료 주장에 포함하지 않는다. "
              "특정 타일 q(s)의 1% 수렴이나 실험적 정확도도 이 결과로 인증하지 않는다.", "",
              "## 재현 자료", "",
              "큰 원시 event·후보·코드 스냅샷은 D의 원래 study에 둔다. 이 폴더의 summary.json에는 독립 재집계값, "
              "각 checkpoint의 쌍별 비교·벡터 범위·검열 상태와 입력 해시를 보존한다. "
              "anytime_curves.pdf는 같은 그림의 벡터 출력이다. 사후 패널/수치 설정 비용은 실시간 600초 경쟁 밖이며 summary.json에 별도 기록했다.", "",
              f"- 연구 자료: `{summary['study']}`",
              f"- 보고서 생성 UTC: {summary['generated_utc']}",
              f"- 생성 코드 SHA-256: `{summary['source_sha256']}`", ""]
    return "\n".join(lines)


def render_curves(summary, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(3, 3, figsize=(11.8, 9.3), sharex=True, sharey=True, layout="constrained")
    lookup = {(r["case_id"], r["seed"], r["policy"]): r for r in summary["live"]["arms"]}
    for row, case in enumerate(summary["frozen_case_ids"]):
        for col, seed in enumerate(SEEDS):
            ax = axes[row, col]
            any_curve = False
            for policy in POLICIES:
                points = lookup[(case, seed, policy)]["curve"]
                if not points:
                    continue
                any_curve = True
                ax.step([p["t_s"] for p in points], [p["peak_ratio"] for p in points], where="post",
                        color=COLORS[policy], linestyle=STYLES[policy], linewidth=1.6)
            if not any_curve:
                ax.text(.5, .5, "No certified baseline", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"{case}  |  seed {seed}", loc="left", fontsize=9)
            ax.set_xlim(0, 600)
            ax.set_ylim(0, 1.06)
            ax.set_xticks([0, 200, 400, 600])
            ax.grid(axis="y", color="#D6D8DB", linewidth=.5)
            if row == 2:
                ax.set_xlabel("Seconds from process launch")
            if col == 0:
                ax.set_ylabel("Fine peak / original-wall peak")
    handles = [Line2D([0], [0], color=COLORS[p], linestyle=STYLES[p], linewidth=1.8, label=p) for p in POLICIES]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Best completed eligible fine result, including the original wall", fontsize=13)
    for suffix in ("png", "pdf"):
        fig.savefig(out/("anytime_curves."+suffix), dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    study, out = args.study.resolve(), args.out.resolve()
    require(study.drive.lower() == "d:" and out.drive.lower() == "d:", "study and report must remain on D")
    require(out.is_relative_to(CODE_ROOT/"results") and not out.is_relative_to(study),
            "report output must be a separate repository results directory")
    # The main campaign audit is the completion gate; no report is built from
    # partial arms. These controls write only under the chosen report folder.
    summary = make_summary(study)
    out.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(out/".matplotlib")
    render_curves(summary, out)
    for name, value in summary["input_hashes"].items():
        require(digest(study/name) == value, "input evidence changed while rendering: "+name)
    (out/"summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    (out/"REPORT.ko.md").write_text(korean_report(summary), encoding="utf-8")
    print("REPORT", out/"REPORT.ko.md", flush=True)


if __name__ == "__main__":
    main()
