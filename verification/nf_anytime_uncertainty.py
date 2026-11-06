"""Locked post-trial candidates and finite numerical-setting sensitivity.

This reports paired changes and finite-setting ranges within the unchanged
transport model. The ranges are neither continuum error bounds nor confidence
intervals, and they do not certify one-percent convergence of local tile q(s).
"""
import argparse
from itertools import combinations
import math
from pathlib import Path
from time import perf_counter

from verification.nf_anytime_reference import (
    CHECKPOINTS, POLICIES, bind_evidence, configure_environment, cost_summary, evaluate_score,
    finite_positive, freeze_post, load_audited_arms, measurement_record,
    require, start_invocation, validate_candidate, validate_frozen,
    validate_measurement, verify_external_hashes,
)


BASELINE_ID = "original_wall"
NATIVE = "cubic1280"
NATIVE_RTOL = 1e-10
SETTINGS = {
    "cubic640": ("cubic", 640, .00025, 321, .0005),
    "cubic1280": ("cubic", 1280, .000125, 641, .00025),
    "quintic640": ("quintic", 640, .00025, 321, .0005),
    "quintic1280": ("quintic", 1280, .000125, 641, .00025),
}
SCOPE = ("observed on-time eligible fine-scored candidates in this equilibrium only; "
         "nearest native peak does not mean every geometrically neighboring design")
RANGE_SCOPE = ("finite numerical-setting ranges within one fixed transport model; "
               "not probabilistic confidence intervals, continuum error bounds or experimental accuracy")
NATIVE_REPLAY_SCOPE = ("all available genuine native observations: eligibility, positive finite peak/A2 and "
                       "observable presence; no previous native observation means no replay comparison; "
                       "geometry-only eligibility evidence does not demonstrate heat-observable agreement")


def close_native(a, b):
    return abs(a-b) <= NATIVE_RTOL * max(abs(a), abs(b))


def has_native_measurement(event):
    """Separate actual fine outcomes from unexecuted interruption placeholders.

    A complete geometry rejection need not have thermal observables. Conversely,
    an interrupted CPU parity guard can leave an actual GPU measurement without
    an objective. The runner's synthetic Deadline record has neither evidence.
    """
    state = event.get("measurement")
    if event.get("stage") != "fine" or event.get("cached") or event.get("end_s") is None or state is None:
        return False
    values = (state.get("peak"), state.get("A2"), event.get("objective"))
    if (any(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
            for value in values) or state.get("ledger")):
        return True
    if any(str(error).partition(":")[0].strip() == "Deadline" for error in state.get("errors", [])):
        return False
    return event.get("status") in ("completed", "late")


def lock_candidate_set(arms):
    """Pure selection: lock designs before any additional numerical calculations."""
    require(bool(arms), "audited arms are required")
    case_ids = {arm["case_id"] for arm in arms}
    require(len(case_ids) == 1, "candidate selection must stay within one equilibrium")
    arm_keys = [(arm["seed"], arm["policy"]) for arm in arms]
    require(len(set(arm_keys)) == len(arm_keys), "duplicate arm in candidate selection")
    registry, observations, pool_observations, pool = {}, {}, {}, {}
    locked = {BASELINE_ID: dict(id=BASELINE_ID, candidate=None,
                               reasons=[dict(kind="original_wall_reference")], native_observations=[])}
    checkpoint_choices, selected_ids, nearest, unavailable = [], set(), [], []

    def register(identity, candidate):
        validate_candidate(identity, candidate)
        if identity in registry:
            require(registry[identity] == candidate, "candidate identity has inconsistent coefficients")
        registry[identity] = candidate

    def include(identity, reason):
        if identity not in locked:
            locked[identity] = dict(id=identity, candidate=registry[identity], reasons=[], native_observations=[])
        locked[identity]["reasons"].append(reason)

    # Replay every genuine native observation, including late and ineligible
    # outcomes. Nearest-peer admission remains completed, on-time and eligible.
    for arm in arms:
        require(arm["budget_s"] == 600 and arm["policy"] in POLICIES, "unregistered arm budget or policy")
        require(set(arm["checkpoints"]) == {str(t) for t in CHECKPOINTS}, "missing or added audited checkpoints")
        if arm.get("audit_status") == "process_failed":
            unavailable.append(dict(seed=arm["seed"], policy=arm["policy"], reason="process_failed; no certified checkpoint"))
        baseline = arm.get("baseline")
        if baseline is not None:
            validate_measurement(baseline)
            locked[BASELINE_ID]["native_observations"].append(dict(seed=arm["seed"], policy=arm["policy"], measurement=baseline))
        for index, event in enumerate(arm["events"]):
            if event.get("stage") not in ("coarse", "fine"):
                continue
            register(event["id"], event["candidate"])
            end = event.get("end_s")
            if end is None:
                continue  # An unfinished evaluation has no measured outcome.
            require(isinstance(end, (int, float)) and math.isfinite(end) and end >= 0,
                    "invalid evaluation completion clock")
            state = event.get("measurement")
            if state is not None:
                validate_measurement(state)
            identity = event["id"]
            evidence = dict(seed=arm["seed"], policy=arm["policy"], event_index=index,
                            end_s=end, status=event.get("status"), late=end > arm["budget_s"],
                            measurement=state)
            if has_native_measurement(event):
                observations.setdefault(identity, []).append(evidence)
            if (not arm.get("native_pool_eligible", True) or event["stage"] != "fine"
                    or event.get("status") != "completed" or end > arm["budget_s"]
                    or event.get("cached") or state is None or not state["eligible"]):
                continue
            pool_observations.setdefault(identity, []).append(evidence)
            if identity in pool:
                require(close_native(pool[identity]["peak"], state["peak"])
                        and close_native(pool[identity]["A2"], state["A2"]),
                        "conflicting stored native observations for the same candidate")
            else:
                pool[identity] = state

    for arm in arms:
        for cutoff in CHECKPOINTS:
            choice = arm["checkpoints"][str(cutoff)]
            item = dict(seed=arm["seed"], policy=arm["policy"], cutoff_s=cutoff,
                        available=choice["available"], id=None)
            if choice["available"]:
                require(choice["kind"] in ("original_wall", "new_design"), "unknown checkpoint design kind")
                if choice["kind"] == "original_wall":
                    require(choice.get("id") is None and choice.get("candidate") is None,
                            "original-wall checkpoint has new-design coefficients")
                    item["id"] = BASELINE_ID
                else:
                    identity = choice["id"]
                    register(identity, choice["candidate"])
                    evidence = [record for record in pool_observations.get(identity, [])
                                if record["seed"] == arm["seed"] and record["policy"] == arm["policy"]
                                and record["end_s"] <= cutoff]
                    require(evidence and all(finite_positive(choice.get(key)) for key in ("peak", "A2")),
                            "chosen candidate lacks eligible native evidence before its checkpoint")
                    require(any(close_native(record["measurement"]["peak"], choice["peak"])
                                and close_native(record["measurement"]["A2"], choice["A2"]) for record in evidence),
                            "chosen checkpoint observables differ from its native evidence")
                    include(identity, dict(kind="checkpoint", seed=arm["seed"], policy=arm["policy"], cutoff_s=cutoff))
                    selected_ids.add(identity)
                    item["id"] = identity
                item.update(peak=choice["peak"], A2=choice["A2"])
            checkpoint_choices.append(item)
        overruns = [(index, event) for index, event in enumerate(arm["events"])
                    if event.get("stage") in ("coarse", "fine") and event.get("end_s") is not None
                    and event["end_s"] > arm["budget_s"]]
        if overruns:
            index, event = overruns[-1]
            include(event["id"], dict(kind="last_overrun", seed=arm["seed"], policy=arm["policy"],
                                      event_index=index, stage=event["stage"], end_s=event["end_s"],
                                      status=event.get("status")))
    for identity in sorted(selected_ids):
        others = [other for other in pool if other != identity]
        other = min(others, key=lambda oid: (abs(pool[oid]["peak"]-pool[identity]["peak"]), oid)) if others else None
        record = dict(selected_id=identity, neighbor_id=other,
                      native_peak=pool[identity]["peak"],
                      neighbor_native_peak=pool[other]["peak"] if other else None,
                      absolute_native_peak_distance=abs(pool[other]["peak"]-pool[identity]["peak"]) if other else None,
                      scope=SCOPE)
        nearest.append(record)
        if other is not None:
            include(other, dict(kind="nearest_observed_native_peak", selected_id=identity,
                                absolute_native_peak_distance=record["absolute_native_peak_distance"]))
    for identity, row in locked.items():
        if identity != BASELINE_ID:
            row["native_observations"] = observations.get(identity, [])
    return dict(case_id=next(iter(case_ids)), candidates=[locked[BASELINE_ID]] +
                [locked[key] for key in sorted(locked) if key != BASELINE_ID],
                checkpoint_choices=checkpoint_choices, nearest_observed=nearest,
                observed_native_pool_ids=sorted(pool), observed_native_pool_size=len(pool),
                unavailable_arms=unavailable, selection_scope=SCOPE,
                native_observation_scope="genuine fine measurements, including late results and completed geometry ineligibility; no cached requests or unexecuted Deadline placeholders",
                additional_numerical_results_used_for_selection=False)


def missing_measurement(reason="setting_not_measured"):
    return dict(eligible=False, errors=[reason], ledger={}, peak=None, A2=None)


def checked_settings(measurements):
    require(isinstance(measurements, dict) and set(measurements).issubset(SETTINGS), "unknown numerical setting")
    return {name: validate_measurement(measurements.get(name, missing_measurement())) for name in SETTINGS}


def observable_range(measurements, name):
    values = [row[name] for row in measurements.values() if row["eligible"]]
    return [min(values), max(values)] if values else None


def assess_candidate(baseline, candidate):
    """Separate paired improvements from cross-setting range separation."""
    baseline, candidate = checked_settings(baseline), checked_settings(candidate)
    settings = {}
    for name in SETTINGS:
        original, design = baseline[name], candidate[name]
        valid = original["eligible"] and design["eligible"]
        row = dict(baseline=original, candidate=design, paired_eligible=valid,
                   peak_ratio=None, peak_reduction_pct=None, A2_ratio=None,
                   A2_change_pct=None, A2_decreased=None, paired_improvement=False)
        if valid:
            row.update(peak_ratio=design["peak"]/original["peak"],
                       peak_reduction_pct=100*(1-design["peak"]/original["peak"]),
                       A2_ratio=design["A2"]/original["A2"],
                       A2_change_pct=100*(design["A2"]/original["A2"]-1),
                       A2_decreased=design["A2"] < original["A2"],
                       paired_improvement=design["peak"] < original["peak"])
        settings[name] = row
    complete = all(row["paired_eligible"] for row in settings.values())
    base_range = observable_range(baseline, "peak")
    candidate_range = observable_range(candidate, "peak")
    ratios = [row["peak_ratio"] for row in settings.values() if row["paired_eligible"]]
    return dict(settings=settings, all_settings_eligible=complete,
                all_settings_paired_improvement=complete and all(r["paired_improvement"] for r in settings.values()),
                candidate_max_below_baseline_min=bool(complete and candidate_range[1] < base_range[0]),
                baseline_peak_range=base_range, candidate_peak_range=candidate_range,
                baseline_A2_range=observable_range(baseline, "A2"),
                candidate_A2_range=observable_range(candidate, "A2"),
                paired_peak_ratio_range=[min(ratios), max(ratios)] if ratios else None,
                A2_decreased_in_any_setting=any(r["A2_decreased"] is True for r in settings.values()),
                ranges_complete=complete, range_scope=RANGE_SCOPE,
                continuum_error_bound=False, confidence_interval=False,
                local_q_1pct_certified=False, experimental_accuracy_claim=False)


def compare_candidate_ranges(first, second):
    first, second = checked_settings(first), checked_settings(second)
    complete = all(first[name]["eligible"] and second[name]["eligible"] for name in SETTINGS)
    a, b = observable_range(first, "peak"), observable_range(second, "peak")
    overlap = max(a[0], b[0]) <= min(a[1], b[1]) if complete else None
    return dict(first_peak_range=a, second_peak_range=b, all_settings_eligible=complete,
                possible_overlap=overlap,
                first_range_strictly_below_second=bool(complete and a[1] < b[0]),
                second_range_strictly_below_first=bool(complete and b[1] < a[0]),
                fine_ordering_resolved_by_ranges=bool(complete and not overlap),
                range_scope=RANGE_SCOPE)


def compare_native(measured, observations):
    validate_measurement(measured)
    comparisons = []
    for observation in observations:
        previous = validate_measurement(observation["measurement"])
        eligible_same = measured["eligible"] == previous["eligible"]
        presence = {key: dict(stored=finite_positive(previous.get(key)),
                              replayed=finite_positive(measured.get(key))) for key in ("peak", "A2")}
        mismatches = [key for key, value in presence.items() if value["stored"] != value["replayed"]]
        differences = {key: abs(measured[key]-previous[key])/max(abs(measured[key]), abs(previous[key]))
                       for key, value in presence.items() if value["stored"] and value["replayed"]}
        passed = eligible_same and not mismatches and all(value <= NATIVE_RTOL for value in differences.values())
        comparisons.append(dict(source={key: value for key, value in observation.items() if key != "measurement"},
                                stored_measurement=previous, eligibility_matches=eligible_same,
                                observable_presence=presence, observable_presence_matches=not mismatches,
                                observable_presence_mismatches=mismatches,
                                relative_differences=differences, observable_compared_count=len(differences),
                                eligibility_only=not differences and not mismatches, passed=passed))
    matched = all(row["passed"] for row in comparisons)
    return dict(passed=matched, all_available_observations_matched=matched,
                compared_count=len(comparisons),
                previously_scored=bool(comparisons), relative_tolerance=NATIVE_RTOL,
                observable_compared_count=sum(row["observable_compared_count"] for row in comparisons),
                eligibility_only_count=sum(row["eligibility_only"] for row in comparisons),
                absolute_tolerance=0., comparison_scope=NATIVE_REPLAY_SCOPE, comparisons=comparisons)


def method_pairs(selection, measurements):
    pairs = []
    groups = {}
    for choice in selection["checkpoint_choices"]:
        groups.setdefault((choice["seed"], choice["cutoff_s"]), []).append(choice)
    for (seed, cutoff), choices in sorted(groups.items()):
        for first, second in combinations(sorted(choices, key=lambda c: POLICIES.index(c["policy"])), 2):
            available = first["available"] and second["available"]
            result = dict(seed=seed, cutoff_s=cutoff, first_policy=first["policy"],
                          second_policy=second["policy"], first_id=first["id"], second_id=second["id"],
                          both_checkpoint_choices_available=available, same_design=available and first["id"] == second["id"])
            if available:
                result.update(compare_candidate_ranges(measurements[first["id"]], measurements[second["id"]]))
            else:
                result.update(possible_overlap=None, fine_ordering_resolved_by_ranges=False)
            pairs.append(result)
    return pairs


def assess_case(selection, measurements):
    identities = {row["id"] for row in selection["candidates"]}
    require(set(measurements) == identities and BASELINE_ID in measurements, "measurements differ from locked candidate set")
    baseline = measurements[BASELINE_ID]
    candidates = []
    for row in selection["candidates"]:
        identity = row["id"]
        candidate = dict(id=identity, candidate=row["candidate"], reasons=row["reasons"],
                         native_reproducibility=compare_native(measurements[identity][NATIVE], row["native_observations"]))
        candidate.update(assess_candidate(baseline, measurements[identity]))
        candidates.append(candidate)
    nearest_pairs = []
    for selected in selection["nearest_observed"]:
        pair = dict(selected)
        if selected["neighbor_id"] is None:
            pair.update(possible_overlap=None, fine_ordering_resolved_by_ranges=False,
                        diagnostic_available=False)
        else:
            pair.update(compare_candidate_ranges(measurements[selected["selected_id"]],
                                                 measurements[selected["neighbor_id"]]), diagnostic_available=True)
        nearest_pairs.append(pair)
    replay = [row["native_reproducibility"] for row in candidates]
    matched = all(row["all_available_observations_matched"] for row in replay)
    coverage = dict(locked_candidate_count=len(replay),
                    previously_scored_candidate_count=sum(row["previously_scored"] for row in replay),
                    unobserved_candidate_count=sum(not row["previously_scored"] for row in replay),
                    compared_observation_count=sum(row["compared_count"] for row in replay),
                    compared_observable_count=sum(row["observable_compared_count"] for row in replay),
                    eligibility_only_observation_count=sum(row["eligibility_only_count"] for row in replay))
    return dict(case_id=selection["case_id"], candidates=candidates,
                native_reproducibility_passed=matched, all_available_native_observations_matched=matched,
                native_replay_coverage=coverage, native_replay_scope=NATIVE_REPLAY_SCOPE,
                method_choice_pairs=method_pairs(selection, measurements),
                nearest_observed_pairs=nearest_pairs,
                range_scope=RANGE_SCOPE, candidate_selection_scope=SCOPE,
                local_q_1pct_certified=False, global_optimum_claim=False)


def uncertainty_details():
    return dict(settings={key: list(value) for key, value in SETTINGS.items()},
                native_replay="every locked candidate and the original wall are re-evaluated at native cubic1280",
                native_comparison="all stored genuine fine measurements for locked candidates, including late and ineligible outcomes; compare each positive finite observable regardless of eligibility and reject observable-presence changes; exclude unexecuted Deadline placeholders; nearest-peer pool remains on-time eligible only",
                native_replay_scope=NATIVE_REPLAY_SCOPE,
                native_relative_tolerance=NATIVE_RTOL, native_absolute_tolerance=0.,
                selection="all audited 200/400/600 checkpoint new designs, each arm's last overrun coarse or fine design, original wall, and closest distinct observed on-time native-peak peer per chosen design; tie ID",
                selected_set_locked_before_calculation=True, nearest_scope=SCOPE,
                decisions="all-settings paired improvement, candidate max below baseline min, A2 changes, ranges and chosen-method-pair overlap",
                range_scope=RANGE_SCOPE, local_q_1pct_certified=False, tuning=False)


def lock_study_selection(study, frozen, store):
    arms, hashes = load_audited_arms(study, frozen)
    bind_evidence(store, "live-evidence.json", hashes)
    selection = dict(schema="nf-anytime-uncertainty-selection.v1", protocol_sha256=store.sha256("protocol.json"),
                     live_evidence_hashes=hashes,
                     cases=[lock_candidate_set([arm for arm in arms if arm["case_id"] == case["case_id"]])
                            for case in frozen["cases"]],
                     additional_numerical_results_used_for_selection=False)
    bind_evidence(store, "selection.json", selection)
    return selection, hashes


def run_setting(store, case, selected, setting, reference, invocation):
    """One shared preparation per case/setting; individual outcomes are immutable."""
    from d3d_optimization.boundary_evaluator import prepare_boundary
    from d3d_optimization.campaign_geometry import Candidate, evaluator
    from verification.nf_base_budget import measurement
    prefix = case["case_id"] + "/scores/"
    protocol_hash, selection_hash = store.sha256("protocol.json"), store.sha256("selection.json")

    def name(identity):
        return prefix + identity + "." + setting + ".json"

    def provenance(row):
        return dict(id=row["id"], candidate=row["candidate"], case_id=case["case_id"], setting=setting,
                    setting_tuple=list(SETTINGS[setting]), protocol_sha256=protocol_hash, selection_sha256=selection_hash)

    def load(row):
        record = store.read(name(row["id"]))
        require(all(record[key] == value for key, value in provenance(row).items()), "reused uncertainty score provenance differs")
        validate_measurement(record["measurement"])
        return record

    rows = selected["candidates"]
    missing = [row for row in rows if not store.exists(name(row["id"]))]
    if missing:
        base_row = next(row for row in rows if row["id"] == BASELINE_ID)
        previous_base = load(base_row) if store.exists(name(BASELINE_ID)) else None
        prepared_at = perf_counter()
        engine, prep_failure = None, None
        if previous_base and previous_base.get("preparation_failed"):
            # A failed request is not replaced by a later successful attempt.
            prep_failure = previous_base["error"]
            baseline = previous_base["measurement"]
        else:
            try:
                method, count, step, nodes, wall_step = SETTINGS[setting]
                original = prepare_boundary(Path(case["root"]), "OMP", reference,
                                            method=method, count=count, step=step, nodes=nodes, wall_step=wall_step)
                baseline = measurement_record(measurement(original.baseline, reference))
                engine = evaluator(original, 20, "gpu")
                del original
            except Exception as exc:
                prep_failure = "preparation_failed: " + type(exc).__name__ + ": " + str(exc)
                baseline = missing_measurement(prep_failure)
            preparation = dict(case_id=case["case_id"], setting=setting, setting_tuple=list(SETTINGS[setting]),
                               protocol_sha256=protocol_hash, selection_sha256=selection_hash,
                               status="failed" if prep_failure else "completed", error=prep_failure,
                               baseline=baseline, seconds=perf_counter()-prepared_at)
            store.write("preparation/" + case["case_id"] + "/" + invocation.split("/")[-1] + "/" + setting + ".json", preparation)
        base_score = dict(provenance(base_row), measurement=baseline,
                          status="failed" if prep_failure else "completed", error=prep_failure,
                          preparation_failed=prep_failure is not None, seconds=0.,
                          cost_note="original-wall transport is included in the shared preparation cost")
        if previous_base:
            require(compare_native(baseline, [dict(measurement=previous_base["measurement"])])["passed"],
                    "resumed original-wall preparation differs from sealed measurement")
        else:
            store.write(name(BASELINE_ID), base_score)
        for row in rows:
            if row["id"] == BASELINE_ID or store.exists(name(row["id"])):
                continue
            if prep_failure:
                score = dict(measurement=missing_measurement(prep_failure),
                             status="not_evaluated_preparation_failed", error=prep_failure, seconds=0.)
            else:
                score = evaluate_score(engine, Candidate.from_record(row["candidate"]))
            store.write(name(row["id"]), dict(provenance(row), **score))
        del engine
    return ({row["id"]: load(row)["measurement"] for row in rows},
            {name(row["id"]): store.sha256(name(row["id"])) for row in rows})


def prepare_protocol(study):
    from verification import nf_anytime_reference
    return freeze_post(study, "uncertainty", uncertainty_details(),
                       [Path(__file__), Path(nf_anytime_reference.__file__)])


def run_uncertainty(study, *, select_only=False):
    configure_environment(study, "uncertainty")
    from dataclasses import asdict
    from d3d_optimization.campaign_state import keep_awake, process_lock
    from d3d_optimization.reference_transport import ReferencePoint
    from verification.nf_anytime_campaign import versions
    start = perf_counter()
    with process_lock(Path(study) / "campaign.lock"), process_lock(Path(study) / "posttrial.lock"), keep_awake():
        frozen, store = prepare_protocol(study)
        require(versions() == frozen["packages"], "runtime package versions changed")
        selection, hashes = lock_study_selection(study, frozen, store)
        if select_only:
            return selection
        if store.exists("summary.json"):
            return store.read("summary.json")
        invocation = start_invocation(store)
        reference = ReferencePoint(**frozen["reference"])
        require(asdict(reference) == frozen["reference"], "reference settings changed")
        try:
            for case, selected in zip(frozen["cases"], selection["cases"]):
                require(case["case_id"] == selected["case_id"], "locked equilibrium order changed")
                quality_name = case["case_id"] + "/quality.json"
                if store.exists(quality_name):
                    continue
                validate_frozen(study)
                measured = {row["id"]: {} for row in selected["candidates"]}
                evidence = {}
                for setting in SETTINGS:
                    values, recorded = run_setting(store, case, selected, setting, reference, invocation)
                    for identity, state in values.items():
                        measured[identity][setting] = state
                    evidence.update(recorded)
                    print("UNCERTAINTY", case["case_id"], setting, len(values), "designs", flush=True)
                quality = dict(schema="nf-anytime-uncertainty-quality.v1", protocol_sha256=store.sha256("protocol.json"),
                               selection_sha256=store.sha256("selection.json"), score_hashes=evidence,
                               **assess_case(selected, measured))
                store.write(quality_name, quality)
            validate_frozen(study)
            verify_external_hashes(study, hashes)
            store.write(invocation + ".finished.json", dict(status="completed", seconds=perf_counter()-start))
        except Exception as exc:
            if not store.exists(invocation + ".finished.json"):
                store.write(invocation + ".finished.json", dict(status="failed", seconds=perf_counter()-start,
                                                               error=type(exc).__name__ + ": " + str(exc)))
            raise
        cases = []
        for case in frozen["cases"]:
            name = case["case_id"] + "/quality.json"
            quality = store.read(name)
            cases.append(dict(case_id=case["case_id"], quality_sha256=store.sha256(name),
                              candidate_count=len(quality["candidates"]),
                              native_reproducibility_passed=quality["native_reproducibility_passed"],
                              all_available_native_observations_matched=quality["all_available_native_observations_matched"],
                              native_replay_coverage=quality["native_replay_coverage"],
                              all_settings_eligible_count=sum(r["all_settings_eligible"] for r in quality["candidates"]),
                              paired_improvement_count=sum(r["all_settings_paired_improvement"] for r in quality["candidates"]),
                              separated_from_baseline_count=sum(r["candidate_max_below_baseline_min"] for r in quality["candidates"]),
                              method_choice_pairs=quality["method_choice_pairs"]))
        coverage = {key: sum(case["native_replay_coverage"][key] for case in cases)
                    for key in cases[0]["native_replay_coverage"]}
        matched = all(case["all_available_native_observations_matched"] for case in cases)
        summary = dict(schema="nf-anytime-uncertainty-summary.v1", protocol_sha256=store.sha256("protocol.json"),
                       selection_sha256=store.sha256("selection.json"), cases=cases, costs=cost_summary(store),
                       native_reproducibility_passed=matched, all_available_native_observations_matched=matched,
                       native_replay_coverage=coverage, native_replay_scope=NATIVE_REPLAY_SCOPE,
                       range_scope=RANGE_SCOPE, local_q_1pct_certified=False, global_optimum_claim=False)
        store.write("summary.json", summary)
        store.verify_all()
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--freeze-only", action="store_true", help="seal code and protocol without any physical calculation")
    parser.add_argument("--select-only", action="store_true", help="seal all audited candidates without additional physical calculation")
    args = parser.parse_args()
    require(not (args.freeze_only and args.select_only), "choose only one preparation action")
    study = args.study.resolve()
    require(study.drive.lower() == "d:", "study must remain on D as requested")
    if args.freeze_only:
        from d3d_optimization.campaign_state import process_lock
        with process_lock(study / "posttrial.lock"):
            prepare_protocol(study)
        print("Frozen", study / "uncertainty" / "protocol.json", flush=True)
    else:
        result = run_uncertainty(study, select_only=args.select_only)
        print("UNCERTAINTY", "SELECTION LOCKED" if args.select_only else "COMPLETE", len(result["cases"]), "equilibria", flush=True)


if __name__ == "__main__":
    main()
