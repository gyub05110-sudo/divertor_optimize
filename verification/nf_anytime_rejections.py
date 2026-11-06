"""Targeted post-hoc native audit of observed coarse thermal rejections.

This additional diagnostic is separate from the preregistered experiment. It
does not change its winners, selection rules, or wall-time comparisons. Every
fresh, completed, on-time coarse thermal return that was jointly ineligible is
included, without filtering by heat value, error type, or prior fine outcome.
"""
import argparse
from copy import deepcopy
from dataclasses import asdict
import math
from pathlib import Path
from time import perf_counter

from verification import nf_anytime_reference as evidence
from verification import nf_anytime_uncertainty as numerical
from verification.nf_anytime_uncertainty import BASELINE_ID, NATIVE, SETTINGS, run_setting


KIND = "screen-rejections"
SCOPE = ("targeted post-hoc audit of observed coarse thermal eligibility rejections, evaluated only at "
         "native cubic1280; not a preregistered endpoint, a global false-negative rate, or an equal-wall-time "
         "comparison; original primary winners are unchanged")


def _qualifies(event, budget):
    if (event.get("stage") != "coarse" or event.get("status") != "completed"
            or event.get("cached") or event.get("physics_evaluated") is not True):
        return False
    end = event.get("end_s")
    evidence.require(isinstance(end, (int, float)) and not isinstance(end, bool)
                     and math.isfinite(end) and end >= 0, "invalid coarse completion clock")
    if end > budget:
        return False
    state = evidence.validate_measurement(event.get("measurement"))
    return not state["eligible"]


def select_rejections(arms, evidence_hashes=None):
    """Derive the entire target set from audited events, never the review memo."""
    observations, targets, uncertified = {}, {}, []
    keys = [(arm["case_id"], arm["seed"], arm["policy"]) for arm in arms]
    evidence.require(len(set(keys)) == len(keys), "duplicate live arm")
    for arm in arms:
        evidence.require(arm["budget_s"] == 600, "registered live deadline changed")
        certified = arm.get("native_pool_eligible", True)
        if not certified:
            uncertified.append(dict(case_id=arm["case_id"], seed=arm["seed"], policy=arm["policy"]))
        prefix = "live/" + arm["case_id"] + "/" + str(arm["seed"]) + "/" + arm["policy"]
        for index, event in enumerate(arm["events"]):
            if event.get("stage") != "coarse":
                continue
            evidence.require(arm["policy"] == "coarse_adaptive", "coarse event outside the coarse policy")
            evidence.validate_candidate(event["id"], event["candidate"])
            key = (arm["case_id"], event["id"])
            observation = dict(seed=arm["seed"], policy=arm["policy"], event_index=index,
                               event=deepcopy(event), event_stream=prefix + "/events.jsonl",
                               result_file=prefix + "/result.json", certified=certified)
            if evidence_hashes is not None:
                observation.update(event_stream_sha256=evidence_hashes.get(observation["event_stream"]),
                                   result_sha256=evidence_hashes.get(observation["result_file"]))
            observations.setdefault(key, []).append(observation)
            if not certified or not _qualifies(event, arm["budget_s"]):
                continue
            if evidence_hashes is not None:
                evidence.require(observation["event_stream_sha256"] and observation["result_sha256"],
                                 "qualifying observation lacks raw live evidence hashes")
            observation["incumbent_600"] = deepcopy(arm["checkpoints"]["600"])
            if key not in targets:
                targets[key] = dict(id=event["id"], candidate=deepcopy(event["candidate"]),
                                    qualifying_observations=[])
            targets[key]["qualifying_observations"].append(observation)
    cases = []
    for case in sorted({key[0] for key in targets}):
        candidates = [dict(id=BASELINE_ID, candidate=None, qualifying_observations=[], coarse_observations=[])]
        for key in sorted(key for key in targets if key[0] == case):
            candidates.append(dict(targets[key], coarse_observations=observations[key]))
        cases.append(dict(case_id=case, candidates=candidates))
    return dict(cases=cases, examined_case_ids=sorted({arm["case_id"] for arm in arms}), examined_arm_count=len(arms),
                selected_case_count=len(cases), selected_candidate_count=len(targets),
                qualifying_observation_count=sum(len(row["qualifying_observations"]) for row in targets.values()),
                uncertified_arms=uncertified, complete_live_coverage=not uncertified,
                post_hoc=True, scope=SCOPE, primary_winners_changed=False)


def _comparison(scored, reference, available):
    state = scored["measurement"]
    peak = available and evidence.finite_positive(state.get("peak")) and evidence.finite_positive(reference.get("peak"))
    area = available and evidence.finite_positive(state.get("A2")) and evidence.finite_positive(reference.get("A2"))
    jointly_eligible = available and reference.get("eligible", True) and scored["status"] == "completed" and state["eligible"]
    return dict(reference=reference, reference_available=available,
                peak_ratio=state["peak"]/reference["peak"] if peak else None,
                peak_reduction_pct=100*(1-state["peak"]/reference["peak"]) if peak else None,
                A2_change_pct=100*(state["A2"]/reference["A2"]-1) if area else None,
                strictly_lower_observed_peak=state["peak"] < reference["peak"] if peak else None,
                eligible_peak_improvement=bool(jointly_eligible and peak and state["peak"] < reference["peak"]),
                diagnostic_only=True)


def assess_rejections(selected, scores):
    identities = {row["id"] for row in selected["candidates"]}
    evidence.require(set(scores) == identities and BASELINE_ID in scores, "native scores differ from locked rejection set")
    for scored in scores.values():
        evidence.validate_measurement(scored["measurement"])
        evidence.require(scored["status"] in ("completed", "failed", "not_evaluated_preparation_failed"),
                         "unknown native score status")
    baseline = scores[BASELINE_ID]
    rows = []
    for selected_row in selected["candidates"]:
        if selected_row["id"] == BASELINE_ID:
            continue
        scored = scores[selected_row["id"]]
        eligible = scored["measurement"]["eligible"] if scored["status"] == "completed" else None
        comparisons, seen = [], set()
        for observation in selected_row["qualifying_observations"]:
            arm = (observation["seed"], observation["policy"])
            if arm in seen:
                continue
            seen.add(arm)
            incumbent = observation["incumbent_600"]
            comparisons.append(dict(seed=arm[0], policy=arm[1],
                                    **_comparison(scored, incumbent, incumbent["available"])))
        rows.append(dict(selected_row, fine=scored, fine_eligible=eligible,
                         coarse_eligibility_false_negative=eligible,
                         diagnostic_status="fine_evaluation_failed" if eligible is None else
                                           "fine_eligible_false_negative" if eligible else "fine_ineligible",
                         vs_baseline=_comparison(scored, baseline["measurement"],
                                                 baseline["status"] == "completed" and baseline["measurement"]["eligible"]),
                         vs_live_coarse_incumbents=comparisons))
    return dict(case_id=selected["case_id"], baseline=baseline, rows=rows,
                false_negative_count=sum(row["coarse_eligibility_false_negative"] is True for row in rows),
                fine_ineligible_count=sum(row["coarse_eligibility_false_negative"] is False for row in rows),
                fine_evaluation_failed_count=sum(row["coarse_eligibility_false_negative"] is None for row in rows),
                post_hoc=True, scope=SCOPE, primary_winners_changed=False)


def source_paths():
    return (Path(__file__).resolve(), Path(evidence.__file__).resolve(), Path(numerical.__file__).resolve())


def protocol_details():
    runner_source = Path(run_setting.__code__.co_filename).resolve()
    evidence.require(runner_source == Path(numerical.__file__).resolve(), "unexpected imported native runner")
    return dict(post_hoc=True, scope=SCOPE, primary_winners_changed=False,
                selection="all certified fresh completed on-time coarse events with physics_evaluated=true and eligible=false; case-deduplicated; no peak/error-type/fine-outcome filter",
                setting=list(SETTINGS[NATIVE]), original_wall_included=True, radius_mm=20,
                imported_runner=dict(qualified_name=run_setting.__module__ + "." + run_setting.__name__,
                                     source_file=runner_source.relative_to(evidence.CODE_ROOT).as_posix(),
                                     source_sha256=evidence.digest(runner_source)),
                membership_locked_before_additional_native_calculations=True,
                diagnostics="fine eligibility/false-negative status and paired peak/A2 changes versus original wall and corresponding frozen 600-second coarse-adaptive incumbent",
                same_wall_time_comparison=False, continuum_error_bound=False)


def prepare_protocol(study):
    return evidence.freeze_post(study, KIND, protocol_details(), source_paths())


def lock_selection(study, frozen, store):
    arms, hashes = evidence.load_audited_arms(study, frozen)
    evidence.bind_evidence(store, "live-evidence.json", hashes)
    selection = dict(schema="nf-anytime-screen-rejections-selection.v1",
                     protocol_sha256=store.sha256("protocol.json"), live_evidence_hashes=hashes,
                     **select_rejections(arms, hashes))
    evidence.bind_evidence(store, "selection.json", selection)
    return selection, hashes


def run_selected_case(store, case, selected, reference, invocation, *, runner=None):
    """Score only the sealed set at native resolution; preserve every failed request."""
    locked = store.read("selection.json")
    evidence.require(next(row for row in locked["cases"] if row["case_id"] == case["case_id"]) == selected,
                     "rejection selection must be sealed before native scoring")
    prefix, result_name = case["case_id"] + "/scores/", case["case_id"] + "/result.json"
    metadata = dict(protocol_sha256=store.sha256("protocol.json"), selection_sha256=store.sha256("selection.json"))
    previous = store.read(result_name) if store.exists(result_name) else None
    returned, returned_hashes = None, None
    if previous is None:
        returned, returned_hashes = (runner or run_setting)(store, case, selected, NATIVE, reference, invocation)
    scores, hashes = {}, {}
    for row in selected["candidates"]:
        name = prefix + row["id"] + "." + NATIVE + ".json"
        scored = store.read(name)
        identity = dict(id=row["id"], candidate=row["candidate"], case_id=case["case_id"], setting=NATIVE, **metadata)
        evidence.require(all(scored[key] == value for key, value in identity.items()), "native rejection score provenance differs")
        scores[row["id"]], hashes[name] = scored, store.sha256(name)
    if returned is not None:
        evidence.require(returned_hashes == hashes and returned == {key: value["measurement"] for key, value in scores.items()},
                         "native runner return differs from sealed score evidence")
    result = dict(schema="nf-anytime-screen-rejections-case.v1", score_hashes=hashes, **metadata,
                  **assess_rejections(selected, scores))
    if previous is None:
        store.write(result_name, result)
    else:
        evidence.require(previous == result, "recomputed rejection diagnostic differs from sealed result")
    return result


def run_audit(study, *, select_only=False):
    evidence.configure_environment(study, KIND)
    from d3d_optimization.campaign_state import keep_awake, process_lock
    from d3d_optimization.reference_transport import ReferencePoint
    from verification.nf_anytime_campaign import versions
    started = perf_counter()
    with process_lock(Path(study)/"campaign.lock"), process_lock(Path(study)/"posttrial.lock"), keep_awake():
        frozen, store = prepare_protocol(study)
        evidence.require(versions() == frozen["packages"], "runtime package versions changed")
        selection, live_hashes = lock_selection(study, frozen, store)
        if select_only:
            return selection
        if store.exists("summary.json"):
            return store.read("summary.json")
        invocation = evidence.start_invocation(store)
        reference = ReferencePoint(**frozen["reference"])
        evidence.require(asdict(reference) == frozen["reference"], "original operating inputs changed")
        cases = {row["case_id"]: row for row in frozen["cases"]}
        summaries = []
        try:
            for selected in selection["cases"]:
                evidence.validate_frozen(study)
                result = run_selected_case(store, cases[selected["case_id"]], selected, reference, invocation)
                summaries.append(dict(case_id=selected["case_id"],
                                      result_sha256=store.sha256(selected["case_id"] + "/result.json"),
                                      **{key: result[key] for key in ("false_negative_count", "fine_ineligible_count", "fine_evaluation_failed_count")}))
                print("SCREEN REJECTIONS", selected["case_id"], result["false_negative_count"], "fine-eligible", flush=True)
            evidence.validate_frozen(study)
            evidence.verify_external_hashes(study, live_hashes)
            store.write(invocation + ".finished.json", dict(status="completed", seconds=perf_counter()-started))
        except Exception as exc:
            if not store.exists(invocation + ".finished.json"):
                store.write(invocation + ".finished.json", dict(status="failed", seconds=perf_counter()-started,
                                                               error=type(exc).__name__ + ": " + str(exc)))
            raise
        summary = dict(schema="nf-anytime-screen-rejections-summary.v1", protocol_sha256=store.sha256("protocol.json"),
                       selection_sha256=store.sha256("selection.json"), cases=summaries,
                       selected_candidate_count=selection["selected_candidate_count"],
                       qualifying_observation_count=selection["qualifying_observation_count"],
                       complete_live_coverage=selection["complete_live_coverage"],
                       **{key: sum(row[key] for row in summaries) for key in
                          ("false_negative_count", "fine_ineligible_count", "fine_evaluation_failed_count")},
                       costs=evidence.cost_summary(store), post_hoc=True, scope=SCOPE, primary_winners_changed=False)
        store.write("summary.json", summary)
        store.verify_all()
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--freeze-only", action="store_true")
    options.add_argument("--select-only", action="store_true")
    args = parser.parse_args()
    study = args.study.resolve()
    evidence.require(study.drive.lower() == "d:", "study must remain on D")
    if args.freeze_only:
        from d3d_optimization.campaign_state import process_lock
        with process_lock(study/"campaign.lock"), process_lock(study/"posttrial.lock"):
            prepare_protocol(study)
        print("SCREEN REJECTIONS PROTOCOL FROZEN", flush=True)
    else:
        result = run_audit(study, select_only=args.select_only)
        print("SCREEN REJECTIONS", "SELECTION LOCKED" if args.select_only else "COMPLETE",
              result["selected_candidate_count"], "candidates", flush=True)


if __name__ == "__main__":
    main()
