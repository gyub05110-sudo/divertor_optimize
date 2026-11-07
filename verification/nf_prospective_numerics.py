"""Seal all primary and long-budget designs before four-setting transport replay."""
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
import platform
from time import perf_counter

from verification import nf_prospective_report as report
from verification.nf_anytime_reference import (
    EvidenceStore, bind_evidence, configure_environment, cost_summary, evaluate_score, measurement_record,
    require, start_invocation, validate_candidate, validate_measurement, verify_external_hashes,
)
from verification.nf_anytime_uncertainty import (
    NATIVE, NATIVE_RTOL, RANGE_SCOPE, SETTINGS, assess_candidate, compare_native, has_native_measurement,
    missing_measurement,
)


def lock_case(case_id, arms):
    candidates = {report.BASE: dict(id=report.BASE, candidate=None, reasons=[dict(kind="original_wall_reference")],
                                    native_observations=[])}
    choices = []
    for arm in arms:
        require(arm["case_id"] == case_id, "cannot pool candidates across equilibria")
        for cutoff in report.cutoffs(arm["budget_s"]):
            if cutoff not in (600, report.LONG):
                continue
            e = report.endpoint(arm, cutoff)
            choice = dict(seed=arm["seed"], policy=arm["policy"], cutoff_s=cutoff,
                          available=e["available"], id=e["selected_id"], status=e["status"], peak=e["peak"], A2=e["A2"])
            choices.append(choice)
            if not e["available"]:
                continue
            identity, cand = e["selected_id"], e["candidate"]
            if identity != report.BASE:
                validate_candidate(identity, cand)
            row = candidates.setdefault(identity, dict(id=identity, candidate=cand, reasons=[], native_observations=[]))
            require(row["candidate"] == cand, "inconsistent coefficients in selected geometry registry")
            row["reasons"].append(dict(kind="checkpoint", **choice))
    for arm in arms:
        if not arm["raw_certified"]:
            continue
        for index, e in enumerate(arm["events"]):
            if e["stage"] == "baseline" and e.get("measurement") is not None:
                identity = report.BASE
            elif e.get("id") in candidates and has_native_measurement(e):
                identity = e["id"]
                require(candidates[identity]["candidate"] == e["candidate"], "native event coefficients changed")
            else:
                continue
            state = validate_measurement(e["measurement"])
            candidates[identity]["native_observations"].append(dict(seed=arm["seed"], policy=arm["policy"],
                event_index=index, stage=e["stage"], status=e["status"], end_s=e["end_s"], measurement=state))
    for choice in choices:
        if choice["available"]:
            require(any(x["seed"] == choice["seed"] and x["policy"] == choice["policy"] and
                        x["status"] == "completed" and x["end_s"] <= choice["cutoff_s"] and x["measurement"]["eligible"] and
                        (x["measurement"]["peak"],x["measurement"]["A2"]) == (choice["peak"],choice["A2"])
                        for x in candidates[choice["id"]]["native_observations"]), "selection lacks its own on-time native measurement")
    return dict(case_id=case_id, candidates=list(candidates.values()), choices=choices)


def compare_designs(fast, direct, same_design=False):
    """Compute every measured finite pair, and qualify claims with joint eligibility."""
    rows, gaps, signs, fast_peaks, direct_peaks, finite_gaps = {}, [], [], [], [], []
    by_threshold = defaultdict(list)
    for setting in SETTINGS:
        a, b = (validate_measurement(v.get(setting, missing_measurement())) for v in (fast,direct))
        finite = report.audit.positive(a.get("peak")) and report.audit.positive(b.get("peak"))
        finite_a2 = report.audit.positive(a.get("A2")) and report.audit.positive(b.get("A2"))
        gap = report.margin(a["peak"],b["peak"]) if finite else None
        eligible = bool(a["eligible"] and b["eligible"])
        rows[setting] = dict(fast=a, direct=b, both_eligible=eligible, finite_margin=gap,
                             finite_A2_fast_over_direct=a["A2"]/b["A2"] if finite_a2 else None,
                             eligible_margin=gap if eligible else None)
        if finite:
            finite_gaps.append(gap["gap_pct"])
        if eligible:
            gaps.append(gap["gap_pct"])
            signs.append(gap["sign"])
            fast_peaks.append(a["peak"])
            direct_peaks.append(b["peak"])
            for key, classification in gap["practical"].items():
                by_threshold[key].append(classification)
    complete = len(gaps) == len(SETTINGS)
    relation, cross_range = "incomplete", None
    if complete:
        if max(fast_peaks) < min(direct_peaks):
            relation = "fast_below_direct"
        elif max(direct_peaks) < min(fast_peaks):
            relation = "direct_below_fast"
        else:
            relation = "overlap_or_touch"
        cross_range = [report.margin(max(fast_peaks),min(direct_peaks))["gap_pct"],
                       report.margin(min(fast_peaks),max(direct_peaks))["gap_pct"]]
    practical = {}
    for threshold in report.THRESHOLDS:
        key = f"{threshold:g}"
        kinds = by_threshold[key]
        practical[key] = "incomplete" if not complete else kinds[0] if len(set(kinds)) == 1 else "mixed"
    return dict(settings=rows, same_design=same_design, finite_setting_count=len(finite_gaps),
                eligible_setting_count=len(gaps), finite_gap_range_pct=report.span(finite_gaps),
                eligible_gap_range_pct=report.span(gaps), all_settings_eligible=complete,
                all_settings_sign=("incomplete" if not complete else signs[0] if len(set(signs)) == 1 else "mixed"),
                all_settings_practical=practical, range_relation=relation,
                fast_peak_range=report.span(fast_peaks), direct_peak_range=report.span(direct_peaks),
                conservative_cross_setting_gap_range_pct=cross_range)


def native_replay(score, observations):
    actual = compare_native(score["measurement"], observations)
    observed = actual["compared_count"] > 0
    completed = score["status"] == "completed"
    thermally_observed = any(report.audit.positive(o["measurement"].get(k)) for o in observations for k in ("peak","A2"))
    matched = bool(completed and actual["passed"]) if observed else None
    actual.update(passed=matched, all_available_observations_matched=matched,
                  actual_score_status=score["status"], genuine_prior_observations=observed,
                  observables_reproduced=bool(matched and actual["observable_compared_count"]) if thermally_observed else None)
    return actual


def assess_case(selection, scores):
    expected = {c["id"] for c in selection["candidates"]}
    require(set(scores) == expected and all(set(s) == set(SETTINGS) for s in scores.values()), "scores differ from locked case/setting registry")
    measures = {i:{s:r["measurement"] for s,r in vals.items()} for i,vals in scores.items()}
    base = measures[report.BASE]
    candidates = []
    for c in selection["candidates"]:
        identity = c["id"]
        comparison = assess_candidate(base, measures[identity])
        eligible_pairs = [s for s in comparison["settings"].values() if s["paired_eligible"]]
        comparison["A2_decreased_settings"] = sum(s["A2_decreased"] for s in eligible_pairs)
        comparison["joint_improvement_settings"] = sum(s["peak_ratio"] < 1 and s["A2_ratio"] > 1 for s in eligible_pairs)
        comparison["finite_pairs_including_ineligible"] = {
            setting:dict(peak_ratio=(m["peak"]/base[setting]["peak"] if report.audit.positive(m.get("peak")) and
                        report.audit.positive(base[setting].get("peak")) else None),
                        A2_ratio=(m["A2"]/base[setting]["A2"] if report.audit.positive(m.get("A2")) and
                        report.audit.positive(base[setting].get("A2")) else None))
            for setting,m in measures[identity].items()}
        candidates.append(dict(id=identity, candidate=c["candidate"], reasons=c["reasons"],
            settings={s:{k:r.get(k) for k in ("status","seconds","error","measurement")} for s,r in scores[identity].items()},
            native_replay=native_replay(scores[identity][NATIVE],c["native_observations"]), baseline_comparison=comparison))
    paired_choices = defaultdict(dict)
    for c in selection["choices"]:
        group = paired_choices[(c["seed"],c["cutoff_s"])]
        require(c["policy"] not in group, "duplicate locked method choice")
        group[c["policy"]] = c
    pairs = []
    for (seed,cutoff), group in paired_choices.items():
        require(set(group) == set(report.POLICIES), "missing registered method in numerical pair")
        a, b = (group[p] for p in report.POLICIES)
        available = a["available"] and b["available"]
        a_settings = measures[a["id"]] if a["available"] else {}
        b_settings = measures[b["id"]] if b["available"] else {}
        pairs.append(dict(case_id=selection["case_id"], seed=seed, cutoff_s=cutoff, fast_id=a["id"], direct_id=b["id"],
                          both_choices_available=available, native_fast=a["peak"], native_direct=b["peak"],
                          **compare_designs(a_settings,b_settings, same_design=bool(available and a["id"] == b["id"]))))
    return dict(case_id=selection["case_id"], candidates=candidates, pairs=pairs)


def require_locked(study):
    require((Path(study)/"numerics"/"selection.json").is_file(), "run --lock-only to seal the full study selection before --score")


def protocol_and_selection(study, frozen, arms, hashes, store):
    from verification.nf_anytime_campaign import local_source_closure, versions
    files = local_source_closure([Path(__file__).resolve(), Path(report.__file__).resolve()])
    sources = {p.relative_to(report.CODE_ROOT).as_posix():report.audit.digest(p) for p in files}
    protocol = dict(schema="nf-prospective-setting-protocol.v1", registration_sha256=hashes["frozen.json"],
        study_protocol_sha256=hashes["protocol.ko.md"], development_choice_sha256=hashes["development-choice.json"],
        settings={k:list(v) for k,v in SETTINGS.items()}, source="OMP", radius_mm=20, reference=frozen["reference"],
        backend="gpu", native=NATIVE, native_relative_tolerance=NATIVE_RTOL, native_absolute_tolerance=0.,
        selection="original wall and both validated method choices at 600 s for each seed and 1800 s for the first seed per equilibrium",
        denominator="G = 100*(q_direct-q_fast)/q_direct; both signs use the direct comparator",
        thresholds_pct=list(report.THRESHOLDS), source_sha256=sources, packages=versions(), python_version=platform.python_version(),
        live_packages=frozen["packages"], live_python_version=frozen["python_version"], range_scope=RANGE_SCOPE)
    bind_evidence(store, "protocol.json", protocol)
    for p in files:
        name = p.relative_to(report.CODE_ROOT).as_posix()
        snapshot = "source/"+name
        if not store.exists(snapshot):
            store.write_bytes(snapshot, p.read_bytes())
        require(store.sha256(snapshot) == sources[name], "numerical source snapshot differs")
    selected = [dict(lock_case(c["case_id"], [a for a in arms if a["case_id"] == c["case_id"]]),
                     input_hashes=c["hashes"], time_ms=c["time_ms"]) for c in frozen["cases"]]
    selection = dict(schema="nf-prospective-setting-selection.v1", protocol_sha256=store.sha256("protocol.json"),
        primary_evidence_hashes=hashes, seeds=frozen["seeds"], cases=selected, new_heat_used_for_selection=False)
    bind_evidence(store, "selection.json", selection)
    return selection


def score_name(case, identity, setting):
    return f"{case}/scores/{identity}.{setting}.json"


def validate_score(row):
    measured = validate_measurement(row["measurement"])
    require(row["status"] in ("completed","failed","not_evaluated_preparation_failed") and
            report.audit.finite(row["seconds"]) and row["seconds"] >= 0, "numerical score disposition or clock invalid")
    if measured["eligible"]:
        require(row["status"] == "completed", "failed score has eligible measurement")
        report.audit.check_measurement(measured, report.audit.Checks())


def score_setting(store, case, selected, setting, reference, invocation):
    """Preserve each request, including an actual baseline if GPU construction fails."""
    from d3d_optimization.boundary_evaluator import prepare_boundary
    from d3d_optimization.campaign_geometry import Candidate, evaluator
    from verification.nf_base_budget import measurement
    cid = case["case_id"]
    require(cid == selected["case_id"] and setting in SETTINGS, "setting request not registered")
    shared = dict(case_id=cid, setting=setting, setting_tuple=list(SETTINGS[setting]),
                  protocol_sha256=store.sha256("protocol.json"), selection_sha256=store.sha256("selection.json"))
    candidates = selected["candidates"]
    names = {c["id"]:score_name(cid,c["id"],setting) for c in candidates}
    pending = [c for c in candidates if not store.exists(names[c["id"]])]

    def record(c, outcome):
        return dict(shared, id=c["id"], candidate=c["candidate"], **outcome)

    if pending:
        existing_base = store.read(names[report.BASE]) if store.exists(names[report.BASE]) else None
        base_row = next(c for c in candidates if c["id"] == report.BASE)
        setup_at = perf_counter()
        engine, original, failure = None, None, None
        setup_file = f"{cid}/setup/{setting}.json"
        old_setup = store.read(setup_file) if store.exists(setup_file) else None
        if old_setup and old_setup["error"]:
            # Resume an interrupted write after a terminal preparation failure.
            failure = old_setup["error"]
            base_out = old_setup["base_outcome"]
        else:
            method, count, step, nodes, wall_step = SETTINGS[setting]
            try:
                original = prepare_boundary(Path(case["root"]), "OMP", reference, method=method, count=count,
                                            step=step, nodes=nodes, wall_step=wall_step)
                base_state = measurement_record(measurement(original.baseline, reference))
                base_out = dict(status="completed", measurement=base_state, error=None, seconds=0.)
            except Exception as exc:
                failure = f"original-wall preparation: {type(exc).__name__}: {exc}"
                base_out = dict(status="failed", measurement=missing_measurement(failure), error=failure, seconds=0.)
            if original is not None and failure is None and any(c["id"] != report.BASE for c in pending):
                try:
                    engine = evaluator(original, 20, "gpu")
                except Exception as exc:
                    failure = f"candidate evaluator: {type(exc).__name__}: {exc}"
            acquisition = dict(shared, base_outcome=base_out, error=failure,
                               status="failed" if failure else "completed", seconds=perf_counter()-setup_at)
            store.write(f"preparation/{cid}/{invocation.rsplit('/',1)[-1]}/{setting}.json", acquisition)
            if old_setup is None:
                store.write(setup_file, acquisition)
            elif base_out["status"] == "completed":
                require(compare_native(base_out["measurement"], [dict(measurement=old_setup["base_outcome"]["measurement"])])["passed"],
                        "resumed same-setting original wall differs")
        if existing_base is None:
            store.write(names[report.BASE], record(base_row,base_out))
        for c in pending:
            if c["id"] == report.BASE:
                continue
            if failure:
                outcome = dict(status="not_evaluated_preparation_failed", error=failure,
                               measurement=missing_measurement(failure), seconds=0.)
            else:
                outcome = evaluate_score(engine,Candidate.from_record(c["candidate"]))
            store.write(names[c["id"]], record(c,outcome))
    result = {}
    for c in candidates:
        row = store.read(names[c["id"]])
        require(all(row.get(k) == v for k,v in record(c,{}).items()), "numerical score provenance differs")
        validate_score(row)
        result[c["id"]] = row
    return result


def power_ledger_summary(scores):
    rows = []
    for score in scores:
        m = score["measurement"]
        ledger = m["ledger"]
        if (score["status"] != "completed" or not report.audit.positive(ledger.get("prescribed")) or
                not all(report.audit.finite(ledger.get(k)) for k in ("deposited","unhit"))):
            continue
        p = ledger["prescribed"]
        rows.append(dict(eligible=m["eligible"],closure=abs(ledger["deposited"]+ledger["unhit"]-p)/p,
                         unhit=ledger["unhit"]/p, deposited=ledger["deposited"]/p))
    eligible = [r for r in rows if r["eligible"]]
    return dict(registered=len(scores),measured_ledgers=len(rows),eligible_ledgers=len(eligible),
        max_relative_closure_all=max((r["closure"] for r in rows),default=None),
        max_relative_closure_eligible=max((r["closure"] for r in eligible),default=None),
        max_unhit_fraction_all=max((r["unhit"] for r in rows),default=None),
        max_unhit_fraction_eligible=max((r["unhit"] for r in eligible),default=None),
        eligible_deposited_over_prescribed_range=report.span([r["deposited"] for r in eligible]))


def summarize(store, selection):
    """Re-evaluate the report arithmetic from sealed native and non-native raw scores."""
    cases, evidence = [], {}
    for c in selection["cases"]:
        values = {}
        for design in c["candidates"]:
            identity = design["id"]
            values[identity] = {}
            for setting in SETTINGS:
                name = score_name(c["case_id"],identity,setting)
                row = store.read(name)
                require(row["selection_sha256"] == store.sha256("selection.json") and
                        row["protocol_sha256"] == store.sha256("protocol.json") and
                        (row["case_id"],row["id"],row["candidate"],row["setting"],row["setting_tuple"]) ==
                        (c["case_id"],identity,design["candidate"],setting,list(SETTINGS[setting])), "rehydrated score provenance differs")
                validate_score(row)
                values[identity][setting] = row
                evidence[name] = store.sha256(name)
        cases.append(assess_case(c,values))
    designs = [d for c in cases for d in c["candidates"]]
    returned = [r for d in designs for r in d["settings"].values()]
    status_counts = dict(Counter(r["status"] for r in returned))
    replay = [d["native_replay"] for d in designs]
    nonmatching = sum(r["passed"] is False for r in replay)
    missing_prior = sum(not r["genuine_prior_observations"] for r in replay)
    ineligible = sum(r["status"] == "completed" and not r["measurement"]["eligible"] for r in returned)
    failed = sum(r["status"] != "completed" for r in returned)
    preparations = [store.read(name) for name in store.artifacts("*.json") if name.startswith("preparation/")]
    setting_costs = {s:dict(preparation_seconds=sum(r["seconds"] for r in preparations if r["setting"] == s),
                            score_seconds=sum(d["settings"][s]["seconds"] for d in designs)) for s in SETTINGS}
    status = ("assessed_with_score_failures" if failed else "assessed_with_native_mismatch" if nonmatching else
              "assessed_with_ineligible" if ineligible else "assessed_with_unobserved_native" if missing_prior else "assessed")
    return dict(schema="nf-prospective-numerics-summary.v1", status=status,
                protocol_sha256=store.sha256("protocol.json"), selection_sha256=store.sha256("selection.json"),
                primary_evidence_hashes=selection["primary_evidence_hashes"], score_hashes=evidence,
                range_scope=RANGE_SCOPE, cases=cases, costs=cost_summary(store), costs_by_setting=setting_costs,
                power_ledgers_by_setting={s:power_ledger_summary([d["settings"][s] for d in designs]) for s in SETTINGS},
                registered_design_case_count=len(designs), registered_setting_requests=len(returned),
                actual_status_counts=status_counts, failed_scores=failed, completed_ineligible_scores=ineligible,
                native_coverage=dict(previously_observed_design_cases=len(designs)-missing_prior, without_prior=missing_prior,
                    compared_observations=sum(r["compared_count"] for r in replay),
                    compared_observables=sum(r["observable_compared_count"] for r in replay),
                    matched_observable_design_cases=sum(r["observables_reproduced"] is True for r in replay),
                    mismatched_design_cases=nonmatching),
                primary_hashes_verified_after_scoring=True)


def run(study, *, score):
    study = Path(study).resolve()
    require(study.drive.lower() == "d:", "study and all numerical evidence must remain on D")
    require((study/"frozen.json").is_file(), "holdout registration must exist before numerical verification")
    if score:
        require_locked(study)
    configure_environment(study,"prospective-numerics")
    from d3d_optimization.campaign_state import keep_awake, process_lock
    from d3d_optimization.reference_transport import ReferencePoint
    with process_lock(study/"campaign.lock"), process_lock(study/"posttrial.lock"):
        frozen, arms, hashes, errors = report.read_study(study)
        require(not errors, "live evidence failed audit: " + str(errors))
        require("audit.json" in hashes and all(a["status"] in (*report.CERTIFIED,"process_failed") for a in arms),
                "all registered attempts and the final supervisor audit must be present before selection")
        store = EvidenceStore(study/"numerics")
        store.verify_all()
        if not store.exists("selection.json"):
            require(not store.artifacts("*.started.json") and not store.artifacts("*.cubic1280.json"),
                    "numerical runs preceded full selection lock")
        locked = protocol_and_selection(study,frozen,arms,hashes,store)
        if not score:
            verify_external_hashes(study,hashes)
            return locked
        reference = ReferencePoint(**frozen["reference"])
        require(asdict(reference) == frozen["reference"], "numerical operating input differs")
        if not store.exists("summary.json"):
            with keep_awake():
                invocation = start_invocation(store)
                start = perf_counter()
                try:
                    for case,selected in zip(frozen["cases"],locked["cases"],strict=True):
                        for setting in SETTINGS:
                            values = score_setting(store,case,selected,setting,reference,invocation)
                            print("NUMERICS",case["case_id"],setting,len(values),flush=True)
                    verify_external_hashes(study,hashes)
                    # Inputs are hashed again after the last setting, independently
                    # of the primary-file hashes (which describe live artifacts).
                    for case in frozen["cases"]:
                        for name,wanted in case["hashes"].items():
                            require(report.audit.digest(Path(case["root"])/name) == wanted, "staged numerical input changed")
                except BaseException as exc:
                    store.write(invocation+".finished.json",dict(status="failed",seconds=perf_counter()-start,error=f"{type(exc).__name__}: {exc}"))
                    raise
                store.write(invocation+".finished.json",dict(status="completed",seconds=perf_counter()-start))
        verify_external_hashes(study,hashes)
        result = summarize(store,locked)
        bind_evidence(store,"summary.json",result)
        store.verify_all()
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--lock-only", action="store_true", help="seal all cases and choices without starting heat calculations")
    mode.add_argument("--score", action="store_true", help="compute four settings using a previously sealed selection")
    args = parser.parse_args()
    result = run(args.study,score=args.score)
    print(result.get("status","selection_locked"), len(result["cases"]), "equilibria", flush=True)
    if result.get("failed_scores") or result.get("native_coverage",{}).get("mismatched_design_cases"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
