"""Read-only reconstruction of the November 7 paired cold-search experiment."""
import argparse
from collections import Counter, defaultdict
from itertools import combinations
import json
import math
from pathlib import Path
from statistics import median

from verification import nf_anytime_final_audit as audit
from verification.nf_anytime_reference import require


SEEDS = (20711, 20729, 20747)
FAST, DIRECT = POLICIES = ("coarse_adaptive", "fine_adaptive")
SHORT = (200, 400, 600)
LONG = 1800
SETTINGS = ("cubic640", "cubic1280", "quintic640", "quintic1280")
BASE = "original_wall"
THRESHOLDS = (.5, 1., 2.)
PHASES = ("population_initialization", "de_trial", "coordinate_polish", "external_validation")
GEOMETRY_ERRORS = {"continuous_shape_constraint_failed","complete_wall_invalid","material_length_exceeds_limit",
                   "geometry_constraint_failed"}
REFERENCE = dict(power_W=1e6, density_m3=1e19, chi_perp_m2_s=.3, lambda_q_m=.003,
                 gamma=7., kappa0=2000.)
CERTIFIED = ("completed", "arm_failed", "early_stopped")
CODE_ROOT = Path(__file__).resolve().parents[1]


def cutoffs(budget):
    require(budget in (600, LONG), "unregistered trial budget")
    return (*SHORT, LONG) if budget == LONG else SHORT


def span(values):
    return [min(values), max(values)] if values else None


def margin(fast, direct):
    require(audit.positive(fast) and audit.positive(direct), "two positive finite peaks required")
    gap = 100.*(direct-fast)/direct
    return dict(gap_pct=gap, sign="fast" if fast < direct else "direct" if direct < fast else "tie",
                practical={f"{t:g}": "fast" if gap >= t else "direct" if gap <= -t else "within"
                           for t in THRESHOLDS})


def check_endpoints(record, raw):
    """Recompute every worker checkpoint from the flush-ordered raw event stream."""
    require(record["events"] == raw, "raw event stream differs from worker result")
    budget = record["budget_s"]
    expected = {str(t): audit.checkpoint(raw, t, budget) for t in cutoffs(budget)}
    saved = record["checkpoints"]
    # Legacy files store the long-budget endpoint under final. Some wrappers also
    # expose it by cutoff; in both formats it must equal the raw reconstruction.
    require(set(saved) in (set(map(str, SHORT)), set(expected)), "unknown worker checkpoint schedule")
    for name, value in saved.items():
        require(value == expected[name], f"raw checkpoint {name} differs from worker")
    require(record["final"] == expected[str(int(budget))], "raw final checkpoint differs from worker")
    bases = [e for e in raw if e["stage"] == "baseline" and e["status"] == "completed"]
    require(len(bases) <= 1, "multiple original-wall measurements")
    base = bases[0] if bases else {}
    require(record.get("baseline") == base.get("measurement") and
            record.get("baseline_ready_s") == base.get("end_s"), "raw baseline differs from worker")


def phase_counts(events, level, kind, cutoff):
    scored = [e for e in events if e["stage"] == level and e.get("proposal_kind") == kind]
    done = [e for e in scored if e.get("status") == "completed" and e["end_s"] <= cutoff]
    fresh = [e for e in done if not e.get("cached")]
    deposition = [e for e in fresh if e.get("physics_evaluated") is True and audit.returned_deposition(e)]
    geom = [e for e in fresh if not audit.returned_deposition(e) and GEOMETRY_ERRORS.intersection(e["measurement"]["errors"])]
    return dict(requests=len(done), started=sum(e["start_s"] < cutoff for e in scored),
                fresh=len(fresh), cached=len(done)-len(fresh), deposition=len(deposition),
                eligible_deposition=sum(e["measurement"]["eligible"] for e in deposition),
                geometry_rejections=len(geom), returned_ineligible=sum(not e["measurement"]["eligible"] for e in deposition),
                completed_seconds=sum(e["end_s"]-e["start_s"] for e in done),
                exposed_seconds=sum(max(0., min(cutoff, e["end_s"])-e["start_s"]) for e in scored))


def endpoint(arm, cutoff):
    """Retain the observed prefix on failure, certifying only fully observed budgets."""
    require(cutoff in cutoffs(arm["budget_s"]), "cutoff outside registered budget")
    raw = arm["events"]
    verified = arm["raw_certified"]
    status = arm["status"]
    horizon = float(cutoff)
    if status == "arm_failed":
        last_observed = max((e["end_s"] for e in raw if audit.finite(e.get("end_s"))), default=0.)
        failed_at = next((e["end_s"] for e in raw if e.get("status") == "failed"), last_observed)
        horizon = min(cutoff, failed_at)
    elif not verified:
        horizon = min(cutoff, arm["elapsed_s"]) if arm.get("elapsed_s") is not None else None
    elif status == "early_stopped":
        horizon = min(cutoff, arm["elapsed_s"])
    snap = audit.checkpoint(raw, horizon, arm["budget_s"]) if verified else None
    base_event = next((e for e in raw if verified and e["stage"] == "baseline" and
                       e.get("status") == "completed" and e["end_s"] <= horizon), None)
    base = base_event["measurement"] if base_event else None
    observed = snap["incumbent"] if snap else None
    available = bool(observed and (status not in ("arm_failed","early_stopped") or cutoff <= horizon))
    selected = observed if available else None
    first = snap["first_5pct_s"] if snap else None
    if first is not None:
        reason = "observed"
    elif status not in CERTIFIED:
        reason = status
    elif status == "arm_failed" and cutoff > horizon:
        reason = "arm_failed"
    elif not base or not base["eligible"]:
        reason = "baseline_unavailable"
    elif status == "early_stopped" and horizon < cutoff:
        reason = "early_stopped"
    else:
        reason = "deadline_without_5pct"
    active = "coarse" if arm["policy"] == FAST else "fine"
    phases = {k: phase_counts(raw, "fine" if k == "external_validation" else active, k, horizon)
              for k in PHASES} if verified else {}
    censored = first is None and verified and base is not None and base["eligible"]
    row = dict(case_id=arm["case_id"], seed=arm["seed"], policy=arm["policy"], cutoff_s=cutoff,
               status=status, available=available, observed_until_s=horizon, raw_certified=verified,
               baseline_status=("eligible" if base["eligible"] else "ineligible") if base else "unavailable",
               baseline=base, baseline_ready_s=base_event["end_s"] if base_event else None,
               observed_before_stop=observed if verified and not available else None,
               selected_id=None, candidate=None, kind=None, peak=None, A2=None, q_ratio=None, A2_ratio=None,
               peak_reduction_pct=None, A2_change_pct=None, joint_improvement=None, A2_decreased=None,
               power_relative_error=None, selected_ledger=None, deposited_over_prescribed=None,
               deposited_over_baseline=None, unhit_fraction=None,
               counts=snap["counts"] if snap else None, phases=phases, de=phases.get("de_trial"),
               first_5pct=dict(event_observed=first is not None, time_s=first, reason=reason,
                               right_censored=censored, exposure_end_s=horizon, censor_time_s=horizon if censored else None))
    if available:
        if selected["kind"] == BASE:
            chosen = base
        else:
            chosen = next(e["measurement"] for e in raw if e["stage"] == "fine" and
                          e.get("id") == selected["id"] and e.get("status") == "completed" and
                          e["end_s"] <= horizon and not e.get("cached") and e["measurement"]["eligible"] and
                          (e["measurement"]["peak"], e["measurement"]["A2"]) == (selected["peak"], selected["A2"]))
        q, a = chosen["peak"], chosen["A2"]
        ledger = chosen["ledger"]
        row.update(kind=selected["kind"], selected_id=selected["id"] or BASE, candidate=selected["candidate"],
                   peak=q, A2=a, q_ratio=q/base["peak"], A2_ratio=a/base["A2"],
                   peak_reduction_pct=100.*(base["peak"]-q)/base["peak"],
                   A2_change_pct=100.*(a-base["A2"])/base["A2"], A2_decreased=a < base["A2"],
                   joint_improvement=q < base["peak"] and a > base["A2"], selected_ledger=ledger,
                   deposited_over_prescribed=ledger["deposited"]/ledger["prescribed"],
                   deposited_over_baseline=(ledger["deposited"]/base["ledger"]["deposited"]
                                            if base["ledger"]["deposited"] > 0 else None),
                   unhit_fraction=ledger["unhit"]/ledger["prescribed"],
                   power_relative_error=abs(ledger["deposited"]+ledger["unhit"]-ledger["prescribed"])/ledger["prescribed"])
    return row


def matching_baselines(a,b):
    if not a or not b:
        return None
    if a == b:
        return True
    if not a["eligible"] or not b["eligible"]:
        return False
    # The physical input hashes are checked independently. Share the numerical
    # replay tolerance for thermal observables and the solver tolerance for power.
    if any(not math.isclose(a[k],b[k],rel_tol=1e-10,abs_tol=0.) for k in ("peak","A2")):
        return False
    x,y = a["ledger"],b["ledger"]
    scale = max(x["prescribed"],y["prescribed"])
    return all(abs(x[k]-y[k]) <= 2e-12*scale for k in ("prescribed","deposited","unhit"))


def pair_endpoints(fast, direct):
    keys = ("case_id", "seed", "cutoff_s")
    require(all(fast[k] == direct[k] for k in keys) and (fast["policy"], direct["policy"]) == POLICIES,
            "unmatched case, seed, method or cutoff")
    same_base = matching_baselines(fast["baseline"], direct["baseline"])
    available = fast["available"] and direct["available"] and same_base is True
    gap = margin(fast["peak"], direct["peak"]) if available else None
    eligible_de = direct["de"]["eligible_deposition"] if direct["de"] else None
    hit = [x["first_5pct"]["time_s"] for x in (fast, direct)]
    timing = ("fast" if hit[0] < hit[1] else "direct" if hit[0] > hit[1] else "equal"
              ) if available and None not in hit else "not_both_observed"
    return dict(**{k: fast[k] for k in keys}, fast=fast, direct=direct, available=available,
                baseline_equal=same_base, baseline_exact=(fast["baseline"] == direct["baseline"]) if same_base is not None else None, margin=gap,
                same_design=(fast["selected_id"] == direct["selected_id"]) if available else None,
                direct_eligible_de=eligible_de, first_5pct_order=timing,
                fast_practical_and_direct_de=bool(gap and gap["practical"]["1"] == "fast" and eligible_de))


def pair_counts(pairs):
    valid = [p for p in pairs if p["available"]]
    thresholds = {}
    for t in THRESHOLDS:
        key = f"{t:g}"
        counts = Counter(p["margin"]["practical"][key] for p in valid)
        thresholds[key] = {label: counts[label] for label in ("fast", "within", "direct")}
    signed = Counter(p["margin"]["sign"] for p in valid)
    gaps = [p["margin"]["gap_pct"] for p in valid]
    return dict(registered=len(pairs), comparable=len(valid), unavailable=len(pairs)-len(valid),
                thresholds=thresholds, signs={k:signed[k] for k in ("fast", "tie", "direct")},
                direct_eligible_de_pairs=sum(bool(p["direct_eligible_de"]) for p in valid),
                fast_practical_and_direct_de=sum(p["fast_practical_and_direct_de"] for p in pairs),
                gap_seed_median_pct=median(gaps) if gaps else None, gap_seed_range_pct=span(gaps),
                first_5pct_order_counts=dict(Counter(p["first_5pct_order"] for p in pairs)))


def geometry_clusters(pairs):
    designs = defaultdict(list)
    vectors = {}
    for p in pairs:
        for key in ("fast", "direct"):
            e = p[key]
            if e["available"]:
                identity = e["selected_id"]
                designs[identity].append(dict(policy=e["policy"], seed=e["seed"]))
                vector = e["candidate"]["vector"] if e["candidate"] else None
                require(vectors.setdefault(identity, vector) == vector, "inconsistent selected geometry")
    distances = []
    for left, right in combinations([key for key in designs if vectors[key] is not None], 2):
        distances.append(dict(ids=[left, right],
                              max_coefficient_difference_m=max(abs(x-y) for x, y in zip(vectors[left], vectors[right]))))
    return dict(selection_count=sum(map(len, designs.values())), unique_identities=len(designs),
                identities=dict(designs), distinct_vector_distances=distances,
                distance_meaning="maximum absolute difference in the sixteen OMP displacement coefficients")


def summarize_live(arms):
    groups = defaultdict(dict)
    for arm in arms:
        key = (arm["case_id"], arm["seed"])
        require(arm["policy"] not in groups[key], "duplicate arm")
        groups[key][arm["policy"]] = arm
    pairs = []
    for methods in groups.values():
        require(set(methods) == set(POLICIES), "registered pair member omitted")
        fast, direct = (methods[k] for k in POLICIES)
        require(fast["budget_s"] == direct["budget_s"], "unmatched budgets")
        for t in cutoffs(fast["budget_s"]):
            pairs.append(pair_endpoints(endpoint(fast, t), endpoint(direct, t)))
    case_groups = []
    case_ids = list(dict.fromkeys(p["case_id"] for p in pairs))
    for cid in case_ids:
        for t in (*SHORT, LONG):
            current = [p for p in pairs if p["case_id"] == cid and p["cutoff_s"] == t]
            if current:
                counts = pair_counts(current)
                complete = counts["comparable"] == counts["registered"]
                med = counts["gap_seed_median_pct"]
                case_groups.append(dict(case_id=cid, cutoff_s=t, **counts, geometry=geometry_clusters(current),
                    complete_seed_set=complete,
                    median_at_least_1pct_fast=complete and med >= 1.,
                    median_at_least_1pct_direct=complete and med <= -1.))
    by_cutoff = {}
    for t in (*SHORT, LONG):
        current = [p for p in pairs if p["cutoff_s"] == t]
        if current:
            counts = pair_counts(current)
            # Across cases, seeds are repeated searches and neither medians nor
            # intervals of the pooled margins are the registered case estimand.
            for key in ("gap_seed_median_pct", "gap_seed_range_pct"):
                counts.pop(key)
            groups_at_t = [g for g in case_groups if g["cutoff_s"] == t]
            by_cutoff[str(t)] = dict(distinct_equilibria=len(groups_at_t),
                complete_equilibria=sum(g["complete_seed_set"] for g in groups_at_t),
                equilibria_fast_median_1pct=sum(g["median_at_least_1pct_fast"] for g in groups_at_t),
                equilibria_direct_median_1pct=sum(g["median_at_least_1pct_direct"] for g in groups_at_t), **counts)
    return dict(pairs=pairs, case_groups=case_groups, by_cutoff=by_cutoff)


def summarize_cost(arm):
    raw = [e for e in arm["events"] if isinstance(e,dict) and isinstance(e.get("stage"),str)]
    stages = Counter()
    for e in raw:
        if audit.finite(e.get("start_s")) and audit.finite(e.get("end_s")):
            stages[e["stage"]] += e["end_s"]-e["start_s"]
    late = [e for e in raw if e.get("stage") in ("coarse","fine") and e.get("status") == "late"]
    base = next((e for e in raw if e.get("stage") == "baseline" and e.get("status") == "completed"), {})
    proc = arm.get("process") or {}
    return dict(case_id=arm["case_id"], seed=arm["seed"], policy=arm["policy"], status=arm["status"],
                error=arm.get("error"), budget_s=arm["budget_s"], parent_elapsed_s=proc.get("elapsed_s"),
                worker_elapsed_s=arm.get("worker_elapsed_s", arm.get("elapsed_s")) if arm["raw_certified"] else None,
                startup_s=arm.get("startup_s"), baseline_ready_s=base.get("end_s"),
                stage_seconds=dict(stages), returned_events=len(raw), late_scores=len(late),
                late_score_events=[{k:e.get(k) for k in ("id","stage","start_s","end_s","measurement","cached")} for e in late],
                failure_events=[{k:e.get(k) for k in ("stage","start_s","end_s","error")} for e in raw if e.get("status") == "failed"])


def load_arm(study, job, registration_sha256, hashes, *, sealed=False):
    study = Path(study).resolve()
    folder = (study/job["relative_dir"]).resolve()
    expected = study/"live"/job["case_id"]/str(job["seed"])/job["policy"]
    require(folder == expected.resolve() and folder.is_relative_to(study/"live"), "unregistered output path")
    arm = dict(**{k:job[k] for k in ("case_id","seed","policy","budget_s","population","projected_seeds")},
               status="missing", events=[], process=None, elapsed_s=None, raw_certified=False, error=None)
    if not folder.exists():
        return arm
    try:
        names = ("job.json","events.jsonl","result.json","process.json","stdout.log")
        byte_hashes = {}
        for name in names:
            path = folder/name
            if path.is_file():
                byte_hashes[name] = audit.digest(path)
                hashes[path.relative_to(study).as_posix()] = byte_hashes[name]
        if "process.json" not in byte_hashes:
            arm["status"] = "unfinished"
            return arm
        proc = audit.read(folder/"process.json")
        arm["process"] = proc
        job_file = audit.read(folder/"job.json")
        if sealed:
            from verification.nf_anytime_reference import EvidenceStore
            receipt = EvidenceStore(folder)
            require(receipt.sha256("job.json") == byte_hashes["job.json"] and
                    receipt.sha256("process.json") == byte_hashes["process.json"], "parent evidence seal differs")
            env = proc["environment"]
            for variable,leaf in (("CUPY_CACHE_DIR","cupy"),("CUDA_CACHE_PATH","cuda"),("TEMP","tmp"),
                                  ("TMP","tmp"),("MPLCONFIGDIR","matplotlib"),("HF_HOME","huggingface")):
                require(Path(env[variable]).resolve() == folder/"cache"/leaf, "nonprivate worker cache")
            require(env["NF_PROSPECTIVE_JOB_SHA256"] == byte_hashes["job.json"] and
                    all(env[k] == "1" for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")),
                    "worker environment differs")
        require(job_file == dict(job, registration_sha256=registration_sha256), "job differs from frozen registration")
        require(proc["job_sha256"] == byte_hashes.get("job.json") and
                proc["registration_sha256"] == registration_sha256, "process job or registration hash differs")
        start, stop = proc["launch_ns"], proc["exit_observed_ns"]
        require(type(start) is int and type(stop) is int and 0 <= start <= stop and audit.finite(proc["elapsed_s"]) and
                abs(proc["elapsed_s"]-(stop-start)/1e9) <= 1e-8, "invalid parent process clock")
        arm["elapsed_s"] = proc["elapsed_s"]
        for file, field in (("events.jsonl","events_sha256"),("result.json","result_sha256"),("stdout.log","stdout_sha256")):
            require(proc[field] == byte_hashes.get(file), "process hash differs for " + file)
        failed_process = proc["returncode"] not in (0, 1) or (proc["returncode"] == 1 and "result.json" not in byte_hashes)
        raw, partial = (audit.read_events(folder/"events.jsonl", allow_partial_tail=failed_process)
                        if "events.jsonl" in byte_hashes else ([], False))
        arm.update(events=raw, partial_event_tail=partial)
        if failed_process:
            arm.update(status="process_failed", error=proc.get("error") or f"worker exit {proc['returncode']}")
            return arm
        require("events.jsonl" in byte_hashes, "missing raw stream for finalized worker")
        rec = audit.read(folder/"result.json")
        for key in ("case_id","policy","seed","budget_s","input_root"):
            require(rec[key] == job[key], "worker identity differs for " + key)
        require(rec["launch_ns"] == start and rec["job_sha256"] == byte_hashes["job.json"] and
                rec["registration_sha256"] == registration_sha256, "worker provenance differs")
        require(rec["elapsed_s"] <= proc["elapsed_s"]+1e-8, "worker elapsed exceeds process time")
        check_endpoints(rec, raw)
        from verification.nf_prospective_worker import audit_arm
        verified = audit_arm(rec, raw, **{k:job[k] for k in ("policy","seed","budget_s","population","projected_seeds")})
        require(proc["returncode"] == int(not verified["normal"]), "worker exit differs from verified disposition")
        status = "completed" if verified["normal"] else "early_stopped" if rec["stop"] == "search_stopped" else "arm_failed"
        arm.update(status=status, raw_certified=True, elapsed_s=rec["elapsed_s"], worker_elapsed_s=rec["elapsed_s"],
                   startup_s=rec["startup_s"], error=rec.get("error") or " | ".join(verified["reasons"]) or None, parity=rec.get("parity"),
                   startup_full_rank=verified.get("startup_full_rank"))
        for e in raw:
            if e["stage"] in ("baseline","coarse","fine") and e.get("measurement"):
                audit.check_measurement(e["measurement"], audit.Checks())
        if rec.get("reference") is not None:
            require(rec["reference"] == REFERENCE, "worker prescribed transport differs")
        return arm
    except (ValueError, OSError, KeyError, TypeError, IndexError, AssertionError) as exc:
        arm.update(status="audit_failed", raw_certified=False, error=f"{type(exc).__name__}: {exc}")
        return arm


def read_study(study):
    from verification.nf_anytime_reference import EvidenceStore
    from verification.nf_prospective_campaign import verify_registration
    from verification.nf_unused_live import load_cases
    study = Path(study).resolve()
    store = EvidenceStore(study)
    frozen = store.read("frozen.json")
    require(frozen["phase"] == "holdout" and frozen["source"] == "OMP" and frozen["radius_mm"] == 20,
            "wrong physical scope or campaign phase")
    verify_registration(store, frozen)
    cases = frozen["cases"]
    require(len(cases) == 3 and cases == load_cases(frozen["source_study"]), "input preflight/frozen equilibrium linkage differs")
    require(frozen["reference"] == REFERENCE and frozen["levels"] == audit.LEVELS, "registered operating point or scoring levels differ")
    require(frozen["seeds"] == list(SEEDS), "registered seeds changed")
    names = ("frozen.json","protocol.ko.md","development-choice.json","development-plan.json","development-audit.json")
    hashes = {name:audit.digest(study/name) for name in names}
    decision = store.read("development-choice.json")
    require(decision["qualified"] and decision["selected"] == frozen["development_choice"] and
            decision["development_audit_sha256"] == hashes["development-audit.json"] and
            decision["development_plan_sha256"] == hashes["development-plan.json"], "development choice differs from freeze")
    for field, name in (("development_choice_sha256","development-choice.json"),("protocol_sha256","protocol.ko.md"),
                        ("study_protocol_sha256","protocol.ko.md")):
        if field in frozen:
            require(frozen[field] == hashes[name], "frozen metadata hash differs for " + name)
    registration = hashes["frozen.json"]
    jobs = frozen["jobs"]
    expected_keys = {(c["case_id"], s, p) for c in cases for s in SEEDS for p in POLICIES}
    keys = [(j["case_id"],j["seed"],j["policy"]) for j in jobs]
    require(len(keys) == len(expected_keys) and set(keys) == expected_keys, "frozen jobs differ from case/seed pairs")
    direct_configs = {(j["population"],j["projected_seeds"]) for j in jobs if j["policy"] == DIRECT}
    require(len(direct_configs) == 1 and direct_configs <= {(24,8),(32,12)}, "inconsistent tuned direct comparator")
    require(direct_configs == {(decision["selected"]["population"],decision["selected"]["projected_seeds"])},
            "held-out direct comparator differs from development selection")
    roots = {c["case_id"]:c["root"] for c in cases}
    for job in jobs:
        require(job["phase"] == "holdout" and job["input_root"] == roots[job["case_id"]] and
                job["budget_s"] == (LONG if job["seed"] == SEEDS[0] else 600) and
                (job["policy"] != FAST or (job["population"],job["projected_seeds"]) == (64,12)), "invalid frozen job geometry, population or budget")
    arms = [load_arm(study,j,registration,hashes,sealed=True) for j in jobs]
    errors = [f"{a['case_id']}/{a['seed']}/{a['policy']} {a['error']}" for a in arms if a["status"] == "audit_failed"]
    previous_exit = None
    for a in arms:
        proc = a["process"]
        if proc and type(proc.get("launch_ns")) is int and type(proc.get("exit_observed_ns")) is int:
            if previous_exit is not None and proc["launch_ns"] < previous_exit:
                errors.append(f"frozen-order workers overlap at {a['case_id']}/{a['seed']}/{a['policy']}")
            previous_exit = proc["exit_observed_ns"]
    dirs = {j["relative_dir"].replace("\\","/") for j in jobs}
    extras = {p.parent.relative_to(study).as_posix() for p in (study/"live").glob("*/*/*/process.json")} - dirs
    if extras:
        errors.append("extra workers outside registration: " + str(sorted(extras)))
    for case in cases:
        bases = [e["measurement"] for a in arms if a["case_id"] == case["case_id"] and a["raw_certified"]
                 for e in a["events"] if e["stage"] == "baseline" and e["status"] == "completed"]
        if bases and any(not matching_baselines(b,bases[0]) for b in bases):
            errors.append("native original-wall baseline differs within " + case["case_id"])
    if store.exists("audit.json"):
        parent = store.read("audit.json")
        hashes["audit.json"] = store.sha256("audit.json")
        require(parent["registration_sha256"] == registration and len(parent["rows"]) == parent["planned_arm_count"] == len(jobs),
                "supervisor audit denominator or study identity differs")
        for j, a, prior in zip(jobs,arms,parent["rows"]):
            require(all(j[k] == prior[k] for k in ("case_id","seed","policy","relative_dir")), "supervisor audit order differs")
            if a["raw_certified"]:
                for cutoff in cutoffs(j["budget_s"]):
                    evaluated = endpoint(a,cutoff)
                    published = prior["checkpoints"][str(cutoff)]
                    if evaluated["available"]:
                        require((evaluated["peak"], evaluated["A2"]) == (published["peak"],published["A2"]) and
                                evaluated["de"]["eligible_deposition"] == published["eligible_thermal_de_trials"],
                                "supervisor and report endpoints differ")
    return frozen, arms, hashes, errors


def baseline_totals(endpoints):
    valid = [e for e in endpoints if e["available"]]
    hits = [e["first_5pct"]["time_s"] for e in endpoints if e["first_5pct"]["event_observed"]]
    return dict(registered=len(endpoints), available=len(valid), improved_q=sum(e["q_ratio"] < 1 for e in valid),
        q_down_A2_up=sum(e["joint_improvement"] for e in valid), A2_down=sum(e["A2_decreased"] for e in valid),
        fallback_original_wall=sum(e["kind"] == BASE for e in valid),
        no_eligible_new_fine=sum(e["counts"]["eligible_fine"] == 0 for e in valid),
        peak_reduction_range_pct=span([e["peak_reduction_pct"] for e in valid]),
        A2_change_range_pct=span([e["A2_change_pct"] for e in valid]),
        first_5pct_observed=len(hits), first_5pct_observed_time_range_s=span(hits),
        first_5pct_reason_counts=dict(Counter(e["first_5pct"]["reason"] for e in endpoints)),
        max_power_closure=max((e["power_relative_error"] for e in valid),default=None),
        max_unhit_fraction=max((e["unhit_fraction"] for e in valid),default=None),
        deposited_over_baseline_range=span([e["deposited_over_baseline"] for e in valid if e["deposited_over_baseline"] is not None]))


def development(study):
    from verification.nf_anytime_reference import EvidenceStore
    from verification.nf_prospective_campaign import verify_process
    from verification.nf_prospective_worker import choose_fine
    store = EvidenceStore(study)
    require(store.exists("development-plan.json"), "development registration missing")
    plan = store.read("development-plan.json")
    plan_hash = store.sha256("development-plan.json")
    parent = store.read("development-audit.json") if store.exists("development-audit.json") else None
    choice = store.read("development-choice.json") if store.exists("development-choice.json") else None
    registered = plan["jobs"]
    require(len(registered) == 6, "development expects all six registered arms")
    if parent:
        require(parent["registration_sha256"] == plan_hash and len(parent["rows"]) == len(registered), "development denominator or phase differs")
    if choice:
        require(parent and choice["development_plan_sha256"] == plan_hash and
                choice["development_audit_sha256"] == store.sha256("development-audit.json"), "development pilot audit linkage differs")
        reproduced = choose_fine(parent["rows"])
        require(all(choice[k] == v for k,v in reproduced.items()), "development selector differs on replay")
    rows = []
    config_columns = ("population","projected_seeds","de_cases","geometric_mean_peak_ratio","qualified","reasons")
    for i,j in enumerate(registered):
        folder = study/j["relative_dir"]
        raw = parent["rows"][i] if parent else {}
        if raw:
            require(all(raw[k] == j[k] for k in ("case_id","seed","policy","population","projected_seeds")), "development worker identity differs")
        p = verify_process(folder,dict(j,registration_sha256=plan_hash)) if (folder/"process.json").exists() else {}
        if raw:
            require(p and raw["process_sha256"] == audit.digest(folder/"process.json"), "development process bytes differ")
        rows.append(dict(**{k:j[k] for k in ("case_id","seed","population","projected_seeds","budget_s","relative_dir")},
            status=raw.get("status","audit_pending" if p else "unfinished" if folder.exists() else "missing"),
            peak=raw.get("peak"), baseline_peak=raw.get("baseline_peak"), peak_ratio=raw.get("peak_ratio"),
            de=raw.get("eligible_thermal_de_trials") if raw.get("status") != "process_failed" else None,
            startup_full_rank=raw.get("startup_full_rank") if raw.get("status") != "process_failed" else None,
            parent_elapsed_s=p.get("elapsed_s"), error=raw.get("error") or p.get("error"),
            process_sha256=audit.digest(folder/"process.json") if p else None))
    cases = [dict(case_id=c["case_id"], time_ms=c["time_ms"], hashes=c["hashes"]) for c in plan["cases"]]
    return dict(registered=len(registered), cases=cases, rows=rows, completed_receipts=sum(r["parent_elapsed_s"] is not None for r in rows),
        parent_seconds=sum(r["parent_elapsed_s"] for r in rows if r["parent_elapsed_s"] is not None),
        qualified=choice["qualified"] if choice else None, selected=choice["selected"] if choice else None,
        reasons=choice["reasons"] if choice else None,
        configurations=[{k:c[k] for k in config_columns} for c in choice["configurations"]] if choice else [],
        plan_sha256=plan_hash, audit_sha256=store.sha256("development-audit.json") if parent else None,
        decision_sha256=store.sha256("development-choice.json") if choice else None)


def load_numerics(study, frozen, arms, hashes):
    from verification.nf_anytime_reference import EvidenceStore, verify_external_hashes
    from verification import nf_prospective_numerics as numeric
    folder = study/"numerics"
    if not (folder/"selection.json").exists():
        return dict(status="not_locked")
    store = EvidenceStore(folder)
    store.verify_all()
    selection,protocol = store.read("selection.json"),store.read("protocol.json")
    require(selection["protocol_sha256"] == store.sha256("protocol.json") and selection["primary_evidence_hashes"] == hashes and
            protocol["registration_sha256"] == hashes["frozen.json"] and selection["seeds"] == list(SEEDS),
            "post-trial selection hash or live evidence differs")
    require(protocol["settings"] == {k:list(v) for k,v in numeric.SETTINGS.items()} and
            protocol["native"] == numeric.NATIVE and protocol["reference"] == frozen["reference"] and
            protocol["native_relative_tolerance"] == numeric.NATIVE_RTOL, "post-trial model/settings differ")
    for name,wanted in protocol["source_sha256"].items():
        require(audit.digest(CODE_ROOT/name) == wanted and store.sha256("source/"+name) == wanted, "numerical implementation changed: "+name)
    for c,locked in zip(frozen["cases"],selection["cases"],strict=True):
        rebuilt = numeric.lock_case(c["case_id"],[a for a in arms if a["case_id"] == c["case_id"]])
        require(locked == dict(rebuilt,input_hashes=c["hashes"],time_ms=c["time_ms"]), "numerical registry differs from independent raw selection")
    if not store.exists("summary.json"):
        return dict(status="locked_scores_incomplete", selection_sha256=store.sha256("selection.json"),
            locked_design_cases=sum(len(c["candidates"]) for c in selection["cases"]),
            completed_score_files=len([n for n in store.artifacts("*.json") if "/scores/" in n]))
    result = numeric.summarize(store,selection)
    require(store.read("summary.json") == result, "numerical report differs from immutable score files")
    verify_external_hashes(study,hashes)
    return dict(result, summary_sha256=store.sha256("summary.json"), original_live_choices_reconstructed_unchanged=True)


def build_summary(study):
    from verification.nf_anytime_reference import verify_external_hashes
    study = Path(study).resolve()
    dev = development(study)
    summary = dict(schema="nf-prospective-report.v1", study=str(study), development=dev,
        metric="G = 100*(q_direct-q_fast)/q_direct in percent, positive means lower fast peak",
        thresholds_pct=list(THRESHOLDS), equality_at_threshold_counts_as_practical=True,
        scope="fixed prescribed transport, anonymous DIII-D equilibria, bounded cold-search budgets and 16 OMP variables")
    if not (study/"frozen.json").is_file():
        return dict(summary,status="development_unqualified" if dev["qualified"] is False else "holdout_not_registered",
                    live=None, numerical=dict(status="not_locked"), errors=[])
    frozen,arms,hashes,errors = read_study(study)
    live = summarize_live(arms)
    by_policy = {}
    for t in (*SHORT,LONG):
        selected = [p for p in live["pairs"] if p["cutoff_s"] == t]
        by_policy[str(t)] = {key:baseline_totals([p[key] for p in selected]) for key in ("fast","direct")}
    costs = [summarize_cost(a) for a in arms]
    try:
        numerical = load_numerics(study,frozen,arms,hashes)
    except (ValueError,OSError,KeyError,TypeError,IndexError,AssertionError) as exc:
        numerical = dict(status="numerical_audit_failed", error=f"{type(exc).__name__}: {exc}")
        errors.append("numerical evidence: " + numerical["error"])
    statuses = dict(Counter(a["status"] for a in arms))
    closed = sum(a["status"] in (*CERTIFIED,"process_failed") for a in arms)
    status = ("audit_failed" if errors else "holdout_incomplete" if closed < len(arms) or "audit.json" not in hashes else
              "numerical_incomplete" if numerical["status"] in ("not_locked","locked_scores_incomplete") else
              "closed_with_exceptions" if numerical["status"] != "assessed" or any(a["status"] != "completed" for a in arms) else "closed")
    input_study = Path(frozen["source_study"])
    screening = audit.read(input_study/"preflight.json")["cases"]
    cohorts = [dict(case_id=c["case_id"],time_ms=c["time_ms"],hashes=c["hashes"],input_root=c["root"]) for c in frozen["cases"]]
    live.update(method_baselines=by_policy, parent_costs=costs, parent_seconds=sum(c["parent_elapsed_s"] for c in costs if c["parent_elapsed_s"] is not None),
                registered_arms=len(arms), with_process_receipt=sum(a["process"] is not None for a in arms), closed_arms=closed,
                status_counts=statuses, audit_errors=errors)
    verify_external_hashes(study,hashes)
    return dict(summary,status=status, errors=errors, reference=frozen["reference"], cases=cohorts,
                input_screening=dict(examined=len(screening),accepted=len(cohorts), source_study=str(input_study),
                    exclusions=[dict(case_id=c["case_id"],reason=c.get("reason"),rank=c.get("rank")) for c in screening if c.get("status") != "passed"]),
                live_evidence_sha256=hashes, source_metadata_sha256=frozen["source_records"], hardware=frozen["hardware"],
                live=live, numerical=numerical)


def fmt(value,digits=3):
    return "—" if value is None else f"{value:.{digits}f}"


def interval(pair,digits=2):
    return "—" if pair is None else f"{fmt(pair[0],digits)} … {fmt(pair[1],digits)}"


def table(headers,rows):
    return "\n".join("| "+" | ".join(map(str,r))+" |" for r in [headers,["---"]*len(headers),*rows])


def write_once(path,payload):
    """Published files can be regenerated exactly, never replaced by a different study."""
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes() == payload, f"different published bytes already exist at {path}")
        return
    with path.open("xb") as stream:
        stream.write(payload)


def hit_text(e):
    hit = e["first_5pct"]
    if hit["event_observed"]:
        return fmt(hit["time_s"],1)
    return f"미관측({fmt(hit['censor_time_s'],1)})" if hit["censor_time_s"] is not None else "확인불가"


def numerical_text(n):
    if "cases" not in n:
        msg = {"not_locked":"추가 수치 계산을 위한 전체 형상 목록이 아직 봉인되지 않았다.",
               "locked_scores_incomplete":"형상 목록을 봉인했고 일부 또는 전체 설정의 재평가 결과를 기다리고 있다.",
               "numerical_audit_failed":"수치 결과의 무결성 확인에서 오류가 발생했다."}[n["status"]]
        return [msg + (" "+n["error"] if n.get("error") else "")]
    coverage = n["native_coverage"]
    text = [
        f"각 평형의 원벽과 두 방법의 600초·1800초 선택을 추가 계산 전에 봉인했다. 중복을 합친 형상·평형 {n['registered_design_case_count']}개, "
        f"cubic640/1280 및 quintic640/1280의 요청 {n['registered_setting_requests']}건에서 계산 실패·준비 실패 {n['failed_scores']}건, "
        f"계산 완료 후 부적격 {n['completed_ineligible_scores']}건이다. 원래의 기록과 선택을 재구성해 봉인 결과와 대조했다.",
        f"native cubic1280에서 기존 실제 관측 {coverage['compared_observations']}건을 대조했다. 양쪽에 모두 양의 유한값이 있는 피크·A₂ 대조는 "
        f"{coverage['compared_observables']}건이다. 관측값이 재현된 형상·평형은 {coverage['matched_observable_design_cases']}개, "
        f"적격성·관측값 유무·상대 허용오차 10⁻¹⁰ 대조 중 불일치가 있는 형상·평형은 {coverage['mismatched_design_cases']}개다. "
        f"기존 native 관측이 없는 형상·평형은 {coverage['without_prior']}개다."]
    labels = {"fast":"고속","direct":"직접","within":"폭 안","mixed":"설정별 상이","incomplete":"미확정"}
    relations = {"fast_below_direct":"고속 범위가 낮음","direct_below_fast":"직접 범위가 낮음",
                 "overlap_or_touch":"중첩·접촉","incomplete":"미확정"}
    pairs = [p for c in n["cases"] for p in c["pairs"]]
    for cutoff in (600,LONG):
        subset = [p for p in pairs if p["cutoff_s"] == cutoff]
        text += [f"### {cutoff}초 선택의 수치 민감도",
            table(["평형","시드","유한/적격 설정","같은 설정 G 범위 (%)","네 설정의 1% 판정","두 피크 범위 비교"],
                  [[p["case_id"],p["seed"],f"{p['finite_setting_count']}/{p['eligible_setting_count']}",
                    interval(p["eligible_gap_range_pct"]),labels[p["all_settings_practical"]["1"]],
                    ("동일 형상 · " if p["same_design"] else "")+relations[p["range_relation"]]] for p in subset])]
    text.append("G 범위는 같은 설정끼리 짝지은 적격 피크로 계산했다. 피크 범위 비교는 한 방법의 네 설정 최댓값이 다른 방법의 최솟값보다 낮은지 검사한다. 두 판정은 서로 다른 질문이다. 각 설정의 부적격 유한 관측도 기계 판독 결과에 남겼다. 이 네 설정의 범위에서 연속체 오차상한이나 통계적 신뢰구간은 추정할 수 없다.")
    rows = []
    for cutoff in (600,LONG):
        for key,label in (("fast","고속"),("direct","직접")):
            assessed, scheduled = [], 0
            for c in n["cases"]:
                designs = {d["id"]:d for d in c["candidates"]}
                for p in c["pairs"]:
                    if p["cutoff_s"] == cutoff:
                        scheduled += 1
                        d = designs.get(p[key+"_id"])
                        if d is not None:
                            assessed.append(d["baseline_comparison"])
            good = [d for d in assessed if d["all_settings_eligible"]]
            gains = [s["peak_reduction_pct"] for d in assessed for s in d["settings"].values() if s["paired_eligible"]]
            a2 = [s["A2_change_pct"] for d in assessed for s in d["settings"].values() if s["paired_eligible"]]
            rows.append([cutoff,label,f"{len(good)}/{scheduled}",sum(d["all_settings_paired_improvement"] for d in good),
                         sum(d["joint_improvement_settings"] == 4 for d in good),sum(d["A2_decreased_settings"] > 0 for d in assessed),
                         interval(span(gains)),interval(span(a2))])
    text += [table(["시간(s)","방법","원벽과 네 설정 적격/등록","모든 설정 피크↓","모든 설정 피크↓·A₂↑",
                    "한 설정 이상 A₂↓","피크 감소 범위(%)","A₂ 변화 범위(%)"],rows)]
    ledger_rows = []
    for name,l in n["power_ledgers_by_setting"].items():
        costs = n["costs_by_setting"][name]
        residual = lambda v: "—" if v is None else f"{v:.3g}"
        ledger_rows.append([name,f"{l['measured_ledgers']}/{l['registered']}", l["eligible_ledgers"],
                            residual(l["max_relative_closure_all"]),residual(l["max_relative_closure_eligible"]),
                            residual(l["max_unhit_fraction_eligible"]),fmt(costs["preparation_seconds"],1),fmt(costs["score_seconds"],1)])
    text.append(table(["설정","전력수지 측정/등록","적격","최대 상대 잔차(전체)","최대 상대 잔차(적격)",
                       "최대 미도달 비율(적격)","공통 준비(s)","후보 점수(s)"],ledger_rows))
    costs = n["costs"]
    text.append(f"사후 수치 계산에서 공통 준비에 {fmt(costs['preparation_seconds'],1)}초, 후보 점수에 "
                f"{fmt(costs['candidate_scoring_seconds'],1)}초가 들었다. 종료 기록이 있는 호출의 내부 경과시간은 "
                f"{fmt(costs['completed_invocation_seconds'],1)}초이며 종료 기록이 없는 호출은 {costs['unfinished_invocation_count']}개다. "
                "원벽 계산은 공통 준비 비용에 포함된다.")
    return text


def korean_report(s):
    dev = s["development"]
    text = ["# 새 DIII-D 평형에서의 제한시간 형상 최적화",
        "DIII-D 원벽에서 출발해 주어진 열원과 수송값을 유지하면서 제한시간 내에 확보한 최저 적격 정밀 피크를 비교했다. "
        "원벽과 후보는 같은 평형·같은 수치 설정에서 평가한다. 유효면적 A₂는 총 침적전력의 제곱을 면적가중 열유속 제곱 적분으로 나눈 값이다."]
    if s["live"] is not None:
        main = s["live"]["by_cutoff"]["600"]
        p = main["thresholds"]["1"]
        text.append(f"**600초 주 비교**에서 등록 {main['registered']}쌍 중 {main['comparable']}쌍에 비교 가능한 피크가 있었다. "
                    f"1% 이상 고속 우세 {p['fast']}쌍, 판정 폭 안 {p['within']}쌍, 1% 이상 직접 정밀 우세 {p['direct']}쌍이다. "
                    f"세 시드가 모두 비교 가능하고 중앙값이 고속 쪽 1% 이상인 평형은 {main['equilibria_fast_median_1pct']}/{main['distinct_equilibria']}개다. "
                    f"고속 쪽 1% 이상이면서 직접 정밀에서 실제 적격 DE trial도 완료한 경우는 {main['fast_practical_and_direct_de']}/{main['registered']}쌍이다.")
        if s["errors"]:
            text.append("**무결성 확인에서 오류가 발생해 위 집계의 판정을 보류한다.** " + " / ".join(s["errors"]))
    text += ["## 직접 정밀 대조군의 개발 선택",
        f"기존 개발 평형 세 개에서 시드 20611로 24/8·32/12(초기 집단/투영 시드)를 600초씩 비교했다. "
        f"등록 {dev['registered']}회 중 종료 기록 {dev['completed_receipts']}회의 부모 경과시간은 {fmt(dev['parent_seconds'],1)}초다."]
    if dev["configurations"]:
        text.append(table(["집단/투영","적격 DE를 완료한 개발 평형","원벽 정규화 피크 기하평균","사전 기준 충족"],
            [[f"{c['population']}/{c['projected_seeds']}",f"{c['de_cases']}/3",fmt(c["geometric_mean_peak_ratio"],5),
              "예" if c["qualified"] else "아니오"] for c in dev["configurations"]]))
    if dev["selected"]:
        chosen = dev["selected"]
        text.append(f"봉인한 직접 정밀 설정은 {chosen['population']}/{chosen['projected_seeds']}이다. 고속 선별은 생산 설정 64/12를 사용한다.")
    if not s["live"]:
        text += ["## 새 평형 전향 시험 미등록",
                 "등록된 새 평형 열계산 및 방법 간 비교 결과가 없다. 개발 단계 상태는 "+s["status"]+"이다."]
        if dev["reasons"]:
            text.extend(dev["reasons"])
        return "\n\n".join(text)+"\n"
    live = s["live"]
    counts = live["status_counts"]
    text += ["## 입력, 시간과 판정",
        f"기하 사전검사에서 {s['input_screening']['examined']}개를 조사해 {len(s['cases'])}개를 채택했다. OMP 1 MW, "
        "n = 10¹⁹ m⁻³, χ⊥ = 0.3 m²/s, λq = 3 mm, γ = 7, κ₀ = 2000과 최소 곡률반경 20 mm를 고정했다. "
        "각 평형의 편집 가능한 OMP 계수는 16개다. 고속은 저비용 cubic320에서 선별한 형상을 cubic1280으로 확인하고 직접 방법은 cubic1280을 쓴다.",
        f"등록 {live['registered_arms']}회, 부모 종료 기록 {live['with_process_receipt']}회, 감사 가능한 종결 상태 {live['closed_arms']}회다. "
        "상태별 실행 수는 "+", ".join(f"{k} {v}" for k,v in counts.items())+"다. "
        "동일 시드의 두 방법을 짝지었고, 시드 20711의 1800초 실행에서 600초 결과를 함께 재구성했다. 시점과 시드는 같은 평형 안의 반복 관측이다.",
        "G = 100 × (q직접 − q고속) / q직접이다. G ≥ 1%면 고속 쪽, G ≤ −1%면 직접 쪽의 실용적 차이다. "
        "1%는 이 시험에서 미리 정한 공학적 판정 폭이다. 표의 판정에는 반올림 전 값을 사용했고 원벽 대비 피크·A₂ 증감은 엄격한 부호로 집계했다. "
        "준비가 끝난 원벽과 그 시점까지 완료한 새 적격 정밀 후보 중 피크가 가장 낮은 형상을 선택했다."]
    base_rows = []
    for c in s["cases"]:
        baseline = next((p[key]["baseline"] for p in live["pairs"] if p["case_id"] == c["case_id"] for key in ("fast","direct")
                         if p[key]["baseline"]),None)
        base_rows.append([c["case_id"],fmt(c["time_ms"],1),fmt(baseline["peak"],5) if baseline else "—",
                          fmt(baseline["A2"],5) if baseline else "—"])
    text.append(table(["평형","시각(ms)","native 원벽 피크(MW/m²)","native 원벽 A₂(m²)"],base_rows))
    for cutoff,title in ((600,"주 평가"),(LONG,"첫 시드 장시간 평가")):
        selected = [p for p in live["pairs"] if p["cutoff_s"] == cutoff]
        text.append(f"## {cutoff}초 {title}")
        rows = []
        for p in selected:
            a,b = p["fast"],p["direct"]
            rows.append([p["case_id"],p["seed"],fmt(a["q_ratio"],4),fmt(b["q_ratio"],4),
                fmt(p["margin"]["gap_pct"],3) if p["margin"] else "—",
                f"{fmt(a['A2_change_pct'],2)} / {fmt(b['A2_change_pct'],2)}",
                f"{hit_text(a)} / {hit_text(b)}",
                "—" if b["de"] is None else f"{b['de']['fresh']} / {b['de']['deposition']} / {b['de']['eligible_deposition']}"])
        text.append(table(["평형","시드","고속 q/q₀","직접 q/q₀","G (%)","A₂ 변화 고속/직접 (%)","최초 5% 고속/직접 (s)",
                           "직접 DE 새 점수/침적/적격"],rows))
        text.append("미관측 옆 괄호는 그 실행을 확인할 수 있었던 시각(s)이다. DE는 초기 집단 평가 이후의 trial이며, 새 침적 횟수에는 캐시 호출이 들어가지 않는다.")
        method_table = []
        for key,label in (("fast","고속"),("direct","직접")):
            m = live["method_baselines"][str(cutoff)][key]
            method_table.append([label,f"{m['available']}/{m['registered']}",m["improved_q"],m["q_down_A2_up"],m["A2_down"],
                                 m["fallback_original_wall"],m["no_eligible_new_fine"],m["first_5pct_observed"],interval(m["peak_reduction_range_pct"])])
        text.append(table(["방법","채택 가능/등록","원벽 대비 피크↓","피크↓·A₂↑","A₂↓","원벽 유지","새 적격 정밀 0","5% 도달","피크 감소 범위 (%)"],method_table))
    text += ["## 시점별 방향과 판정 폭", table(["시간(s)","비교 가능/등록","부호 고속/동률/직접","0.5% 고속/폭 안/직접",
        "1% 고속/폭 안/직접","2% 고속/폭 안/직접","직접 적격 DE가 있는 쌍"],
        [[t,f"{c['comparable']}/{c['registered']}", " / ".join(str(c["signs"][k]) for k in ("fast","tie","direct")),
          *[" / ".join(str(c["thresholds"][th][k]) for k in ("fast","within","direct")) for th in ("0.5","1","2")],
          c["direct_eligible_de_pairs"]] for t,c in live["by_cutoff"].items()])]
    table_rows = []
    for g in live["case_groups"]:
        geom = g["geometry"]
        table_rows.append([g["case_id"],g["cutoff_s"],f"{g['comparable']}/{g['registered']}",fmt(g["gap_seed_median_pct"]),
            interval(g["gap_seed_range_pct"]),f"{geom['unique_identities']}/{geom['selection_count']}"])
    text += [table(["평형","시간(s)","비교 가능/등록","시드 중앙 G (%)","시드별 G 범위 (%)","고유 형상/채택 수"],table_rows),
        "고유 형상은 같은 평형·같은 시점에서 동일한 계수 기록 또는 동일한 원벽 선택을 묶은 수다. 서로 다른 계수 기록 사이의 최대 계수 차이(m)도 저장했다. 초기 집단의 16차원 rank는 초기화 검사에 해당한다. 실제 적격 후보가 16개 방향을 고르게 탐색했다는 결론은 여기서 얻을 수 없다."]
    if any(a["status"] != "completed" for a in live["parent_costs"]):
        text += ["### 종결 예외와 미완료",table(["평형","시드","방법","상태","관측된 사유"],
            [[a["case_id"],a["seed"],a["policy"],a["status"],(a["error"] or "—").replace("|","/")]
             for a in live["parent_costs"] if a["status"] != "completed"]),
            "인증된 기록에서 실패 전에 충분히 관측한 시점의 값은 남겼다. 실패 이후 시점과 최종 기록 검증이 불가능한 실행은 비교 불가로 표시했다."]
    text.append("## 실행 비용과 전력수지")
    ready = [a["baseline_ready_s"] for a in live["parent_costs"] if a["baseline_ready_s"] is not None]
    text.append(f"새 평형의 부모 실행시간 합계는 {fmt(live['parent_seconds'],1)}초다. 원벽 준비 완료는 프로세스 시작 후 {interval(span(ready),1)}초였으며 "
                f"등록 시한을 넘겨 끝난 coarse·fine 점수 이벤트는 {sum(r['late_scores'] for r in live['parent_costs'])}건이다. "
                "부모가 실행을 시작한 시점부터 최적화 비용을 측정했다. 원자료 취득, 기하 사전 검사, 부모의 실행 전 해시 확인은 이 시간 합계에 포함되지 않는다.")
    overhead = []
    for policy,label in ((FAST,"고속"),(DIRECT,"직접")):
        available_costs = [a for a in live["parent_costs"] if a["policy"] == policy]
        def stage_sum(name):
            return sum(a["stage_seconds"].get(name,0.) for a in available_costs)
        overhead.append([label,fmt(sum(a["parent_elapsed_s"] for a in available_costs if a["parent_elapsed_s"] is not None),1),
                         fmt(sum(a["startup_s"] for a in available_costs if a["startup_s"] is not None),1),
                         fmt(stage_sum("prepare_fine"),1),fmt(stage_sum("prepare_coarse"),1),fmt(stage_sum("initialize"),1)])
    text.append(table(["방법","부모 경과 합계(s)","기록된 시작 지연(s)","정밀 준비(s)","저비용 준비(s)","집단 생성(s)"],overhead))
    stages = []
    for key,label in (("fast","고속"),("direct","직접")):
        # One final registered horizon per real process.
        endpoints = [p[key] for p in live["pairs"] if p["cutoff_s"] == (LONG if p["seed"] == SEEDS[0] else 600)]
        for phase in PHASES:
            phase_rows = [e["phases"][phase] for e in endpoints if e["phases"]]
            stages.append([label,phase,sum(v["fresh"] for v in phase_rows),sum(v["geometry_rejections"] for v in phase_rows),
                           sum(v["deposition"] for v in phase_rows),
                           sum(v["eligible_deposition"] for v in phase_rows),fmt(sum(v["completed_seconds"] for v in phase_rows),1)])
    text.append(table(["방법","점수 단계","새 점수","침적 전 기하 탈락","반환 침적","적격 침적","시한 내 완료 이벤트(s)"],stages))
    power_rows = []
    for t in (600,LONG):
        for key,label in (("fast","고속"),("direct","직접")):
            m = live["method_baselines"][str(t)][key]
            power_rows.append([t,label,"—" if m["max_power_closure"] is None else f"{m['max_power_closure']:.3g}",
                               "—" if m["max_unhit_fraction"] is None else f"{m['max_unhit_fraction']:.3g}",
                               interval(m["deposited_over_baseline_range"],9)])
    text += [table(["시간(s)","방법","채택 형상 최대 상대 전력수지 잔차","최대 미도달 전력 비율","채택/원벽 침적전력 범위"],power_rows),
        "전력수지 잔차는 |침적 + 미도달 − 배정| / 배정이다. 보존 점검만으로 국소 도달면을 검증할 수 없다.",
        "## 봉인한 형상의 네 수치 설정 재평가", *numerical_text(s["numerical"]),
        "이 비교의 범위는 세 익명 평형에서 사전 지정한 OMP 운전점과 두 탐색 설정의 시간 내 성능이다. "
        "독립 FIREFLY 수송식에 대한 이전 검증과 이번 고정 수송식의 수치 민감도는 별도 증거다. 이번에 새 FIREFLY 결과나 실험 열유속 예측 정확도를 측정하지 않았다. "
        "[사전 규칙](protocol.ko.md) · [시점별 원값·형상·해시](summary.json)"]
    return "\n\n".join(text)+"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True,help="exclusive publication folder on D; use another folder for interim reports")
    args = parser.parse_args()
    study, out = args.study.resolve(), args.output.resolve()
    require(study.drive.lower() == out.drive.lower() == "d:", "study and publication must stay on D")
    require((study/"development-plan.json").is_file(), "prospective development registration missing")
    require(not out.is_relative_to(study) and not study.is_relative_to(out),
            "publish outside raw study and its parent")
    from d3d_optimization.campaign_state import process_lock
    with process_lock(study/"campaign.lock"), process_lock(study/"posttrial.lock"):
        summary = build_summary(study)
    raw = (json.dumps(summary,sort_keys=True,ensure_ascii=False,allow_nan=False,indent=2)+"\n").encode("utf8")
    publication = {"summary.json":raw, "REPORT.ko.md":korean_report(summary).encode("utf8")}
    for name,data in publication.items():
        require(not (out/name).exists() or (out/name).read_bytes() == data, f"different publication exists at {out/name}")
    for name,data in publication.items():
        write_once(out/name,data)
    print(summary["status"],out/"REPORT.ko.md",flush=True)
    if summary["status"] not in ("closed","closed_with_exceptions","development_unqualified"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
