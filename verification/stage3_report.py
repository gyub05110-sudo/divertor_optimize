"""Audit prospectively scheduled Stage 3 searches and reconstruct their paired outcomes."""
import argparse
from collections import Counter, defaultdict
import csv
from itertools import product
import json
import math
from pathlib import Path
from statistics import median

import numpy as np

from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.campaign_state import public_save
from d3d_optimization.user_output import file_hash
from optimize_divertor import COARSE, FINE
from verification.user_shape_replay import INPUT_KEYS


SEEDS = (20711, 20729, 20747)
ARMS = ("floor-screened", "full-screened", "full-direct")
CONFIG = dict(zip(ARMS, (("floor-only","screened"), ("full","screened"), ("full","direct"))))
CUTOFFS = (200, 400, 600)
NONFLOOR = ("omp_entry", "shelf", "omp_connector")
REQUIRED_FILES = {"geometry.csv", "profile.csv.gz", "donors.csv.gz", "boundary.npz"}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def close(a, b):
    return math.isclose(a, b, rel_tol=2e-14, abs_tol=1e-11)


def positive(x):
    return isinstance(x, (int,float)) and math.isfinite(x) and x > 0


def checkpoint(admitted, cutoff, reference_q, power):
    scores = [e for e in admitted if e["end_s"] <= cutoff]
    best = min(scores, key=lambda e:e["peak_q"]) if scores else None
    q = best["peak_q"] if best else None
    return dict(id=best["id"] if best else None, q_native=q,
                peak_MW_m2=q*power/1e6 if best else None, completion_s=best["end_s"] if best else None,
                eligible_fine_count=len(scores),
                reduction_pct=100*(reference_q-q)/reference_q if best and reference_q is not None else None)


def time_to_threshold(admitted, cutoff, reference_q, censorable, horizon, reason):
    first = next((e for e in admitted if e["end_s"] <= cutoff and reference_q is not None
                  and e["peak_q"] <= .95 * reference_q), None)
    if first:
        return dict(observed_s=first["end_s"], right_censored=False, censor_s=None, reason="observed")
    return dict(observed_s=None, right_censored=censorable,
                censor_s=min(horizon,cutoff) if censorable else None, reason=reason)


def no_run(status, detail):
    first=time_to_threshold([],600,None,False,None,status)
    return dict(status=status, reason=detail, candidate=None, has_nonfloor_coefficients=False,
                nonfloor_max_displacement_m=None, baseline_peak_MW_m2=None, baseline_measurable=False,
                checkpoints={str(t):dict(checkpoint([],t,None,None),availability=status,first_5pct=first) for t in CUTOFFS},
                first_5pct=first)


def inspect_events(events, data):
    args, timing = data["arguments"], data["timings"]
    method, radius, budget = args["method"], args["radius_mm"], args["budget_seconds"]
    proposals, coarse, fine_seen = {}, {}, set()
    counts = Counter()
    admitted = []
    previous_end = 0.
    prep = []
    phases = Counter()
    errors = Counter()
    for e in events:
        stage, status = e["stage"], e.get("status")
        if "end_s" in e:
            start, end = e["start_s"], e["end_s"]
            require(math.isfinite(end) and 0 <= previous_end <= start <= end,
                    f"out-of-order clock for {stage}")
            require(close(e["duration_s"], end-start), f"{stage} duration changed")
            timely = end <= budget
            require(e["within_budget"] is timely, f"{stage} deadline flag contradicts clock")
            previous_end = end
        elif "elapsed_s" in e:
            require(previous_end <= e["elapsed_s"] <= timing["search_finished_s"], f"invalid {stage} elapsed time")
            previous_end = e["elapsed_s"]
        if status == "finished" and stage not in ("proposal","coarse","fine","initialize"):
            prep.append({k:v for k,v in e.items() if k!="status"})
        if stage=="initialize":
            require(method=="direct", "unexpected direct population")
            counts["late"] += int(not timely)
            counts["failed"] += int(status=="failed")
            if status=="completed":
                pop=np.asarray(e["vectors"])
                require(pop.shape==(32,16) and np.isfinite(pop).all() and
                        e["rank"]==np.linalg.matrix_rank((pop-pop.mean(axis=0))/.2, tol=1e-12),
                        "invalid direct population rank")
                require(e["restart"]!=0 or e["rank"]==16, "initial direct population lacks full rank")
        if stage not in ("proposal","coarse","fine"):
            continue
        require(e["start_s"] < budget and status in ("completed","failed","rejected","late"),
                f"invalid {stage} admission")
        require(status=="failed" or (status=="late") != timely, f"inconsistent {stage} completion status")
        if stage in ("coarse","fine") and not timely:
            counts["late_evaluations"] += 1
        counts["late"] += int(not timely)
        counts["failed"] += int(status=="failed")
        errors.update(e.get("errors",[]))
        index, ident = e["index"], e.get("id")
        if e.get("candidate") is not None:
            cand = Candidate.from_record(e["candidate"])
            require(cand.identity==ident and cand.source=="OMP" and cand.radius_mm==radius, "candidate identity mismatch")
            require(args["family"]!="floor-only" or not any(cand.vector[:12]), "floor-only candidate changed another surface")
        if stage=="proposal":
            require(index==counts["proposals"],"nonsequential proposal index")
            counts["proposals"] += 1
            if timely and status=="completed":
                require(ident is not None and e.get("candidate") is not None, "completed proposal without identity")
                proposals[index] = ident
            continue
        require(ident is not None and e.get("candidate") is not None and proposals.get(index)==ident,
                "score not bound to a completed proposal")
        counts[stage] += 1
        valid = timely and status=="completed" and e.get("eligible") is True
        if valid:
            require(positive(e["peak_q"]), "invalid admitted peak")
            require(not any(e.get(k) is False for k in ("geometry_valid","shape_passed","whole_wall_passed")),
                    "eligible score failed a recorded geometry constraint")
            counts["eligible_"+stage] += 1
        if stage=="coarse":
            require(method=="screened" and ident not in coarse, "unexpected coarse evaluation")
            coarse[ident] = e["peak_q"] if valid else None
        else:
            require(ident not in fine_seen, "duplicate fine score")
            fine_seen.add(ident)
            if method=="screened":
                require(coarse.get(ident) is not None and e.get("coarse_peak_q",coarse[ident])==coarse[ident],
                        "fine score lacks eligible coarse evidence")
            else:
                require(not e.get("cached",False), "cached fine score admitted")
            phases[(e.get("proposal_kind","unknown"),"attempted")] += 1
            if valid:
                phases[(e.get("proposal_kind","unknown"),"admitted")] += 1
            improves = valid and (not admitted or e["peak_q"] < min(r["peak_q"] for r in admitted))
            require(e.get("incumbent_updated",False)==improves, "incumbent flag disagrees with raw fine scores")
            if valid:
                admitted.append(e)
    require(previous_end <= timing["search_finished_s"], "search finished before final event")
    require(prep==timing.get("preparation",[]), "preparation clocks disagree")
    counts["eligible_fine"] += 0
    counts["eligible_coarse"] += 0
    for k in ("proposals","coarse","fine","eligible_fine","eligible_coarse","failed","late"):
        if "search" in data:
            require(counts[k]==data["search"]["counts"][k], f"search count mismatch for {k}")
    if "search" in data:
        require(data["search"]["first_feasible_elapsed_s"]==(admitted[0]["end_s"] if admitted else None),
                "first feasible time disagrees with fine journal")
    return admitted, counts, {p:{s:phases[p,s] for s in ("attempted","admitted")} for p,_ in phases}, dict(errors)


def audit_run(data, events):
    args = data["arguments"]
    budget, timing, status = args["budget_seconds"], data["timings"], data["status"]
    require(budget==600 and args["method"] in ("screened","direct"), "unregistered run settings")
    require(tuple(data["levels"]["fine"])==FINE, "native fine level differs")
    if args["method"]=="screened":
        require(tuple(data["levels"]["coarse"])==COARSE, "screen level differs")
    else:
        require(not any(e["stage"]=="prepare_coarse_search" for e in events), "direct arm prepared coarse search")
    duration, write_s = timing["search_finished_s"], timing["elapsed_before_result_write_s"]
    require(0 <= duration <= write_s and close(timing["deadline_overrun_s"], max(0.,duration-budget)), "invalid process clock")
    admitted, counts, phases, errors = inspect_events(events,data)
    power = data.get("reference_SI",{}).get("power_W")
    base = data.get("baseline") or {}
    ready = timing.get("baseline_ready_s")
    base_q = None
    if base.get("physics_measurable"):
        require(positive(power) and positive(base.get("peak_MW_m2")) and base.get("reference_only") is True
                and ready is not None and 0 <= ready <= duration and base["within_budget"] is (ready <= budget),
                "invalid original-wall measurement or baseline time")
        if ready<=budget:
            base_q = base["peak_MW_m2"]/(power/1e6)
    require(not admitted or (base_q is not None and ready<=admitted[0]["start_s"]), "fine score without on-time reference")
    candidate, selected = data.get("candidate"), data.get("selected")
    require(bool(admitted)==bool(candidate)==bool(selected), "selected geometry lacks matching admitted fine record")
    require(status!="qualified" or candidate is not None, "qualified run without a winner")
    if admitted:
        best = min(admitted,key=lambda r:r["peak_q"])
        winner = Candidate.from_record(candidate)
        require(winner.identity==candidate["id"]==best["id"] and winner.record()==best["candidate"], "wrong strict minimum")
        require(selected["candidate_eligible"] is True and selected["physics_measurable"] is True
                and close(selected["peak_MW_m2"],best["peak_q"]*power/1e6), "selected export contradicts fine score")
        require(status in ("qualified","partial_result"), "winner in failed result")
    horizon = min(duration, budget)
    censorable = base_q is not None and data.get("fixed_shelf",{}).get("admissible") is True
    reason = ("deadline_without_5pct" if censorable and horizon==budget else
              "early_stop_without_5pct" if censorable else status if status!="qualified" else "reference_unavailable")
    cuts = {str(t):checkpoint(admitted,t,base_q,power) for t in CUTOFFS}
    for t in CUTOFFS:
        snap=cuts[str(t)]
        if snap["id"]:
            state="eligible_fine"
        elif duration<=t and status in ("constraint_excluded","calculation_error","input_error","partial_result"):
            state=status
        elif ready is not None and ready>t:
            state="reference_pending"
        elif base_q is None:
            state="reference_unavailable"
        else:
            state="reference_only"
        suffix="checkpoint_without_5pct" if t<budget and t<=horizon and censorable else reason
        snap.update(availability=state,first_5pct=time_to_threshold(admitted,t,base_q,censorable,horizon,suffix))
    if selected is not None and "reduction_vs_reference_pct" in selected:
        require(close(cuts["600"]["reduction_pct"],selected["reduction_vs_reference_pct"]), "reduction mismatch")
    return dict(status=status, reason=data.get("reason",status), candidate=candidate,
                has_nonfloor_coefficients=bool(candidate and any(candidate["vector"][:12])),
                nonfloor_max_displacement_m=None, baseline_peak_MW_m2=base.get("peak_MW_m2"),
                baseline_measurable=bool(base.get("physics_measurable")),
                baseline_ready_s=ready, baseline_A2_m2=base.get("A2_m2"), selected_A2_m2=(selected or {}).get("A2_m2"),
                baseline_ledger=base.get("ledger"), selected_ledger=(selected or {}).get("ledger"),
                search_finished_s=duration, elapsed_before_result_write_s=write_s, deadline_overrun_s=timing["deadline_overrun_s"],
                fixed_shelf_admissible=data.get("fixed_shelf",{}).get("admissible"),
                checkpoints=cuts, first_5pct=time_to_threshold(admitted,600,base_q,censorable,horizon,reason),
                late_evaluations=counts["late_evaluations"], journal_counts=dict(counts), fine_phases=phases,
                evaluation_errors=errors)


def paired(rows):
    groups = defaultdict(dict)
    for row in rows:
        group=groups[(row["case_id"],row["seed"])]
        require(row["arm"] not in group, "duplicate arm in matched set")
        group[row["arm"]]=row
    output=[]
    for (case,seed),group in groups.items():
        require(set(group)==set(ARMS), "incomplete arm inventory")
        for arm,comparator in (("full-screened","floor-screened"), ("full-direct","floor-screened"),
                               ("full-screened","full-direct")):
            a,b=group[arm],group[comparator]
            if a["baseline_peak_MW_m2"] is not None and b["baseline_peak_MW_m2"] is not None:
                require(close(a["baseline_peak_MW_m2"],b["baseline_peak_MW_m2"]), "unmatched paired references")
            for t in CUTOFFS:
                x,y=a["checkpoints"][str(t)]["q_native"],b["checkpoints"][str(t)]["q_native"]
                assessed=x is not None and y is not None and a["status"]=="qualified" and b["status"]=="qualified"
                outcome=("unavailable" if not assessed else "arm_lower" if x<y else "comparator_lower" if x>y else "tie")
                moved = a.get("nonfloor_max_displacement_m") or 0.
                nonfloor_changed = bool(a["has_nonfloor_coefficients"] and moved>1e-12)
                trigger = comparator=="floor-screened" and t==600 and outcome=="arm_lower"
                geometry_contrast = comparator=="floor-screened"
                gap = (100*(y-x)/(x if geometry_contrast else y)) if assessed else None
                practical = ("unavailable" if gap is None else "within_1pct" if abs(gap)<=1. else
                             "arm_lower_over_1pct" if gap>1. else "comparator_lower_over_1pct")
                output.append(dict(case_id=case,seed=seed,cutoff_s=t,arm=arm,comparator=comparator,
                                   arm_status=a["status"],comparator_status=b["status"],arm_q_native=x,comparator_q_native=y,
                                   signed_outcome=outcome, registered_gap_pct=gap, practical_outcome=practical,
                                   contrast_scope=("primary_family" if arm=="full-screened" and geometry_contrast else
                                                   "secondary_joint_family_search" if geometry_contrast else "search_method"),
                                   nonfloor_changed=nonfloor_changed if t==600 and x is not None else None,
                                   replay_trigger=trigger))
    return output


def per_shot(rows):
    groups=defaultdict(list)
    for r in rows:
        groups[(r["case_id"],r["arm"])].append(r)
    summary=[]
    for (case,arm), group in groups.items():
        qualified=[r for r in group if r["status"]=="qualified"]
        reductions=[r["checkpoints"]["600"]["reduction_pct"] for r in qualified]
        hits=[r["first_5pct"]["observed_s"] for r in qualified if r["first_5pct"]["observed_s"] is not None]
        summary.append(dict(case_id=case,arm=arm,scheduled_seeds=len(group),qualified_seeds=len(qualified),
                            reached_5pct_seeds=len(hits), improved_seeds=sum(x>0 for x in reductions),
                            nonimproving_qualified_seeds=sum(x<=0 for x in reductions),
                            right_censored_seeds=sum(r["first_5pct"]["right_censored"] for r in group),
                            qualified_median_reduction_pct=median(reductions) if reductions else None,
                            observed_only_median_first_5pct_s=median(hits) if hits else None,
                            statuses=dict(Counter(r["status"] for r in group))))
    return summary



def by_arm(shot_rows, rows):
    output={}
    for arm in ARMS:
        cases=[s for s in shot_rows if s["arm"]==arm]
        runs=[r for r in rows if r["arm"]==arm]
        n=len(runs)
        improved=sum(s["improved_seeds"] for s in cases)
        threshold=sum(s["reached_5pct_seeds"] for s in cases)
        output[arm]=dict(scheduled_shots=len(cases), scheduled_shot_seed_pairs=n,
                         qualified_pairs=sum(s["qualified_seeds"] for s in cases),
                         improved_pairs=improved, end_to_end_improved_pct=100*improved/n,
                         reached_5pct_pairs=threshold, end_to_end_5pct_pct=100*threshold/n,
                         shots_improving_any_seed=sum(s["improved_seeds"]>0 for s in cases),
                         shots_improving_all_seeds=sum(s["improved_seeds"]==s["scheduled_seeds"] for s in cases),
                         shots_5pct_any_seed=sum(s["reached_5pct_seeds"]>0 for s in cases),
                         shots_5pct_all_seeds=sum(s["reached_5pct_seeds"]==s["scheduled_seeds"] for s in cases),
                         right_censored_pairs=sum(s["right_censored_seeds"] for s in cases),
                         outcome_counts=dict(Counter(r["status"] for r in runs)))
    return output


def check_exports(folder, data, row):
    """Check material coordinates, energy/area/peak reconstructions, and the selected vector."""
    from verification.user_shape_results import audit_measurement, verify_geometry
    power=data["reference_SI"]["power_W"] if data.get("reference_SI") else None
    profiles={}
    for source in ("baseline","selected"):
        measurement=data.get(source)
        if measurement and "sha256" in measurement:
            require(set(measurement["sha256"])==REQUIRED_FILES, "missing measurement artifacts")
            if measurement.get("physics_measurable") is False:
                require(source=="baseline", "nonmeasurable selected candidate")
                for name,digest in measurement["sha256"].items():
                    require(file_hash(folder/source/name)==digest, "failed baseline artifact digest mismatch")
            else:
                profiles[source]=audit_measurement(folder/source, measurement,power)
        elif measurement and measurement.get("physics_measurable"):
            raise ValueError(f"missing {source} measurement inventory")
    row["nonfloor_max_displacement_m"]=None
    if "selected" in profiles:
        original,wall=profiles["baseline"],profiles["selected"]
        verify_geometry(original,wall,Candidate.from_record(data["candidate"]))
        shifted=np.linalg.norm(wall.wall-original.wall,axis=1)
        offsets={}
        for name in NONFLOOR:
            part=original.regions==name
            vertices=np.r_[part,False] | np.r_[False,part]
            offsets[name]=float(shifted[vertices].max())
        row["nonfloor_displacements_m"]=offsets
        row["nonfloor_max_displacement_m"]=max(offsets.values())


def verify_source(result, task, case, plan):
    args=result["arguments"]
    require(all(args[k]==v for k,v in plan["settings"].items()), "CLI settings differ from prospective seal")
    require((args["family"],args["method"])==CONFIG[task["arm"]] and
            all(args[k]==task[k] for k in ("family","method","seed","budget_seconds","time_ms")),
            "CLI arguments differ from scheduled task")
    require(Path(args["equilibrium"]).resolve()==Path(case["path"]).resolve(), "different equilibrium path")
    if result.get("inputs"):
        inputs=result["inputs"]
        require(inputs["equilibrium_sha256"]==case["sha256"] and inputs["time_ms"]==case["time_ms"],
                "wrong registered discharge or time")
        require(result["reference_SI"]==dict(power_W=args["power_mw"]*1e6,density_m3=args["density_m3"],
                                             chi_perp_m2_s=args["chi_m2_s"],lambda_q_m=args["lambda_q_mm"]/1000,
                                             gamma=args["gamma"],kappa0=args["kappa0"]),
                "operating point differs from registered arguments")
        for key,path in (("geometry_sha256","d3d_optimization/data/DIIID_geom.json"),
                         ("design_policy_sha256","baseline_review/optimization_region_policy.json"),
                         ("design_reference_sha256","baseline_review/baseline_review.json")):
            require(inputs[key]==plan["files"][path], "immutable input mismatch for "+key)
    for module,digest in result["software"]["code_sha256"].items():
        location=module if module in plan["files"] else "d3d_optimization/"+module
        require(plan["files"].get(location)==digest,"search code differs from seal for "+module)
    expected_modules={Path(k).name for k in plan["files"] if k.startswith("d3d_optimization/") and k.endswith(".py")}
    expected_modules.update(("optimize_divertor.py","run_tangent_curvature_study.py"))
    require(set(result["software"]["code_sha256"])==expected_modules, "incomplete solver code inventory")


def pair_by_shot(pairs):
    groups=defaultdict(list)
    for p in pairs:
        groups[(p["case_id"],p["arm"],p["comparator"],p["cutoff_s"])].append(p)
    summary=[]
    for (case,arm,cmp,cut),group in groups.items():
        gaps=[p["registered_gap_pct"] for p in group if p["registered_gap_pct"] is not None]
        summary.append(dict(case_id=case,arm=arm,comparator=cmp,cutoff_s=cut,scheduled_pairs=len(group),
                            comparable_pairs=len(gaps),median_comparable_gap_pct=median(gaps) if gaps else None,
                            signed_outcomes=dict(Counter(p["signed_outcome"] for p in group)),
                            practical_outcomes=dict(Counter(p["practical_outcome"] for p in group))))
    return summary


def flush_csv(path,rows):
    target=path.with_suffix(".tmp")
    with target.open("w",newline="",encoding="utf-8") as fp:
        writer=csv.DictWriter(fp,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    target.replace(path)


def case_flow(cases,rows):
    summary=[]
    for case in cases.values():
        runs=[r for r in rows if r["case_id"]==case["id"]]
        qualified=[r for r in runs if r["status"]=="qualified"]
        summary.append(dict(case_id=case["id"],rank=case["rank"],source_sha256=case["sha256"],
                            verified_source=True,registration_status=case["status"],registration_reason=case.get("reason"),
                            unique_registered_frame=case["status"]=="ready",
                            measured_baseline_any=any(r["baseline_measurable"] for r in runs),
                            eligible_replacement_any=bool(qualified),
                            improved_replacement_any=any(r["checkpoints"]["600"]["reduction_pct"]>0 for r in qualified),
                            statuses=dict(Counter(r["status"] for r in runs))))
    return summary


def check_process(receipt, metrics):
    if receipt.get("timing_source")=="recovered_cli_result":
        require(receipt["elapsed_s"] is None and receipt["returncode"] is None,
                "recovered CLI supplied unobserved parent timing or exit code")
        return
    elapsed=receipt["elapsed_s"]
    require(positive(elapsed) and elapsed>=metrics["elapsed_before_result_write_s"], "invalid parent process time")
    require(receipt["returncode"]=={"partial_result":1,"calculation_error":1,"input_error":2}.get(metrics["status"],0),
            "return code differs from optimizer outcome")


def run(study, allow_incomplete=False):
    from verification.stage3_campaign import load_seal, verified_receipt
    study=Path(study).resolve()
    plan=load_seal(study)
    cases={c["id"]:c for c in read(study/"inputs"/"registration.json")["cases"]}
    require(len(cases)==12 and len({c["sha256"] for c in cases.values()})==12, "duplicate prospective sources")
    for case in cases.values():
        require(file_hash(case["path"])==case["sha256"], "changed registered equilibrium")
    scheduled={(t["case_id"],t["seed"],t["arm"]) for t in plan["tasks"]}
    require(len(plan["tasks"])==108 and scheduled==set(product(cases,SEEDS,ARMS)), "incomplete prospective inventory")
    rows=[]
    signatures={}
    for task in plan["tasks"]:
        case,seed,arm=task["case_id"],task["seed"],task["arm"]
        folder=study/task["relative_output"]
        common=dict(case_id=case,seed=seed,arm=arm,relative_output=task["relative_output"],
                    registration_status=cases[case]["status"],receipt_status=None,process_elapsed_s=None)
        if not (folder/"receipt.json").exists():
            require(allow_incomplete, "missing terminal receipt at "+str(folder))
            rows.append(dict(common,**no_run("pending","No terminal receipt")))
            continue
        receipt=verified_receipt(folder,task,plan)
        common.update(receipt_status=receipt["status"],result_sha256=receipt["result_sha256"],
                      events_sha256=receipt["events_sha256"],receipt_sha256=file_hash(folder/"receipt.json"),
                      process_elapsed_s=receipt["elapsed_s"],returncode=receipt["returncode"],
                      timing_source=receipt.get("timing_source"))
        if receipt["status"]!="completed":
            metrics=no_run(receipt["status"],receipt["detail"])
        else:
            result=read(folder/"result.json")
            verify_source(result,task,cases[case],plan)
            events=[json.loads(line) for line in (folder/"events.jsonl").read_text(encoding="utf8").splitlines() if line.strip()]
            metrics=audit_run(result,events)
            check_process(receipt,metrics)
            check_exports(folder,result,metrics)
            if result.get("inputs"):
                signature=([result["inputs"][k] for k in INPUT_KEYS],result["reference_SI"])
                require(signatures.setdefault(case,signature)==signature, "unmatched source, geometry or operating point")
        rows.append(dict(common,**metrics))
    pairs=paired(rows)
    grouped=per_shot(rows)
    complete=all(r["status"]!="pending" for r in rows)
    summary=dict(schema="nf-stage3-report.v1",complete=complete,
                 scope="DIII-D public source files at 1400 ms, fixed reduced transport, 600 s cold search",
                 expected_sha256=file_hash(study/"expected.json"),registration_sha256=plan["registration_sha256"],
                 report_driver_sha256=file_hash(__file__),runs=rows,case_flow=case_flow(cases,rows),
                 per_shot=grouped,by_arm=by_arm(grouped,rows),paired=pairs,paired_per_shot=pair_by_shot(pairs),
                 definitions=dict(reduction_pct="100*(reference-native_fine)/reference; signed and uncapped",
                                  family_gap_pct="100*(q_floor/q_full-1); full-direct vs floor is secondary",
                                  method_gap_pct="100*(1-q_full_screened/q_full_direct)",
                                  threshold="native eligible q <= 0.95*reference, from process start",
                                  right_censor="search observed with measurable on-time baseline and admissible shelf",
                                  nonfloor_changed="nonzero first twelve coefficients and measured vertex displacement >1e-12 m",
                                  practical_margin_pct=1,unit="one anonymized source file; seeds paired within each file"))
    flat=[]
    for r in rows:
        for cut,snap in r["checkpoints"].items():
            flat.append(dict(case_id=r["case_id"],seed=r["seed"],arm=r["arm"],cutoff_s=cut,status=r["status"],
                             reason=r["reason"],process_elapsed_s=r["process_elapsed_s"],baseline_ready_s=r.get("baseline_ready_s"),
                             baseline_peak_MW_m2=r["baseline_peak_MW_m2"],**{k:v for k,v in snap.items() if k!="first_5pct"},
                             first_5pct_observed_s=snap["first_5pct"]["observed_s"],
                             first_5pct_right_censored=snap["first_5pct"]["right_censored"],
                             first_5pct_censor_s=snap["first_5pct"]["censor_s"],first_5pct_reason=snap["first_5pct"]["reason"]))
    output=study/"interim" if allow_incomplete else study
    output.mkdir(parents=True,exist_ok=True)
    flush_csv(output/"stage3_runs.csv",flat)
    flush_csv(output/"stage3_pairs.csv",pairs)
    summary["companion_sha256"]={name:file_hash(output/name) for name in ("stage3_runs.csv","stage3_pairs.csv")}
    public_save(output/"stage3_report.json",summary)
    return summary


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path)
    parser.add_argument("--allow-incomplete",action="store_true",help="write diagnostic outputs under interim/")
    args=parser.parse_args(argv)
    result=run(**vars(args))
    print(f"{len(result['runs'])} registered slots; complete={result['complete']}")


if __name__=="__main__":
    main()
