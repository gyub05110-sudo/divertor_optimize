"""Audit the prospectively registered product comparison and the independent first-seed controls."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from statistics import median

import numpy as np

from verification import stage3_report as audit
from verification import stage3_campaign as runner
from verification import stage4_campaign as design


REGIONS = ("floor", "omp_connector", "shelf", "omp_entry")
FAMILY = dict(zip(REGIONS, ("floor-only","floor-connector","floor-shelf","floor-entry")))
BLOCKS = ("floor-entry","floor-shelf","floor-connector","floor-only")
CONFIG = {a:(a,"screened") for a in design.ALL}


def qualified_peak(row, cutoff=600):
    if row["status"] in ("qualified","partial_result"):
        return row["checkpoints"][str(cutoff)]["q_native"]
    return None


def reached_threshold(row, cutoff):
    if qualified_peak(row,cutoff) is None:
        return False
    return row["checkpoints"][str(cutoff)]["peak_MW_m2"] <= .95 * row["baseline_peak_MW_m2"]


def product_pair(auto, full, cutoff):
    a, b = qualified_peak(auto,cutoff), qualified_peak(full,cutoff)
    a5, b5 = reached_threshold(auto,cutoff), reached_threshold(full,cutoff)
    gap = None if a is None or b is None else 100*(a-b)/b
    availability = ("both" if a is not None and b is not None else "full_only" if b is not None else
                    "auto_only" if a is not None else "neither")
    return dict(case_id=auto["case_id"],seed=auto["seed"],cutoff=cutoff,availability=availability,
                auto_status=auto["status"],full_status=full["status"],auto_q=a,full_q=b,gap_pct=gap,
                practical=(None if gap is None else "within_1pct" if .99*b<=a<=1.01*b else
                           "auto_lower" if gap<0 else "full_lower"), auto_5pct=a5,full_5pct=b5,
                missed_full_5pct=b5 and not a5,gained_5pct=a5 and not b5)


def prediction_case(group):
    first = {r["arm"]:r for r in group if r["seed"]==design.SEEDS[0]}
    auto = first["hotspot"]
    predicted = auto.get("predicted_family")
    peaks = {a:qualified_peak(first[a]) for a in design.FOCUSED}
    inventory = all(q is not None for q in peaks.values())
    best = min(peaks.values()) if inventory else None
    best_names = [a for a,q in peaks.items() if q==best] if inventory else []
    competitive = [a for a,q in peaks.items() if q<=1.01*best] if inventory else []
    observed = ("floor-only" in competitive) if inventory else None
    scorable = inventory and predicted in design.FOCUSED
    regret = 100*(peaks[predicted]-best)/best if scorable else None
    forecast = auto.get("predicted_floor_sufficient")
    confusion = None
    if scorable:
        audit.require(forecast is (predicted=="floor-only"), "binary prediction disagrees with family")
        confusion = ("true_" if forecast==observed else "false_")+("positive" if forecast else "negative")
    return dict(case_id=auto["case_id"],seed=auto["seed"],predicted_family=predicted,
                predicted_floor_sufficient=forecast,abstained=auto.get("abstained"),
                forced_qualified_count=sum(q is not None for q in peaks.values()),all_controls=inventory,
                forced_peaks=peaks,forced_status={a:first[a]["status"] for a in design.FOCUSED},
                best_families=best_names,competitive_families=competitive,best_peak=best,
                scorable=scorable,actual_floor_sufficient=observed,confusion=confusion,
                selected_competitive=(predicted in competitive) if scorable else None,
                always_floor_competitive=observed if scorable else None,selected_forced_regret_pct=regret)


def verify_regional_samples(wall, names, points, parent, q, surfaces):
    audit.require(parent.shape == q.shape == (len(points),) and parent.dtype.kind in "iu"
                  and np.all((0 <= parent) & (parent < len(names))), "invalid refined sample mapping")
    delta = wall[parent+1] - wall[parent]
    length_sq = np.einsum("ij,ij->i",delta,delta)
    audit.require(np.all(length_sq > 0), "degenerate material segment")
    offset = points-wall[parent]
    fraction = np.clip(np.einsum("ij,ij->i",offset,delta)/length_sq,0.,1.)
    audit.require(np.all(np.linalg.norm(offset - fraction[:,None]*delta,axis=1) <= 1e-9),
                  "refined sample coordinate differs from material parent")
    regions = names[parent]
    for name, row in surfaces.items():
        values = q[regions == name]
        recorded = row["peak_MW_m2"]
        audit.require((recorded is None and not len(values)) or
                      (len(values) and recorded is not None and audit.close(recorded,float(np.max(values)))),
                      f"regional peak differs from refined samples: {name}")


def check_regional_exports(folder, data):
    for kind in ("baseline","selected"):
        record = data.get(kind) or {}
        if not record.get("physics_measurable"):
            continue
        with (folder/kind/"geometry.csv").open(newline="",encoding="utf-8") as f:
            names = np.asarray([line["outgoing_region"] for line in csv.DictReader(f)][:-1])
        with np.load(folder/kind/"boundary.npz",allow_pickle=False) as export:
            verify_regional_samples(export["wall"],names,export["points"],export["sample_parent"],
                                    export["q_MW_m2"],record["surfaces"])


def audit_decision(data, events):
    policy=data.get("surface_policy")
    selected=None
    if policy is not None:
        base=data["baseline"]
        audit.require(base and base["physics_measurable"] and base["within_budget"], "policy without timely baseline")
        surfaces=base["surfaces"]
        peak=max(v["peak_MW_m2"] for v in surfaces.values() if v["peak_MW_m2"] is not None)
        ties=[n for n,v in surfaces.items() if v["peak_MW_m2"] is not None and
              peak-v["peak_MW_m2"]<=1e-10*peak]
        region=None if any(n not in REGIONS for n in ties) else max(ties,key=lambda n:(
            surfaces[n]["deposited_W"],-REGIONS.index(n)))
        predicted=FAMILY.get(region,"full")
        selected = predicted if data["arguments"]["family"]=="hotspot" else data["arguments"]["family"]
        audit.require(set(policy["tied_regions"])==set(ties) and audit.close(peak,policy["global_peak_MW_m2"])
                      and policy["predicted_family"]==predicted and policy["predicted_region"]==region
                      and policy["predicted_floor_sufficient"] is (None if region is None else region=="floor")
                      and policy["selected_family"]==selected and policy["requested_family"]==data["arguments"]["family"],
                      "independent surface prediction disagrees with recorded decision")
        audit.require(policy["policy_version"]=="regional-hotspot-v1" and
                      policy["source"]=="original_native_fine_baseline.surface_metrics" and
                      policy["tie_tolerance_MW_m2"]==1e-10*peak and
                      policy["fallback_reason"]==("unadjustable_region_in_peak_tie" if region is None else None),
                      "policy identity, abstention or tie tolerance differs from freeze")
        audit.require(set(policy["regional_metrics"])==set(surfaces),"missing regional predictors")
        for n,v in surfaces.items():
            audit.require(all(policy["regional_metrics"][n][k]==v[k] for k in ("peak_MW_m2","deposited_W")),
                          "different predictor and baseline regional data")
        journal=[(i,e) for i,e in enumerate(events) if e["stage"]=="surface_policy_decision"]
        audit.require(len(journal)==1 and journal[0][1]["policy"]==policy,
                      "independent decision lacks presearch event")
        completed = [e for e in events if e["stage"]=="surface_policy" and e["status"]=="finished"]
        audit.require(len(completed)==1 and completed[0]["within_budget"] and
                      data["timings"]["baseline_ready_s"]<=completed[0]["end_s"]==journal[0][1]["elapsed_s"]
                      <=data["arguments"]["budget_seconds"],"decision clock differs from measured presearch phase")
        audit.require(all(i>journal[0][0] for i,e in enumerate(events) if e["stage"] in
                          ("prepare_fine_design","prepare_coarse_search","proposal","coarse","fine")),
                      "prediction recorded after design preparation or proposal")
    else:
        audit.require(not any(e["stage"] in ("surface_policy_decision","prepare_fine_design","prepare_coarse_search",
                              "proposal","coarse","fine") for e in events),"decision missing before search")
        if any(e["stage"]=="surface_policy" and e["status"]=="finished" and e["within_budget"] for e in events):
            audit.require(data["status"]=="calculation_error" and any(e["stage"]=="search_failure" for e in events),
                          "completed prediction phase lacks a decision or failure")
    if selected in design.FOCUSED:
        for e in events:
            if e.get("candidate") is not None:
                vec=e["candidate"]["vector"]
                audit.require(all(not any(vec[4*j:4*j+4]) for j,arm in enumerate(BLOCKS)
                                  if arm not in (selected,"floor-only")),"forbidden coefficient varied")
    witness=data.get("fixed_wall")
    if data["status"]=="geometry_excluded":
        audit.require(witness and witness["excluded"] and not data.get("baseline") and not data.get("candidate")
                      and not any(e["stage"].startswith("prepare_") for e in events),"invalid geometric early exit")
    return policy


def audit_jobs(study, plan, cases, partial):
    rows=[]
    signature={}
    for task in plan["tasks"]:
        cid=task["case_id"]
        folder=study/task["relative_output"]
        row=dict(case_id=cid,arm=task["arm"],seed=task["seed"],relative_output=task["relative_output"],
                 predicted_family=None,predicted_floor_sufficient=None,abstained=None,
                 executed_family=None,hotspot=None,decision_s=None)
        if not (folder/"receipt.json").exists():
            audit.require(partial, f"missing terminal receipt: {folder}")
            row.update(audit.no_run("pending","pending"),receipt_status="pending")
        else:
            receipt=runner.verified_receipt(folder,task,plan)
            row.update(receipt_status=receipt["status"],receipt_sha256=runner.sha(folder/"receipt.json"),
                       result_sha256=receipt["result_sha256"],events_sha256=receipt["events_sha256"],
                       cli_elapsed_s=receipt["elapsed_s"])
            if receipt["status"]!="completed":
                row.update(audit.no_run(receipt["status"],receipt["detail"]))
            else:
                data=runner.read(folder/"result.json")
                audit.verify_source(data,task,cases[cid],plan,CONFIG)
                events=[json.loads(line) for line in (folder/"events.jsonl").read_text(encoding="utf8").splitlines()]
                audit.require(data["software"]["git_head"] == plan["revision"], "worker revision differs from seal")
                metrics=audit.audit_run(data,events)
                audit.check_process(receipt,metrics)
                audit.check_exports(folder,data,metrics)
                check_regional_exports(folder,data)
                row.update(metrics)
                policy=audit_decision(data,events)
                row.update(fixed_wall=data.get("fixed_wall"),fixed_shelf=data.get("fixed_shelf"),
                           equilibrium_topology_diagnostic=data.get("equilibrium_topology_diagnostic"),
                           hotspot=(data.get("baseline") or {}).get("peak_region"))
                if policy:
                    pred={k:policy[k] for k in ("predicted_family","predicted_region","predicted_floor_sufficient",
                                                "regional_metrics","tied_regions","fallback_reason")}
                    audit.require(signature.setdefault(cid,pred)==pred,"baseline prediction changed across arms/seeds")
                    row.update(predicted_family=pred["predicted_family"],predicted_region=pred["predicted_region"],
                               predicted_floor_sufficient=pred["predicted_floor_sufficient"],
                               abstained=pred["predicted_family"]=="full",executed_family=policy["selected_family"],
                               decision_s=next(e["elapsed_s"] for e in events if e["stage"]=="surface_policy_decision"))
        rows.append(row)
    return rows


def write_csv(path, rows):
    columns=list(dict.fromkeys(k for r in rows for k in r))
    target=path.with_suffix(".tmp")
    with target.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,columns,lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in row.items()})
    target.replace(path)


def summarize(rows, cases):
    lookup={(r["case_id"],r["seed"],r["arm"]):r for r in rows}
    pairs=[product_pair(lookup[c,s,"hotspot"],lookup[c,s,"full"],t) for c in cases for s in design.SEEDS
           for t in audit.CUTOFFS]
    predictions=[prediction_case([r for r in rows if r["case_id"]==c]) for c in cases]
    by_case=[]
    for c in cases:
        these=[r for r in pairs if r["case_id"]==c and r["cutoff"]==600]
        gaps=[r["gap_pct"] for r in these if r["gap_pct"] is not None]
        by_case.append(dict(case_id=c,comparable_seeds=len(gaps),
                            median_gap_pct=median(gaps) if len(gaps)==3 else None,
                            worst_gap_pct=max(gaps) if len(gaps)==3 else None,
                            auto_5pct_seeds=sum(r["auto_5pct"] for r in these),
                            full_5pct_seeds=sum(r["full_5pct"] for r in these),
                            missed_full_seeds=[r["seed"] for r in these if r["missed_full_5pct"]]))
    counts={}
    for a in design.PRODUCT:
        group=[r for r in rows if r["arm"]==a]
        counts[a]=dict(slots=len(group),statuses=dict(Counter(r["status"] for r in group)),
                       eligible=sum(qualified_peak(r) is not None for r in group),
                       checkpoints_5pct={str(t):sum(reached_threshold(r,t) for r in group) for t in audit.CUTOFFS},
                       reached_cases=sum(any(v["case_id"]==c and reached_threshold(v,600) for v in group) for c in cases),
                       right_censored=sum(bool(r["first_5pct"]["right_censored"]) for r in group))
    scored=[p for p in predictions if p["scorable"]]
    confusion=Counter(p["confusion"] for p in scored)
    def recall(good, missed):
        n=confusion[good]+confusion[missed]
        return confusion[good]/n if n else None
    return dict(pairs=pairs,predictions=predictions,per_case=by_case,product=counts,
                accuracy=dict(scorable=len(scored),abstentions=sum(p["abstained"] is True for p in predictions),
                              selected_hits=sum(p["selected_competitive"] for p in scored),
                              always_floor_hits=sum(p["always_floor_competitive"] for p in scored),
                              confusion={k:confusion[k] for k in ("true_positive","true_negative","false_positive","false_negative")},
                              floor_sufficient_recall=recall("true_positive","false_negative"),
                              other_needed_recall=recall("true_negative","false_positive")))


def run(study,partial=False):
    plan=design.load_seal(study)
    cases={c["id"]:c for c in runner.read(study/"inputs/registration.json")["cases"]}
    for case in cases.values():
        if case["status"] != "acquisition_failed":
            audit.require(runner.sha(case["path"])==case["sha256"],"source bytes changed")
    rows=audit_jobs(study,plan,cases,partial)
    report=dict(schema="nf-stage4-forecast-report.v1",complete=all(r["status"]!="pending" for r in rows),
                plan_sha256=runner.sha(study/"expected.json"),report_driver_sha256=runner.sha(__file__),
                registration_sha256=plan["registration_sha256"],runs=rows,**summarize(rows,cases))
    out=study/"interim" if partial else study
    out.mkdir(exist_ok=True)
    write_csv(out/"policy-predictions.csv",report["predictions"])
    write_csv(out/"product-pairs.csv",report["pairs"])
    write_csv(out/"registered-runs.csv",
              [{k:v for k,v in r.items() if k in ("case_id","seed","arm","status","reason","predicted_family",
                "executed_family","hotspot","decision_s","baseline_ready_s","baseline_peak_MW_m2","relative_output")} |
               dict(first_5pct_s=r["first_5pct"]["observed_s"],right_censored=r["first_5pct"]["right_censored"],
                    censor_s=r["first_5pct"]["censor_s"],**{f"q_{t}":qualified_peak(r,t) for t in audit.CUTOFFS})
               for r in rows])
    report["companion_sha256"]={p:runner.sha(out/p) for p in ("policy-predictions.csv","product-pairs.csv","registered-runs.csv")}
    runner.public_save(out/"stage4_report.json",report)
    print(len(rows),"scheduled;",report["complete"],report["accuracy"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",required=True,type=Path)
    p.add_argument("--allow-incomplete",action="store_true")
    args=p.parse_args()
    run(args.study.resolve(),args.allow_incomplete)
