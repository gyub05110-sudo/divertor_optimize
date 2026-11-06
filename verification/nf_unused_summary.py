"""Summarize completed held-out runs and independently check observed selection decisions."""
import argparse
from pathlib import Path
import numpy as np
from d3d_optimization.broadening_study import digest, load_json, save
from verification.nf_unused_live import audit_campaign, SEEDS


def metrics(arm):
    baseline, selected = arm["baseline"], arm["incumbent"]
    if baseline is None or selected is None:
        return dict(verified_result=False, stop=arm["stop"], error=arm.get("error"))
    return dict(verified_result=True, stop=arm["stop"], kind=selected["kind"],
                selected_i=selected.get("i"), selected_id=selected.get("id"),
                baseline_peak_MW_m2=baseline["peak"], peak_MW_m2=selected["peak"],
                peak_reduction_pct=100*(1-selected["peak"]/baseline["peak"]),
                baseline_A2_m2=baseline["A2"], A2_m2=selected["A2"],
                A2_change_pct=100*(selected["A2"]/baseline["A2"]-1),
                joint_improvement=selected["peak"]<baseline["peak"] and selected["A2"]>baseline["A2"],
                fine_count=len(arm["fine_ids"]), coarse_count=len(arm["coarse_ids"]),
                first_5pct_s=arm["first_5pct_s"])


def policy_trace(record):
    events=[x for x in record["events"] if x["status"]=="completed"]
    fine=[x for x in events if x["stage"]=="fine"]
    if record["policy"]=="direct":
        if [x["i"] for x in fine] != list(range(len(fine))):
            raise ValueError("direct fine evaluations differ from proposal order")
    else:
        pending=[]
        for event in events:
            if event["stage"]=="coarse":
                pending.append(event)
            if event["stage"]=="fine":
                eligible=[x for x in pending if x["measurement"]["eligible"]]
                expected=min(eligible,key=lambda x:(x["measurement"]["peak"],x["i"]))
                if event["i"]!=expected["i"]:
                    raise ValueError("screen did not select minimum observed eligible coarse peak")
                pending=[]
            # An entirely ineligible full group causes no fine evaluation.
            if len(pending)==4 and not any(x["measurement"]["eligible"] for x in pending):
                pending=[]
    for event in events:
        if event["stage"] in ("fine","coarse"):
            state=event["measurement"]
            if state["eligible"]:
                if state["errors"] or state["ledger"]["balance_error"]>1e-10:
                    raise ValueError("eligible measurement has errors or power residual")


def summarize(study):
    audit=audit_campaign(study)
    rows=[]
    for pair in audit["pairs"]:
        row=dict(case_id=pair["case_id"],seed=pair["seed"],status=pair["status"])
        folder=study/"live"/row["case_id"]/str(row["seed"])
        if pair["status"]!="completed":
            row["errors"]=pair.get("errors",{})
        for policy, arm in pair.get("arms",{}).items():
            row[policy]=metrics(arm)
            path=folder/f"{policy}.json"
            if path.exists():
                raw=load_json(path)
                policy_trace(raw)
                row[policy]["sha256"]=digest(path)
                row[policy]["new_geometry_count"]=len(raw["generated"])
                row[policy]["phase_seconds"]=pair.get("audits",{}).get(policy,{}).get("phase_seconds")
        if pair["status"]=="completed" and all(row[p]["verified_result"] for p in ("screen","direct")):
            direct,screen=row["direct"],row["screen"]
            delta=100*(screen["peak_MW_m2"]/direct["peak_MW_m2"]-1)
            row.update(screen_vs_direct_peak_pct=delta,screen_vs_direct_A2_pct=100*(screen["A2_m2"]/direct["A2_m2"]-1),
                       exact_peak_equal=screen["peak_MW_m2"]==direct["peak_MW_m2"],
                       within_point1pct=abs(delta)<=.1,within_1pct=abs(delta)<=1.,
                       shared_prefix=pair["shared_prefix"])
        rows.append(row)
    comparisons=[r for r in rows if "screen_vs_direct_peak_pct" in r]
    values=[r["screen_vs_direct_peak_pct"] for r in comparisons]
    aggregate=dict(comparable=len(comparisons),strict_screen_lower=sum(x<0 for x in values),
                   strict_direct_lower=sum(x>0 for x in values),exact_ties=sum(x==0 for x in values),
                   within_point1pct=sum(abs(x)<=.1 for x in values),
                   within_1pct=sum(abs(x)<=1 for x in values),
                   median_screen_vs_direct_peak_pct=float(np.median(values)) if values else None)
    return dict(rows=rows,aggregate=aggregate,case_count=len({r["case_id"] for r in rows}),
                seeds=list(SEEDS),protocol_sha256=digest(study/"live/protocol.json"),
                audit_sha256=digest(study/"live/audit.json"),
                summary_driver_sha256=digest(__file__))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path)
    args=parser.parse_args()
    save(args.study/"summary.json",summarize(args.study),exclusive=True)

