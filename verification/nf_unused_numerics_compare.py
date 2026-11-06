"""Recompute the scalar sensitivity comparison from preserved native and alternative scores."""
import argparse
from pathlib import Path

from d3d_optimization.broadening_study import load_json, save


def compare(study):
    summary=load_json(study/"summary.json")
    alternatives={}
    parity=[]
    for path in (study/"numerics").glob("*.json"):
        record=load_json(path)
        if not record["baseline"]["eligible"]:
            raise ValueError("alternative baseline was ineligible")
        alternatives[record["case_id"]]=record
        for candidate in record["candidates"].values():
            if not candidate["alternative"]["eligible"] or not candidate["cpu_check"]["eligible"]:
                raise ValueError("alternative candidate was ineligible")
            parity.extend(candidate["parity"].values())
    rows=[]
    for native in summary["rows"]:
        case=alternatives[native["case_id"]]
        baseline=case["baseline"]
        row=dict(case_id=native["case_id"],seed=native["seed"])
        for policy in ("direct","screen"):
            original=native[policy]
            identity=original["selected_id"]
            score=case["candidates"][identity]["alternative"] if identity else baseline
            values=dict(peak_reduction_pct=100*(1-score["peak"]/baseline["peak"]),
                        A2_change_pct=100*(score["A2"]/baseline["A2"]-1),peak_MW_m2=score["peak"])
            values["peak_reduction_delta_pp"]=values["peak_reduction_pct"]-original["peak_reduction_pct"]
            values["A2_change_delta_pp"]=values["A2_change_pct"]-original["A2_change_pct"]
            row[policy]=values
        delta=100*(row["screen"]["peak_MW_m2"]/row["direct"]["peak_MW_m2"]-1)
        if (delta<0)!=(native["screen_vs_direct_peak_pct"]<0) or (delta==0)!=native["exact_peak_equal"]:
            raise ValueError("a paired decision changed")
        row["screen_vs_direct_peak_pct"]=delta
        rows.append(row)
    return dict(rows=rows,all_nine_paired_orders_retained=True,max_CPU_GPU_relative_difference=max(parity),
        max_peak_reduction_change_pp=max(abs(x[k]["peak_reduction_delta_pp"]) for x in rows for k in ("direct","screen")),
        max_A2_change_pp=max(abs(x[k]["A2_change_delta_pp"]) for x in rows for k in ("direct","screen")),
        total_check_seconds=sum(c["elapsed_s"] for c in alternatives.values()))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,required=True)
    args=parser.parse_args()
    result=compare(args.study)
    path=args.study/"numerics-comparison.json"
    if path.exists():
        if load_json(path)!=result:
            raise ValueError("saved sensitivity summary differs")
        print("Preserved numerical comparison verified.")
    else:
        save(path,result,exclusive=True)

