"""Decompose each fixed full winner's peak advantage at the same source strength and grid."""
import argparse
import csv
import json
import math
from pathlib import Path


LEVELS = ("cubic1280", "cubic640", "quintic1280")
FIELDS = ("case_id", "seed", "arm", "comparator", "level", "native_reproduced", "full_id", "floor_id",
          "ablated_id", "nonfloor_coefficients_all_zero", "full_peak_MW_m2", "floor_peak_MW_m2",
          "ablated_peak_MW_m2", "floor_minus_full_pct_of_full", "ablated_minus_full_pct_of_full",
          "floor_minus_ablated_pct_of_full", "comparison_available")


def table(summary, output):
    data = json.loads(summary.read_bytes())
    if not data["complete"] or len(data["cases"]) != 12:
        raise ValueError("complete 12-file replay required")
    rows = []
    for case in data["cases"]:
        for pair in case["comparisons"]:
            if not pair["replay_trigger"]:
                continue
            for name in LEVELS:
                level = pair["by_level"][name]
                full, floor, ablated = (level[k] for k in ("full_q", "floor_q", "ablation_q"))
                available = all(x is not None for x in (full, floor, ablated))
                if available:
                    family = 100*(floor-full)/full
                    effect = 100*(ablated-full)/full
                    search = 100*(floor-ablated)/full
                    if not math.isclose(family, effect+search, abs_tol=1e-12):
                        raise ValueError("ablation decomposition is inconsistent")
                else:
                    family = 100*(floor-full)/full if full and floor else None
                    effect = search = None
                rows.append(dict(case_id=case["case_id"],seed=pair["seed"],arm=pair["arm"],
                                 comparator=pair["comparator"],level=name,
                                 native_reproduced=pair["native_reproduction_passed"],
                                 full_id=pair["full_id"],floor_id=pair["floor_id"],ablated_id=pair["ablation_id"],
                                 nonfloor_coefficients_all_zero=pair["zero_nonfloor_full"],
                                 full_peak_MW_m2=full,floor_peak_MW_m2=floor,ablated_peak_MW_m2=ablated,
                                 floor_minus_full_pct_of_full=family,ablated_minus_full_pct_of_full=effect,
                                 floor_minus_ablated_pct_of_full=search,comparison_available=available))
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer=csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(len(rows), "fixed pair/grid comparisons")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",required=True,type=Path)
    parser.add_argument("--output",type=Path,
                        default=Path(__file__).resolve().parent/"ablation.csv")
    args=parser.parse_args()
    table(args.summary,args.output)
