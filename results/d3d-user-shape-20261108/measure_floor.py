"""Power accounting on the material floor of the original and both returned walls."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np


def calculate(folder, case, arm, kind):
    report = json.loads((folder/"result.json").read_text(encoding="utf-8"))
    src = folder/kind
    with gzip.open(src/"profile.csv.gz", "rt") as stream:
        cells = [r for r in csv.DictReader(stream) if r["region"] == "floor"]
    with gzip.open(src/"donors.csv.gz", "rt") as stream:
        donors = [r for r in csv.DictReader(stream) if r["first_region"] == "floor"]
    edges = np.r_[float(cells[0]["s_start_m"]), [float(r["s_end_m"]) for r in cells]]
    cumulative = np.r_[0., np.cumsum([float(r["deposited_W"]) for r in cells])]
    s10, s90 = np.interp([.1, .9], cumulative/cumulative[-1], edges)
    watts = [float(r["assigned_W"]) for r in donors]
    means = {k: float(np.average([float(r[k]) for r in donors], weights=watts))
             for k in ("physical_width_m", "parallel_length_m", "chord_incidence_deg")}
    return dict(case=case, family=arm, wall=kind, peak_MW_m2=report[kind]["peak_MW_m2"],
                A2_m2=report[kind]["A2_m2"], first_floor_W=sum(watts), deposited_floor_W=cumulative[-1],
                floor_central_80pct_mm=1000*(s90-s10), **means)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = [calculate(args.study/case/family, case, family, kind)
            for case in ("e806f47778", "1a3155370c") for family, kind in
            (("floor-only", "baseline"), ("floor-only", "selected"), ("full", "selected"))]
    with args.output.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
