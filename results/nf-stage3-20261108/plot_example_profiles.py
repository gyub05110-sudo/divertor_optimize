"""Plot the first qualified source for each reference hotspot, using the first paired seed."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from d3d_optimization.user_output import file_hash
from verification.user_shape_results import audit_measurement, plots


def run(study, output):
    report = json.loads((study/"stage3_report.json").read_bytes())
    if not report["complete"]:
        raise ValueError("complete prospective report required")
    seed = 20711
    by_case = {(r["case_id"],r["arm"]):r for r in report["runs"] if r["seed"] == seed and r["status"] == "qualified"}
    seen = set()
    for case in report["case_flow"]:
        cid = case["case_id"]
        if not all((cid, arm) in by_case for arm in ("floor-screened","full-screened")):
            continue
        source = {}
        for role, arm in (("floor","floor-screened"), ("full","full-screened")):
            row = by_case[cid, arm]
            folder = study/row["relative_output"]
            path = folder/"result.json"
            if file_hash(path) != row["result_sha256"]:
                raise ValueError(f"result changed: {path}")
            source[role] = (folder, json.loads(path.read_bytes()))
        floor_dir, floor = source["floor"]
        hotspot = floor["baseline"]["peak_region"]
        if hotspot in seen:
            continue
        seen.add(hotspot)
        power = floor["reference_SI"]["power_W"]
        series = {"baseline": audit_measurement(floor_dir/"baseline",floor["baseline"],power)}
        for role, (folder, data) in source.items():
            series[role] = audit_measurement(folder/"selected",data["selected"],power)
        out = output/cid
        out.mkdir(parents=True, exist_ok=True)
        plots(out,series,labels=dict(baseline="Reference wall",floor="Floor, screened",full="All surfaces, screened"))
        print(cid, hotspot, seed, out)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,required=True)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parent/"examples")
    args=parser.parse_args()
    run(args.study,args.output)
