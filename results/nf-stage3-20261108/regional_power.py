"""Export paired regional power and whole-wall diagnostics for every qualified timed design."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read(path):
    return json.loads(path.read_bytes())


def write(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(study, output):
    report = read(study / "stage3_report.json")
    if not report["complete"] or len(report["runs"]) != 108:
        raise ValueError("complete timed report required")
    regions, overall = [], []
    identity = ["case_id","seed","arm","candidate"]
    columns = ["reference_peak_MW_m2","selected_peak_MW_m2","reference_hotspot","selected_hotspot",
               "reference_A2_m2","selected_A2_m2","reference_mean_sigma_mm","selected_mean_sigma_mm",
               "reference_mean_connection_m","selected_mean_connection_m","reference_deposited_fraction",
               "selected_deposited_fraction","gap_margin_mm","nonfloor_max_displacement_mm"]
    cols_by_region = ["region","reference_first_hit_kW","selected_first_hit_kW","reference_deposited_kW",
                      "selected_deposited_kW","reference_regional_peak_MW_m2","selected_regional_peak_MW_m2"]
    for row in report["runs"]:
        if row["status"] != "qualified":
            continue
        source = study / row["relative_output"] / "result.json"
        if hashlib.sha256(source.read_bytes()).hexdigest() != row["result_sha256"]:
            raise ValueError(f"result digest changed at {source}")
        data = read(source)
        base, candidate = data["baseline"], data["selected"]
        if set(base["surfaces"]) != set(candidate["surfaces"]):
            raise ValueError(f"different material regions at {source}")
        common = {k: row[k] for k in identity[:-1]}
        common["candidate"] = row["candidate"]["id"]
        values = dict(common)
        for tag, obs in (("reference",base),("selected",candidate)):
            deposited = sum(s["deposited_W"] for s in obs["surfaces"].values())
            if not math.isclose(deposited, data["reference_SI"]["power_W"] * obs["ledger"]["deposited"],
                                abs_tol=1e-5, rel_tol=1e-10):
                raise ValueError(f"regional watts contradict ledger at {source}")
            values.update({tag+"_peak_MW_m2":obs["peak_MW_m2"], tag+"_hotspot":obs["peak_region"],
                           tag+"_A2_m2":obs["A2_m2"], tag+"_deposited_fraction":obs["ledger"]["deposited"],
                           tag+"_mean_sigma_mm":1000*obs["thermal"]["sigma_m"]["power_weighted_mean"],
                           tag+"_mean_connection_m":obs["thermal"]["length_m"]["power_weighted_mean"]})
        overall.append(dict(values, gap_margin_mm=candidate["gap_margin_m"]*1000,
                            nonfloor_max_displacement_mm=row["nonfloor_max_displacement_m"]*1000))
        for name in base["surfaces"]:
            b, c = base["surfaces"][name], candidate["surfaces"][name]
            regions.append(dict(common, region=name, reference_first_hit_kW=b["first_hit_W"]/1000,
                                selected_first_hit_kW=c["first_hit_W"]/1000,
                                reference_deposited_kW=b["deposited_W"]/1000,
                                selected_deposited_kW=c["deposited_W"]/1000,
                                reference_regional_peak_MW_m2=b["peak_MW_m2"],
                                selected_regional_peak_MW_m2=c["peak_MW_m2"]))
    output.mkdir(parents=True, exist_ok=True)
    write(output / "winning-walls.csv", identity+columns, overall)
    write(output / "regional-power.csv", identity+cols_by_region, regions)
    print(len(overall), "qualified designs and", len(regions), "paired region records")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    run(args.study, args.output)
