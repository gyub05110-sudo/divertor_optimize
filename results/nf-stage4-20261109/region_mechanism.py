"""Summarize first-seed native regional measurements from a completed Stage 4 study; print JSON to stdout."""
import argparse
from collections import Counter
import hashlib
import json
from math import isclose
from pathlib import Path
import sys


SEED = 30713
PRODUCT = ("hotspot", "full")
FOCUSED = ("floor-only", "floor-entry", "floor-shelf", "floor-connector")
ARMS = PRODUCT + FOCUSED
EXPORTS = {"geometry.csv", "profile.csv.gz", "donors.csv.gz", "boundary.npz"}
REGIONAL = ("peak_MW_m2", "first_hit_W", "deposited_W")


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def measurement(folder, name, raw):
    if not raw or not raw["physics_measurable"]:
        return None
    hashes = raw["sha256"]
    require(raw["directory"] == name and set(hashes) == EXPORTS, "native export inventory changed")
    for file, expected in hashes.items():
        require(sha(folder / name / file) == expected, f"native export digest differs: {folder/name/file}")
    surfaces = {region:{k:row[k] for k in REGIONAL} for region,row in raw["surfaces"].items()}
    peak = raw["peak_MW_m2"]
    region = raw["peak_region"]
    require(peak > 0 and surfaces[region]["peak_MW_m2"] == peak and
            max(s["peak_MW_m2"] for s in surfaces.values() if s["peak_MW_m2"] is not None) == peak,
            "global peak and native regional peak disagree")
    ties = sorted(n for n,s in surfaces.items() if s["peak_MW_m2"] is not None and peak-s["peak_MW_m2"] <= 1e-10*peak)
    return dict(global_peak_MW_m2=peak, peak_region=region, peak_RZ_m=raw.get("peak_RZ_m"), A2_m2=raw["A2_m2"],
                peak_region_power=surfaces[region], tied_peak_regions=ties, tied_peak_region_values={r:surfaces[r] for r in ties},
                surfaces=surfaces, exports_sha256=hashes)


def difference(left, right):
    """Signed left minus right for global scalars and power on the same original material regions."""
    if left is None or right is None:
        return None
    delta = {}
    for name,r in right["baseline_hotspot_regions"].items():
        l = left["baseline_hotspot_regions"][name]
        delta[name] = {k:(None if l[k] is None or r[k] is None else l[k]-r[k]) for k in REGIONAL}
    return dict(global_peak_MW_m2=left["global_peak_MW_m2"]-right["global_peak_MW_m2"],
                A2_m2=left["A2_m2"]-right["A2_m2"], baseline_hotspot_regions=delta)


def case_mechanism(study, case, prediction, runs):
    cid = case["id"]
    require(set(runs) == set(ARMS) and prediction["case_id"] == cid and prediction["seed"] == SEED,
            "first-seed arm or prediction inventory differs")
    baseline = None
    winners, artifacts = {}, {}
    for arm in ARMS:
        run = runs[arm]
        rel = f"runs/{cid}/{SEED}/{arm}"
        require(run["relative_output"] == rel and (run["case_id"],run["seed"],run["arm"]) == (cid,SEED,arm),
                "reported file, seed or arm differs")
        artifacts[arm] = dict(relative_output=rel, result_sha256=run["result_sha256"], terminal_status=run["status"],
                              receipt_status=run["receipt_status"])
        winners[arm] = None
        if run["receipt_status"] != "completed":
            continue
        folder = study / rel
        path = folder / "result.json"
        require(sha(path) == run["result_sha256"], f"timed result digest differs: {path}")
        data = read(path)
        require((data["arguments"]["seed"],data["arguments"]["family"],data["status"]) == (SEED,arm,run["status"]),
                "native result identity differs")
        if data.get("inputs"):
            require(data["inputs"]["equilibrium_sha256"] == case["sha256"] and data["inputs"]["time_ms"] == case["time_ms"],
                    "native result references a different equilibrium or time")
        base = measurement(folder, "baseline", data.get("baseline"))
        if base:
            if baseline is None:
                baseline = dict(base, source_arm=arm, result_sha256=run["result_sha256"],
                                within_budget=data["baseline"]["within_budget"])
            else:
                require(all(base[k] == baseline[k] for k in ("global_peak_MW_m2","A2_m2","surfaces")),
                        "unmatched original-wall measurements")
            artifacts[arm]["baseline_within_budget"] = data["baseline"]["within_budget"]
        cut = run["checkpoints"]["600"]
        if run["status"] not in ("qualified","partial_result") or cut["id"] is None:
            continue
        selected = data["selected"]
        require(base and data["baseline"]["within_budget"] and selected["candidate_eligible"] and selected["physics_measurable"]
                and data["candidate"] == run["candidate"] and data["candidate"]["id"] == cut["id"]
                and isclose(selected["peak_MW_m2"],cut["peak_MW_m2"],rel_tol=1e-12), "timed winner differs from native export")
        winners[arm] = dict(measurement(folder, "selected", selected), candidate_id=cut["id"],
                            result_sha256=run["result_sha256"], native_fine_completion_s=cut.get("completion_s"))
    if baseline is not None:
        # Track every original tie throughout, including unadjustable material.
        original = baseline["tied_peak_regions"]
        baseline["baseline_hotspot_regions"] = {r:baseline["surfaces"][r] for r in original}
        for value in winners.values():
            if value is not None:
                value["baseline_hotspot_regions"] = {r:value["surfaces"][r] for r in original}
                value["change_from_baseline"] = difference(value, baseline)
    predicted = prediction["predicted_family"]
    best = prediction["best_families"]
    require(prediction["all_controls"] == all(winners[a] is not None for a in FOCUSED),
            "native control availability differs from audited prediction")
    if prediction["scorable"]:
        require(predicted in FOCUSED and best and set(best) <= set(FOCUSED), "invalid scored control comparison")
        label = "hit" if prediction["selected_competitive"] else "miss"
        gaps = {a:difference(winners[predicted],winners[a]) for a in best}
    else:
        label = "abstained" if prediction["abstained"] else "incomplete_controls" if not prediction["all_controls"] else "no_prediction"
        gaps = {}
    return dict(case_id=cid, rank=case["rank"], source_file=case["dataset_path"], source_sha256=case["sha256"],
                registration_status=case["status"], time_ms=case["time_ms"], seed=SEED, classification=label,
                prediction=prediction, native_runs=artifacts, baseline=baseline, winners=winners,
                product_minus_full=difference(winners["hotspot"],winners["full"]), predicted_minus_strict_best=gaps)


def analyze(study):
    report_path = study / "stage4_report.json"
    report = read(report_path)
    require(report.get("complete") is True and report.get("schema") == "nf-stage4-forecast-report.v1",
            "complete final Stage 4 report required")
    reg_path = study / "inputs/registration.json"
    require(report["registration_sha256"] == sha(reg_path) and report["plan_sha256"] == sha(study / "expected.json"),
            "final report is not bound to the supplied registration and plan")
    cases = read(reg_path)["cases"]
    predictions = {r["case_id"]:r for r in report["predictions"]}
    runs = {(r["case_id"],r["arm"]):r for r in report["runs"] if r["seed"] == SEED}
    ids = {case["id"] for case in cases}
    require(len(cases) == len(ids) == len(report["predictions"]) == 24 and set(predictions) == ids and
            len(report["runs"]) == 240 and set(runs) == {(cid,a) for cid in ids for a in ARMS} and
            all(r["status"] != "pending" for r in report["runs"]), "incomplete file/seed inventory")
    output = [case_mechanism(study, c, predictions[c["id"]], {a:runs[c["id"],a] for a in ARMS}) for c in cases]
    return dict(schema="nf-stage4-regional-measurements.v1", seed=SEED, native_level="cubic1280",
                report_sha256=sha(report_path), registration_sha256=report["registration_sha256"], plan_sha256=report["plan_sha256"],
                method="All signed differences are first named wall minus second. Original hotspot regions use the registered relative 1e-10 tie tolerance. Powers describe integrated deposition and first hits on the same named material regions; A2 is the whole-wall effective area. Comparisons describe these fixed, timed winners.",
                classification_counts=dict(Counter(row["classification"] for row in output)), cases=output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path,help="complete study directory containing native measurement exports")
    values = analyze(parser.parse_args().study)
    json.dump(values, sys.stdout, indent=2, ensure_ascii=False, allow_nan=False)
    print()
