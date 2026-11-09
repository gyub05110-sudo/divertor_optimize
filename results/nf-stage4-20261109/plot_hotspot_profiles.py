"""Plot observed first-seed heat-flux samples on the original hotspot material face."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.user_output import file_hash
from verification.user_shape_results import arclength, audit_measurement, match, material_coordinate, require, verify_geometry


SEED = 30713
SERIES = {
    "original": ("Original", "#303940", "o"),
    "hotspot": ("Hotspot-guided", "#166886", "D"),
    "full": ("Full search", "#ae491b", "^"),
    "floor-only": ("Floor only", "#53616a", "s"),
}


def native_face_samples(reference_wall, wall, regions, parent, points, q, face):
    ids = np.flatnonzero(regions == face)
    require(len(ids) and np.all(np.diff(ids) == 1), "example requires one contiguous original material face")
    require(parent.dtype.kind in "iu" and parent.shape == q.shape == (len(points),)
            and np.all((0 <= parent) & (parent < len(regions))), "invalid refined flux sample parents")
    on_face = regions[parent] == face
    require(np.any(on_face), "missing native samples on hotspot face")
    arc = arclength(reference_wall)
    x = material_coordinate(reference_wall, wall, parent[on_face], points[on_face])
    return 1000*(x-arc[ids[0]]), q[on_face], 1000*(arc[ids[-1]+1]-arc[ids[0]])


def checked_run(study, row):
    folder = study/row["relative_output"]
    require(file_hash(folder/"receipt.json") == row["receipt_sha256"], "receipt digest differs from audited report")
    receipt = json.loads((folder/"receipt.json").read_bytes())
    require(receipt["status"] == "completed" and receipt["case_id"] == row["case_id"] and
            receipt["arm"] == row["arm"] and receipt["seed"] == row["seed"] == SEED, "different terminal run")
    for filename, key in (("result.json", "result_sha256"), ("events.jsonl", "events_sha256")):
        require(file_hash(folder/filename) == receipt[key] == row[key], "changed source "+filename)
    result = json.loads((folder/"result.json").read_bytes())
    candidate = Candidate.from_record(result["candidate"])
    require(candidate.identity == result["candidate"]["id"] == row["checkpoints"]["600"]["id"] and
            row["status"] == result["status"] == "qualified", "stored geometry differs from qualified deadline winner")
    match(result["selected"]["peak_MW_m2"], row["checkpoints"]["600"]["peak_MW_m2"], "deadline winner peak")
    return folder, result, candidate


def profile(path, reference, measurement, record, face):
    require(np.array_equal(reference.regions, measurement.regions), "original surface labels differ from export")
    # audit_measurement has verified this NPZ plus source CSVs against the result digests.
    with np.load(path/"boundary.npz", allow_pickle=False) as raw:
        x, q, length = native_face_samples(reference.wall, measurement.wall, measurement.regions,
                                           raw["sample_parent"], raw["points"], raw["q_MW_m2"], face)
    peak = int(np.argmax(q))
    match(q[peak], record["surfaces"][face]["peak_MW_m2"], "regional native sampled peak")
    return dict(x=x, q=q, peak=float(q[peak]), peak_x=float(x[peak]), length=length,
                kW=measurement.metrics["regions"][face]["deposited_W"]/1000.)


def examples(study, report):
    require(report["schema"] == "nf-stage4-forecast-report.v1" and report["complete"], "complete Stage 4 audit required")
    runs = {(r["case_id"], r["arm"]): r for r in report["runs"] if r["seed"] == SEED}
    chosen = {}
    for p in report["predictions"]:  # Audited registration order.
        cid, observed = p["case_id"], p["actual_floor_sufficient"]
        if p["scorable"] and observed not in chosen and all(runs[cid, arm]["status"] == "qualified"
                                                           for arm in ("hotspot", "full", "floor-only")):
            chosen[observed] = cid
    require(chosen == {False: "3264fdd71f", True: "70481becc0"}, "different first eligible inputs by observed class")
    outputs = []
    for observed, face in ((False, "shelf"), (True, "floor")):
        cid = chosen[observed]
        arm_rows = ("hotspot", "full") if observed else ("hotspot", "full", "floor-only")
        profiles = {}
        for arm in arm_rows:
            folder, result, candidate = checked_run(study, runs[cid, arm])
            require(result["baseline"]["peak_region"] == face, "unexpected measured original hotspot")
            power = result["reference_SI"]["power_W"]
            base = audit_measurement(folder/"baseline", result["baseline"], power)
            if arm == "hotspot":
                reference = base
                profiles["original"] = profile(folder/"baseline", reference, base, result["baseline"], face)
            else:
                require(np.array_equal(base.wall, reference.wall) and np.array_equal(base.regions, reference.regions),
                        "separate run used another original material wall")
            selected = audit_measurement(folder/"selected", result["selected"], power)
            verify_geometry(base, selected, candidate)
            profiles[arm] = profile(folder/"selected", reference, selected, result["selected"], face)
        outputs.append((cid, face, profiles))
    return outputs


def plot(study, output):
    raw = (study/"stage4_report.json").read_bytes()
    digest = sha256(raw).hexdigest()
    panels = examples(study, json.loads(raw))
    output.mkdir(parents=True, exist_ok=True)
    style = {"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": .65,
             "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": "white"}
    with plt.rc_context(style):
        fig = plt.figure(figsize=(183/25.4, 106/25.4))
        for index, (cid, face, profiles) in enumerate(panels):
            left = .080 + index*.50
            ax = fig.add_axes([left, .405, .390, .465])
            fig.text(left-.052, .935, "ab"[index], size=11, weight="bold")
            fig.text(left, .935, f"{face.capitalize()}  ·  {cid}", size=9.4, weight="bold")
            ax.set_xlim(0., profiles["original"]["length"])
            ax.set_ylim(-.13, 8.15)
            ax.set_xticks(np.arange(0, profiles["original"]["length"], 50))
            ax.set_yticks([0, 2, 4, 6, 8])
            ax.tick_params(length=2.8, width=.65, pad=2)
            ax.grid(axis="y", color="#dce1e3", lw=.5, zorder=0)
            ax.set_xlabel("Original material distance (mm)", labelpad=2)
            if index == 0:
                ax.set_ylabel(r"Surface-normal heat flux (MW m$^{-2}$)", labelpad=2)
            for role, series in profiles.items():
                _, color, symbol = SERIES[role]
                x, q = series["x"], series["q"]
                if role == "floor-only":
                    # Sparse actual samples make its coincidence with the original visible.
                    visible = np.flatnonzero(q >= .03*series["peak"])
                    targets = np.linspace(min(x[visible]), max(x[visible]), 8)
                    ids = np.unique([visible[np.argmin(abs(x[visible]-t))] for t in targets])
                    ax.plot(x[ids], q[ids], marker=symbol, ms=2.8, mfc="white", mec=color, mew=.65, ls="none", zorder=4)
                else:
                    # Draw every saved value, including coincident endpoints and local extrema.
                    # Unordered points are never connected across material or sampling gaps.
                    ax.plot(x, q, marker="o", ms=.67, mfc=color, mew=0, ls="none", zorder=2)
                    ax.plot(series["peak_x"], series["peak"], marker=symbol, ms=4.0, color=color,
                            mec="white", mew=.50, ls="none", zorder=5)
            table = fig.add_axes([left-.010, .09, .405, .195])
            table.set_axis_off()
            table.set(xlim=(0, 1), ylim=(0, 1))
            for xx, col, alignment in ((0., "Geometry", "left"), (.65, r"Peak (MW m$^{-2}$)", "right"),
                                       (1., "Power (kW)", "right")):
                table.text(xx, 1, col, ha=alignment, va="center", fontsize=7.1, weight="bold")
            for j, (role, series) in enumerate(profiles.items()):
                label, color, marker = SERIES[role]
                y = .76-j*.236
                table.plot(.023, y, marker=marker, ms=4.0, mec=color, mew=.7,
                           mfc="white" if role == "floor-only" else color, ls="none")
                table.text(.062, y, label, va="center", fontsize=7.3)
                table.text(.65, y, f"{series['peak']:.3f}", ha="right", va="center", fontsize=7.3)
                table.text(1., y, f"{series['kW']:.3f}", ha="right", va="center", fontsize=7.3)
            print(cid, face, [(role, round(p["peak"], 6), round(p["kW"], 6)) for role, p in profiles.items()])
        fig.text(.03, .025, "Distance starts at the outboard face endpoint. Peaks use native samples. "
                 "Face power sums deposited cell watts.", fontsize=7.2)
        provenance = f"Stage 4 complete report SHA-256 {digest}. First seed {SEED}."
        try:
            fig.savefig(output/"hotspot-profiles.png", dpi=450, metadata={"Description": provenance})
            fig.savefig(output/"hotspot-profiles.pdf", metadata={"Subject": provenance})
        finally:
            plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    plot(args.study, args.output)
