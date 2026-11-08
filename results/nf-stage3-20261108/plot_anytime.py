"""Plot best completed fine scores from the audited journals and retain the plotted events."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SEEDS = (20711, 20729, 20747)
ARMS = ("floor-screened", "full-screened", "full-direct")
NAMES = ("Floor, screened", "All surfaces, screened", "All surfaces, direct")
COLORS = ("#2166ac", "#c05519", "#267462")
STYLES = ("-", "--", "-.")
MARKERS = ("o", "s", "^")


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def best_events(study, row):
    """Return improvements, including the first admissible score even if worse than the reference."""
    if row["receipt_status"] != "completed":
        return []
    path = study / row["relative_output"] / "events.jsonl"
    if digest(path) != row["events_sha256"]:
        raise ValueError(f"event journal changed: {path}")
    scores = []
    baseline = row["baseline_peak_MW_m2"]
    best = math.inf
    for line in path.read_text(encoding="utf-8").splitlines():
        event = json.loads(line)
        if (event["stage"] != "fine" or event["status"] != "completed" or not event.get("eligible")
                or not event["within_budget"] or event["end_s"] > 600):
            continue
        q = event["peak_q"]
        if q < best:
            best = q
            # The registered power is 1 MW, so q_native has the same numerical value in MW/m².
            scores.append(dict(case_id=row["case_id"], seed=row["seed"], arm=row["arm"], candidate=event["id"],
                               elapsed_s=event["end_s"], peak_MW_m2=q, reduction_pct=100*(baseline-q)/baseline))
    for cutoff, snap in row["checkpoints"].items():
        earlier = [s for s in scores if s["elapsed_s"] <= int(cutoff)]
        if (earlier[-1]["candidate"] if earlier else None) != snap["id"]:
            raise ValueError(f"plot and audited checkpoint disagree: {row['relative_output']} at {cutoff}")
        if earlier and not math.isclose(earlier[-1]["reduction_pct"], snap["reduction_pct"], abs_tol=1e-10):
            raise ValueError("different peak reduction in plotted series")
    return scores


def plot(study, output):
    report = read(study / "stage3_report.json")
    if not report["complete"] or len(report["runs"]) != 108:
        raise ValueError("a complete audited stage 3 report is required")
    rows = report["runs"]
    series = {(r["case_id"], r["arm"], r["seed"]): best_events(study, r) for r in rows}
    cases = [c["id"] for c in read(study / "inputs" / "registration.json")["cases"]]
    visible = [c for c in cases if any(series[c, a, s] for a in ARMS for s in SEEDS)]
    output.mkdir(parents=True, exist_ok=True)
    values = [point for events in series.values() for point in events]
    keys = ("case_id", "seed", "arm", "candidate", "elapsed_s", "peak_MW_m2", "reduction_pct")
    with (output / "anytime-events.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(values)
    if not visible:
        print("No eligible fine score in the 12 registered sources. No heat-flux curve.")
        return

    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                         "font.family":"DejaVu Sans", "pdf.fonttype":42})
    fig, axes = plt.subplots(len(visible), 3, squeeze=False, sharex=True, sharey="row",
                             figsize=(10.2, 1.5+2.35*len(visible)))
    global_min = min(0, min(v["reduction_pct"] for v in values))
    global_max = max(5, max(v["reduction_pct"] for v in values))
    for i, case in enumerate(visible):
        local = [v["reduction_pct"] for v in values if v["case_id"] == case]
        expanded = max(map(abs, local)) < 1
        lower = min(0, min(local)) if expanded else global_min
        upper = max(.1, max(local)) if expanded else global_max
        margin = .045 * (upper - lower)
        for j, (arm, color) in enumerate(zip(ARMS, COLORS)):
            ax = axes[i, j]
            ax.axhline(0, lw=.7, color="#777777")
            if not expanded:
                ax.axhline(5, lw=.8, color="#666666", ls=":")
            ax.set_ylim(lower-margin, upper+margin)
            ax.grid(axis="y", color="#d5d5d5", lw=.45)
            for k, seed in enumerate(SEEDS):
                points = series[case, arm, seed]
                if not points:
                    continue
                times = [p["elapsed_s"] for p in points]
                reductions = [p["reduction_pct"] for p in points]
                ax.step(times+[600], reductions+[reductions[-1]], where="post", color=color,
                        ls=STYLES[k], lw=1.2, marker=MARKERS[k], ms=4,
                        markevery=range(len(times)), mfc="white", mew=.85, alpha=.95)
            if i == 0:
                ax.set_title(NAMES[j], color=color, weight="bold", fontsize=11)
            if j == 0:
                suffix = "\n(expanded scale)" if expanded else ""
                ax.set_ylabel(case+"\nPeak reduction (%)"+suffix)
            ax.set_xticks([0, 200, 400, 600])
            ax.set_xlim(0, 604)
            if i == len(visible)-1:
                ax.set_xlabel("Elapsed time (s)")
    handles = [Line2D([], [], ls=style, color="#303030", marker=marker, mfc="white", ms=4,
                      label=str(seed)) for seed, style, marker in zip(SEEDS, STYLES, MARKERS)]
    handles.append(Line2D([], [], ls=":", color="#666666", label="5% threshold"))
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(.53,.985))
    fig.tight_layout(rect=(0,0,1,1-.49/fig.get_figheight()), h_pad=1.0)
    for ext in ("pdf","png"):
        fig.savefig(output / ("anytime."+ext), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("Rendered", len(visible), "of", len(cases), "sources with eligible scores;", len(values), "records")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    plot(args.study, args.output)
