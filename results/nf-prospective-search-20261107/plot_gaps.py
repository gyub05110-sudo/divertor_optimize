"""Plot the preregistered peak comparison and four-setting numerical range."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def draw(study, output):
    numerical = json.loads((study / "numerics/summary.json").read_text(encoding="utf-8"))
    rows = [row for case in numerical["cases"] for row in case["pairs"]]
    counts = {t: sum(row["cutoff_s"] == t for row in rows) for t in (600, 1800)}
    if counts != {600: 9, 1800: 3}:
        raise ValueError(f"Missing registered paired choices: {counts}")
    edge, accent, band = "#386582", "#12394f", "#eef0ed"
    fig, axes = plt.subplots(2, 1, figsize=(7.8, 6.6), sharex=True,
                             gridspec_kw={"height_ratios": [9, 3]})
    fig.subplots_adjust(left=.29, right=.96, bottom=.17, top=.94, hspace=.4)
    for ax, cutoff in zip(axes, (600, 1800)):
        pairs = [r for r in rows if r["cutoff_s"] == cutoff]
        ax.axvspan(-1, 1, color=band, zorder=0)
        ax.axvline(0, color="#67737b", linewidth=.75)
        labels = []
        for index, row in enumerate(pairs):
            labels.append(f"{row['case_id']}   {row['seed']}")
            if row["all_settings_eligible"]:
                low, high = row["eligible_gap_range_pct"]
                ax.plot([low, high], [index, index], color=edge, linewidth=2, solid_capstyle="round")
                ax.plot([low, high], [index, index], "|", color=edge, markersize=7)
            if row["both_choices_available"]:
                qc, qf = row["native_fast"], row["native_direct"]
                ax.plot(100*(qf-qc)/qf, index, "o", markersize=4.5, color=accent, zorder=3)
            else:
                ax.text(.99, index, "pair unavailable", transform=ax.get_yaxis_transform(),
                        ha="right", va="center", fontsize=8, color="#60666a")
        ax.set_yticks(range(len(pairs)), labels, fontsize=9)
        ax.set_ylim(len(pairs)-.45, -.55)
        ax.set_title(f"{cutoff:,} s" + ("   All registered seeds" if cutoff == 600 else "   Seed 20711"),
                     loc="left", fontsize=11, weight="bold")
        ax.grid(axis="x", color="#c2c8cc", linewidth=.45, alpha=.65)
        ax.tick_params(axis="both", length=0, pad=5)
        ax.spines[["top", "left", "right"]].set_visible(False)
        ax.spines["bottom"].set_color("#8d989e")
    axes[1].set_xlabel("Peak reduction relative to direct fine search, G (%)", labelpad=7)
    handles = [Line2D([], [], color=accent, marker="o", linestyle="none", markersize=5,
                      label="Native endpoint"),
               Line2D([], [], color=edge, linewidth=2, marker="|", label="Four-setting range"),
               Patch(facecolor=band, edgecolor="#d8dcd6", label="±1% comparison band")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
    output.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(output / f"paired-peak-gaps.{ext}", dpi=225, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    draw(options.study, options.output)
