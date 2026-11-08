"""Plot the same fixed geometries on each common grid, normalized to that grid's original wall."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASES = ("e806f47778", "1a3155370c")
LEVELS = ("cubic640", "cubic1280", "quintic1280")


PANEL = ("g1_zero", "one_amplitude_+20mm", "one_amplitude_+40mm", "floor_winner", "full_winner")


def plot(study, output):
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), sharey=True)
    colors = ("#176777", "#192c53", "#b04b2f")
    for ax, case in zip(axes, CASES):
        folder = study/case/"replay"
        protocol = json.loads((folder/"protocol.json").read_text(encoding="utf-8"))
        report = json.loads((folder/"report.json").read_text(encoding="utf-8"))
        if not report["complete"]:
            raise ValueError("unfinished replay")
        for j, level in enumerate(LEVELS):
            grid = report["results"][level]
            original = grid["baseline"]["peak_MW_m2"]
            y = [1.]
            for role in PANEL:
                score = grid["candidates"][protocol["roles"][role]["id"]]
                y.append(score["peak_MW_m2"]/original if score.get("candidate_eligible") else np.nan)
            ax.scatter(np.arange(len(y)) + .12*(j-1), y, color=colors[j], s=33,
                       marker=("o", "s", "^")[j], label=level, zorder=3)
        ax.axhline(1, color="#8c929a", lw=.8, ls="--")
        ax.set_xticks(np.arange(6), ["Original", "G1 only\n0 mm", "Bend\n+20 mm", "Bend\n+40 mm",
                                    "Floor\n4 coeff.", "Four faces\n16 coeff."])
        ax.set(title=f"DIII-D  {case}   |   1400 ms", xlim=(-.44, 5.47))
        ax.grid(axis="y", alpha=.17)
        ax.spines[["right", "top"]].set_visible(False)
        ax.text(.03, .04, "Bends of −20 and −40 mm failed the shape constraints.",
                transform=ax.transAxes, fontsize=8, color="#4c5864")
    axes[0].set_ylabel("Whole-wall peak / original-wall peak")
    axes[0].set_ylim(bottom=0)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(.5, 0))
    fig.subplots_adjust(left=.09, right=.98, bottom=.23, top=.89, wspace=.08)
    fig.savefig(output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    plot(args.study, args.output)
