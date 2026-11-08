"""Replot two completed candidate streams directly from the published evidence ZIP."""
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
STUDY = "nf-prospective-search-20261107"
SOURCE = HERE.parent / STUDY
CASE = "e806f47778"
SEED = 20711
ARMS = (("coarse_adaptive", "Screen + verify", "#14746e"),
        ("fine_adaptive", "Direct fine", "#965135"))


def read_streams():
    manifest = json.loads((SOURCE / "evidence-manifest.json").read_text(encoding="utf8"))["members"]
    data = {}
    with ZipFile(SOURCE / "evidence.zip") as archive:
        for arm, _, _ in ARMS:
            name = f"studies/{STUDY}/live/{CASE}/{SEED}/{arm}/events.jsonl"
            raw = archive.read(name)
            if hashlib.sha256(raw).hexdigest() != manifest[name]["sha256"]:
                raise ValueError(f"event digest mismatch for {arm}")
            rows = [json.loads(line) for line in raw.splitlines()]
            events = [r for r in rows if r["stage"] == "fine" and r["status"] == "completed"
                      and r["end_s"] <= 1800 and not r["cached"] and r["measurement"]["eligible"]]
            data[arm] = [r["measurement"] for r in events]
    return data


def plot():
    data = read_streams()
    fig, (all_ax, zoom_ax) = plt.subplots(1, 2, figsize=(10, 4.15), layout="constrained")
    handles = []
    for arm, label, color in ARMS:
        points = data[arm]
        best = min(points, key=lambda x: x["peak"])
        within = [x for x in points if x["peak"] <= 1.01 * best["peak"]]
        wider = max(within, key=lambda x: x["A2"])
        handles.append(Line2D([], [], marker="o", linestyle="none", color=color, label=label))
        for ax in (all_ax, zoom_ax):
            ax.scatter([x["peak"] for x in points], [x["A2"] for x in points],
                       s=21, alpha=.52, color=color, edgecolors="none")
            ax.scatter(best["peak"], best["A2"], marker="*", s=140, color=color,
                       edgecolors="black", linewidths=.7, zorder=4)
        if arm == "coarse_adaptive":
            zoom_ax.axvspan(best["peak"], best["peak"] * 1.01, color=color, alpha=.10, zorder=0)
            zoom_ax.scatter(wider["peak"], wider["A2"], marker="D", facecolors="white",
                            edgecolors=color, linewidths=1.6, s=65, zorder=5)
        print(arm, len(points), best["peak"], best["A2"], wider["peak"], wider["A2"])
    handles.extend([Line2D([], [], marker="*", color="#444", markersize=11, linestyle="none",
                           label="Minimum peak"), Line2D([], [], marker="D", color=ARMS[0][2],
                           markerfacecolor="white", linestyle="none", label="Largest area within screen +1%")])
    all_ax.set_title("All completed fine evaluations", loc="left", fontsize=11)
    zoom_ax.set_title("Low-peak designs", loc="left", fontsize=11)
    zoom_ax.set_xlim(2.258, 2.327)
    zoom_ax.set_ylim(.53, .76)
    for ax in (all_ax, zoom_ax):
        ax.grid(axis="y", linewidth=.5, alpha=.32)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Whole-wall peak (MW/m²)")
    all_ax.set_ylabel("Effective wetted area A₂ (m²)")
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"DIII-D equilibrium {CASE}   •   seed {SEED}   •   1800 s", fontsize=13)
    fig.savefig(HERE / "peak-area.png", dpi=210)
    plt.close(fig)


if __name__ == "__main__":
    plot()
