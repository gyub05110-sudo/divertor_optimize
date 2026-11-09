"""Plot paired product gaps and independent forced-family comparisons from a complete Stage 4 audit."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import isclose, isfinite
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator


SEEDS = (30713, 30731, 30757)
FAMILIES = ("floor-only", "floor-entry", "floor-shelf", "floor-connector")
HEADERS = ("Floor\nonly", "Floor +\nentry", "Floor +\nshelf", "Floor +\nconnector")
INK, BLUE, LIGHT, AMBER = "#242c32", "#245b79", "#dceaf0", "#993f0a"
STYLE = {"font.family": "DejaVu Sans", "font.size": 7.5, "axes.labelcolor": INK,
         "text.color": INK, "xtick.color": INK, "ytick.color": INK, "axes.linewidth": .6,
         "pdf.fonttype": 42, "savefig.facecolor": "white", "figure.facecolor": "white"}


def figure_data(report):
    if report.get("complete") is not True or report.get("schema") != "nf-stage4-forecast-report.v1":
        raise ValueError("a complete audited Stage 4 report is required")
    paired = defaultdict(dict)
    for pair in report["pairs"]:
        if pair["cutoff"] != 600:
            continue
        case, seed = pair["case_id"], pair["seed"]
        if seed in paired[case]:
            raise ValueError("duplicate 600 s seed comparison")
        paired[case][seed] = pair["gap_pct"]
    rows = []
    scorable = forecast_hits = floor_hits = 0
    predictions = report["predictions"]
    ids = [p["case_id"] for p in predictions]
    if len(set(ids)) != len(ids) or set(ids) != set(paired):
        raise ValueError("product and prediction input IDs differ")
    for p in predictions:
        case = p["case_id"]
        if p["seed"] != SEEDS[0] or set(paired[case]) != set(SEEDS):
            raise ValueError("unexpected seed inventory")
        gaps = list(paired[case].values())
        seed_range = None
        if all(v is not None for v in gaps):
            if not all(isfinite(v) for v in gaps):
                raise ValueError("nonfinite product gap")
            seed_range = (min(gaps), median(gaps), max(gaps))
        forced_gaps, minima, competitive, forecast = None, set(), set(), None
        if p["all_controls"]:
            peaks = p["forced_peaks"]
            if set(peaks) != set(FAMILIES) or not all(isfinite(q) and q > 0 for q in peaks.values()):
                raise ValueError("incomplete or invalid forced-family peaks")
            best = min(peaks.values())
            minima = {f for f in FAMILIES if peaks[f] == best}
            competitive = {f for f in FAMILIES if peaks[f] <= 1.01*best}
            if minima != set(p["best_families"]) or competitive != set(p["competitive_families"]):
                raise ValueError("forced-family display disagrees with the audited comparison")
            forced_gaps = tuple(100*(peaks[f]-best)/best for f in FAMILIES)
            if p["scorable"]:
                forecast = p["predicted_family"]
                if forecast not in FAMILIES:
                    raise ValueError("scored prediction requires a focused family")
                scorable += 1
                forecast_hits += forecast in competitive
                floor_hits += "floor-only" in competitive
        if seed_range is not None or forced_gaps is not None:
            rows.append(dict(case_id=case, seed_range=seed_range, forced_gaps=forced_gaps,
                             minima=minima, competitive=competitive, forecast=forecast))
    return dict(rows=rows, total=len(predictions), scorable=scorable, forecast_hits=forecast_hits, floor_hits=floor_hits)


def gap_text(value, competitive):
    if value == 0:
        return "0"
    if value < .01:
        return "<0.01"
    if not competitive and round(value, 2) <= 1:
        return ">1"
    if competitive and value < 1 and round(value, 2) >= 1:
        return "≤1"
    if isclose(value, 1, rel_tol=1e-12) and competitive:
        return "1"
    return f"{value:.2f}".rstrip("0").rstrip(".") if value < 100 else f"{value:.3g}"


@plt.rc_context(STYLE)
def make_figure(data):
    rows, total = data["rows"], data["total"]
    n = len(rows)
    width = 183.
    row_height = 5.7
    plot_height = row_height*max(n, 2)
    bottom, top = 22., 22.
    height = bottom+plot_height+top
    figure = plt.figure(figsize=(width/25.4, height/25.4))
    left, right = 21., 95.
    def label(x, y, message, **kwargs):
        return figure.text(x/width, y/height, message, va="center", **kwargs)
    def panel(x, panel_width):
        return figure.add_axes([x/width, bottom/height, panel_width/width, plot_height/height])
    product, forced = panel(left, 62.), panel(right, 85.)
    label(1., height-3.5, "a", fontsize=11, weight="bold")
    label(5., height-3.5, "Hotspot-guided versus full search", fontsize=9, weight="bold")
    label(right, height-3.5, "b", fontsize=11, weight="bold")
    label(right+4, height-3.5, "Forced search families", fontsize=9, weight="bold")
    paired = [r["seed_range"] for r in rows if r["seed_range"] is not None]
    inventories = sum(r["forced_gaps"] is not None for r in rows)
    label(1., height-8.5, f"600 s · median and range of 3 paired seeds · {len(paired)}/{total} inputs", fontsize=7.2)
    label(right, height-8.5, f"600 s · first seed · gap to observed best (%) · {inventories}/{total} inputs", fontsize=7.2)
    for ax in (product, forced):
        ax.set_ylim(max(n, 1)-.5, -.5)
        ax.spines[:].set_visible(False)
        ax.tick_params(length=0, pad=2)
    product.set_yticks(range(n), [r["case_id"] for r in rows], fontsize=7.3)
    product.text(0, 1, "Input", transform=product.transAxes, ha="right", va="bottom", fontsize=7.3)
    product.spines["bottom"].set_visible(True)
    product.spines["bottom"].set_color("#717981")
    product.tick_params(axis="x", length=2, width=.55)
    product.axvspan(-1., 1., color=LIGHT, ec="none", zorder=0)
    product.axvline(0, color="#596872", lw=.65, zorder=1)
    lower = min([-1., *[p[0] for p in paired]])
    upper = max([1., *[p[2] for p in paired]])
    margin = .08*(upper-lower)
    ticks = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]).tick_values(lower-margin, upper+margin)
    product.set_xticks(ticks)
    product.set_xlim(ticks[0], ticks[-1])
    product.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:+g}".replace("-", "−") if v else "0"))
    product.set_xlabel(r"Peak difference, $100(q_{\rm auto}/q_{\rm full}-1)$ (%)", labelpad=2.5)
    forced.set_xlim(-.5, 3.5)
    forced.set_yticks([])
    forced.set_xticks(range(4), HEADERS, fontsize=7.2)
    forced.xaxis.tick_top()
    for i, row in enumerate(rows):
        if row["seed_range"] is None:
            product.text(.5, i, "Incomplete seed pairs", transform=product.get_yaxis_transform(), ha="center", va="center",
                         fontsize=6.6, color="#5d6267", bbox=dict(facecolor="white", edgecolor="none", pad=.8))
        else:
            lo, mid, hi = row["seed_range"]
            product.hlines(i, lo, hi, color=BLUE, lw=1.05, zorder=3)
            product.plot([lo, hi], [i, i], "|", ms=4, mew=.85, color=BLUE, zorder=3)
            product.plot(mid, i, "o", color=BLUE, ms=4.2, zorder=4)
        if row["forced_gaps"] is None:
            forced.text(1.5, i, "Incomplete family comparison", ha="center", va="center", color="#5d6267", fontsize=7)
            continue
        for j, (family, value) in enumerate(zip(FAMILIES, row["forced_gaps"])):
            exact, close = family in row["minima"], family in row["competitive"]
            fill = BLUE if exact else LIGHT if close else "white"
            forced.add_patch(Rectangle((j-.49, i-.47), .98, .94, fc=fill, ec="#b6c1c7", lw=.45))
            if family == row["forecast"]:
                forced.add_patch(Rectangle((j-.46, i-.41), .92, .82, fill=False, ec="white", lw=2.55))
                forced.add_patch(Rectangle((j-.46, i-.41), .92, .82, fill=False, ec=AMBER, lw=1.45))
            forced.text(j, i, gap_text(value, close), ha="center", va="center", color="white" if exact else INK,
                        fontsize=7.3, weight="bold" if exact else "normal")
    if not n:
        product.text(.5, .5, "No complete seed pairs", transform=product.transAxes, ha="center", va="center")
        forced.text(.5, .5, "No complete family comparisons", transform=forced.transAxes, ha="center", va="center")
    label(1., 2.5, "Shading marks the ±1% practical band.", fontsize=7)
    handles = [Patch(facecolor=BLUE, label="Best"), Patch(facecolor=LIGHT, edgecolor="#8c9aa1", lw=.5, label="Within 1%"),
               Patch(facecolor="white", edgecolor=AMBER, lw=1.45, label="Forecast")]
    figure.legend(handles=handles, ncol=3, loc="center left", bbox_to_anchor=(right/width, 13/height),
                  frameon=False, fontsize=7, borderaxespad=0, handlelength=1.05, handletextpad=.4, columnspacing=1.05)
    denominator = data["scorable"]
    summary = (f"Within 1% · forecast {data['forecast_hits']}/{denominator} · "
               f"always floor {data['floor_hits']}/{denominator}" if denominator else "No scorable family forecasts.")
    label(right, 7.5, summary, fontsize=7)
    if denominator and data["floor_hits"] in (0, denominator):
        observed = "Floor sufficient for all scored inputs." if data["floor_hits"] else "Additional face needed for all scored inputs."
        label(right, 3., observed, fontsize=7)
    return figure


def plot(study, output):
    raw = (study/"stage4_report.json").read_bytes()
    data = figure_data(json.loads(raw))
    digest = sha256(raw).hexdigest()
    output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        figure = make_figure(data)
        try:
            provenance = f"Source: stage4_report.json SHA-256 {digest}"
            figure.savefig(output/"policy-validation.png", dpi=450, metadata={"Description": provenance})
            figure.savefig(output/"policy-validation.pdf", metadata={"Subject": provenance})
        finally:
            plt.close(figure)
    print(f"Rendered {len(data['rows'])}/{data['total']} registered input rows from report SHA-256 {digest}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True, help="directory containing the complete stage4_report.json")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    plot(args.study, args.output)
