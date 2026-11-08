"""Plot the local magnetic-pitch bound and the original shelf incidence."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from d3d_optimization.case import load_case
from d3d_optimization.independent_sources import split_design_wall
from verification.nf_equilibrium_times import REVIEW, SHOT
from verification.nf_shelf_corner_check import first_passing_distance
from .measure_shelf import BUNDLE, FAILED, checked_root


def profiles(root, case, method):
    review = json.loads((root / REVIEW).read_text(encoding="utf8"))
    shelf = next(s for s in split_design_wall(review)[1] if s["name"] == "shelf")
    start, end = np.asarray(shelf["endpoints_RZ_m"])
    length = np.linalg.norm(start - end)
    tangent = (start - end) / length
    normal = [-tangent[1], tangent[0]]
    field = load_case(root / SHOT, case["time_ms"], method=method)[0]

    def angles(distance):
        xyz = end + np.atleast_1d(distance)[:, None] * tangent
        b = field.magnetic(xyz)
        bound = np.degrees(np.arctan2(np.linalg.norm(b[:, :2], axis=1), abs(b[:, 2])))
        actual = field.incidence(xyz, np.broadcast_to(normal, xyz.shape))
        if np.ndim(distance) == 0:
            return actual[0], bound[0]
        return actual, bound

    floor = shelf["angle_floor_deg"]
    first = [first_passing_distance(lambda d: angles(d)[i], floor, length) for i in (0, 1)]
    return angles, floor, first


def plot(input_study):
    cases = json.loads((BUNDLE / "summary.json").read_text(encoding="utf8"))["cases"]
    fig, axes = plt.subplots(1, 2, figsize=(9.1, 3.8), layout="constrained", sharey=True)
    for ax, case_id in zip(axes, FAILED, strict=True):
        case = next(c for c in cases if c["case_id"] == case_id)
        root = checked_root(case, input_study)
        for method in ("cubic", "quintic"):
            angles, floor, first = profiles(root, case, method)
            print(case_id, method, "first passing distance (mm), wall and pitch",
                  [None if x is None else round(x * 1000, 5) for x in first])
            if method != "cubic":
                continue
            mm = np.linspace(0, 90, 451)
            actual, bound = angles(mm / 1000)
            ax.plot(mm, actual, color="#126661", lw=2.5, label="Original shelf incidence")
            ax.plot(mm, bound, color="#955033", lw=2.3, linestyle="--", label="Maximum for an axisymmetric normal")
            ax.axhline(floor, color="#444", lw=1.4, linestyle=":", label="1.5° constraint")
            if first[1] is not None and first[1] > 0:
                ax.axvspan(0, first[1] * 1000, color="#955033", alpha=.1)
            ax.set_xlim(0, 90)
        ax.set_title(case_id, loc="left", fontsize=12)
        ax.set_xlabel("Distance from fixed end on original shelf (mm)")
        ax.grid(axis="y", linewidth=.55, alpha=.28)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Angle (degrees)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", frameon=False, ncol=3, fontsize=8)
    fig.suptitle("Incidence and magnetic-pitch limits at the DIII-D shelf", fontsize=14)
    fig.savefig(Path(__file__).with_name("shelf-limits.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-study", type=Path)
    plot(parser.parse_args().input_study)
