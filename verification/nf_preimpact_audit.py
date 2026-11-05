"""Verify saved ring profiles and bounds, and export a paired comparison on physical wall arclength."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


CASES = ("baseline", "OMP50", "OMP20")
PRIMARY = "leg-1.305"
GRID = "2mm_0"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def recompute(power, area):
    q = power / area
    return {"peak": float(q.max()), "A2": float(power.sum() ** 2 / np.dot(power, q))}


def row_for(case, mode, grid, fine, coarse, profile):
    power = profile["resolved_MW"]
    old = profile["original_MW"]
    area = profile["area"]
    unknown = fine["unknown"]
    assert np.all(area > 0) and np.all(power >= 0) and old.min() > -1e-12
    assert abs(old.sum() - 1) < 1e-10 and abs(power.sum() + unknown - 1) < 1e-10
    assert np.all(profile["unknown_bin_cap_MW"] <= unknown + 1e-12)
    row = dict(case=case, mode=mode, grid=grid, unresolved_MW=unknown,
               max_ledger_error_MW=max(abs(old.sum()-1), abs(power.sum()+unknown-1)),
               old_peak_s_m=float(profile["s"][np.argmax(old / area)]),
               routed_peak_s_m=float(profile["s"][np.argmax(power / area)]),
               max_local_delta_over_old_peak=np.max(abs(power-old)/area) / np.max(old/area),
               half_sum_abs_delta_MW=abs(power-old).sum()/2)
    for name, values in (("old", old), ("routed", power)):
        for metric, observed in recompute(values, area).items():
            np.testing.assert_allclose(observed, fine[name][metric], rtol=2e-12)
            row[f"{name}_{metric}"] = observed
    for metric in ("peak", "A2"):
        row[f"coarse_fine_{metric}_relative"] = coarse["routed"][metric] / fine["routed"][metric] - 1
        for bound in ("arbitrary", "flux_conserving"):
            lo, hi = fine[bound][metric]
            assert lo <= fine["routed"][metric] + 1e-12 and fine["routed"][metric] <= hi + 1e-12
            row[f"{bound}_{metric}_lo"] = lo
            row[f"{bound}_{metric}_hi"] = hi
    return row


def strict_order(rows, bound, metric):
    order = CASES if metric == "A2" else CASES[::-1]
    return all(rows[a][f"{bound}_{metric}_hi"] < rows[b][f"{bound}_{metric}_lo"]
               for a, b in zip(order, order[1:]))


def render_profiles(profiles, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 2, figsize=(10, 6.9), sharex=True, layout="constrained",
                             gridspec_kw={"width_ratios": [1.4, 1]})
    for case, (left, right) in zip(CASES, axes):
        s, a, old, new = (profiles[case][k] for k in ("s", "area", "original_MW", "resolved_MW"))
        old, new = old / a, new / a
        left.plot(s, old, color="#17496C", linewidth=1.5, label="Original")
        left.plot(s, new, color="#C8521E", linewidth=1.25, linestyle="--", label="First contact")
        left.set(ylim=(-0.06, 2.9), ylabel=f"{case}\n$q$ (MW/m²)")
        right.axhline(0, color="#7c8387", linewidth=.7)
        right.plot(s, new-old, color="#713C73", linewidth=1.2)
        right.set(ylim=(-.355, .355), ylabel="$\\Delta q$ (MW/m²)")
        for ax in (left, right):
            ax.grid(axis="y", alpha=.20)
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False, loc="upper left")
    axes[0, 0].set_title("Heat flux, 2 mm rings", loc="left")
    axes[0, 1].set_title("First contact minus original", loc="left")
    for ax in axes[-1]:
        ax.set(xlim=(4.42, 4.98), xlabel="Wall arclength from original vertex (m)")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def audit(study, evidence, source):
    primary = {}
    cases = {}
    rows = []
    profiles = {}
    direct = read(study / "raycheck/result.json")
    inputs = read(evidence)["inputs"]
    table = np.genfromtxt(source, names=True, delimiter=",")
    weights = table["power_fraction"]
    np.testing.assert_allclose(weights, table["power_W"]/1e6, rtol=2e-14)
    for case in CASES:
        entry = inputs[case]
        assert sha256(entry["path"]) == entry["sha256"]
        with np.load(entry["path"]) as arc:
            np.testing.assert_allclose(arc["donor_region_power"].sum(axis=1), weights,
                                       rtol=1e-12, atol=1e-17)
            ids = np.argsort(arc["thermal_sigma_m"])
            median = arc["thermal_sigma_m"][ids[np.searchsorted(np.cumsum(weights[ids]), .5)]]
        fine_path = study / f"{case}_fine"
        fine = read(fine_path / "result.json")
        coarse = read(study / f"{case}_coarse/result.json")
        for record in (coarse, fine):
            assert record["archived_hit_arc_matches"] == record["donors"] == len(weights)
            assert not record["mismatch_intervals"]
        cases[case] = {"donors": len(weights), "max_hit_distance_m": fine["archived_hit_max_distance_m"],
                       "power_weighted_sigma_median_m": float(median), "input_sha256": entry["sha256"]}
        for mode, grids in fine["measures"].items():
            for grid, result in grids.items():
                with np.load(fine_path / f"{mode}_{grid}.npz") as profile:
                    row = row_for(case, mode, grid, result, coarse["measures"][mode][grid], profile)
                    rows.append(row)
                    if (mode, grid) == (PRIMARY, GRID):
                        profiles[case] = dict(profile)
                        primary[case] = row
        for n, level in direct["levels"].items():
            with np.load(study / f"raycheck/{case}_{n}.npz") as ray:
                np.testing.assert_array_equal(ray["area"], profiles[case]["area"])
                measured = recompute(ray["first_hit_MW"], ray["area"])
                for key, value in measured.items():
                    np.testing.assert_allclose(level["rows"][case]["direct"][key], value, rtol=2e-13)
    orders = {}
    for mode in fine["measures"]:
        orders[mode] = {}
        for grid in fine["measures"][mode]:
            selected = {r["case"]: r for r in rows if (r["mode"], r["grid"]) == (mode, grid)}
            orders[mode][grid] = {b: {m: strict_order(selected, b, m) for m in ("peak", "A2")}
                                  for b in ("arbitrary", "flux_conserving")}
    results = dict(primary=primary, cases=cases, bound_orders=orders,
                   max_ledger_error_MW=max(r["max_ledger_error_MW"] for r in rows),
                   max_coarse_fine_relative={m: max(abs(r[f"coarse_fine_{m}_relative"]) for r in rows)
                                             for m in ("peak", "A2")},
                   inputs={"source_sha256": sha256(source), "audit_sha256": sha256(evidence),
                           "shot_sha256": sha256(Path(__file__).resolve().parents[1] /
                                                 "input/D3D/d3d/d3d_shot_1a3155370c.parquet")},
                   direct_check=direct)
    with (study / "comparisons.csv").open("w", newline="") as dest:
        writer = csv.DictWriter(dest, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    with (study / "profiles_2mm.csv").open("w", newline="") as dest:
        writer = csv.writer(dest)
        writer.writerow(("case", "s_m", "R_m", "Z_m", "area_m2", "old_MW", "first_hit_MW",
                         "old_q_MW_per_m2", "first_hit_q_MW_per_m2"))
        for case, data in profiles.items():
            for s, rz, area, old, new in zip(data["s"], data["rz"], data["area"], data["original_MW"],
                                              data["resolved_MW"]):
                writer.writerow((case, s, *rz, area, old, new, old/area, new/area))
    render_profiles(profiles, study / "profiles_2mm.png")
    (study / "verified.json").write_text(json.dumps(results, indent=2), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--source-csv", required=True, type=Path)
    args = parser.parse_args()
    audit(args.study, args.evidence, args.source_csv)
