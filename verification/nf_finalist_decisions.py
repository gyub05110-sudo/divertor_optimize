"""Separate fixed-design peak decisions from local profile qualification."""
import argparse
from pathlib import Path
import sys

import numpy as np

from d3d_optimization.broadening_study import digest, load_json, save
from verification.nf_finalist_resolution import FIXTURE, INPUTS, ROOT

sys.path.insert(0, str(ROOT / "viewers/d3d-progression-20261104"))
from distribution_policy_v3_final import distribution_compare, from_npz

CONFIGS = ("cubic_640", "quintic_640", "cubic_1280", "quintic_1280")
EDGES = (("cubic_640", "quintic_640"), ("cubic_1280", "quintic_1280"),
         ("cubic_640", "cubic_1280"), ("quintic_640", "quintic_1280"))


def span(values):
    return [float(min(values)), float(max(values))]


def hot_intervals(u, q, fraction):
    """Intervals of the saved piecewise-linear profile above a fraction of its own peak.

    Keep the fixed wall seam at 0/1; a patch spanning it is displayed in two pieces.
    """
    u, q = np.asarray(u), np.asarray(q)
    if not 0 < fraction < 1 or u.shape != q.shape or len(u) < 2 or np.any(np.diff(u) <= 0):
        raise ValueError("need ordered samples and a threshold strictly between zero and one")
    if not np.all(np.isfinite(q)) or np.max(q) <= 0 or np.min(q) < 0:
        raise ValueError("invalid heat flux")
    threshold = fraction * np.max(q)
    hot = q >= threshold
    crossing = np.flatnonzero(hot[:-1] != hot[1:])
    roots = u[crossing] + np.diff(u)[crossing] * (threshold-q[crossing]) / np.diff(q)[crossing]
    edges = np.r_[u[0] if hot[0] else [], roots, u[-1] if hot[-1] else []]
    return [(float(a), float(b)) for a,b in edges.reshape(-1,2) if b > a]


def hot_region(profiles, fraction=.99):
    canonical = profiles["cubic_1280"]
    ranges = []
    peaks = {}
    for setting, p in profiles.items():
        u, q = p["continuous"]
        ranges.extend(hot_intervals(u, q, fraction))
        x = u[np.argmax(q)]
        peaks[setting] = dict(s_m=float(p["length"]*x),
                              RZ_m=[float(np.interp(x,p["u"],c)) for c in p["wall"].T])
    merged = []
    for a, b in sorted(ranges):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a,b])
    return dict(peak_by_configuration=peaks, fraction=fraction, canonical="cubic_1280",
                units="metres from wall origin, union across the four numerical cases",
                s_intervals_m=[[a*canonical["length"], b*canonical["length"]] for a,b in merged])


def compare_peaks(first, second):
    if not first or first.keys() != second.keys():
        raise ValueError("need the same nonempty numerical configurations")
    if not all(np.isfinite(v) and v > 0 for v in (*first.values(), *second.values())):
        raise ValueError("peaks must be positive and finite")
    ratios = {key: first[key] / second[key] - 1 for key in first}
    delta = span(ratios.values())
    paired = "changes_or_ties"
    if delta[1] < 0:
        paired = "first_lower"
    elif delta[0] > 0:
        paired = "second_lower"
    intervals = "overlap"
    if max(first.values()) < min(second.values()):
        intervals = "first_lower"
    elif max(second.values()) < min(first.values()):
        intervals = "second_lower"
    return dict(paired_order=paired, interval_order=intervals, paired_delta_fraction=delta,
                paired_by_configuration=ratios, first_peak=span(first.values()), second_peak=span(second.values()))


def read_checked_profiles(root, indexes, key, identity):
    evidence = {}
    for setting, index in indexes.items():
        row = index["cases"][key]
        if row["identity"] != identity:
            raise ValueError(f"candidate identity mismatch: {setting}/{key}")
        if not row["eligible"]:
            continue
        path = root / setting / row["profile"]
        if digest(path) != row["sha256"]:
            raise ValueError(f"source hash mismatch: {path}")
        reference = row["reference"]
        power = reference["power_W"]/1e6
        with np.load(path) as f:
            np.testing.assert_allclose(float(f["assigned_power_W"]), reference["power_W"], rtol=0, atol=0)
            cell = f["cell_power"] * power
            area = np.diff(f["area_edges"])
            if np.any(area <= 0) or np.any(cell < 0):
                raise ValueError("invalid cell power or area")
            q = f["q_MW_m2"]
            np.testing.assert_allclose(cell.sum(), power, rtol=0, atol=2e-12)
            np.testing.assert_allclose(q.max(), row["peak"], rtol=1e-12)
            np.testing.assert_allclose(cell.sum()**2 / np.sum(cell**2/area), row["A2"], rtol=2e-9)
        evidence[setting] = from_npz(path, power_MW=power)
    return evidence


def run(root):
    cases = load_json(FIXTURE)["candidates"]
    inputs = {name:digest(ROOT/name) for name in INPUTS}
    indexes = {s:load_json(root/s/"index.json") for s in CONFIGS}
    for setting, idx in indexes.items():
        if (f"{idx['method']}_{idx['count']}" != setting or idx["fixture_sha256"] != digest(FIXTURE)
                or idx["input_sha256"] != inputs):
            raise ValueError(f"inconsistent experiment inputs at {setting}")
    report = dict(fixture_sha256=digest(FIXTURE), inputs=inputs,
                  index_sha256={s:digest(root/s/"index.json") for s in CONFIGS}, operating={})
    for chi in (.3, .6):
        groups = {}
        for name in ("baseline", *cases):
            key = f"chi{chi}_{name}"
            rows = {s:indexes[s]["cases"][key] for s in CONFIGS}
            identity = cases[name]["identity"] if name in cases else None
            profiles = read_checked_profiles(root, indexes, key, identity)
            comparisons = {}
            for a, b in EDGES:
                if a in profiles and b in profiles:
                    cmp = distribution_compare(profiles[a], profiles[b])
                    comparisons[f"{a}_vs_{b}"] = {k:cmp[k] for k in (
                        "passed", "failures", "linf_relative", "l1_relative", "a2_relative",
                        "power_relative", "wall_distance_m", "wall_length_difference_m")}
            groups[name] = dict(identity=identity,
                eligible_all=len(profiles) == len(CONFIGS),
                failed_settings={s:r["errors"] for s,r in rows.items() if not r["eligible"]},
                peak={s:rows[s]["peak"] for s in profiles}, A2={s:rows[s]["A2"] for s in profiles},
                profile_comparisons=comparisons,
                all_four_profile_pairs_pass=(len(comparisons)==4 and all(c["passed"] for c in comparisons.values())))
            if len(profiles) == len(CONFIGS):
                groups[name]["near_peak_region"] = hot_region(profiles)
        base = groups["baseline"]
        for name in cases:
            g = groups[name]
            if base["eligible_all"] and g["eligible_all"]:
                g["vs_baseline"] = compare_peaks(g["peak"], base["peak"])
                g["A2_gain_fraction"] = span([g["A2"][s]/base["A2"][s]-1 for s in CONFIGS])
        pairs = {}
        for radius in (20,50):
            same_radius = [n for n,c in cases.items() if c["candidate"]["radius_mm"]==radius]
            first, = [n for n in same_radius if cases[n]["role"] == "published"]
            other, = [n for n in same_radius if n != first]
            g1, g2 = groups[first], groups[other]
            pairs[str(radius)] = dict(first=first, second=other, status="insufficient_eligible_settings")
            if g1["eligible_all"] and g2["eligible_all"]:
                pairs[str(radius)] = dict(first=first, second=other, status="compared",
                                         **compare_peaks(g1["peak"],g2["peak"]))
        report["operating"][str(chi)] = dict(geometry=groups, nearest_neighbors=pairs)
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profiles", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    save(args.output, run(args.profiles), exclusive=True)
