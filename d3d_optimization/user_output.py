"""Portable whole-wall records for a measured reference or accepted candidate."""
import csv
import gzip
import hashlib
from pathlib import Path

import numpy as np

from .boundary_evaluator import export_boundary
from .reference_transport import thermal_summary


def file_hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def segment_regions(result):
    ids = np.zeros(len(result["wall"]) - 1, int)
    for j, name in enumerate(result["region_names"][1:], 1):
        ids[result["surface_segments"][name]] = j
    return ids


def reference_junction(wall, floor_segment):
    a, b, c = np.asarray(wall)[floor_segment:floor_segment+3]
    left, right = b-a, c-b
    cross = left[0]*right[1] - left[1]*right[0]
    angle = float(np.degrees(np.arctan2(abs(cross), left @ right)))
    return dict(tangent_error_deg=angle, required_below_deg=1e-8, passed=angle < 1e-8)


def surface_metrics(result, weights, power_W):
    region = segment_regions(result)
    dep = result["deposition"]
    hit = result["hits"].segment
    hit_region = np.full(len(hit), -1)
    hit_region[hit >= 0] = region[hit[hit >= 0]]
    rows = {}
    for j, name in enumerate(result["region_names"]):
        cell = region[dep["cell_parent"]] == j
        sampled = region[dep["sample_parent"]] == j
        rows[name] = dict(first_hit_W=float(np.sum(weights[hit_region == j]) * power_W),
                          deposited_W=float(np.sum(dep["cell_power"][cell]) * power_W),
                          peak_MW_m2=float(np.max(dep["q"][sampled]) * power_W / 1e6)
                          if np.any(sampled) else None)
    return rows


def export_measurement(directory, evaluator, result):
    """Write immutable material vertex IDs, cell averages, refined peak and donor paths."""
    directory.mkdir()
    dep = result["deposition"]
    power = evaluator.reference.power_W
    region = segment_regions(result)
    names = np.asarray(result["region_names"])
    wall = result["wall"]
    with (directory / "geometry.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["vertex", "R_m", "Z_m", "outgoing_region"])
        writer.writerows((i, *point, names[region[i]] if i < len(region) else "")
                         for i, point in enumerate(wall))
    receiver = dep["wall"]
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(receiver, axis=0), axis=1))]
    with gzip.open(directory / "profile.csv.gz", "wt", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["wall_segment", "region", "s_start_m", "s_end_m", "R_start_m",
                         "Z_start_m", "R_end_m", "Z_end_m", "area_m2", "deposited_W", "cell_mean_MW_m2"])
        for i, (parent, area, fraction) in enumerate(zip(dep["cell_parent"], dep["areas"], dep["cell_power"])):
            writer.writerow((parent, names[region[parent]], arc[i], arc[i+1], *receiver[i],
                             *receiver[i+1], area, fraction*power, fraction*power/(1e6*area)))
    export_boundary(directory / "boundary.npz", evaluator, result)
    hit = result["hits"].segment
    with gzip.open(directory / "donors.csv.gz", "wt", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["source_x_m", "assigned_W", "first_region", "first_R_m", "first_Z_m",
                         "chord_incidence_deg", "analytic_shelf_incidence_deg", "parallel_length_m", "physical_width_m"])
        shelf_angles = result.get("shelf_direct_incidence_deg", np.full(len(hit), np.nan))
        for i, parent in enumerate(hit):
            writer.writerow((evaluator.source.x[i], evaluator.source.weights[i]*power,
                             names[region[parent]] if parent >= 0 else "unhit", *result["hits"].points[i],
                             result["maps"][0].angles[i], shelf_angles[i], result["thermal"]["length_m"][i],
                             dep["sigma_physical_m"][i]))
    cell_power = np.asarray(dep["cell_power"])
    squared = np.sum(cell_power**2 / dep["areas"])
    peak = int(np.argmax(dep["q"]))
    row = dict(physics_measurable=bool(result["physics_measurable"]), candidate_eligible=bool(result["feasible"]),
               peak_MW_m2=float(dep["q"][peak]*power/1e6),
               peak_RZ_m=dep["points"][peak].tolist(),
               peak_region=names[region[dep["sample_parent"][peak]]],
               A2_m2=float(cell_power.sum()**2/squared) if squared > 0 else None,
               surfaces=surface_metrics(result, evaluator.source.weights, power),
               ledger=dep["ledger"], thermal=thermal_summary(result["thermal"], evaluator.source.weights, evaluator.reference),
               receiving_angles=result.get("receiving_angles"), shelf_exposure=result.get("shelf_exposure"),
               separatrix_witness=result.get("boundary_witness"),
               errors=result.get("errors", []), precise_length_ratios=result.get("precise_length_ratios"),
               geometry_screen=result.get("geometry_screen"), shape_constraints=result.get("shape_constraints"),
               directory=directory.name,
               sha256={p.name: file_hash(p) for p in sorted(directory.iterdir())})
    if result.get("reference_only"):
        floor = next(s for s in evaluator.surfaces if s["name"] == "floor")
        row["original_floor_junction"] = reference_junction(evaluator.wall, floor["segment_index"])
    constraint = result.get("shape_constraints")
    if constraint:
        row["gap_margin_m"] = constraint["gap_m"] - constraint["gap_minimum_m"]
        row["curvature_margin_per_m"] = {name: 1/constraint["radius_minimum_m"] - value
                                         for name, value in constraint["curvature_max_1_m"].items()}
    return row
