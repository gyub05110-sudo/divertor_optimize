"""Check immutable material segments against the supplied closed LCFS in terminal no-design cases."""
import argparse
import csv
import hashlib
from pathlib import Path
import sys

import numpy as np
from shapely.geometry import LineString, Point, Polygon, mapping

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from d3d_optimization.independent_sources import split_design_wall
from d3d_optimization.tangent_geometry import Curve
from d3d_optimization.user_input import OMP_SURFACES, REFERENCE, _parquet
from d3d_optimization.user_physics import _engineering_regions
from optimize_divertor import serial
from verification.stage3_campaign import canonical_hash, public_save, read, require, sha


def coverage(shape, plasma):
    """Keep positive-length intersections and point contacts under the native closed-polygon predicate."""
    if not shape.intersects(plasma):
        return None
    overlap = shape.intersection(plasma)
    pieces = list(overlap.geoms) if hasattr(overlap, "geoms") else [overlap]
    samples = []
    for part in pieces:
        if part.geom_type == "Point":
            samples.append(part)
        elif part.geom_type == "LineString":
            samples.extend([Point(part.coords[0]), part.interpolate(.5, normalized=True), Point(part.coords[-1])])
    sample = max(samples, key=lambda p: p.distance(plasma.boundary))
    interior = plasma.contains(sample)
    return dict(intersection=mapping(overlap), overlap_length_m=float(overlap.length),
                witness_RZ_m=list(sample.coords[0]), witness_strictly_interior=interior,
                witness_distance_to_LCFS_boundary_m=float(sample.distance(plasma.boundary)))


def check_case(study, case, task, result):
    folder = study/task["relative_output"]
    inputs = result["inputs"]
    raw = Path(case["path"]).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == case["sha256"] == inputs["equilibrium_sha256"], "raw source changed")
    require(case["time_ms"] == inputs["time_ms"] and inputs["design_reference_sha256"] == sha(REFERENCE),
            "fixed equilibrium time or engineering reference changed")
    field, lcfs, _, _ = _parquet(raw, case["time_ms"])
    plasma = Polygon(lcfs)
    require(plasma.is_valid and not plasma.is_empty, "invalid LCFS polygon")
    review = dict(read(REFERENCE), legs=inputs["baseline_lower_strike_legs"])
    require(hashlib.sha256(np.asarray(review["wall"],dtype="<f8").tobytes()).hexdigest() == inputs["wall_sha256"],
            "original material wall vertices changed")
    xpoint = result["equilibrium_topology"]["xpoint_RZ_m"]
    review, basis = _engineering_regions(review, field, float(field.values([xpoint])[0]))
    wall, surfaces, original = split_design_wall(review)
    named = {s["segment_index"]: s for s in surfaces}
    active = {s["name"] for s in surfaces if s["name"] in OMP_SURFACES}
    require(active == set(OMP_SURFACES), "incomplete adjustable geometry")

    measurement = folder/"baseline"/"geometry.csv"
    require(sha(measurement) == result["baseline"]["sha256"]["geometry.csv"], "baseline geometry export changed")
    with measurement.open(newline="",encoding="utf8") as stream:
        data = list(csv.DictReader(stream))
    exported = np.array([[float(r["R_m"]),float(r["Z_m"])] for r in data])
    names = np.array([r["outgoing_region"] for r in data[:-1]])
    fine_nodes = result["levels"]["fine"][3]
    intersections, endpoint_hits = [], []
    node = 0
    n_fixed = 0
    for i, (a,b) in enumerate(zip(wall[:-1],wall[1:])):
        surface = named.get(i)
        name = surface["name"] if surface else "fixed"
        stop = node+(fine_nodes-1 if surface else 1)
        require(np.allclose(exported[[node,stop]],[a,b],rtol=0,atol=1e-11) and
                np.all(names[node:stop]==name), "split/material IDs disagree with saved baseline geometry")
        meta = dict(original_material_segment=int(original[i]), split_segment=i,
                    exported_material_segment_first=node, exported_material_segment_last=stop-1,
                    region=name, segment_RZ_m=[a,b])
        if name in OMP_SURFACES:
            for side, vertex in enumerate(surface["endpoints_RZ_m"]):
                hit = coverage(Point(vertex),plasma)
                if hit is not None:
                    endpoint_hits.append(dict(meta,endpoint=side,**hit))
        else:
            n_fixed += 1
            if surface:
                # OMP vectors append eight zero IMP coefficients. The native
                # 1e-7 m chord screen samples a zero-displacement IMP spline at 161 points.
                ca,cb = surface["endpoints_RZ_m"]
                coords = Curve(ca,cb,np.zeros(4)).point(np.linspace(0,1,161))
                coords[0],coords[-1] = ca,cb
            else:
                coords = [a,b]
            hit = coverage(LineString(coords),plasma)
            if hit is not None:
                intersections.append(dict(meta,**hit))
        node = stop
    require(node == len(exported)-1, "incomplete baseline material mapping")
    return dict(case_id=case["id"], time_ms=case["time_ms"], raw_sha256=case["sha256"],
                baseline_geometry_csv_sha256=sha(measurement), original_wall_sha256=inputs["wall_sha256"],
                LCFS_RZ_little_endian_float64_sha256=hashlib.sha256(lcfs.astype("<f8").tobytes()).hexdigest(),
                LCFS_vertices=len(lcfs), engineering_region_basis=basis, split_segments=len(wall)-1,
                fixed_segments_tested=n_fixed, adjustable_endpoints_tested=2*len(active),
                fixed_segment_intersections=intersections, fixed_adjustable_endpoint_intersections=endpoint_hits,
                all_registered_candidate_geometries_rejected_by_lcfs=bool(intersections or endpoint_hits))


def run(study):
    study = Path(study).resolve()
    registry = study/"inputs"/"registration.json"
    expected = study/"expected.json"
    plan = read(expected)
    require(sha(expected) == (study/"expected.sha256").read_text().strip() and
            sha(registry) == plan["registration_sha256"], "source registration or plan changed")
    plan_hash = canonical_hash(plan)
    for name in ("d3d_optimization/campaign_geometry.py", "d3d_optimization/tangent_evaluator.py",
                 "d3d_optimization/tangent_geometry.py", "d3d_optimization/tangent_direct.py",
                 "d3d_optimization/independent_sources.py", "d3d_optimization/user_physics.py",
                 "d3d_optimization/user_input.py", "baseline_review/baseline_review.json"):
        require(sha(ROOT/name) == plan["files"][name], "sealed wall or geometry predicate changed")
    cases = {c["id"]:c for c in read(registry)["cases"]}
    output = {}
    for task in plan["tasks"]:
        folder = study/task["relative_output"]
        terminal = folder/"receipt.json"
        if not terminal.is_file():
            continue
        receipt = read(terminal)
        if receipt["status"] != "completed" or receipt["detail"] != "no_design":
            continue
        case_id = task["case_id"]
        result = read(folder/"result.json")
        require(result["status"]=="no_design" and not result.get("candidate") and
                receipt["result_sha256"]==sha(folder/"result.json") and
                receipt["job_sha256"]==canonical_hash(task) and receipt["plan_sha256"]==plan_hash,
                "unverified no-design result")
        require(case_id==receipt["case_id"] and result["arguments"]["seed"]==task["seed"] and
                result["inputs"]["equilibrium_sha256"]==cases[case_id]["sha256"], "wrong evidence identity")
        if case_id not in output:
            output[case_id] = check_case(study,cases[case_id],task,result)
            output[case_id]["terminal_no_design_runs"] = []
        output[case_id]["terminal_no_design_runs"].append(dict(seed=task["seed"],arm=task["arm"],
                  receipt_sha256=sha(terminal),result_sha256=receipt["result_sha256"],
                  journal_sha256_from_sealed_receipt=receipt["events_sha256"],search_counts=result["search"]["counts"],
                  search_errors=result["search"]["errors"],baseline_measurable=result["baseline"]["physics_measurable"]))
    report = dict(schema="nf-stage3-fixed-wall-domain.v1", expected_sha256=sha(expected), script_sha256=sha(__file__),
                  index_convention="Zero-based original and split wall segment IDs, plus segment rows in baseline geometry.csv",
                  predicate="LineString(assembled_candidate_wall).intersects(Polygon(supplied_LCFS))",
                  implication="The material wall outside the four OMP portions and both ends of every adjustable portion "
                              "are fixed for every permitted candidate. Any point shared by these invariants and the closed "
                              "LCFS polygon forces the complete-wall intersection predicate to reject every permitted candidate "
                              "at this registered frame. This proof assumes the supplied wall, LCFS polygon, fixed regions "
                              "and closed-geometry rule; a negative witness check leaves feasibility undecided.",
                  depth="Distance from a reported interior sample to the LCFS boundary, in metres; no maximum-depth claim.",
                  cases=list(output.values()))
    public_save(study/"fixed-wall-diagnostics.json",serial(report))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    report = run(parser.parse_args().study)
    for row in report["cases"]:
        print(row["case_id"], len(row["fixed_segment_intersections"]), "fixed segment intersections;",
              len(row["fixed_adjustable_endpoint_intersections"]), "active endpoint contacts")

