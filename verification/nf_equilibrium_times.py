"""Test two preregistered time slices without changing the supplied DIII-D wall or OMP design region."""
import argparse
from pathlib import Path
import shutil

import numpy as np
from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.campaign_geometry import Candidate, evaluator, evaluate_candidate
from d3d_optimization.case import load_case, separatrix
from d3d_optimization.continuous_boundary import extended_leg, continuous_design
from d3d_optimization.frozen_mapping import first_hits
from d3d_optimization.independent_sources import split_design_wall
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_search_benchmark import ROOT, prep, measurement

TIMES = (1200, 1600)
SHOT = Path("input/D3D/d3d/d3d_shot_1a3155370c.parquet")
REVIEW = Path("baseline_review/baseline_review.json")
GEOMETRY = Path("input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json")
MANIFEST = GEOMETRY.with_name("manifest.json")
PANEL = ROOT/"results/d3d-base-budget-20261105/paired/OMP/5161/panel.json"


def geometry_at(time_ms):
    old = load_json(ROOT/MANIFEST)
    original = load_json(ROOT/REVIEW)
    if digest(ROOT/SHOT) != old["shot_sha256"] or digest(ROOT/GEOMETRY) != old["geometry"]["sha256"]:
        raise ValueError("source hash mismatch")
    wall = np.asarray(original["wall"])
    field, lcfs, axis = load_case(ROOT/SHOT, time_ms)
    rsep, xp, psi = separatrix(field, lcfs, axis)
    legs = {}
    for source, key in (("IMP", "inner"), ("OMP", "outer")):
        hits = []
        for step in (.0005,.00025):
            path = extended_leg(field, xp, psi, source, step)
            hit = first_hits([path], wall)
            if hit.segment[0] < 0:
                raise ValueError(f"no first hit on {key} separatrix")
            hits.append((int(hit.segment[0]), hit.points[0].tolist()))
        dist = float(np.linalg.norm(np.subtract(hits[0][1], hits[1][1])))
        if hits[0][0] != hits[1][0] or dist > 1e-5:
            raise ValueError(f"unresolved {key} separatrix first hit")
        legs[key] = dict(segment=hits[1][0], strike_RZ_m=hits[1][1], refinement_m=dist)
    info = dict(time_ms=time_ms, axis_RZ_m=axis.tolist(), X_RZ_m=xp.tolist(),
                omp_sep_R_m=float(rsep), legs=legs)
    if legs["inner"]["segment"] != 69:
        return info, None, f"inner strike on {legs['inner']['segment']}, current engineering rule requires 69"
    # Carry only the fields that define the material wall and original permitted segments.
    # The production preparation computes pressure/field/source paths afresh for this time.
    review = {k:original[k] for k in ("wall","optimization_footprints","optimization_policy")}
    review.update(legs=legs, sources={"OMP":{}, "IMP":{}}, time_ms=time_ms)
    checked, _ = continuous_design(review, field, psi)
    _, new, _ = split_design_wall(checked)
    _, old_surfaces, _ = split_design_wall(original)
    active = {"omp_entry","shelf","omp_connector","floor"}
    get = lambda surfaces: {s["name"]:s["endpoints_RZ_m"] for s in surfaces if s["name"] in active}
    if get(new) != get(old_surfaces):
        raise ValueError("OMP adjustable endpoints changed")
    info["active_OMP_endpoints_RZ_m"] = get(new)
    return info, review, None


def preflight(out):
    out.mkdir(parents=True, exist_ok=False)
    for t in TIMES:
        folder = out / str(t)
        try:
            info, review, reason = geometry_at(t)
        except ValueError as exc:
            info, review, reason = dict(time_ms=t), None, str(exc)
        info.update(reason=reason, passed=reason is None, shot_sha256=digest(ROOT/SHOT),
                    review_sha256=digest(ROOT/REVIEW), geometry_sha256=digest(ROOT/GEOMETRY),
                    panel_sha256=digest(PANEL), driver_sha256=digest(__file__))
        if review:
            staged = folder/"input-root"
            save(staged/REVIEW, review, exclusive=True)
            save(staged/MANIFEST, dict(time_ms=t, shot_sha256=info["shot_sha256"]), exclusive=True)
            for item in (SHOT, GEOMETRY):
                (staged/item).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT/item, staged/item)
        save(folder/"preflight.json", info, exclusive=True)
        print(t, info["passed"], reason, info.get("legs"), flush=True)


def measure_time(folder):
    info = load_json(folder/"preflight.json")
    if info["panel_sha256"] != digest(PANEL) or info["driver_sha256"] != digest(__file__):
        raise ValueError("preflight or candidate panel changed")
    report = dict(time_ms=info["time_ms"], preflight_sha256=digest(folder/"preflight.json"),
                  status="outside_implementation_domain", reason=info["reason"])
    if info["passed"]:
        candidates = [Candidate.from_record(row["candidate"]) for row in load_json(PANEL)["rows"]]
        report["status"] = "measured"
        for level in (0,2):
            key = ("coarse" if level == 0 else "fine")
            try:
                original, seconds = prep(folder/"input-root", "OMP", ReferencePoint(), level)
            except ValueError as exc:
                report.update(status="baseline_preparation_failed", reason=str(exc), failed_level=key)
                break
            e = evaluator(original, 20, "gpu")
            rows = [dict(i=i, id=c.identity, **measurement(evaluate_candidate(e,c), original.reference))
                    for i,c in enumerate(candidates)]
            report[key] = dict(baseline=measurement(original.baseline,original.reference),
                               prep_s=seconds, rows=rows)
            if level == 0:
                ids = []
                for i in range(0, len(rows), 4):
                    valid = [r for r in rows[i:i+4] if r["eligible"]]
                    if valid:
                        ids.append(min(valid, key=lambda r: (r["peak"],r["i"]))["i"])
                save(folder/"shortlist.json", dict(ids=ids, data=report[key]), exclusive=True)
            else:
                ids = load_json(folder/"shortlist.json")["ids"]
                baseline = report[key]["baseline"]["peak"]
                best = min([baseline]+[r["peak"] for r in rows if r["eligible"]])
                picked = min([baseline]+[rows[i]["peak"] for i in ids if rows[i]["eligible"]])
                report.update(best=best, selected_peak=picked, baseline_peak=baseline,
                              selected_excess_fraction=picked/best-1, shortlisted=ids)
            del e, original
    save(folder/"result.json", report, exclusive=True)
    print(info["time_ms"], report["status"], report["reason"], flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--phase", required=True, choices=("geometry","heat"))
    args = p.parse_args()
    if args.phase == "geometry":
        preflight(args.output)
    else:
        for t in TIMES:
            measure_time(args.output/str(t))
