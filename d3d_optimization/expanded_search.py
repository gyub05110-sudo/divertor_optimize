"""Expanded, reproducible shape search; fixed EFIT and target-only 1.5 deg."""
import csv
import json
from pathlib import Path
import shutil
import time
import numpy as np
from scipy.optimize import differential_evolution
from scipy.stats import qmc
from d3d_optimization.case import (
    load_case, separatrix, pretrace, ShapeEvaluator, save_json, sha256,
    select_target, rows_for, plots, write_angle_scope, WIDTHS_MM, ANGLE_SCOPE)

import argparse
parser = argparse.ArgumentParser(description="Expanded fixed-equilibrium D3D search")
parser.add_argument("--output",required=True,type=Path)
parser.add_argument("--reference-case",required=True,type=Path)
parser.add_argument("--min-cases",type=int,default=1000)
parser.add_argument("--generations",type=int,default=60)
parser.add_argument("--bound-mm",type=float,default=100.)
args = parser.parse_args()
if args.min_cases < 1000 or args.generations < 0 or not 0 < args.bound_mm <= 1000:
    parser.error("min-cases >=1000, generations >=0, and 0<bound-mm<=1000 required")
output = args.output.resolve()
if (output/"manifest.json").exists():
    raise FileExistsError("refusing to overwrite an existing expanded result")
output.mkdir(parents=True,exist_ok=True)
previous = args.reference_case.resolve()
bound = args.bound_mm/1000
manifest = json.loads((previous / "manifest.json").read_text(encoding="utf-8"))
assert sha256(manifest["shot_path"]) == manifest["shot_sha256"]
for name in ["DIIID_geom.json", "geometry_provenance.json"]:
    shutil.copyfile(previous / name, output / name)
assert sha256(output / "DIIID_geom.json") == manifest["geometry"]["sha256"]
geometry = json.loads((output / "DIIID_geom.json").read_text())
wall = np.asarray(geometry["limiter"], dtype=float)
wall = np.vstack([wall, wall[0]])
target_index = select_target(wall, manifest["geometry"]["sha256"], manifest["target_configuration"])
field, lcfs, axis = load_case(manifest["shot_path"], manifest["time_ms"])
r_sep, xp, psi_sep = separatrix(field, lcfs, axis)
edges, paths, inside, _, drift = pretrace(field, r_sep, axis[1], wall, 320, .002)
e = ShapeEvaluator(field, lcfs, wall, paths, inside, edges, r_sep, axis[1], target_index)
base_target = np.array([float(max(m.q[np.isin(m.segment, e.baseline["target_segments"])],default=0.)) for m in e.baseline["maps"]])
if np.any(base_target<=0):
    manifest.update(numerically_verified_feasible=False,certified_heat_optimum=False,
                    optimization_status="baseline_target_unilluminated")
    save_json(output/"manifest.json",manifest)
    raise SystemExit("target_unilluminated: no positive baseline target peak to normalize the search")
start = time.monotonic()
seen, records, pool = {}, [], []
counts = {"evaluated": 0, "geometry_valid": 0, "feasible": 0}
best = None
stage = "baseline"
log = (output / "evaluations.jsonl").open("w", encoding="utf-8")


def evaluate(coefficients):
    global best
    key = tuple(np.round(coefficients, 12))
    if key in seen:
        return seen[key]
    result = e.evaluate(coefficients)
    counts["evaluated"] += 1
    row = {"case": counts["evaluated"], "stage": stage,
           "coefficients_m": np.asarray(coefficients).tolist(),
           "geometry_valid": result["geometry_valid"], "feasible": result["feasible"]}
    if result["geometry_valid"]:
        counts["geometry_valid"] += 1
        global_ratios = np.array([max(m.q)/b for m, b in zip(result["maps"], e.baseline_peaks)])
        target_peaks = np.array([max(m.q[np.isin(m.segment, result["target_segments"])], default=0.) for m in result["maps"]])
        target_ratios = target_peaks / base_target
        primary = float(max(global_ratios))
        secondary = float(max(target_ratios))
        mean_secondary = float(np.mean(target_ratios))
        row.update(worst_all_wall_peak_ratio=primary, worst_target_peak_ratio=secondary,
                   mean_target_peak_ratio=mean_secondary, target_peak_ratios=target_ratios.tolist(),
                   all_wall_peak_ratios=global_ratios.tolist(), min_angle_deg=result["angle"],
                   angle_application=result["angle_application"], target_power=result["target_power"],
                   max_displacement_m=float(np.max(np.linalg.norm(result["target"]-e.target, axis=1))),
                   target_length_m=result["length"], unmapped_power=result["missing"],
                   constraint_penalties=result["constraint_penalties"])
        if result["feasible"]:
            counts["feasible"] += 1
            # Primary objective remains all-wall worst peak. Secondary chooses
            # among numerically equal primary scores; never retain a worsened
            # all-wall score just to improve target-only heat.
            if best is None or primary < best[0]-1e-10 or (
                    abs(primary-best[0]) <= 1e-10 and (secondary, mean_secondary) < best[1:3]):
                best = (primary, secondary, mean_secondary, np.asarray(coefficients).copy(), result)
                save_json(output / "best_checkpoint.json", row)
            pool.append((primary, secondary, mean_secondary, np.asarray(coefficients).tolist()))
        objective = 1e4*(primary-1)+secondary+.01*mean_secondary+1000*result["penalty"]
    else:
        objective = 1e6+result["penalty"]
    row["search_objective"] = float(objective)
    records.append(row)
    log.write(json.dumps(row, ensure_ascii=False, allow_nan=False)+"\n")
    seen[key] = float(objective)
    if counts["evaluated"] % 128 == 0:
        log.flush()
        summary = "no feasible candidate yet" if best is None else f"global={best[0]:.6f}, target={best[1]:.6f}"
        print(f"{stage}: {counts}; best {summary}; elapsed={time.monotonic()-start:.1f}s", flush=True)
    return float(objective)


evaluate(np.zeros(4))
for i, limit in enumerate([min(.02,bound), min(.05,bound), bound]):
    stage = f"Sobol_{limit*1000:g}mm"
    points = qmc.Sobol(4, scramble=True, seed=20261003+i).random_base2(10)
    for coefficients in limit*(2*points-1):
        evaluate(coefficients)
# Count distinct geometrically valid designs, not rejected shapes or cache hits.
while counts["geometry_valid"] < args.min_cases:
    stage = "extra_valid_20mm"
    points = qmc.Sobol(4, scramble=True, seed=20261003+counts["evaluated"]).random_base2(9)
    for coefficients in min(.02,bound)*(2*points-1):
        evaluate(coefficients)

optimizer_runs = []
for seed in [17, 42, 93]:
    stage = f"DE_seed_{seed}"
    ranked = sorted(pool)
    rng = np.random.default_rng(seed)
    init = np.vstack([np.asarray([p[3] for p in ranked[:24]]).reshape(-1,4), rng.uniform(-bound, bound, (40, 4))])
    if best is not None:
        init[0] = best[3]
    generations = [0]
    def callback(x, convergence):
        generations[0] += 1
        if generations[0] % 10 == 0:
            print(f"{stage} generation {generations[0]}; {counts}; best={best[:3] if best is not None else None}", flush=True)
    opt = differential_evolution(evaluate, [(-bound, bound)]*4, init=init, seed=seed,
                                 maxiter=args.generations, tol=0, atol=0, polish=False, callback=callback)
    optimizer_runs.append({"seed": seed, "generations": int(opt.nit), "nfev": int(opt.nfev),
                           "message": str(opt.message), "search_objective": float(opt.fun)})
log.close()
assert counts["geometry_valid"] >= args.min_cases
if best is None:
    manifest.update(counts=counts,numerically_verified_feasible=False,certified_heat_optimum=False,
                    optimization_status="no_feasible_candidate")
    save_json(output/"manifest.json",manifest)
    raise SystemExit("No candidate satisfied the fixed constraints; source and angle floor were not changed.")
result, coefficients = best[4], best[3]
rows = rows_for(e, e.baseline, "baseline")+rows_for(e, result, "candidate")
with (output / "width_scan.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
np.savez_compressed(output / "case_arrays.npz", R=field.r, Z=field.z, psi=field.psi, lcfs=lcfs,
                    baseline_wall=wall, candidate_wall=result["wall"], target=result["target"],
                    coefficients=coefficients, x_sources=e.sources[0].x)
manifest.update(coefficients_m=coefficients.tolist(), results=rows, feasible_candidate=result["feasible"],
                minimum_target_angle_deg=result["angle"], angle_application=result["angle_application"],
                angle_constraint_scope=ANGLE_SCOPE, constraint_penalties=result["constraint_penalties"],
                target_original_segment=target_index, tubes=320, step_m=.002, preflight_only=False,
                numerically_verified_feasible=False, certified_heat_optimum=False,
                coefficient_bound_m=bound, max_length_ratio=1.2, optimizer=optimizer_runs,
                elapsed_seconds=time.monotonic()-start, counts=counts,
                objective_policy="minimize worst all-wall peak ratio; break equal-primary ties with worst target peak ratio then mean target peak ratio",
                optimization_outcome={"worst_peak_ratio": best[0], "worst_target_peak_ratio": best[1],
                                      "mean_target_peak_ratio": best[2], "peak_reduction_percent": 100*(1-best[0])})
manifest.pop("angle_application_diagnostic", None)
save_json(output / "manifest.json", manifest)
save_json(output / "numerical_checks.json", [])
write_angle_scope(output, e, result)
plots(output, field, lcfs, xp, e, result)
print("DONE", counts, "best", best[:3], "coefficients", coefficients.tolist(), flush=True)
