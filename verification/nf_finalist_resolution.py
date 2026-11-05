"""Fixed-design OMP finalist checks across magnetic interpolation and numerical resolution."""
import argparse
from dataclasses import asdict
from pathlib import Path
from time import perf_counter

from d3d_optimization.boundary_evaluator import prepare_boundary, export_boundary
from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.campaign_geometry import Candidate, evaluator, evaluate_candidate
from d3d_optimization.operating_sweep import reuse_boundary
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_search_benchmark import ROOT, measurement
from verification.validate_boundary_transport import LEVELS

FIXTURE = ROOT / "verification/fixtures/nf_robust_pairs.json"
INPUTS = ("input/D3D/d3d/d3d_shot_1a3155370c.parquet", "baseline_review/baseline_review.json",
          "input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json")


def run(output, count, method, backend):
    folder = output / f"{method}_{count}"
    folder.mkdir(parents=True, exist_ok=False)
    _, _, step, nodes, wall_step = LEVELS[{640: 1, 1280: 2}[count]]
    cases = load_json(FIXTURE)["candidates"]
    reference = ReferencePoint()
    start = perf_counter()
    original = prepare_boundary(ROOT, "OMP", reference, count=count, step=step, nodes=nodes,
                                wall_step=wall_step, method=method)
    report = dict(method=method, count=count, step=step, nodes=nodes, wall_step=wall_step,
                  backend=backend, initial_prep_s=perf_counter()-start, fixture_sha256=digest(FIXTURE),
                  input_sha256={name:digest(ROOT/name) for name in INPUTS},
                  driver_sha256=digest(__file__), cases={})
    for chi in (.3, .6):
        ref = ReferencePoint(chi_perp_m2_s=chi)
        base = original if ref == reference else reuse_boundary(original, ref)
        if not base.baseline["feasible"]:
            raise ValueError("original wall failed its fixed baseline constraints")
        solvers = {radius:evaluator(base, radius, backend) for radius in (20,50)}
        for name in ("baseline", *cases):
            if name == "baseline":
                calc = base.baseline
                identity = None
            else:
                candidate = Candidate.from_record(cases[name]["candidate"])
                identity = candidate.identity
                if identity != cases[name]["identity"]:
                    raise ValueError("candidate fixture changed")
                calc = evaluate_candidate(solvers[candidate.radius_mm], candidate)
            row = measurement(calc, ref)
            path = folder / f"chi{chi}_{name}.npz"
            if "deposition" in calc:
                export_boundary(path, base, calc)
                row.update(profile=path.name, sha256=digest(path))
            if backend == "gpu" and chi == .3 and name == "R20_stage18":
                cpu = evaluator(base, candidate.radius_mm, "cpu")
                check = measurement(evaluate_candidate(cpu, candidate), ref)
                parity = abs(row["peak"]/check["peak"]-1)
                if parity > 1e-10 or row["eligible"] != check["eligible"]:
                    raise ValueError("finalist backend mismatch")
                row["cpu_peak_relative_difference"] = parity
            row.update(identity=identity, reference=asdict(ref))
            report["cases"][path.stem] = row
            print(method, count, path.stem, row.get("peak"), row["eligible"], flush=True)
        del solvers, base
    report["total_s"] = perf_counter()-start
    save(folder/"index.json", report, exclusive=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--method", required=True, choices=("cubic", "quintic"))
    p.add_argument("--count", type=int, required=True, choices=(640,1280))
    p.add_argument("--backend", default="gpu", choices=("gpu","cpu"))
    args = p.parse_args()
    run(args.output.resolve(), args.count, args.method, args.backend)
