"""Find a low-peak OMP divertor shape from an explicitly supplied DIII-D equilibrium."""
from time import perf_counter

if __name__ == "__main__":
    _started = perf_counter()

import argparse
from contextlib import contextmanager
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
COARSE = ("cubic", 320, .0005, 161, .001)
FINE = ("cubic", 1280, .000125, 641, .00025)


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("expected a finite positive number")
    return number


def arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equilibrium", type=Path, required=True, help="DIII-D single-discharge parquet")
    parser.add_argument("--time-ms", type=float, required=True, help="exact EFIT time in ms")
    parser.add_argument("--geometry", type=Path, help="DIII-D material geometry JSON")
    parser.add_argument("--design-policy", type=Path, help="DIII-D region and fixed endpoint policy JSON")
    parser.add_argument("--output", type=Path, required=True, help="new directory, must not exist")
    parser.add_argument("--budget-seconds", type=positive, default=600.)
    parser.add_argument("--power-mw", type=positive, default=1.)
    parser.add_argument("--density-m3", type=positive, default=1e19)
    parser.add_argument("--chi-m2-s", type=positive, default=.3)
    parser.add_argument("--lambda-q-mm", type=positive, default=3.)
    parser.add_argument("--gamma", type=positive, default=7.)
    parser.add_argument("--kappa0", type=positive, default=2000.)
    parser.add_argument("--radius-mm", type=int, choices=(20, 50, 100), default=20)
    parser.add_argument("--seed", type=int, default=91001)
    parser.add_argument("--backend", choices=("gpu", "cpu"), default="cpu")
    parser.add_argument("--family", choices=("full", "floor-only"), default="full")
    args = parser.parse_args(argv)
    if not math.isfinite(args.time_ms) or args.seed < 0:
        parser.error("expected finite time and nonnegative seed")
    return args


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    if hasattr(value, "tolist"):
        return serial(value.tolist())
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def software():
    modules = sorted((ROOT / "d3d_optimization").glob("*.py")) + [Path(__file__), ROOT / "run_tangent_curvature_study.py"]
    versions = {}
    for name in ("numpy", "scipy", "shapely", "pyarrow"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    return dict(python=sys.version, platform=platform.platform(), git_head=revision,
                packages=versions, code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                               for p in modules})


def gpu_environment():
    gpu = None
    if "cupy" in sys.modules:
        cp = sys.modules["cupy"]
        try:
            runtime = cp.cuda.runtime
            props = runtime.getDeviceProperties(cp.cuda.Device().id)
            gpu = dict(cupy=cp.__version__, device=props["name"].decode(),
                       cuda_driver=runtime.driverGetVersion(), cuda_runtime=runtime.runtimeGetVersion())
        except (AttributeError, RuntimeError) as exc:
            gpu = dict(metadata_error=str(exc))
    return gpu


def run(args, started, report, observe):
    deadline = started + args.budget_seconds
    report.update(status="budget_exhausted", reason="time budget exhausted during preparation")
    if perf_counter() >= deadline:
        return
    @contextmanager
    def stage(name):
        begin = perf_counter()
        observe(dict(stage=name, status="started", elapsed_s=begin-started))
        try:
            yield
        finally:
            end = perf_counter()
            record = dict(stage=name, start_s=begin-started, end_s=end-started,
                          duration_s=end-begin, within_budget=end <= deadline)
            report["timings"].setdefault("preparation", []).append(record)
            observe(dict(status="finished", **record))

    if not args.equilibrium.is_file():
        report.update(status="input_error", reason=f"equilibrium file missing: {args.equilibrium}")
        return

    from d3d_optimization.user_input import load_user_input
    from d3d_optimization.user_physics import prepare_user, shelf_limits
    from d3d_optimization.user_search import search_peak
    from d3d_optimization.reference_transport import ReferencePoint
    from d3d_optimization.campaign_geometry import evaluator
    from d3d_optimization.user_output import export_measurement

    if perf_counter() >= deadline:
        return
    try:
        with stage("load_input"):
            user = load_user_input(args.equilibrium, args.time_ms, args.geometry, args.design_policy)
        reference = ReferencePoint(power_W=args.power_mw*1e6, density_m3=args.density_m3,
                                   chi_perp_m2_s=args.chi_m2_s, lambda_q_m=args.lambda_q_mm/1000.,
                                   gamma=args.gamma, kappa0=args.kappa0)
    except (ValueError, OSError) as exc:
        report.update(status="input_error", reason=str(exc))
        return
    report.update(inputs=user.provenance, reference_SI=asdict(reference),
                  equilibrium_topology={k: user.review[k] for k in ("axis_RZ_m", "xpoint_RZ_m", "legs")})
    if perf_counter() >= deadline:
        return
    with stage("shelf_limits"):
        limits = shelf_limits(user, FINE[0])
    report["fixed_shelf"] = limits
    if perf_counter() >= deadline:
        return
    with stage("prepare_fine_reference"):
        original = prepare_user(user, reference, FINE)
    ready = perf_counter()
    base = original.baseline
    report["timings"]["baseline_ready_s"] = ready - started
    selected = None
    owner = None
    if not base.get("physics_measurable", False):
        report.update(status="calculation_error", reason="reference heat calculation failed",
                      baseline=dict(physics_measurable=False, errors=base.get("errors")))
    elif ready > deadline:
        report["reason"] = "reference calculation completed after the search deadline"
    elif not limits["admissible"]:
        report.update(status="constraint_excluded",
                      reason="fixed shelf tangent or local magnetic pitch violates receiving-angle floor under positive kernel support")
    else:
        with stage("prepare_fine_design"):
            fine = evaluator(original, args.radius_mm, args.backend)
        coarse = None
        if perf_counter() < deadline:
            with stage("prepare_coarse_search"):
                coarse = evaluator(prepare_user(user, reference, COARSE), args.radius_mm, args.backend)
        report["timings"]["search_ready_s"] = perf_counter()-started
        if coarse is not None and perf_counter() < deadline:
            found = search_peak(coarse, fine, radius_mm=args.radius_mm, seed=args.seed,
                                started=started, deadline=deadline, baseline=base,
                                family=args.family, observe=observe)
            selected = found["selected"]
            owner = fine
            candidate = found["winner"]
            report["search"] = {k: v for k, v in found.items() if k not in ("selected", "winner")}
            report.update(status=found["status"], reason=found.get("reason", found["status"]),
                          candidate=dict(id=candidate.identity, **candidate.record()) if candidate else None)
            if found.get("termination") == "failed":
                report.update(status="partial_result" if selected is not None else "calculation_error",
                              reason="search stopped on a calculation failure; see search.errors and events.jsonl")
    report["timings"]["search_finished_s"] = perf_counter()-started
    if "deposition" in base:
        report["baseline"] = dict(export_measurement(args.output / "baseline", original, base),
                                  within_budget=ready <= deadline, reference_only=True)
    if selected is not None:
        report["selected"] = export_measurement(args.output / "selected", owner, selected)
        ratio = report["selected"]["peak_MW_m2"] / report["baseline"]["peak_MW_m2"]
        report["selected"]["reduction_vs_reference_pct"] = 100*(1-ratio)


def main(argv=None, started=None):
    started = started if started is not None else perf_counter()
    args = arguments(argv)
    try:
        args.output.mkdir(exist_ok=False, parents=True)
    except OSError as exc:
        print(exc, file=sys.stderr)
        return 2
    report = dict(schema_version=1, objective="minimum eligible fine whole-wall peak_MW_m2",
                  arguments=vars(args), levels=dict(coarse=COARSE, fine=FINE), timings={},
                  baseline=None, selected=None, candidate=None)
    with (args.output / "events.jsonl").open("x", encoding="utf-8") as stream:
        def observe(event):
            stream.write(json.dumps(serial(event), allow_nan=False) + "\n")
            stream.flush()
        try:
            report["software"] = software()
            run(args, started, report, observe)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            report.update(status="calculation_error", reason=f"{type(exc).__name__}: {exc}")
    timing = report["timings"]
    timing.setdefault("search_finished_s", perf_counter()-started)
    if "software" in report:
        report["software"]["gpu"] = gpu_environment()
    elapsed = perf_counter()-started
    timing["deadline_overrun_s"] = max(0., timing["search_finished_s"]-args.budget_seconds)
    timing["elapsed_before_result_write_s"] = elapsed
    target = args.output / "result.json"
    temp = target.with_suffix(".tmp")
    temp.write_text(json.dumps(serial(report), indent=2, allow_nan=False), encoding="utf-8")
    temp.replace(target)
    print(f"{report['status']} | {target.resolve()}")
    return {"input_error": 2, "calculation_error": 1, "partial_result": 1}.get(report["status"], 0)


if __name__ == "__main__":
    raise SystemExit(main(started=_started))
