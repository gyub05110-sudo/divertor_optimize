"""Frozen serial cold-start experiments; read-only event reconstruction at fixed deadlines."""
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
from dataclasses import asdict
from importlib.metadata import PackageNotFoundError, distributions, version
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from time import perf_counter_ns
import numpy as np

from d3d_optimization.broadening_study import digest, load_json as read, save
from d3d_optimization.campaign_state import keep_awake, process_lock
from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_unused_live import load_cases, input_record, INPUTS
from verification.validate_boundary_transport import LEVELS
from verification.nf_anytime_search import POPULATION, PROJECTED_SEEDS

CODE_ROOT = Path(__file__).resolve().parents[1]
SEEDS = (20611, 20627, 20639)
POLICIES = ("coarse_adaptive", "fine_adaptive", "random_fine")
CHECKPOINTS = (200, 400, 600)
SOURCE_RECORDS = ("protocol.ko.md", "selection.json", "downloads.json", "preflight.json")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def arm_order(case_index, seed):
    shift = (case_index + SEEDS.index(seed)) % len(POLICIES)
    return POLICIES[shift:] + POLICIES[:shift]


def source_files():
    paths = list((CODE_ROOT / "d3d_optimization").glob("*.py"))
    names = ("verification/nf_anytime_campaign.py", "verification/nf_anytime_search.py",
             "verification/nf_unused_live.py", "verification/nf_base_live.py",
             "verification/nf_base_budget.py", "verification/nf_search_benchmark.py",
             "verification/validate_boundary_transport.py", "verification/validate_reference_transport.py",
             "optimize_boundary_transport.py", "run_tangent_curvature_study.py")
    return local_source_closure(paths + [CODE_ROOT / name for name in names])


def local_source_closure(initial):
    """Freeze transitive repository imports without capturing unrelated later analysis files."""
    pending, found = list(initial), set()
    while pending:
        path = pending.pop()
        if path in found:
            continue
        found.add(path)
        package = list(path.relative_to(CODE_ROOT).parts[:-1])
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8-sig"))):
            modules = []
            if isinstance(node, ast.Import):
                modules = [item.name.split(".") for item in node.names]
            elif isinstance(node, ast.ImportFrom):
                base = package[:len(package)-node.level+1] if node.level else []
                base += node.module.split(".") if node.module else []
                modules = [base] + [base + item.name.split(".") for item in node.names if item.name != "*"]
            for parts in modules:
                if not parts:
                    continue
                module = CODE_ROOT.joinpath(*parts)
                candidates = [module.with_suffix(".py"), module / "__init__.py"]
                candidates += [CODE_ROOT.joinpath(*parts[:n]) / "__init__.py" for n in range(1,len(parts))]
                pending.extend(p for p in candidates if p.is_file() and p not in found)
    return sorted(found)


def versions():
    found = {}
    for name in ("numpy", "scipy", "pandas", "pyarrow", "shapely", "cupy", "cupy-cuda12x", "cupy-cuda11x"):
        try:
            found[name] = version(name)
        except PackageNotFoundError:
            pass
    runtime = Path(os.environ.get("DIVERTOR_GPU_RUNTIME", "D:/DivertorOptimization/runtime/cupy"))
    for package in distributions(path=[str(runtime)]):
        found[package.metadata["Name"]] = package.version
    return found


def hardware():
    proc = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"],
                          check=True, capture_output=True, text=True)
    return dict(platform=platform.platform(), machine=platform.machine(), cpu=platform.processor(),
                gpu=proc.stdout.strip().splitlines(), serial_workers=True, blas_threads=1)


def freeze(study, source):
    study, source = Path(study).resolve(), Path(source).resolve()
    protocol = dict(schema="nf-anytime-campaign.v1", code_root=str(CODE_ROOT), source_study=str(source),
        seeds=list(SEEDS), policies=list(POLICIES), checkpoints_s=list(CHECKPOINTS), budget_s=max(CHECKPOINTS),
        source="OMP", radius_mm=20, population=POPULATION, projected_seeds=PROJECTED_SEEDS,
        reference=asdict(ReferencePoint()), levels=dict(coarse=list(LEVELS[0]), fine=list(LEVELS[2])),
        hardware=hardware(),
        cases=load_cases(source), method_order="Latin rotation of case index plus seed index",
        source_records={name: digest(source / name) for name in SOURCE_RECORDS},
        original_reference_hashes={name: digest(CODE_ROOT / name) for name in INPUTS},
        study_protocol_sha256=digest(study / "protocol.ko.md"),
        source_hashes={p.relative_to(CODE_ROOT).as_posix(): digest(p) for p in source_files()},
        python_version=platform.python_version(), packages=versions(),
        timing="parent monotonic launch timestamp through completion of each evaluation",
        exclusions=["input acquisition and geometry preflight", "parent code/input hashing",
                    "empty cache creation", "final result serialization and report generation",
                    "post-trial numerical checks and finite-panel reference"],
        scope="same three previously studied public equilibria; no claim of new unseen equilibrium selection",
        cache_policy="fresh private D directories for each child; no simultaneous study workers")
    path = study / "frozen.json"
    if path.exists():
        old = read(path)
        created = old.pop("created_utc")
        require(old == protocol, "frozen protocol, code or input changed; mixed continuation refused")
        for name, expected in (old["source_hashes"] | old["original_reference_hashes"]).items():
            require(digest(study / "code-snapshot" / name) == expected, "frozen code/input snapshot changed")
        return dict(old, created_utc=created)
    protocol["created_utc"] = datetime.now(timezone.utc).isoformat()
    save(path, protocol, exclusive=True)
    for name, expected in protocol["source_hashes"].items():
        target = study / "code-snapshot" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(CODE_ROOT / name, target)
        require(digest(target) == expected, "source snapshot differs")
    for name, expected in protocol["original_reference_hashes"].items():
        target = study / "code-snapshot" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(CODE_ROOT / name, target)
        require(digest(target) == expected, "original reference snapshot differs")
    return protocol


def reconstruct(record, cutoff):
    """Independent of runner snapshot and optimizer state: only finished eligible fine scores count."""
    require(0 <= cutoff <= record["budget_s"], "checkpoint outside registered budget")
    baseline = record.get("baseline")
    ready = record.get("baseline_ready_s")
    selected = None
    if baseline is not None and ready is not None and ready <= cutoff:
        require(baseline["eligible"], "ineligible reference baseline")
        require(all(math.isfinite(baseline[k]) and baseline[k] > 0 for k in ("peak", "A2")),
                "invalid reference observable")
        selected = dict(kind="original_wall", id=None, candidate=None,
                        peak=baseline["peak"], A2=baseline["A2"])
    count = 0
    thermal_count = 0
    first = None
    for event in record["events"]:
        start, end = event["start_s"], event["end_s"]
        require(math.isfinite(start) and math.isfinite(end) and 0 <= start <= end, "invalid event clock")
        if event["stage"] != "fine":
            continue
        state = event.get("measurement")
        if state and state["eligible"]:
            require(all(math.isfinite(state[k]) and state[k] > 0 for k in ("peak", "A2")),
                    "invalid eligible observable")
        if end > cutoff or event["status"] != "completed" or state is None or event.get("cached"):
            continue
        count += 1
        thermal_count += bool(event.get("physics_evaluated", False))
        require(selected is not None, "candidate fine evaluation admitted before original reference ready")
        if not state["eligible"]:
            continue
        if state["peak"] < selected["peak"]:
            selected = dict(kind="new_design", id=event["id"], candidate=event.get("candidate"),
                            peak=state["peak"], A2=state["A2"])
        if first is None and state["peak"] <= .95 * baseline["peak"]:
            first = end
    result = dict(cutoff_s=cutoff, available=selected is not None, fine_count=count,
                  fine_thermal_count=thermal_count, first_5pct_s=first,
                  first_5pct_right_censored=first is None)
    if selected is not None:
        result.update(selected, peak_reduction_pct=100 * (1 - selected["peak"] / baseline["peak"]),
                      A2_change_pct=100 * (selected["A2"] / baseline["A2"] - 1))
    return result


def validate_record(record, raw):
    """Validate causality, physical design identities, and append-only event evidence."""
    require(raw == record["events"], "raw event stream differs from result events")
    budget = record["budget_s"]
    previous = record["startup_s"]
    require(math.isfinite(previous) and previous >= 0, "invalid startup clock")
    known = {}
    for row in record["candidates"]:
        c = Candidate.from_record(row["candidate"])
        require(c.identity == row["id"] and c.source == "OMP" and c.radius_mm == 20,
                "candidate identity or constraints changed")
        require(c.identity not in known, "duplicate candidate registration")
        known[c.identity] = c.record()
    terminal = False
    baseline_events = []
    for event in raw:
        start, end, status = event["start_s"], event["end_s"], event["status"]
        require(not terminal, "event started after terminal event")
        require(math.isfinite(start) and math.isfinite(end) and previous <= start < budget and end >= start,
                "event overlaps, starts after deadline, or has invalid clock")
        require(status in ("completed", "late", "failed"), "unknown event status")
        if status == "completed":
            require(end <= budget, "late event marked completed")
        elif status == "late":
            require(end >= budget, "on-time event marked late")
        terminal = status != "completed"
        previous = end
        if "id" in event:
            require(event["id"] in known and event.get("candidate") == known[event["id"]],
                    "event design differs from candidate registry")
        if event["stage"] in ("fine", "coarse"):
            require(bool(baseline_events) and start >= baseline_events[0]["end_s"],
                    "score started before completed reference baseline")
            require("id" in event, "score lacks design identity")
            state = event.get("measurement")
            require(state is not None, "score lacks measurement, including late/failure")
            if state["eligible"]:
                require(all(math.isfinite(state[k]) and state[k] > 0 for k in ("peak", "A2")),
                        "invalid eligible measurement")
        if event["stage"] == "baseline" and status == "completed":
            baseline_events.append(event)
    require(math.isfinite(record["elapsed_s"]) and record["elapsed_s"] >= previous,
            "final clock precedes last event")
    require(len(baseline_events) <= 1, "multiple baseline measurements")
    if baseline_events:
        require(record["baseline"] == baseline_events[0]["measurement"] and
                record["baseline_ready_s"] == baseline_events[0]["end_s"], "baseline value or completion differs")
    else:
        require(record.get("baseline") is None and record.get("baseline_ready_s") is None,
                "baseline not backed by completed event")
    for cutoff in CHECKPOINTS:
        if cutoff <= budget:
            reconstruct(record, cutoff)
    return True


def validate_decisions(record):
    """Reject hidden promotion choices and inconsistent published checkpoint selections."""
    pending, seen, promotions, used = [], set(), record.get("promotions", []), 0
    for event in record["events"]:
        if event["stage"] == "coarse":
            require(record["policy"] == "coarse_adaptive", "coarse score outside coarse policy")
        if event["stage"] == "fine" and record["policy"] == "coarse_adaptive":
            require(event.get("search_stage") == "coarse_promotion", "coarse policy bypassed promotion")
        if event["stage"] == "coarse" and event["status"] == "completed" and not event.get("cached"):
            require(len(pending) < 4, "new coarse score started before completed group was promoted")
            if event["measurement"]["eligible"] and event["id"] not in seen:
                seen.add(event["id"])
                pending.append(event)
        if event["stage"] == "fine" and event.get("search_stage") == "coarse_promotion":
            require(record["policy"] == "coarse_adaptive" and len(pending) == 4, "incomplete promotion group")
            expected = min(pending, key=lambda e: (e["measurement"]["peak"], e["start_s"]))
            require(event["id"] == expected["id"] and used < len(promotions), "incorrect coarse choice")
            chosen = promotions[used]
            require(chosen["ids"] == [e["id"] for e in pending] and chosen["selected"] == expected["id"],
                    "promotion record differs from observed coarse group")
            require(pending[-1]["end_s"] <= chosen["at_s"] <= event["start_s"], "invalid promotion clock")
            pending, used = [], used + 1
    if len(promotions) > used:
        require(len(promotions) == used + 1 and len(pending) == 4, "unmatched promotion records")
        chosen = promotions[used]
        expected = min(pending, key=lambda e: (e["measurement"]["peak"], e["start_s"]))
        require(chosen["ids"] == [e["id"] for e in pending] and chosen["selected"] == expected["id"],
                "terminal promotion choice differs")
        require(record["stop"] in ("deadline", "failed"), "unexecuted promotion without terminal cause")
    require(record["policy"] == "coarse_adaptive" or not promotions, "other policy has coarse promotions")
    require(set(record["checkpoints"]) == {str(t) for t in CHECKPOINTS}, "missing or added checkpoints")
    for t in CHECKPOINTS:
        actual, expected = record["checkpoints"][str(t)], reconstruct(record, t)
        require(actual["available"] == expected["available"], "checkpoint availability differs")
        if expected["available"]:
            for name in ("kind", "id", "candidate", "peak", "A2"):
                require(actual["incumbent"][name] == expected[name], "checkpoint incumbent differs: " + name)
        else:
            require(actual["incumbent"] is None, "checkpoint leaked unavailable incumbent")
        require(actual["first_5pct_s"] == expected["first_5pct_s"], "checkpoint first improvement differs")
        require(actual["counts"]["fine"] == expected["fine_count"], "checkpoint fine request count differs")
    require(record["final"] == record["checkpoints"][str(max(CHECKPOINTS))], "final checkpoint differs")
    return True


def validate_initialization(record):
    config = record["adaptive_config"]
    require(config["population"] == POPULATION and config["projected_seeds"] == PROJECTED_SEEDS,
            "population settings differ from frozen native defaults")
    initialized = None
    for event in record["events"]:
        if event["stage"] == "initialize":
            require(record["policy"] != "random_fine", "random policy initialized an adaptive population")
            evidence = event.get("initial_population")
            if evidence is None:
                require(event["status"] != "completed", "completed initialization lacks population evidence")
                continue
            x = np.asarray(evidence["vectors"], float)
            require(x.ndim == 2 and x.shape[1] == 16 and len(x) > 0 and np.isfinite(x).all(),
                    "invalid population vectors")
            span = np.ptp(x, axis=0)
            singular = np.linalg.svd((x-x.mean(axis=0)) / .2, compute_uv=False)
            rank = int(np.count_nonzero(singular > 1e-12))
            full = bool(np.all(span > 1e-12) and rank == 16)
            require(evidence["count"] == len(x) and evidence["dimension"] == 16 and evidence["rank"] == rank and
                    evidence["full_space_passed"] == full and np.allclose(evidence["axis_span"], span, rtol=0, atol=1e-15) and
                    np.allclose(evidence["scale"], .2, rtol=0, atol=1e-15) and
                    np.allclose(evidence["singular_values"], singular, rtol=1e-12, atol=1e-15),
                    "population diversity evidence differs from recorded vectors")
            if event["restart"] == 0:
                require(evidence["startup"] and evidence["full_space_required"] and
                        evidence["archive_count"] == 0 and evidence["inherited_count"] == 0,
                        "cold initial population contains inherited designs")
                if event["status"] == "completed":
                    require(full and len(x) == POPULATION, "initial population failed full-space gate")
                    initialized = event["end_s"]
        if event["stage"] in ("fine", "coarse") and record["policy"] != "random_fine":
            require(initialized is not None and initialized <= event["start_s"],
                    "adaptive score precedes passed initialization gate")
    return True


def check_inputs(case):
    require(input_record(case["root"]) == {k: v for k, v in case.items() if k != "case_id"},
            "staged equilibrium, original material wall or constraints changed")


def worker(study, case_id, seed, policy, started_ns):
    from verification.nf_anytime_search import run_arm
    frozen = read(study / "frozen.json")
    case = next(c for c in frozen["cases"] if c["case_id"] == case_id)
    require(seed in SEEDS and policy in POLICIES, "unregistered child arguments")
    out = study / "live" / case_id / str(seed) / policy
    record = run_arm(policy, Path(case["root"]), seed, out, started_ns, budget_s=max(CHECKPOINTS))
    record.update(case_id=case_id, input_root=case["root"], launch_ns=started_ns, pid=os.getpid(),
                  campaign_sha256=digest(study / "frozen.json"))
    save(out / "result.json", record, exclusive=True)
    return int(record["stop"] == "failed")


def run_one(study, case, seed, policy):
    check_inputs(case)
    out = study / "live" / case["case_id"] / str(seed) / policy
    out.mkdir(parents=True, exist_ok=False)
    locations = {name: out / "cache" / leaf for name, leaf in
                 (("CUPY_CACHE_DIR", "cupy"), ("CUDA_CACHE_PATH", "cuda"), ("TEMP", "tmp"),
                  ("TMP", "tmp"), ("MPLCONFIGDIR", "matplotlib"), ("HF_HOME", "huggingface"))}
    for path in locations.values():
        path.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", PYTHONIOENCODING="utf-8", **{k: str(v) for k, v in locations.items()})
    args = [sys.executable, "-m", "verification.nf_anytime_campaign", "--study", str(study),
            "--case", case["case_id"], "--seed", str(seed), "--policy", policy, "--worker", "--started-ns"]
    print(f"START {case['case_id']} seed={seed} {policy} budget={max(CHECKPOINTS)}s", flush=True)
    with (out / "stdout.log").open("x", encoding="utf8") as log:
        start = perf_counter_ns()
        error = None
        try:
            proc = subprocess.run([*args, str(start)], cwd=CODE_ROOT, env=env, stdout=log,
                                  stderr=subprocess.STDOUT, check=False, timeout=max(CHECKPOINTS) + 600)
            returncode = proc.returncode
        except (OSError, subprocess.TimeoutExpired) as exc:
            returncode, error = None, str(exc)
        finish = perf_counter_ns()
    process = dict(launch_ns=start, exit_observed_ns=finish, elapsed_s=(finish-start)/1e9,
                   returncode=returncode, error=error,
                   environment={k: env[k] for k in (*locations, "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                  "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE")})
    save(out / "process.json", process, exclusive=True)
    check_inputs(case)
    print(f"FINISH {case['case_id']} seed={seed} {policy} rc={returncode} elapsed={process['elapsed_s']:.2f}s", flush=True)


def audit(study):
    protocol = read(study / "frozen.json")
    rows = []
    previous_exit = None
    for index, case in enumerate(protocol["cases"]):
        for seed in SEEDS:
            baselines = []
            for policy in arm_order(index, seed):
                out = study / "live" / case["case_id"] / str(seed) / policy
                row = dict(case_id=case["case_id"], seed=seed, policy=policy, status="not_run")
                try:
                    if not (out / "process.json").exists():
                        rows.append(row)
                        continue
                    process = read(out / "process.json")
                    row["process"] = process
                    require(previous_exit is None or previous_exit <= process["launch_ns"],
                            "study workers overlapped or registered ordering changed")
                    previous_exit = process["exit_observed_ns"]
                    if process["returncode"] not in (0, 1):
                        row["status"] = "process_failed"
                        rows.append(row)
                        continue
                    record = read(out / "result.json")
                    require(record["policy"] == policy and record["seed"] == seed, "result case mismatch")
                    require(record["input_root"] == case["root"], "result input mismatch")
                    require(record["launch_ns"] == process["launch_ns"], "launch clock changed")
                    require(record["budget_s"] == max(CHECKPOINTS), "result budget changed")
                    require(record["campaign_sha256"] == digest(study / "frozen.json"), "result protocol mismatch")
                    raw = [json.loads(line) for line in (out / "events.jsonl").read_text(encoding="utf8").splitlines()]
                    validate_record(record, raw)
                    validate_decisions(record)
                    validate_initialization(record)
                    require((process["returncode"] == 1) == (record["stop"] == "failed"), "return code and stop disagree")
                    require(record["elapsed_s"] <= process["elapsed_s"], "worker duration exceeds observed process")
                    row.update(status="arm_failed" if record["stop"] == "failed" else "completed",
                        stop=record["stop"], error=record.get("error"), result_sha256=digest(out / "result.json"),
                        checkpoints={str(t): reconstruct(record, t) for t in CHECKPOINTS},
                        late_events=[e for e in record["events"] if e["end_s"] > max(CHECKPOINTS)],
                        event_count=len(record["events"]))
                    if record.get("baseline") is not None:
                        baselines.append(record["baseline"])
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    row.update(status="audit_failed", error=str(exc))
                rows.append(row)
            require(all(b == baselines[0] for b in baselines), "methods disagree on original reference baseline")
    result = dict(schema="nf-anytime-audit.v1", campaign_sha256=digest(study / "frozen.json"), rows=rows,
                  status_counts=dict(Counter(r["status"] for r in rows)))
    save(study / "audit.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--source-study", type=Path)
    parser.add_argument("--case")
    parser.add_argument("--seed", type=int, choices=SEEDS)
    parser.add_argument("--policy", choices=POLICIES)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--skip-completed", action="store_true")
    parser.add_argument("--freeze-only", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--started-ns", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    study = args.study.resolve()
    require(study.drive.lower() == "d:", "study must be stored on D as requested")
    if args.worker:
        require(all(v is not None for v in (args.case, args.seed, args.policy, args.started_ns)), "missing worker arguments")
        raise SystemExit(worker(study, args.case, args.seed, args.policy, args.started_ns))
    require(args.started_ns is None, "started-ns is worker-only")
    source = args.source_study or Path(read(study / "frozen.json")["source_study"])
    protocol = freeze(study, source)
    if args.case:
        require(args.case in {c["case_id"] for c in protocol["cases"]}, "unregistered case")
    if args.freeze_only:
        print("Frozen", study / "frozen.json", flush=True)
        return
    if not args.audit_only:
        with process_lock(study / "campaign.lock"), keep_awake():
            for index, case in enumerate(protocol["cases"]):
                if args.case and case["case_id"] != args.case:
                    continue
                for seed in SEEDS:
                    if args.seed is not None and seed != args.seed:
                        continue
                    for policy in arm_order(index, seed):
                        if args.policy and policy != args.policy:
                            continue
                        freeze(study, source)
                        out = study / "live" / case["case_id"] / str(seed) / policy
                        if args.skip_completed and (out / "process.json").exists():
                            rows = audit(study)["rows"]
                            old = next(r for r in rows if (r["case_id"], r["seed"], r["policy"]) ==
                                       (case["case_id"], seed, policy))
                            require(old["status"] in ("completed", "arm_failed", "process_failed"), "existing run fails audit")
                            continue
                        run_one(study, case, seed, policy)
                        freeze(study, source)
                        print(audit(study)["status_counts"], flush=True)
    print(audit(study)["status_counts"], flush=True)


if __name__ == "__main__":
    main()
