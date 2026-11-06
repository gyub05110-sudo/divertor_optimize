"""Serial fresh-process paired trials on preregistered, separately staged equilibria.

Input: STUDY/preflight.json with a cases list containing case_id, status and root.
Only status='passed' (or passed=true) cases run. Roots may be study-relative.
The original nf_base_live policy and nf_base_budget generator remain unchanged.
"""
import argparse
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter_ns

from d3d_optimization.broadening_study import digest, load_json as read, save
from d3d_optimization.independent_sources import split_design_wall
from verification import nf_base_live as live

CODE_ROOT = Path(__file__).resolve().parents[1]
SEEDS = (20611, 20627, 20639)
BUDGET = 200
SHOT = "input/D3D/d3d/d3d_shot_1a3155370c.parquet"
REVIEW = "baseline_review/baseline_review.json"
GEOMETRY = "input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json"
MANIFEST = "input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json"
INPUTS = (SHOT, REVIEW, GEOMETRY, MANIFEST)
ENDPOINTS = ("stop", "error", "baseline", "incumbent", "best_eligible_new", "first_5pct_s",
             "coarse_ids", "fine_ids", "eligible_fine", "completed_draws", "rejected",
             "startup_s", "elapsed_s")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def input_record(root):
    """Check the actual physical baseline, independently of staged metadata claims."""
    root = Path(root).resolve()
    hashes = {name: digest(root / name) for name in INPUTS}
    review, original = read(root / REVIEW), read(CODE_ROOT / REVIEW)
    require(review["wall"] == original["wall"], "staged supplied material wall differs")
    require(hashes[GEOMETRY] == digest(CODE_ROOT / GEOMETRY), "supplied geometry bytes differ")
    require(review["optimization_policy"] == original["optimization_policy"], "engineering policy differs")
    active = {"omp_entry", "shelf", "omp_connector", "floor"}

    def endpoints(item):
        return {s["name"]: s["endpoints_RZ_m"] for s in split_design_wall(item)[1]
                if s["name"] in active}

    require(endpoints(review) == endpoints(original), "supplied OMP design endpoints differ")
    manifest = read(root / MANIFEST)
    require(manifest["shot_sha256"] == hashes[SHOT], "staged shot hash mismatch")
    require(math.isfinite(manifest["time_ms"]), "finite equilibrium time required")
    return dict(root=str(root), time_ms=manifest["time_ms"], hashes=hashes,
                supplied_wall_equal=True, supplied_geometry_equal=True, omp_endpoints_equal=True)


def load_cases(study):
    study = Path(study).resolve()
    rows = read(study / "preflight.json")["cases"]
    ids = [r["case_id"] for r in rows]
    require(len(set(ids)) == len(ids), "duplicate preflight case_id")
    cases = []
    downloads = read(study / "downloads.json")["files"]
    for row in rows:
        if row.get("status") != "passed" and row.get("passed") is not True:
            continue
        case_id = row["case_id"]
        require(case_id and all(c.isalnum() or c in "_-" for c in case_id), "unsafe case_id")
        root = Path(row["root"])
        if not root.is_absolute():
            root = study / root
        info = input_record(root)
        expected = {name.replace("\\", "/"): value for name, value in row["inputs"].items()}
        require(expected == info["hashes"], "preflight input hashes differ from staged bytes")
        source = row["source"]
        require(sum(entry == source for entry in downloads) == 1, "preflight source differs from download record")
        raw = Path(source["local"])
        require(raw.stat().st_size == source["bytes"] and digest(raw) == source["sha256"]
                == info["hashes"][SHOT], "raw download or staged shot hash differs from preflight source")
        require(row.get("time_ms", info["time_ms"]) == info["time_ms"], "preflight time mismatch")
        cases.append(dict(case_id=case_id, **info))
    return cases


def arm_order(case_index, seed):
    order = ("direct", "screen")
    return order if (case_index * len(SEEDS) + SEEDS.index(seed)) % 2 == 0 else order[::-1]


def prepare_campaign(study):
    study = Path(study).resolve()
    selection = read(study / "selection.json")
    require(selection["seeds"] == list(SEEDS) and selection["budget_s"] == BUDGET,
            "selection seeds or budget differ from preregistered protocol")
    require(selection["protocol_sha256"] == digest(study / "protocol.ko.md"),
            "selection/protocol hash mismatch")
    preflight = read(study / "preflight.json")
    for key, path in (("protocol_sha256", study / "protocol.ko.md"),
                      ("selection_sha256", study / "selection.json"),
                      ("downloads_sha256", study / "downloads.json"),
                      ("driver_sha256", CODE_ROOT / "verification/nf_unused_inputs.py"),
                      ("geometry_driver_sha256", CODE_ROOT / "verification/nf_equilibrium_times.py")):
        require(preflight[key] == digest(path), f"preflight {key} mismatch")
    for name, expected in preflight.get("original_sha256", {}).items():
        require(digest(CODE_ROOT / name.replace("\\", "/")) == expected,
                "preflight original input changed: " + name)
    require(digest(live.frozen.__file__) == live.FROZEN_SHA, "frozen generator checksum mismatch")
    dependencies = list((CODE_ROOT / "d3d_optimization").glob("*.py"))
    dependencies += [Path(__file__), Path(live.__file__), Path(live.frozen.__file__),
                     CODE_ROOT / "verification/nf_search_benchmark.py",
                     CODE_ROOT / "verification/validate_boundary_transport.py",
                     CODE_ROOT / "verification/validate_reference_transport.py",
                     CODE_ROOT / "optimize_boundary_transport.py",
                     CODE_ROOT / "run_tangent_curvature_study.py"]
    protocol = dict(schema="unused-equilibrium-live.v1", seeds=list(SEEDS), budget_s=BUDGET,
        source="OMP", radius_mm=live.frozen.RADIUS, size=live.frozen.SIZE,
        batch=live.frozen.BATCH, max_draws=live.MAX_DRAWS, backend="gpu",
        levels=dict(coarse=list(live.LEVELS[0]), fine=list(live.LEVELS[2])),
        reference=asdict(live.frozen.ReferencePoint()), cases=load_cases(study),
        source_hashes={str(p.relative_to(CODE_ROOT)).replace("\\", "/"): digest(p)
                       for p in sorted(set(dependencies))},
        study_hashes={name: digest(study / name) for name in
                      ("protocol.ko.md", "selection.json", "downloads.json", "preflight.json")},
        python_version=platform.python_version(),
        arm_order_rule="alternate parity of selected case index * 3 + seed index",
        budget_exclusions=["source acquisition and staging", "geometry admission preflight",
                           "parent hashing and GPU inventory", "empty cache directory creation",
                           "result serialization and cache cleanup"],
        cache_policy="independent empty CUDA and CuPy disk caches for every fresh worker",
        startup_policy="parent monotonic timestamp immediately before subprocess launch")
    path = study / "live" / "protocol.json"
    if path.exists():
        old = read(path)
        old.pop("created_utc", None)
        require(old == protocol, "campaign protocol/code/inputs changed; refusing mixed continuation")
    else:
        save(path, dict(protocol, created_utc=datetime.now(timezone.utc).isoformat()), exclusive=True)
    return protocol


def run_worker(out, input_root, seed, policy, started_ns):
    """Keep launch clock unchanged; all imports and configuration count as startup."""
    require(seed in SEEDS, "unregistered seed")
    live.SEED, live.SOURCE = seed, "OMP"
    live.frozen.ROOT = Path(input_root).resolve()
    record = live.run_arm(policy, started_ns, budget_s=BUDGET)
    record.update(pid=os.getpid(), launch_ns=started_ns, seed=seed,
                  input_root=str(live.frozen.ROOT))
    save(Path(out) / f"{policy}.json", record, exclusive=True)
    return 1 if record["stop"] == "failed" else 0


def audit_arm(record):
    """Reconstruct endpoints from on-time completed events; never trust late values."""
    budget = record["budget_s"]
    require(budget == BUDGET, "arm budget differs")
    require(math.isfinite(record["startup_s"]) and record["startup_s"] >= 0, "invalid startup")
    last = record["startup_s"]
    fine, coarse, draws = [], [], []
    phases = Counter()
    terminal = False
    for event in record["events"]:
        start, end = event["start_s"], event["end_s"]
        require(not terminal, "work started after failed or late event")
        require(math.isfinite(start) and math.isfinite(end) and last <= start < budget and end >= start,
                "invalid or overlapping event clock")
        last = end
        status = event["status"]
        require(status in ("completed", "failed", "late_discarded"), "unknown event status")
        if status == "completed":
            require(end <= budget, "late event recorded as completed")
            if event["stage"] == "fine":
                fine.append(event)
            if event["stage"] == "coarse":
                coarse.append(event)
            if event["stage"] == "draw":
                draws.append(event)
        else:
            terminal = True
            require("measurement" not in event, "failed or late measurement leaked")
            if status == "late_discarded":
                require(end > budget, "on-time event marked late")
        phases[event["stage"]] += end - start
    require(record["elapsed_s"] >= last, "arm elapsed time precedes its events")
    generated = record["generated"]
    accepted = [e for e in draws if "id" in e]
    require(len(accepted) == len(generated), "accepted draws disagree with generated candidates")
    for i, (event, row) in enumerate(zip(accepted, generated)):
        require(row["i"] == i and event["i"] == i and event["id"] == row["id"]
                and event["attempt"] == row["cumulative_draws"]
                and live.frozen.Candidate.from_record(row["candidate"]).identity == row["id"],
                "candidate identity or draw order mismatch")
    require(record["completed_draws"] == len(draws), "completed draw count differs")
    require(dict(Counter(e["rejection"] for e in draws if "rejection" in e)) == record["rejected"],
            "draw rejection counts differ")
    for level, events in (("fine", fine), ("coarse", coarse)):
        require(record[level + "_ids"] == [e["i"] for e in events], level + " IDs differ")
        for event in events:
            require(0 <= event["i"] < len(generated)
                    and generated[event["i"]]["id"] == event["id"], "score candidate mismatch")
    measured_original = any(e["stage"] == "measure_original" and e["status"] == "completed"
                            for e in record["events"])
    baseline = record["baseline"]
    require((baseline is not None) == measured_original, "baseline missing or admitted before completion")
    incumbent = None if baseline is None else dict(kind="original_wall", peak=baseline["peak"], A2=baseline["A2"])
    best, first = None, None
    eligible = 0
    for event in fine:
        state = event["measurement"]
        if not state["eligible"]:
            continue
        eligible += 1
        require(baseline is not None, "fine result precedes baseline")
        require(all(math.isfinite(state[k]) and state[k] > 0 for k in ("peak", "A2")),
                "nonpositive/nonfinite eligible endpoint")
        new = dict(i=event["i"], id=event["id"], peak=state["peak"], A2=state["A2"])
        if best is None or new["peak"] < best["peak"]:
            best = new
        if new["peak"] < incumbent["peak"]:
            incumbent = dict(kind="new_design", **new)
        if first is None and new["peak"] <= .95 * baseline["peak"]:
            first = event["end_s"]
    require(record["incumbent"] == incumbent, "incumbent differs from on-time reconstruction")
    require(record["best_eligible_new"] == best, "best new candidate differs")
    require(record["first_5pct_s"] == first and record["eligible_fine"] == eligible,
            "improvement timestamp or eligibility count differs")
    return dict(valid=True, phase_seconds=dict(phases), late_events=[e for e in record["events"]
                if e["status"] == "late_discarded"], eligible_fine=eligible)


def audit_pair(out):
    out = Path(out)
    protocol = read(out / "protocol.json")
    arms, audits, errors, processes = {}, {}, {}, {}
    for policy in ("direct", "screen"):
        path = out / f"{policy}.json"
        process_path = out / f"{policy}_process.json"
        if not path.exists():
            errors[policy] = "worker produced no arm record"
            continue
        try:
            record = read(path)
            processes[policy] = read(process_path)
            require(record["policy"] == policy and record["seed"] == protocol["seed"]
                    and record["input_root"] == protocol["case"]["root"], "worker configuration differs")
            proc = processes[policy]
            require(record["launch_ns"] == proc["launch_ns"] and proc["elapsed_s"] >= record["elapsed_s"],
                    "launch or process elapsed mismatch")
            require(proc["returncode"] == (1 if record["stop"] == "failed" else 0),
                    "process exit disagrees with saved arm status")
            audits[policy] = audit_arm(record)
            arms[policy] = record
        except (ValueError, KeyError, TypeError, OSError, IndexError) as exc:
            errors[policy] = str(exc)
    if len(processes) == 2:
        try:
            first, second = (processes[p] for p in protocol["arm_order"])
            require(first["exit_observed_ns"] <= second["launch_ns"],
                    "workers overlap or violate registered serial order")
        except (ValueError, KeyError, TypeError) as exc:
            errors["arm_order"] = str(exc)
    prefix = live.prefix_audit(arms["direct"]["generated"], arms["screen"]["generated"]) if len(arms) == 2 else None
    if prefix is not None:
        prefix["informative"] = prefix["length"] > 0
        if not prefix["equal"]:
            errors["shared_prefix"] = "candidate geometry prefix mismatch"
    baselines = [r["baseline"] for r in arms.values() if r["baseline"] is not None]
    baseline_equal = None if len(baselines) != 2 else baselines[0] == baselines[1]
    if baseline_equal is False:
        errors["baseline"] = "paired original baseline differs"
    status = "audit_failed" if errors else "completed"
    if not errors and any(r["stop"] == "failed" for r in arms.values()):
        status = "arm_failed"
    if any(isinstance(p, dict) and p.get("returncode") not in (0, 1) for p in processes.values()):
        status = "process_failed"
    return dict(case_id=protocol["case"]["case_id"], seed=protocol["seed"], status=status,
                arm_order=protocol["arm_order"], shared_prefix=prefix, baseline_equal=baseline_equal,
                arms={p: {k: r.get(k) for k in ENDPOINTS} for p, r in arms.items()},
                audits=audits, errors=errors, processes=processes)


def run_pair(study, case, seed, case_index):
    study = Path(study).resolve()
    require(seed in SEEDS, "unregistered seed")
    require(input_record(case["root"]) == {k: v for k, v in case.items() if k != "case_id"},
            "case inputs changed before pair")
    out = study / "live" / case["case_id"] / str(seed)
    out.mkdir(parents=True, exist_ok=False)
    order = arm_order(case_index, seed)
    try:
        gpu = subprocess.run(["nvidia-smi", "--query-gpu=index,name,driver_version",
                              "--format=csv,noheader"], capture_output=True, text=True, check=False)
        inventory = dict(returncode=gpu.returncode, rows=gpu.stdout.splitlines())
    except OSError as exc:
        inventory = dict(error=str(exc))
    save(out / "protocol.json", dict(case=case, seed=seed, arm_order=list(order), budget_s=BUDGET,
         campaign_protocol_sha256=digest(study / "live/protocol.json"), gpu_inventory=inventory,
         created_utc=datetime.now(timezone.utc).isoformat()), exclusive=True)
    for policy in order:
        with TemporaryDirectory(dir=out, prefix=policy + "-cache-") as temp:
            cache = Path(temp)
            locations = {name: cache / leaf for name, leaf in
                         (("CUPY_CACHE_DIR", "cupy"), ("CUDA_CACHE_PATH", "cuda"),
                          ("TEMP", "tmp"), ("TMP", "tmp"), ("MPLCONFIGDIR", "matplotlib"),
                          ("HF_HOME", "huggingface"))}
            for path in locations.values():
                path.mkdir(parents=True, exist_ok=True)
            env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       **{name: str(path) for name, path in locations.items()})
            args = [sys.executable, "-m", "verification.nf_unused_live", "--study", str(study),
                    "--output", str(out), "--input-root", case["root"], "--seed", str(seed),
                    "--arm", policy, "--started-ns"]
            print(f"{case['case_id']} seed={seed} {policy}: {BUDGET}s from launch", flush=True)
            with (out / f"{policy}.log").open("x", encoding="utf8") as log:
                start = perf_counter_ns()
                error = None
                try:
                    proc = subprocess.run([*args, str(start)], cwd=CODE_ROOT, env=env,
                                          stdout=log, stderr=subprocess.STDOUT, check=False)
                    returncode = proc.returncode
                except OSError as exc:
                    returncode, error = None, str(exc)
                finish = perf_counter_ns()
            save(out / f"{policy}_process.json", dict(returncode=returncode, error=error,
                 launch_ns=start, exit_observed_ns=finish, elapsed_s=(finish-start)/1e9,
                 environment={k: env[k] for k in (*locations, "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE")}), exclusive=True)
        require(input_record(case["root"]) == {k: v for k, v in case.items() if k != "case_id"},
                "input changed during pair; retained existing records")
    save(out / "comparison.json", audit_pair(out), exclusive=True)
    return out


def audit_campaign(study):
    study = Path(study).resolve()
    protocol = prepare_campaign(study)
    pairs = []
    for case in protocol["cases"]:
        for seed in SEEDS:
            out = study / "live" / case["case_id"] / str(seed)
            if (out / "protocol.json").exists():
                record = audit_pair(out)
                if (out / "comparison.json").exists():
                    require(record == read(out / "comparison.json"), "published pair audit differs")
            else:
                record = dict(case_id=case["case_id"], seed=seed, status="not_run")
            pairs.append(record)
    report = dict(schema="unused-equilibrium-live-audit.v1", pairs=pairs,
                  exclusions=[r for r in read(study / "preflight.json")["cases"]
                              if r.get("status") != "passed" and r.get("passed") is not True])
    save(study / "live" / "audit.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--case")
    parser.add_argument("--seed", type=int, choices=SEEDS)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--skip-completed", action="store_true", help="explicitly preserve existing complete pairs")
    parser.add_argument("--output", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--input-root", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--arm", choices=("direct", "screen"), help=argparse.SUPPRESS)
    parser.add_argument("--started-ns", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.arm is not None:
        if args.output is None or args.input_root is None or args.seed is None or args.started_ns is None:
            parser.error("workers require output, input-root, seed and started-ns")
        raise SystemExit(run_worker(args.output, args.input_root, args.seed, args.arm, args.started_ns))
    if any(v is not None for v in (args.output, args.input_root, args.started_ns)):
        parser.error("worker arguments require --arm")
    protocol = prepare_campaign(args.study)
    if args.case and args.case not in {r["case_id"] for r in protocol["cases"]}:
        parser.error("case is not a passed preregistered case")
    if not args.audit_only:
        for index, case in enumerate(protocol["cases"]):
            if args.case and case["case_id"] != args.case:
                continue
            for seed in SEEDS:
                if args.seed is not None and seed != args.seed:
                    continue
                out = args.study.resolve() / "live" / case["case_id"] / str(seed)
                if args.skip_completed and (out / "comparison.json").exists():
                    require(audit_pair(out) == read(out / "comparison.json"), "existing pair audit differs")
                    print(f"Preserved completed pair {case['case_id']} seed={seed}", flush=True)
                    continue
                run_pair(args.study, case, seed, index)
    report = audit_campaign(args.study)
    print(dict(Counter(r["status"] for r in report["pairs"])), flush=True)


if __name__ == "__main__":
    main()
