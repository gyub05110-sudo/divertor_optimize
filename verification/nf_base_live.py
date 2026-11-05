"""Serial wall-clock confirmation of the frozen DIII-D cold-design benchmark."""
import argparse
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
import os
from pathlib import Path
import platform
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter_ns

import numpy as np
from verification import nf_base_budget as frozen
from verification.validate_boundary_transport import LEVELS

SEED = 12547
SOURCE = "OMP"
MAX_DRAWS = 600
FROZEN_SHA = "7940f20406c755b14ca3c52cc159348d9599b2a18f3536eedd727db7e1b01c28"
ARMS = ("direct", "screen")
NS_PER_S = 1_000_000_000


def draw(geometry, rng, accepted):
    """One draw of frozen.panel, so a deadline can interrupt rejected attempts."""
    x = frozen.propose(geometry, rng)
    if x is None:
        return None, "projection"
    candidate = frozen.Candidate(SOURCE, frozen.RADIUS, tuple(x))
    co = candidate.bind(geometry)
    if not geometry.constraints(co)["passed"] or max(geometry.precise_lengths(co)[0].values()) > 1.2:
        return None, "geometry_limits"
    if not frozen.geometry_check(geometry, co):
        return None, "whole_wall"
    if any(np.max(abs(x - np.asarray(r["candidate"]["vector"]))) < .001 for r in accepted):
        return None, "near_duplicate"
    return candidate, None


def candidate_row(candidate, index, attempt):
    return dict(i=index, candidate=candidate.record(), id=candidate.identity, cumulative_draws=attempt)


def prefix_audit(direct, screen):
    pairs = list(zip(direct, screen))
    equal = all(a["id"] == b["id"] and a["cumulative_draws"] == b["cumulative_draws"]
                for a, b in pairs)
    return dict(length=len(pairs), equal=equal)


class Deadline(Exception):
    """No more work can be admitted by this arm."""


def run_arm(policy, started_ns, *, clock_ns=perf_counter_ns, budget_s=frozen.BUDGET):
    if policy not in ARMS:
        raise ValueError(policy)
    deadline_ns = started_ns + round(budget_s * NS_PER_S)

    def elapsed(t):
        return (t - started_ns) / NS_PER_S

    result = dict(policy=policy, budget_s=budget_s, startup_s=elapsed(clock_ns()), events=[],
                  generated=[], rejected=Counter(), completed_draws=0, coarse_ids=[], fine_ids=[],
                  eligible_fine=0, baseline=None, incumbent=None, best_eligible_new=None,
                  first_5pct_s=None)

    def timed(stage, operation, *args, **labels):
        start = clock_ns()
        if start >= deadline_ns:
            raise Deadline
        event = dict(stage=stage, start_s=elapsed(start), **labels)
        result["events"].append(event)
        try:
            value = operation(*args)
        except Exception:
            event.update(end_s=elapsed(clock_ns()), status="failed")
            raise
        end = clock_ns()
        event["end_s"] = elapsed(end)
        if end > deadline_ns:
            event["status"] = "late_discarded"
            raise Deadline
        event["status"] = "completed"
        return value, event

    def next_row(attempt):
        candidate, reason = draw(high, rng, result["generated"])
        if reason:
            return None, None, reason
        row = candidate_row(candidate, len(result["generated"]), attempt)
        return candidate, row, None

    def score_one(level, evaluator, candidate, row):
        state, event = timed(level, frozen.score, evaluator, candidate, i=row["i"], id=row["id"])
        event["measurement"] = state
        result[level + "_ids"].append(row["i"])
        if level == "fine" and state["eligible"]:
            result["eligible_fine"] += 1
            new = dict(i=row["i"], id=row["id"], peak=state["peak"], A2=state["A2"])
            best = result["best_eligible_new"]
            if best is None or new["peak"] < best["peak"]:
                result["best_eligible_new"] = new
            if new["peak"] < result["incumbent"]["peak"]:
                result["incumbent"] = dict(kind="new_design", **new)
            if new["peak"] <= .95 * result["baseline"]["peak"] and result["first_5pct_s"] is None:
                result["first_5pct_s"] = event["end_s"]
        return state

    try:
        ref = frozen.ReferencePoint()
        rng = np.random.default_rng(SEED)
        (original, _), _ = timed("prepare_fine", frozen.prep, frozen.ROOT, SOURCE, ref, 2)
        baseline, _ = timed("measure_original", frozen.measurement, original.baseline, ref)
        result["baseline"] = baseline
        result["incumbent"] = dict(kind="original_wall", peak=baseline["peak"], A2=baseline["A2"])
        high, _ = timed("build_fine", frozen.evaluator, original, frozen.RADIUS, "gpu")
        del original
        if policy == "screen":
            (original, _), _ = timed("prepare_coarse", frozen.prep, frozen.ROOT, SOURCE, ref, 0)
            low, _ = timed("build_coarse", frozen.evaluator, original, frozen.RADIUS, "gpu")
            del original
        group = []
        for attempt in range(1, MAX_DRAWS + 1):
            (candidate, row, reason), event = timed("draw", next_row, attempt, attempt=attempt)
            result["completed_draws"] = attempt
            if reason:
                result["rejected"][reason] += 1
                event["rejection"] = reason
            else:
                result["generated"].append(row)
                event.update(i=row["i"], id=row["id"])
                if policy == "direct":
                    score_one("fine", high, candidate, row)
                else:
                    state = score_one("coarse", low, candidate, row)
                    group.append((state, candidate, row))
            if policy == "screen" and (len(group) == frozen.BATCH or attempt == MAX_DRAWS):
                eligible = [entry for entry in group if entry[0]["eligible"]]
                if eligible:
                    _, candidate, row = min(eligible, key=lambda entry: (entry[0]["peak"], entry[2]["i"]))
                    score_one("fine", high, candidate, row)
                group.clear()
            if len(result["generated"]) == frozen.SIZE:
                break
        result["stop"] = "candidate_limit" if len(result["generated"]) == frozen.SIZE else "draw_limit"
    except Deadline:
        result["stop"] = "deadline"
    except Exception as exc:
        result.update(stop="failed", error=f"{type(exc).__name__}: {exc}")
    result["elapsed_s"] = elapsed(clock_ns())
    return result


def run_pair(out):
    if frozen.digest(frozen.__file__) != FROZEN_SHA:
        raise ValueError("frozen replay driver checksum mismatch")
    out.mkdir(parents=True, exist_ok=False)
    gpu = subprocess.run(["nvidia-smi", "--query-gpu=index,name,driver_version",
                          "--format=csv,noheader"], check=True, capture_output=True, text=True)
    inputs = (frozen.SHOT, "baseline_review/baseline_review.json",
              "input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json")
    environment = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible and not all(device.strip().lstrip("-").isdigit() for device in visible.split(",")):
        visible = "non-ordinal selector"
    environment["CUDA_VISIBLE_DEVICES"] = visible
    protocol = dict(seed=SEED, source=SOURCE, radius_mm=frozen.RADIUS, budget_s=frozen.BUDGET,
                    size=frozen.SIZE, batch=frozen.BATCH, max_draws=MAX_DRAWS, arm_order=ARMS,
                    levels=dict(coarse=LEVELS[0], fine=LEVELS[2]), backend="gpu",
                    reference=asdict(frozen.ReferencePoint()), frozen_sha256=FROZEN_SHA,
                    driver_sha256=frozen.digest(__file__), inputs={p: frozen.digest(frozen.ROOT/p) for p in inputs},
                    python_version=platform.python_version(),
                    gpu_inventory=gpu.stdout.splitlines(), created_utc=datetime.now(timezone.utc).isoformat(),
                    environment=environment)
    frozen.save(out/"protocol.json", protocol, exclusive=True)
    arms = {}
    for policy in ARMS:
        with TemporaryDirectory(dir=out, prefix=policy + "-cache-") as temp:
            cupy_cache, cuda_cache = Path(temp)/"cupy", Path(temp)/"cuda"
            cupy_cache.mkdir()
            cuda_cache.mkdir()
            env = dict(os.environ, CUPY_CACHE_DIR=str(cupy_cache), CUDA_CACHE_PATH=str(cuda_cache))
            args = [sys.executable, "-m", "verification.nf_base_live", "--output", str(out),
                    "--arm", policy, "--started-ns"]
            print(f"Starting {policy} ({frozen.BUDGET} s)", flush=True)
            start = perf_counter_ns()
            proc = subprocess.run([*args, str(start)], env=env, cwd=frozen.ROOT)
            finish = perf_counter_ns()
            frozen.save(out/f"{policy}_process.json", dict(returncode=proc.returncode,
                        launch_ns=start, exit_observed_ns=finish, elapsed_s=(finish-start)/NS_PER_S),
                        exclusive=True)
            proc.check_returncode()
        arms[policy] = frozen.load_json(out/f"{policy}.json")
    audit = prefix_audit(arms["direct"]["generated"], arms["screen"]["generated"])
    endpoints = ("stop", "baseline", "incumbent", "best_eligible_new", "first_5pct_s",
                 "coarse_ids", "fine_ids", "eligible_fine", "completed_draws", "rejected", "elapsed_s")
    comparison = {p: {key: arms[p][key] for key in endpoints} for p in ARMS}
    frozen.save(out/"comparison.json", dict(shared_prefix=audit, arms=comparison), exclusive=True)
    if not audit["equal"]:
        raise RuntimeError("candidate geometry prefix mismatch; see arm records")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--arm", choices=ARMS, help=argparse.SUPPRESS)
    parser.add_argument("--started-ns", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    out = args.output.resolve()
    if (args.arm is None) != (args.started_ns is None):
        parser.error("--arm and --started-ns must be supplied together by the supervisor")
    if args.arm is None:
        run_pair(out)
    else:
        record = run_arm(args.arm, args.started_ns)
        record.update(pid=os.getpid(), launch_ns=args.started_ns)
        frozen.save(out/f"{args.arm}.json", record, exclusive=True)
        if record["stop"] == "failed":
            raise RuntimeError(record["error"])


if __name__ == "__main__":
    main()
