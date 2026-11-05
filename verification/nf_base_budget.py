"""Cold-design matched-budget pilot on the unmodified DIII-D reference wall.

Actual GPU evaluations supply a common table of costs and outcomes. Each policy
sees only its own completed evaluations during the serial-budget replay.
"""
import argparse
from collections import Counter
from dataclasses import asdict
import gc
import platform
from pathlib import Path
from time import perf_counter

import numpy as np
from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.campaign_geometry import Candidate, evaluator, evaluate_candidate, geometry_check
from d3d_optimization.campaign_search import scales
from d3d_optimization.imp_tangent import project_imp
from d3d_optimization.reference_transport import ReferencePoint
from run_tangent_curvature_study import project_floor
from verification.nf_search_benchmark import ROOT, measurement, prep

SEEDS = (5161, 7753, 9923)
SOURCES = ("OMP", "IMP")
SIZE = 12
BUDGET = 200
BATCH = 4
RADIUS = 20
SHOT = "input/D3D/d3d/d3d_shot_1a3155370c.parquet"


def replay(rows, policy, low_setup, high_setup, baseline, budget=BUDGET):
    if policy not in ("screen", "direct"):
        raise ValueError(policy)
    elapsed = high_setup
    result = dict(policy=policy, budget_s=budget, peak=None, A2=None, selected=None,
                  first_5pct_s=None, fine_ids=[], seen_coarse=[], eligible_fine=0,
                  elapsed_last_completed_s=0, best_eligible_new=None)
    if elapsed > budget:
        return result
    result.update(peak=baseline, elapsed_last_completed_s=elapsed)
    if policy == "screen":
        elapsed += low_setup
    if elapsed > budget:
        return result

    def confirm(row):
        nonlocal elapsed
        state = row["fine"]
        elapsed += state["seconds"]
        if elapsed > budget:
            return False
        result["elapsed_last_completed_s"] = elapsed
        result["fine_ids"].append(row["i"])
        if state["eligible"]:
            result["eligible_fine"] += 1
            best = result["best_eligible_new"]
            if best is None or state["peak"] < best["peak"]:
                result["best_eligible_new"] = dict(i=row["i"], peak=state["peak"], A2=state["A2"])
            if state["peak"] < result["peak"]:
                result.update(selected=row["i"], peak=state["peak"], A2=state["A2"])
            if state["peak"] <= baseline * .95 and result["first_5pct_s"] is None:
                result["first_5pct_s"] = elapsed
        return True

    for start in range(0, len(rows), BATCH):
        group = rows[start:start+BATCH]
        seen = []
        for row in group:
            elapsed += row["generation_s"]
            if policy == "direct":
                if not confirm(row):
                    return result
            else:
                elapsed += row["coarse"]["seconds"]
                if elapsed > budget:
                    return result
                result["elapsed_last_completed_s"] = elapsed
                result["seen_coarse"].append(row["i"])
                if row["coarse"]["eligible"]:
                    seen.append(row)
        if policy == "screen" and seen:
            winner = min(seen, key=lambda r: (r["coarse"]["peak"], r["i"]))
            if not confirm(winner):
                break
    return result


def propose(e, rng):
    source = e.metadata["source"]
    if source == "IMP":
        target = rng.uniform(-.024, .024, 8)
        return project_imp(e, target, rng.uniform(.3, .7))
    floor = project_floor(e, rng.uniform(-.045, .045, 4))
    if floor is None:
        return None
    # Other editable OMP faces have different chord lengths and curvature limits.
    small = rng.uniform(-1, 1, 12) * scales(e)[:12] * rng.uniform(.06, .45)
    return np.r_[small, floor]


def panel(e, seed):
    source = e.metadata["source"]
    rng = np.random.default_rng(seed)
    accepted = []
    rejected = Counter()
    last = perf_counter()
    for attempt in range(1, 601):
        x = propose(e, rng)
        if x is None:
            rejected["projection"] += 1
            continue
        c = Candidate(source, RADIUS, tuple(x))
        co = c.bind(e)
        if not e.constraints(co)["passed"] or max(e.precise_lengths(co)[0].values()) > 1.2:
            rejected["geometry_limits"] += 1
            continue
        if not geometry_check(e, co):
            rejected["whole_wall"] += 1
            continue
        # Remove near-duplicate displacement vectors independent of thermal outcome.
        active = x[:8] if source == "IMP" else x
        if any(np.max(abs(active-np.asarray(r["candidate"]["vector"])[:len(active)])) < .001
               for r in accepted):
            rejected["near_duplicate"] += 1
            continue
        now = perf_counter()
        accepted.append(dict(i=len(accepted), candidate=c.record(), id=c.identity,
                             generation_s=now-last, cumulative_draws=attempt))
        last = now
        if len(accepted) == SIZE:
            break
    return dict(rows=accepted, rejected=dict(rejected), attempts=attempt,
                generation_s=sum(row["generation_s"] for row in accepted),
                trailing_rejected_s=perf_counter()-last)


def score(e, candidate):
    t = perf_counter()
    out = evaluate_candidate(e, candidate)
    measured = measurement(out, e.reference)
    return dict(seconds=perf_counter()-t, **measured)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    ref = ReferencePoint()
    save(out/"protocol.json", dict(seeds=SEEDS, sources=SOURCES, size=SIZE, budget_s=BUDGET,
         batch=BATCH, radius_mm=RADIUS, reference=asdict(ref), python=platform.python_version(),
         driver_sha256=digest(__file__), inputs={p:digest(ROOT/p) for p in (SHOT,
         "baseline_review/baseline_review.json", "input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json")}),
         exclusive=True)
    low = {}
    low_setup = {}
    # Freeze ALL source/seed candidate orders before reading any candidate heat scores.
    for source in SOURCES:
        original, sec = prep(ROOT, source, ref, 0)
        t = perf_counter()
        low[source] = evaluator(original, RADIUS, "gpu")
        low_setup[source] = sec + perf_counter()-t
        for seed in SEEDS:
            frozen = panel(low[source], seed)
            save(out/source/str(seed)/"panel.json", frozen, exclusive=True)
            print("frozen", source, seed, len(frozen["rows"]), frozen["attempts"], flush=True)
        del original
    for source in SOURCES:
        original, sec = prep(ROOT, source, ref, 2)
        baseline = measurement(original.baseline, ref)
        t = perf_counter()
        high = evaluator(original, RADIUS, "gpu")
        high_setup = sec + perf_counter()-t
        save(out/source/"setup.json", dict(high_s=high_setup, low_s=low_setup[source],
             baseline=baseline), exclusive=True)
        parity = None
        for seed in SEEDS:
            folder = out/source/str(seed)
            frozen = load_json(folder/"panel.json")
            rows = []
            for item in frozen["rows"]:
                candidate = Candidate.from_record(item["candidate"])
                row = dict(item)
                # Alternate measurement order; replay charges each observed cost once per arm.
                levels = [("coarse", low[source]), ("fine", high)]
                if item["i"] % 2:
                    levels.reverse()
                for name, e in levels:
                    row[name] = score(e, candidate)
                if parity is None and row["coarse"]["eligible"]:
                    cpu = evaluator(low[source], RADIUS, "cpu")
                    check = score(cpu, candidate)
                    if not check["eligible"]:
                        raise ValueError("CPU/GPU eligibility mismatch")
                    diff = abs(check["peak"]/row["coarse"]["peak"]-1)
                    if diff > 1e-10:
                        raise ValueError("CPU/GPU peak mismatch")
                    parity = dict(id=candidate.identity, relative_peak=diff, CPU_s=check["seconds"],
                                  GPU_s=row["coarse"]["seconds"])
                    save(out/source/"parity.json", parity, exclusive=True)
                    del cpu
                rows.append(row)
                save(folder/f"pair_{row['i']:02}.json", row, exclusive=True)
                print("measured", source, seed, row["i"], row["fine"].get("peak"), flush=True)
            policies = {name:replay(rows, name, low_setup[source], high_setup, baseline["peak"])
                        for name in ("screen","direct")}
            for record in policies.values():
                if record["selected"] is None and record["peak"] is not None:
                    record["A2"] = baseline["A2"]
            feasible = [r for r in rows if r["fine"]["eligible"]]
            best = min(feasible, key=lambda r:r["fine"]["peak"], default=None)
            finite = min(baseline["peak"], best["fine"]["peak"]) if best else baseline["peak"]
            save(folder/"result.json", dict(panel_sha256=digest(folder/"panel.json"), policies=policies,
                 finite_best_including_baseline=finite, best_new=best["i"] if best else None,
                 fine_eligible=len(feasible), candidates=len(rows)), exclusive=True)
        del high, original
        gc.collect()


if __name__ == "__main__":
    main()
