"""Enumerate rank-guided and uniform promotions on the archived DIII-D panels."""
import argparse
from collections import Counter
from itertools import batched, product
import json
from math import isclose, prod
from pathlib import Path
from statistics import fmean, median

from verification import nf_base_budget_audit as archived


ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "results/d3d-base-budget-20261105"
OUTPUT = ROOT / "results/nf-finalist-robustness-20261106/budget-value.json"


def fine_peak(row, baseline):
    """A failed verification consumes the slot. The original wall remains available."""
    return min(baseline, row["fine"]["peak"]) if row["fine"]["eligible"] else baseline


def relative_regret(peak, oracle):
    return 100 * (peak / oracle - 1)


def summarize_panel(folder, setup, batch):
    archived_result = archived.read(folder / "result.json")
    rows = [archived.read(path) for path in sorted(folder.glob("pair_*.json"))]
    baseline = setup["baseline"]["peak"]
    oracle = min(fine_peak(row, baseline) for row in rows)
    archived.check(oracle == archived_result["finite_best_including_baseline"], "finite oracle mismatch")
    blocks = list(batched(rows, batch))
    groups = [[r for r in block if r["coarse"]["eligible"]] for block in blocks]
    archived.check(len(groups) == 3 and all(groups), f"{folder}: undefined three-call policy")
    leaders = [min(group, key=lambda r: (r["coarse"]["peak"], r["i"])) for group in groups]
    score_best = min(fine_peak(r, baseline) for r in leaders)
    replay, _ = archived.replay_scores(rows, setup, "screen", batch, float("inf"))
    leader_ids = [r["i"] for r in leaders]
    archived.check(replay["peak"] == score_best and replay["fine_ids"] == leader_ids,
                   f"{folder}: leaders differ from the original policy")

    common_costs = {
        "fine_setup_including_original_wall_s": setup["high_s"],
        "coarse_setup_s": setup["low_s"],
        "generation_s": sum(r["generation_s"] for r in rows),
        "all_coarse_evals_s": sum(r["coarse"]["seconds"] for r in rows),
    }
    common_time = sum(common_costs.values())
    ranked_time = common_time + sum(r["fine"]["seconds"] for r in leaders)
    # The exhaustive reference can skip all coarse work and the coarse setup.
    exhaustive_time = setup["high_s"] + sum(r["generation_s"] + r["fine"]["seconds"] for r in rows)
    block_best = [min((r["fine"]["peak"] for r in block if r["fine"]["eligible"]), default=None)
                  for block in blocks]
    minima, times = [], []
    for picks in product(*groups):
        minima.append(min(fine_peak(r, baseline) for r in picks))
        times.append(common_time + sum(r["fine"]["seconds"] for r in picks))

    # Check the whole CDF by a second calculation, without enumerating tuples.
    # All selected peaks exceed x iff every group supplied one such candidate.
    counts = Counter(minima)
    total = prod(map(len, groups))
    archived.check(total == len(minima), "random outcome denominator")
    for peak in {baseline, *(fine_peak(r, baseline) for group in groups for r in group)}:
        above = prod(sum(fine_peak(r, baseline) > peak for r in group) for group in groups)
        archived.check(sum(n for value, n in counts.items() if value > peak) == above,
                       f"{folder}: random CDF mismatch at {peak}")
    expected_time = common_time + sum(fmean(r["fine"]["seconds"] for r in group) for group in groups)
    archived.check(isclose(expected_time, fmean(times), abs_tol=1e-10), "random time mismatch")

    score_regret = relative_regret(score_best, oracle)
    random_regrets = [relative_regret(peak, oracle) for peak in minima]
    matched_or_better = sum(peak <= score_best for peak in minima)
    return {
        "source": folder.parent.name, "seed": int(folder.name),
        "baseline_MW_m2": baseline, "panel_oracle_MW_m2": oracle,
        "fine_eligible_n": sum(r["fine"]["eligible"] for r in rows),
        "coarse_eligible_by_group": [len(group) for group in groups],
        "ranked": {
            "fine_ids": leader_ids, "best_MW_m2": score_best, "regret_pct": score_regret,
            "reduction_from_baseline_pct": 100 * (1 - score_best / baseline),
            "fine_block_best_hits": [
                leader["fine"]["eligible"] and leader["fine"]["peak"] == minimum
                for leader, minimum in zip(leaders, block_best, strict=True)
            ],
            "total_replay_s": ranked_time,
        },
        "uniform": {
            "denominator": total,
            "best_distribution": [
                {"best_MW_m2": peak, "regret_pct": relative_regret(peak, oracle), "count": counts[peak]}
                for peak in sorted(counts)
            ],
            "mean_best_MW_m2": fmean(minima), "median_best_MW_m2": median(minima),
            "mean_regret_pct": fmean(random_regrets), "median_regret_pct": median(random_regrets),
            "beats_ranked_count": sum(peak < score_best for peak in minima),
            "ties_ranked_count": counts[score_best],
            "matches_or_beats_ranked_probability": matched_or_better / total,
            "attains_panel_oracle_probability": counts[oracle] / total,
            "min_total_replay_s": min(times), "mean_total_replay_s": fmean(times),
            "max_total_replay_s": max(times),
        },
        "mean_regret_gap_pp": fmean(random_regrets) - score_regret,
        "common_costs": common_costs,
        "exhaustive_12_fine_total_replay_s": exhaustive_time,
        "exhaustive_over_ranked_time_ratio": exhaustive_time / ranked_time,
    }


def analyze():
    paired = STUDY / "paired"
    audit, old_rows = archived.audit(paired, STUDY / "historical")
    archived.check(audit == archived.read(paired / "audit/audit.json"), "archived file audit mismatch")
    protocol = archived.read(paired / "protocol.json")
    panels = []
    for source in protocol["sources"]:
        setup = archived.read(paired / source / "setup.json")
        for seed in protocol["seeds"]:
            panels.append(summarize_panel(paired / source / str(seed), setup, protocol["batch"]))
    old_screen = [r for r in old_rows if r["policy"] == "screen"]
    old_direct = [r for r in old_rows if r["policy"] == "direct"]
    ties = sum((s["source"], s["seed"], s["peak"], s["selected"]) ==
               (d["source"], d["seed"], d["peak"], d["selected"])
               for s, d in zip(old_screen, old_direct, strict=True))
    live = archived.read(STUDY / "live/comparison.json")["arms"]
    return {
        "paired_audit_sha256": archived.sha(paired / "audit/audit.json"),
        "live_comparison_sha256": archived.sha(STUDY / "live/comparison.json"),
        "fine_calls_per_policy_per_panel": 3,
        "panels": panels,
        "equal_panel_mean": {
            "ranked_regret_pct": fmean(p["ranked"]["regret_pct"] for p in panels),
            "uniform_expected_regret_pct": fmean(p["uniform"]["mean_regret_pct"] for p in panels),
            "primary_uniform_minus_ranked_regret_pp": fmean(p["mean_regret_gap_pp"] for p in panels),
            "ranked_attains_oracle_panels": sum(p["ranked"]["regret_pct"] == 0 for p in panels),
            "uniform_matches_or_beats_ranked_probability": fmean(
                p["uniform"]["matches_or_beats_ranked_probability"] for p in panels),
        },
        "archived_200s_replay": {
            "same_final_design_and_peak_panels": ties,
            "screen_fine_calls": sum(r["fine_n"] for r in old_screen),
            "direct_fine_calls": sum(r["fine_n"] for r in old_direct),
            "screen_attains_oracle_panels": sum(r["regret_pct"] == 0 for r in old_screen),
            "direct_attains_oracle_panels": sum(r["regret_pct"] == 0 for r in old_direct),
        },
        "archived_200s_live": {
            name: {
                "best_MW_m2": arm["incumbent"]["peak"], "fine_calls": len(arm["fine_ids"]),
                "coarse_calls": len(arm["coarse_ids"]), "observed_stop_s": arm["elapsed_s"],
                "stop": arm["stop"],
            } for name, arm in live.items()
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="create the result once")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = analyze()
    if args.write:
        with args.output.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        archived.check(result == archived.read(args.output), "recomputation differs from the published result")
    print(f"{len(result['panels'])} panels checked. Mean regret gap "
          f"{result['equal_panel_mean']['primary_uniform_minus_ranked_regret_pp']:.6f} pp.")


if __name__ == "__main__":
    main()
