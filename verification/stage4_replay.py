"""Post-budget fixed-wall numerical sensitivity for the preregistered Stage 4 contrasts."""
import argparse
from math import isclose
from pathlib import Path

from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.campaign_state import process_lock
from d3d_optimization.user_output import file_hash
from verification import stage3_replay as fixed
from verification import stage4_campaign as design
from verification import stage4_report as report
from verification.nf_prospective_inputs import publish
from verification.stage3_report import positive, read, require


def plan_case(study, case, source, report_hash):
    seed = design.SEEDS[0]
    first = [r for r in source["runs"] if r["case_id"] == case and r["seed"] == seed]
    runs = {r["arm"]: r for r in first}
    require(len(first) == len(design.ALL) and set(runs) == set(design.ALL), "first-seed run inventory differs")
    product = report.product_pair(runs["hotspot"], runs["full"], 600)
    saved_pairs = [p for p in source["pairs"] if p["case_id"] == case and p["seed"] == seed and p["cutoff"] == 600]
    require(saved_pairs == [product], "reported product comparison differs from timed runs")
    prediction = report.prediction_case(first)
    require([p for p in source["predictions"] if p["case_id"] == case] == [prediction],
            "reported forced-family comparison differs from timed runs")

    gap = product["gap_pct"]
    discordant = product["auto_5pct"] != product["full_5pct"]
    reasons = []
    if product["practical"] in ("auto_lower", "full_lower"):
        reasons.append("absolute_gap_over_1pct")
    if gap is not None and discordant:
        reasons.append("discordant_5pct")
    regret = prediction["selected_forced_regret_pct"]
    best = next((arm for arm in design.FOCUSED if arm in prediction["best_families"]), None)
    pairs = [dict(contrast="product", seed=seed, left="hotspot", right="full", trigger=bool(reasons),
                  reasons=reasons, availability=product["availability"], timed_gap_pct=gap, timed_practical=product["practical"],
                  timed_left_5pct=product["auto_5pct"], timed_right_5pct=product["full_5pct"]),
             dict(contrast="forced_family", seed=seed, left=prediction["predicted_family"], right=best,
                  trigger=prediction["selected_competitive"] is False, timed_gap_pct=regret, scorable=prediction["scorable"],
                  best_families=prediction["best_families"], forced_status=prediction["forced_status"],
                  forced_qualified_count=prediction["forced_qualified_count"], abstained=prediction["abstained"])]
    common = None
    candidates = {}
    for pair in pairs:
        if not pair["trigger"]:
            continue
        native = {}
        for arm in (pair["left"], pair["right"]):
            run = runs[arm]
            path = study / run["relative_output"] / "result.json"
            require(file_hash(path) == run["result_sha256"], "original timed result digest differs")
            raw = read(path)
            require(raw["arguments"]["seed"] == seed and raw["arguments"]["family"] == arm and
                    raw["candidate"] == run["candidate"] and run["status"] in ("qualified","partial_result"),
                    "wrong timed winner")
            candidate = Candidate.from_record(raw["candidate"])
            checkpoint = run["checkpoints"]["600"]
            require(candidate.identity == checkpoint["id"] == raw["candidate"]["id"], "timed winner identity differs")
            cfg = dict(inputs=raw["inputs"], reference_SI=raw["reference_SI"], radius_mm=raw["arguments"]["radius_mm"],
                       backend=raw["arguments"]["backend"], search_code_sha256=raw["software"]["code_sha256"])
            require(common is None or cfg == common, "fixed operating conditions differ across timed runs")
            common = cfg
            scale = cfg["reference_SI"]["power_W"] / 1e6
            winner_q = checkpoint["q_native"]
            require(positive(winner_q) and isclose(winner_q * scale, raw["selected"]["peak_MW_m2"], rel_tol=1e-12),
                    "timed winner peak differs")
            candidates[candidate.identity] = candidate.record()
            native[arm] = dict(id=candidate.identity, q=winner_q, A2_m2=raw["selected"]["A2_m2"],
                               baseline_q=raw["baseline"]["peak_MW_m2"]/scale, baseline_A2_m2=raw["baseline"]["A2_m2"],
                               result_sha256=run["result_sha256"])
        pair["native"] = native
    return dict(schema="nf-stage4-fixed-replay.v1", case_id=case, search_report_sha256=report_hash,
                replay_driver_sha256=file_hash(__file__), fixed_wall_driver_sha256=file_hash(fixed.__file__),
                first_seed_statuses={a:runs[a]["status"] for a in design.ALL}, common=common, pairs=pairs,
                candidates=candidates, levels={k:list(fixed.LEVELS[k]) for k in ("cubic1280","cubic640","quintic1280")})


def adjudicate(pair, grids):
    if not pair["trigger"]:
        return dict(pair, adjudicated=False)
    native = pair["native"]
    comparisons, reproduction = {}, {}
    for level, grid in grids.items():
        shelf = (grid.get("fixed_shelf") or {}).get("admissible") is True
        baseline = grid.get("baseline") or {}
        base_q = baseline.get("q_native") if baseline.get("physics_measurable") else None
        records, values, reductions = {}, {}, {}
        for arm in (pair["left"], pair["right"]):
            expected = native[arm]
            score = grid.get("candidates", {}).get(expected["id"], {})
            q = score.get("q_native")
            usable = bool(shelf and score.get("candidate_eligible") and score.get("physics_measurable") and positive(q))
            values[arm] = q if usable else None
            reductions[arm] = 100*(base_q-q)/base_q if usable and positive(base_q) else None
            records[arm] = dict(id=expected["id"], q_native=values[arm], A2_m2=score.get("A2_m2"),
                                eligible=usable, reduction_pct=reductions[arm])
            if level == "cubic1280":
                reproduction[arm] = dict(candidate=bool(shelf and fixed.native_matches(score, expected["q"], expected["A2_m2"])),
                                         baseline=fixed.native_matches(baseline, expected["baseline_q"],
                                                                       expected["baseline_A2_m2"], reference=True))
        left, right = values[pair["left"]], values[pair["right"]]
        gap = 100*(left-right)/right if left is not None and right is not None else None
        def reaches(arm):
            return None if reductions[arm] is None else values[arm] <= .95*base_q
        comparisons[level] = dict(status=grid["status"], error=grid.get("error"), fixed_shelf_admissible=shelf,
                                  baseline_q=base_q, baseline_A2_m2=baseline.get("A2_m2"), candidates=records,
                                  gap_pct=gap, sign=None if gap is None else
                                  "left_lower" if gap < 0 else "right_lower" if gap > 0 else "equal",
                                  within_1pct=None if gap is None else .99*right <= left <= 1.01*right,
                                  left_5pct=reaches(pair["left"]), right_5pct=reaches(pair["right"]))
    return dict(pair, adjudicated=True, native_reproduction_relative_tolerance=1e-9, native_reproduction=reproduction,
                native_reproduction_passed=bool(reproduction and all(all(v.values()) for v in reproduction.values())),
                by_level=comparisons)


def load_replay_user(common):
    source = common["inputs"]
    paths = fixed._replay_input_paths(source)
    user = fixed.load_user_input(paths["equilibrium"], source["time_ms"], paths["geometry"], paths["design_policy"])
    require(all(source[k] == user.provenance[k] for k in fixed.INPUT_KEYS), "fixed physical input changed")
    require(fixed.software()["code_sha256"] == common["search_code_sha256"], "solver code differs from timed search")
    return user


def checked_level(folder, name, protocol, user, protocol_sha):
    checkpoint = folder / (name + ".json")
    directory = folder / name
    if not checkpoint.exists():
        grid = fixed.run_level(directory, name, protocol, user)
        paths = [directory / "setup.json", *sorted(directory.glob("*/measurement.json"))]
        hashes = {p.relative_to(directory).as_posix():file_hash(p) for p in paths}
        publish(checkpoint, dict(protocol_sha256=protocol_sha, files_sha256=hashes, grid=fixed.serial(grid)))
    record = read(checkpoint)
    require(record["protocol_sha256"] == protocol_sha, "replay grid protocol digest differs")
    for path, sha in record["files_sha256"].items():
        require(file_hash(directory / path) == sha, "replay record digest differs: " + path)
    grid = record["grid"]
    setup = read(directory / "setup.json")
    require(all(grid[k] == v for k,v in setup.items()), "grid setup differs")
    if grid["status"] == "measured":
        for role in ["baseline", *protocol["candidates"]]:
            saved = grid["baseline"] if role == "baseline" else grid["candidates"][role]
            require(fixed.load_measurement(directory / role) == saved, "fixed-wall measurement differs: " + role)
    return grid


def case_replay(folder, protocol):
    folder.mkdir(parents=True, exist_ok=True)
    with process_lock(folder / "replay.lock"):
        protocol_path = folder / "protocol.json"
        publish(protocol_path, protocol)
        protocol_sha = file_hash(protocol_path)
        grids = {}
        if protocol["common"]:
            user = load_replay_user(protocol["common"])
            for name in protocol["levels"]:
                grids[name] = checked_level(folder, name, protocol, user, protocol_sha)
        result = dict(complete=True, case_id=protocol["case_id"], protocol_sha256=protocol_sha,
                      search_report_sha256=protocol["search_report_sha256"],
                      first_seed_statuses=protocol["first_seed_statuses"],
                      trigger_count=sum(p["trigger"] for p in protocol["pairs"]),
                      fixed_unique_geometry_count=len(protocol["candidates"]), results=grids,
                      comparisons=[adjudicate(p, grids) for p in protocol["pairs"]],
                      grid_checkpoints_sha256={name:file_hash(folder / (name+".json")) for name in grids})
        publish(folder / "report.json", result)
        return result


def run(study):
    study = Path(study).resolve()
    plan = design.load_seal(study)
    report_path = study / "stage4_report.json"
    source = read(report_path)
    require(source["complete"] is True and source["plan_sha256"] == file_hash(study / "expected.json") and
            source["registration_sha256"] == plan["registration_sha256"] and
            source["report_driver_sha256"] == file_hash(report.__file__), "complete audited timed campaign required")
    registry = read(study / "inputs/registration.json")
    cases = {row["id"]:row for row in registry["cases"]}
    require(source["runs"] == report.audit_jobs(study, plan, cases, partial=False), "audited timed run records changed")
    digest = file_hash(report_path)
    dest = study / "replays"
    finished = []
    with process_lock(dest / "replay.lock"):
        for case in cases:
            protocol = plan_case(study, case, source, digest)
            folder = dest / case
            row = case_replay(folder, protocol)
            finished.append(dict(case_id=case, report_sha256=file_hash(folder / "report.json"),
                                 trigger_count=row["trigger_count"], fixed_unique_geometry_count=row["fixed_unique_geometry_count"],
                                 comparisons=row["comparisons"], first_seed_statuses=row["first_seed_statuses"]))
            print(case, row["trigger_count"], "triggered comparisons;", row["fixed_unique_geometry_count"],
                  "unique winners", flush=True)
        summary = dest / "summary.json"
        publish(summary, dict(schema="nf-stage4-replay-summary.v1", complete=True, search_report_sha256=digest,
                              cases=finished))
        publish(study / "stage4_replay.json", summary.read_bytes())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    run(parser.parse_args().study)
