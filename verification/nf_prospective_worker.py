"""Configured cold SearchCase workers and independent adjudication of prospective arms.

The supervisor starts one Python process per arm and supplies its launch clock.
Raw observations include work crossing the deadline. The adjudicator replays the
active DE target sequence and selects the best eligible, completed fine return.
"""
import argparse
from collections import OrderedDict
from dataclasses import asdict
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

from d3d_optimization.reference_transport import ReferencePoint
from verification import nf_anytime_final_audit as evidence


CONFIGURATIONS = (("coarse_adaptive", 64, 12), ("fine_adaptive", 24, 8), ("fine_adaptive", 32, 12))
DEVELOPMENT_CASES = ("deb842cf1e", "fad1624516", "4d14a35fcb")
DEVELOPMENT_SEED = 20611
CHECKPOINTS = {600: (200, 400, 600), 1800: (200, 400, 600, 1800)}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _configuration(policy, population, projected_seeds, budget_s):
    _require(type(population) is int and type(projected_seeds) is int and
             (policy, population, projected_seeds) in CONFIGURATIONS, "unregistered policy or population configuration")
    _require(evidence.finite(budget_s) and budget_s in CHECKPOINTS, "unregistered worker budget")
    return dict(policy=policy, population=population, projected_seeds=projected_seeds,
                budget_s=int(budget_s), checkpoints_s=list(CHECKPOINTS[budget_s]))


def _engine(config):
    # Give this invocation its own module globals, including the globals used
    # by run_arm's instrumented classes. Native SearchCase methods are shared
    # read-only; two direct API invocations cannot change one another's sizes.
    from verification import nf_anytime_search

    spec = nf_anytime_search.__spec__
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    engine.POPULATION = config["population"]
    engine.PROJECTED_SEEDS = config["projected_seeds"]
    engine.CHECKPOINTS = tuple(config["checkpoints_s"])

    class ConfiguredArm(engine._Arm):
        def fail(self, event, exc):
            # The native gate enforces POPULATION dynamically but its original
            # diagnostic says "64". Preserve the vectors and record the actual
            # registered size in both the failure event and the returned error.
            p = event.get("initial_population")
            if (isinstance(exc, ValueError) and event["stage"] == "initialize" and
                    event["restart"] == 0 and p is not None and
                    (p["count"] != engine.POPULATION or not p["full_space_passed"])):
                error = ValueError(f"cold startup requires {engine.POPULATION} members, all 16 axis spans > "
                                   f"1e-12 and affine rank 16; observed {p['count']} members and rank {p['rank']}")
                super().fail(event, error)
                raise error from exc
            super().fail(event, exc)

    engine._Arm = ConfiguredArm
    return engine


def run_configured(policy, population, projected_seeds, input_root, seed, out, started_ns, budget_s):
    config = _configuration(policy, population, projected_seeds, budget_s)
    _require(type(seed) is int and seed >= 0 and type(started_ns) is int and started_ns > 0,
             "a nonnegative integer seed and the parent monotonic launch timestamp are required")
    out = Path(out)
    for name in ("events.jsonl", "result.json", "anytime_progress.json"):
        if (out/name).exists():
            raise FileExistsError(out/name)
    engine = _engine(config)
    record = engine.run_arm(policy, Path(input_root), seed, out, started_ns, budget_s=budget_s)
    record["levels"] = {level: list(settings) for level, settings in record["levels"].items()}
    record.update(prospective_config=config, input_root=str(Path(input_root).resolve()),
                  launch_ns=started_ns, pid=os.getpid())
    return record


def _population(event, expected_size, checks):
    """Independently derive the recorded extrema, centroid, spans, spectrum and rank."""
    startup = event["restart"] == 0
    p = event.get("initial_population")
    if p is None:
        checks.expect(event["status"] != "completed", "initial_population_missing", event["restart"])
        return None
    x = np.asarray(p["vectors"], dtype=float)
    checks.expect(x.ndim == 2 and x.shape[1] == 16 and len(x) > 0 and
                  np.isfinite(x).all() and (np.abs(x) <= .1+evidence.TOL).all(),
                  "population_vector_bounds", event["restart"])
    center, span = x.mean(axis=0), np.ptp(x, axis=0)
    singular = np.linalg.svd((x-center)/.2, compute_uv=False)
    rank = int(np.count_nonzero(singular > evidence.TOL))
    full = bool(rank == 16 and (span > evidence.TOL).all())
    for field, value in (("axis_min", x.min(axis=0)), ("axis_max", x.max(axis=0)), ("axis_span", span),
                         ("center", center), ("scale", np.full(16, .2)), ("singular_values", singular)):
        actual = np.asarray(p[field], dtype=float)
        checks.expect(actual.shape == value.shape and np.allclose(actual, value, rtol=1e-12, atol=1e-15),
                      "population_numeric_evidence", field)
    checks.expect(p["rank"] == rank and p["count"] == len(x) and p["dimension"] == 16 and
                  p["rank_method"] == "centered_and_fixed_bound_width_scaled" and
                  p["rank_tolerance"] == p["axis_span_tolerance"] == evidence.TOL and
                  p["full_space_passed"] is full, "population_rank", event["restart"])
    archive, inherited = event["archive_count"], event["inherited_count"]
    checks.expect(all(type(n) is int and n >= 0 for n in (archive, inherited)) and
                  archive == p["archive_count"] and inherited == p["inherited_count"] and
                  p["startup"] is startup and p["full_space_required"] is startup and
                  event["full_space_required"] is startup and
                  p["archive_origin"] == ("current_run" if archive or inherited else "none"),
                  "population_origin", event["restart"])
    if event["status"] == "completed":
        checks.expect(len(x) == expected_size and (not startup or (full and archive == inherited == 0)),
                      "configured_cold_population", (len(x), expected_size, rank, archive, inherited))
    return dict(count=len(x), rank=rank, axis_span=span.tolist(), full=full, vectors=p["vectors"])


def _candidate_registry(record, checks):
    registry = {}
    for i, row in enumerate(record["candidates"]):
        key = evidence.candidate_id(row["candidate"], checks)
        checks.expect(row["i"] == i and row["id"] == key and key not in registry and
                      evidence.finite(row["first_seen_s"]) and 0 <= row["first_seen_s"] <= record["elapsed_s"],
                      "candidate_registry", i)
        registry[key] = row
    return registry


def _parity(parities, raw, checks):
    checks.expect(set(parities) == {"fine", "coarse"}, "parity_levels", tuple(parities))
    for level, p in parities.items():
        eligible = [e for e in raw if e["stage"] == level and e["measurement"]["eligible"]]
        first = next((e for e in eligible if e["status"] == "completed" and not e["cached"]), None)
        if first is not None:
            checks.expect(p is not None and p["candidate"] == first["candidate"], "parity_first_eligible", level)
        if p is not None:
            checks.expect(p.get("passed") is True and all(evidence.finite(p.get(k)) and p[k] >= 0 for k in
                          ("peak_relative_difference", "density_scaled_max_difference", "quadrature_absolute_difference",
                           "cpu_seconds", "gpu_seconds")) and p["peak_relative_difference"] <= 1e-10,
                          "native_cpu_gpu_parity", level)
            key = evidence.candidate_id(p["candidate"], checks)
            checks.expect(any(e["id"] == key and evidence.returned_deposition(e) for e in eligible),
                          "parity_candidate", level)


def _promotion(pending, item, latest_time, checks):
    checks.expect(len(pending) == 4 and isinstance(item, dict), "promotion_quartet", len(pending))
    winner = min(pending, key=lambda e: (e["measurement"]["peak"], e["start_s"]))
    checks.expect(item["ids"] == [e["id"] for e in pending] and item["selected"] == winner["id"] and
                  pending[-1]["end_s"] <= item["at_s"] <= latest_time, "promotion_selection", item)
    return winner["id"]


def _audit_events(record, raw, active, config, checks):
    """Replay fresh/cached scores, the DE target index and four-member coarse promotions."""
    budget = config["budget_s"]
    prelude = ["cold_inputs", "prepare_fine", "baseline", "build_fine"]
    if active == "coarse":
        prelude += ["prepare_coarse", "build_coarse"]
    stages = [e["stage"] for e in raw]
    checks.expect(stages[:len(prelude)] == prelude[:len(stages)] and
                  not set(stages[len(prelude):]).difference({"initialize", "proposal", "coarse", "fine", "progress_io"}),
                  "cold_setup_sequence", stages[:len(prelude)])
    registry = _candidate_registry(record, checks)
    previous, terminal = record["startup_s"], False
    checks.expect(evidence.finite(previous) and previous >= 0, "startup_clock", previous)
    startup = baseline = targets = None
    generations = target_index = total = restarts = promoted = 0
    best_seen = False
    memo = dict(fine=OrderedDict(), coarse=OrderedDict())
    pending, coarse_seen = [], set()
    for position, event in enumerate(raw):
        stage, status = event["stage"], event["status"]
        start, end = event["start_s"], event["end_s"]
        checks.expect(not terminal and evidence.finite(start) and evidence.finite(end) and
                      previous <= start < budget and start <= end, "event_clock", position)
        checks.expect(status in ("completed", "late", "failed") and (status != "completed" or end <= budget) and
                      (status != "late" or end >= budget), "event_status", position)
        previous, terminal = end, status != "completed"
        if "id" in event:
            row = registry.get(event["id"])
            checks.expect(row is not None and row["candidate"] == event.get("candidate") and row["first_seen_s"] <= end,
                          "event_candidate", position)
        if stage == "baseline":
            if "measurement" in event:
                evidence.check_measurement(event["measurement"], checks)
            if status == "completed":
                checks.expect("measurement" in event, "baseline_measurement", position)
                baseline = event
        if stage == "initialize":
            checks.expect(event["level"] == active and type(event["restart"]) is int and
                          event["restart"] == restarts, "initialization_sequence", position)
            pop = _population(event, config["population"], checks)
            if pop is not None and restarts == 0:
                startup = {k: v for k, v in pop.items() if k != "vectors"}
                startup["status"] = status
            if status == "completed":
                targets = pop["vectors"]
                generations = target_index = 0
            restarts += 1
        if stage == "proposal":
            checks.expect(event["level"] == active and targets is not None, "proposal_cold_gate", position)
            if "id" in event:
                underlying = "coordinate_polish" if best_seen and total % 19 == 0 else "differential_evolution"
                kind = "coordinate_polish" if underlying == "coordinate_polish" else (
                    "population_initialization" if generations == 0 else "de_trial")
                checks.expect(event.get("search_stage") == underlying and event.get("proposal_kind") == kind and
                              event.get("index") == target_index and event.get("generation") == generations,
                              "native_proposal_context", position)
                if kind == "population_initialization":
                    checks.expect(event["candidate"]["vector"] == targets[target_index],
                                  "initial_population_request", position)
            else:
                checks.expect(status != "completed", "proposal_identity", position)
        if stage not in ("coarse", "fine"):
            continue
        checks.expect(targets is not None and baseline is not None and baseline["measurement"]["eligible"] and
                      baseline["end_s"] <= start, "score_after_cold_setup", position)
        state = event["measurement"]
        evidence.check_measurement(state, checks)
        if stage == active:
            before = raw[position-1]
            checks.expect(before["stage"] == "proposal" and before["status"] == "completed" and
                          all(before.get(k) == event.get(k) for k in
                              ("id", "candidate", "search_stage", "proposal_kind", "generation", "index")),
                          "score_proposal_origin", position)
            checks.expect(len(pending) < 4, "promotion_before_new_proposal", position)
        else:
            checks.expect(active == "coarse" and stage == "fine" and event.get("search_stage") == "coarse_promotion"
                          and event.get("proposal_kind") == "external_validation" and
                          promoted < len(record["promotions"]), "external_promotion", position)
            winner = _promotion(pending, record["promotions"][promoted], start, checks)
            checks.expect(winner == event["id"], "promoted_fine_identity", position)
            pending, promoted = [], promoted+1
        if status != "completed":
            continue
        cache = memo[stage]
        key = event["id"]
        checks.expect(type(event.get("cached")) is bool and event["cached"] is (key in cache) and
                      type(event.get("physics_evaluated")) is bool and evidence.finite(event.get("objective")),
                      "fresh_score_evidence", position)
        if event["cached"]:
            old = cache[key]
            checks.expect(all(old[k] == event[k] for k in ("measurement", "objective", "physics_evaluated")),
                          "cached_outcome", position)
            cache.move_to_end(key)
        else:
            cache[key] = event
            if len(cache) > 2000:
                cache.popitem(last=False)
        if state["eligible"]:
            checks.expect(event["physics_evaluated"] and evidence.returned_deposition(event),
                          "eligible_returned_thermal_score", position)
            if stage == "coarse" and not event["cached"] and key not in coarse_seen:
                coarse_seen.add(key)
                pending.append(event)
        if stage == active:
            if event["search_stage"] == "differential_evolution":
                target_index = (target_index+1) % len(targets)
                generations += target_index == 0
            total += 1
            best_seen |= state["eligible"]
    checks.expect(evidence.finite(record["elapsed_s"]) and record["elapsed_s"] >= previous and
                  (record["stop"] != "deadline" or record["elapsed_s"] >= budget), "completion_clock", record["elapsed_s"])
    checks.expect(record["stop"] in ("deadline", "failed", "search_stopped") and
                  ((record["stop"] == "failed") == bool(record.get("error"))) and
                  (not any(e["status"] == "failed" for e in raw) or record["stop"] == "failed") and
                  (not any(e["status"] == "late" for e in raw) or record["stop"] == "deadline"),
                  "completion_status", record["stop"])
    checks.expect(record["baseline"] == (baseline["measurement"] if baseline else None) and
                  record["baseline_ready_s"] == (baseline["end_s"] if baseline else None), "baseline_record", budget)
    extra = record["promotions"][promoted:]
    if extra:
        checks.expect(active == "coarse" and len(extra) == 1 and len(pending) == 4 and
                      record["stop"] in ("failed", "deadline"), "terminal_promotion", len(extra))
        _promotion(pending, extra[0], record["elapsed_s"], checks)
    checks.expect(not record["rejected"] and record["completed_draws"] == 0, "adaptive_draw_counts", active)
    _parity(record["parity"], raw, checks)
    return startup, max(0, restarts-1)


def _endpoint(raw, snapshot, active, cutoff, baseline):
    requests = [e for e in raw if e["stage"] == active and e.get("proposal_kind") == "de_trial" and
                e["start_s"] < cutoff]
    completed = [e for e in requests if e["status"] == "completed" and e["end_s"] <= cutoff]
    fresh = [e for e in completed if not e["cached"]]
    thermal = [e for e in fresh if e["physics_evaluated"] and evidence.returned_deposition(e)]
    chosen = snapshot["incumbent"]
    base = baseline if snapshot["available"] else None
    return dict(t_s=cutoff, available=snapshot["available"], selected=chosen,
                peak=chosen["peak"] if chosen else None, A2=chosen["A2"] if chosen else None,
                baseline_peak=base["peak"] if base else None, baseline_A2=base["A2"] if base else None,
                peak_ratio=chosen["peak"]/base["peak"] if chosen else None,
                A2_ratio=chosen["A2"]/base["A2"] if chosen else None,
                first_5pct_s=snapshot["first_5pct_s"], counts=snapshot["counts"],
                de_proposals=sum(e["stage"] == "proposal" and e.get("proposal_kind") == "de_trial" and
                                 e["status"] == "completed" and e["end_s"] <= cutoff for e in raw),
                de_requests=len(requests), completed_de_requests=len(completed), fresh_de_requests=len(fresh),
                returned_thermal_de_trials=len(thermal),
                eligible_thermal_de_trials=sum(e["measurement"]["eligible"] for e in thermal))


def audit_arm(record, raw, *, policy, population, projected_seeds, budget_s, seed):
    """Check a finalized worker against its independently supplied registered arguments.

    Inconsistent evidence raises ValueError. An honestly recorded failure still
    has an auditable causal curve and returns normal=False with the reasons.
    """
    config = _configuration(policy, population, projected_seeds, budget_s)
    raw = list(raw)
    checks = evidence.Checks()
    checks.expect(record["schema"] == "nf-anytime-search.v1" and record["policy"] == policy and
                  type(record["seed"]) is int and record["seed"] == seed and record["budget_s"] == budget_s and
                  record["prospective_config"] == config and record["adaptive_config"] ==
                  dict(population=population, projected_seeds=projected_seeds,
                       initial_population_gate="startup_only_full_affine_dimension"), "registered_configuration", config)
    checks.expect(record["source"] == "OMP" and record["radius_mm"] == 20 and
                  record["levels"] == evidence.LEVELS, "native_physics_settings", policy)
    if "reference" in record or any(e["stage"] == "prepare_fine" for e in raw):
        checks.expect(record.get("reference") == asdict(ReferencePoint()), "fixed_operating_inputs", policy)
    checks.expect(raw == record["events"], "raw_event_equality", "events.jsonl differs from result")
    active = "coarse" if policy == "coarse_adaptive" else "fine"
    startup, restarts = _audit_events(record, raw, active, config, checks)
    endpoints = {}
    for t in config["checkpoints_s"]:
        snapshot = evidence.checkpoint(raw, t, budget=budget_s)
        checks.expect(record["checkpoints"].get(str(t)) == snapshot, "causal_checkpoint", t)
        endpoints[str(t)] = _endpoint(raw, snapshot, active, t, record["baseline"])
    final_key = str(config["budget_s"])
    checks.expect(set(record["checkpoints"]) == set(endpoints) and record["final"] == record["checkpoints"][final_key],
                  "causal_final", budget_s)
    reasons = []
    if record["stop"] != "deadline":
        reasons.append(record.get("error") or record["stop"])
    startup_passed = bool(startup and startup["full"] and startup["status"] == "completed")
    if not startup_passed:
        reasons.append("registered cold startup did not complete")
    if not endpoints[final_key]["available"]:
        reasons.append("eligible fine original wall unavailable by deadline")
    late = [dict(id=e["id"], stage=e["stage"], proposal_kind=e["proposal_kind"], status=e["status"],
                 start_s=e["start_s"], end_s=e["end_s"], measurement=e["measurement"],
                 returned_deposition=evidence.returned_deposition(e))
            for e in raw if e["stage"] in ("coarse", "fine") and e["end_s"] > budget_s]
    summary = dict(**config, seed=seed, normal=not reasons, reasons=reasons, stop=record["stop"],
                   startup_full_rank=startup_passed, startup_population=startup, restarts=restarts, de_level=active,
                   elapsed_s=record["elapsed_s"], checks=dict(checks.counts), checkpoints=endpoints,
                   parity=record["parity"], late_scores=late)
    summary.update(endpoints[final_key])
    for name in ("case_id", "input_root", "job_sha256", "registration_sha256", "pid", "launch_ns"):
        if name in record:
            summary[name] = record[name]
    return summary


def choose_fine(development_observations):
    """Apply the registered rule to both pilot configurations on all three inputs."""
    groups = {n: {} for policy, n, _ in CONFIGURATIONS if policy == "fine_adaptive"}
    for obs in development_observations:
        case, n = obs["case_id"], obs["population"]
        _require((obs["policy"], n, obs["projected_seeds"]) in CONFIGURATIONS[1:] and
                 case in DEVELOPMENT_CASES and obs["seed"] == DEVELOPMENT_SEED and obs["budget_s"] == 600,
                 "observation does not belong to the six registered development arms")
        _require(case not in groups[n], "duplicate development arm")
        _require(type(obs["normal"]) is bool and type(obs["startup_full_rank"]) is bool and
                 type(obs["eligible_thermal_de_trials"]) is int and obs["eligible_thermal_de_trials"] >= 0 and
                 type(obs["de_requests"]) is int and obs["de_requests"] >= obs["eligible_thermal_de_trials"],
                 "invalid development arm diagnostics")
        if obs.get("peak_ratio") is not None:
            _require(evidence.positive(obs["peak"]) and evidence.positive(obs["baseline_peak"]) and
                     obs["peak_ratio"] == obs["peak"]/obs["baseline_peak"], "development peak ratio differs from endpoint")
        groups[n][case] = dict(obs)
    configurations, eligible, incomplete = [], [], []
    for policy, n, projected in CONFIGURATIONS[1:]:
        by_case = groups[n]
        reasons = []
        missing = [case for case in DEVELOPMENT_CASES if case not in by_case]
        if missing:
            reason = f"population {n} lacks development observations for {', '.join(missing)}"
            reasons.append(reason)
            incomplete.append(reason)
        de_cases = sum(row["eligible_thermal_de_trials"] > 0 for row in by_case.values())
        if de_cases < 2:
            reasons.append(f"eligible fresh returned DE deposition in {de_cases} of 3 inputs")
        for case, row in by_case.items():
            if not row["normal"] or not row["startup_full_rank"] or row.get("peak_ratio") is None:
                reasons.append(f"{case} lacks a normal full-rank arm with a final eligible endpoint")
        ratios = [row["peak_ratio"] for row in by_case.values() if evidence.positive(row.get("peak_ratio"))]
        gm = math.exp(math.fsum(math.log(r) for r in ratios)/3) if len(ratios) == 3 else None
        settings = dict(policy=policy, population=n, projected_seeds=projected)
        configurations.append(dict(**settings, by_case=by_case, geometric_mean_peak_ratio=gm,
                                   de_cases=de_cases, qualified=not reasons, reasons=reasons))
        if not reasons:
            # Exact products avoid an arbitrary tolerance or rounding a near
            # tie while ranking the geometric means of three positive ratios.
            product = math.prod(Fraction(row["peak"])/Fraction(row["baseline_peak"]) for row in by_case.values())
            eligible.append(((product, -de_cases, n), settings))
    selected = min(eligible, key=lambda pair: pair[0])[1] if eligible and not incomplete else None
    reasons = incomplete or ([] if selected is not None else ["neither fine configuration met the development rule"])
    return dict(qualified=selected is not None, selected=selected, configurations=configurations, reasons=reasons)


def _bound_job(args):
    payload = (args.out/"job.json").read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    _require(os.environ.get("NF_PROSPECTIVE_JOB_SHA256") == digest, "parent job SHA256 mismatch")
    job = evidence.strict_json(payload)
    _require(Path(job["input_root"]).resolve() == args.input_root.resolve() and
             all(job[k] == getattr(args, k) for k in ("policy", "population", "projected_seeds", "budget_s", "seed")),
             "worker arguments differ from registered job")
    _require(isinstance(job.get("case_id"), str) and job["case_id"] and isinstance(job.get("phase"), str) and
             job["phase"] and isinstance(job.get("registration_sha256"), str) and
             len(job["registration_sha256"]) == 64 and all(c in "0123456789abcdef" for c in job["registration_sha256"]),
             "job lacks its case, phase or registration SHA256")
    return job, digest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", required=True, choices=("coarse_adaptive", "fine_adaptive"))
    for name in ("population", "projected-seeds", "seed", "started-ns"):
        parser.add_argument("--"+name, required=True, type=int)
    parser.add_argument("--budget-s", required=True, type=float)
    for name in ("input-root", "out"):
        parser.add_argument("--"+name, required=True, type=Path)
    args = parser.parse_args(argv)
    job, digest = _bound_job(args)
    params = dict(policy=args.policy, population=args.population, projected_seeds=args.projected_seeds,
                  seed=args.seed, budget_s=args.budget_s)
    record = run_configured(**params, input_root=args.input_root, out=args.out, started_ns=args.started_ns)
    record.update(case_id=job["case_id"], phase=job["phase"], job_sha256=digest,
                  registration_sha256=job["registration_sha256"])
    with (args.out/"result.json").open("x", encoding="utf-8") as stream:
        json.dump(record, stream, allow_nan=False, indent=2)
        stream.flush()
        os.fsync(stream.fileno())
    raw, _ = evidence.read_events(args.out/"events.jsonl")
    decision = audit_arm(record, raw, **params)
    if not decision["normal"]:
        print(" | ".join(decision["reasons"]), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
