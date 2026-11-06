"""Independent, read-only audit of the frozen 27-arm anytime experiment.

Run AFTER all timed workers terminate. This imports no search, evaluator, or
campaign module and performs no physical calculation. Its endpoint arithmetic
uses raw events, never the runner's snapshot/reconstruction functions. Only the
explicit --out artifact is written, exclusively; source evidence is untouched.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np


CODE_ROOT = Path(__file__).resolve().parents[1]
SEEDS = (20611, 20627, 20639)
POLICIES = ("coarse_adaptive", "fine_adaptive", "random_fine")
CHECKPOINTS = (200, 400, 600)
LEVELS = dict(coarse=["cubic", 320, .0005, 161, .001],
              fine=["cubic", 1280, .000125, 641, .00025])
SHOT = "input/D3D/d3d/d3d_shot_1a3155370c.parquet"
REVIEW = "baseline_review/baseline_review.json"
GEOMETRY = "input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json"
MANIFEST = "input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json"
INPUTS = (SHOT, REVIEW, GEOMETRY, MANIFEST)
TOL = 1e-12
SCOPE = ("same fixed transport model and frozen inputs; independent audit of recorded selection, "
         "timing and provenance, not experimental validation or a global-optimality certificate")


class AuditError(ValueError):
    def __init__(self, code, detail):
        self.code = code
        super().__init__(code + ": " + str(detail))


class Checks:
    def __init__(self):
        self.counts = Counter()

    def expect(self, condition, code, detail):
        self.counts[code] += 1
        if not condition:
            raise AuditError(code, detail)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def positive(value):
    return finite(value) and value > 0


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            result.update(block)
    return result.hexdigest()


def strict_json(payload):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise AuditError("json_duplicate_key", key)
            result[key] = value
        return result

    def invalid(value):
        raise AuditError("json_nonfinite", value)

    return json.loads(payload, object_pairs_hook=pairs, parse_constant=invalid)


def read(path):
    return strict_json(Path(path).read_text(encoding="utf-8"))


def read_events(path, allow_partial_tail=False):
    rows, partial = [], False
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        try:
            rows.append(strict_json(line))
        except json.JSONDecodeError:
            if not allow_partial_tail or index != len(lines)-1:
                raise
            partial = True
    return rows, partial


def verify_hash(path, expected, checks):
    checks.expect(digest(path) == expected, "file_hash", str(path))


def safe_child(root, relative, checks):
    root, relative = Path(root).resolve(), Path(relative)
    target = (root / relative).resolve()
    checks.expect(not relative.is_absolute() and ".." not in relative.parts and
                  target.is_relative_to(root) and target != root, "safe_path", str(relative))
    return target


def candidate_id(candidate, checks):
    checks.expect(isinstance(candidate, dict) and set(candidate) == {"source", "radius_mm", "vector"},
                  "candidate_schema", candidate)
    x = candidate["vector"]
    checks.expect(candidate["source"] == "OMP" and type(candidate["radius_mm"]) is int and
                  candidate["radius_mm"] == 20 and isinstance(x, list) and len(x) == 16 and
                  all(finite(v) and abs(v) <= .1+TOL for v in x), "candidate_bounds", candidate)
    return hashlib.sha256(json.dumps(candidate, sort_keys=True, allow_nan=False).encode()).hexdigest()


def check_measurement(value, checks):
    checks.expect(isinstance(value, dict) and type(value.get("eligible")) is bool and
                  isinstance(value.get("errors"), list) and isinstance(value.get("ledger"), dict) and
                  all(key in value for key in ("peak", "A2")), "measurement_schema", value)
    if value["eligible"]:
        checks.expect(positive(value["peak"]) and positive(value["A2"]), "eligible_observables", value)
        ledger = value["ledger"]
        checks.expect(all(finite(ledger.get(k)) for k in ("prescribed", "deposited", "unhit", "balance_error")) and
                      positive(ledger["prescribed"]) and ledger["deposited"] >= 0 and
                      0 <= ledger["unhit"] <= .001 and 0 <= ledger["balance_error"] <= 2e-12,
                      "eligible_power_ledger", ledger)
        error = abs(ledger["deposited"]+ledger["unhit"]-ledger["prescribed"])/ledger["prescribed"]
        checks.expect(error <= 2e-12 and abs(error-ledger["balance_error"]) <= 1e-15,
                      "power_balance", ledger)


def returned_deposition(event):
    state = event.get("measurement") or {}
    return (positive(state.get("peak")) and positive(state.get("A2")) and
            positive(state.get("ledger", {}).get("prescribed")))


def checkpoint(events, cutoff, budget=600.):
    """Independent raw-event endpoint, including all of the published counters."""
    cutoff, admission = float(cutoff), min(float(cutoff), float(budget))
    baseline = next((e for e in events if e["stage"] == "baseline" and
                     e.get("status") == "completed" and e["end_s"] <= admission and
                     e.get("measurement", {}).get("eligible") is True), None)
    chosen = None
    if baseline:
        value = baseline["measurement"]
        chosen = dict(kind="original_wall", id=None, candidate=None, peak=value["peak"], A2=value["A2"])
    counts = dict(fine=0, coarse=0, eligible_fine=0, eligible_coarse=0,
                  cached=0, proposals=0, rejected=0, late=0, failed=0)
    first = None
    for event in events:
        end = event.get("end_s")
        if end is None or end > cutoff:
            continue
        status, stage = event.get("status"), event["stage"]
        if status in ("late", "failed"):
            counts[status] += 1
        if status != "completed" or end > admission:
            continue
        counts["proposals"] += stage in ("proposal", "draw")
        counts["rejected"] += bool(event.get("rejection"))
        if stage not in ("coarse", "fine"):
            continue
        if event.get("cached"):
            counts["cached"] += 1
            continue
        counts[stage] += 1
        value = event["measurement"]
        if not value["eligible"]:
            continue
        counts["eligible_"+stage] += 1
        if stage == "fine" and chosen is not None:
            if first is None and value["peak"] <= .95*baseline["measurement"]["peak"]:
                first = end
            if value["peak"] < chosen["peak"]:
                chosen = dict(kind="new_design", id=event["id"], candidate=event["candidate"],
                              peak=value["peak"], A2=value["A2"])
    return dict(t_s=cutoff, available=baseline is not None, incumbent=chosen,
                first_5pct_s=first, counts=counts)


def vector_coverage(events):
    unique = {e["id"]: e["candidate"]["vector"] for e in events}
    if not unique:
        return dict(requests=len(events), unique_designs=0, rank=0, axis_span=None,
                    singular_values=[], full_affine_dimension_observed=False)
    x = np.asarray(list(unique.values()), float)
    singular = np.linalg.svd((x-x.mean(axis=0))/.2, compute_uv=False)
    span = np.ptp(x, axis=0)
    rank = int(np.count_nonzero(singular > TOL))
    return dict(requests=len(events), unique_designs=len(x), rank=rank,
                axis_span=span.tolist(), singular_values=singular.tolist(),
                full_affine_dimension_observed=bool(rank == 16 and np.all(span > TOL)))


def observed_coverage(events):
    output = {}
    for cutoff in CHECKPOINTS:
        output[str(cutoff)] = {}
        for level in ("coarse", "fine"):
            rows = [e for e in events if e["stage"] == level and e["status"] == "completed"
                    and e["end_s"] <= cutoff and not e.get("cached")]
            output[str(cutoff)][level] = dict(
                fresh_scores=vector_coverage(rows),
                returned_deposition=vector_coverage([e for e in rows if returned_deposition(e)]),
                eligible=vector_coverage([e for e in rows if e["measurement"]["eligible"]]),
                proposal_kind_counts=dict(Counter(e.get("proposal_kind", "missing") for e in rows)))
    return output


def population_evidence(event, checks):
    evidence = event.get("initial_population")
    if evidence is None:
        checks.expect(event["status"] != "completed", "population_missing", event)
        return None
    x = np.asarray(evidence["vectors"], float)
    checks.expect(x.ndim == 2 and x.shape[1] == 16 and len(x) > 0 and
                  np.isfinite(x).all() and np.all(abs(x) <= .1+TOL), "population_bounds", event["restart"])
    span, center = np.ptp(x, axis=0), x.mean(axis=0)
    singular = np.linalg.svd((x-center)/.2, compute_uv=False)
    rank = int(np.count_nonzero(singular > TOL))
    full = bool(rank == 16 and np.all(span > TOL))
    checks.expect(evidence["count"] == len(x) and evidence["dimension"] == 16 and
                  evidence["rank"] == rank and evidence["full_space_passed"] is full and
                  evidence["rank_method"] == "centered_and_fixed_bound_width_scaled" and
                  evidence["rank_tolerance"] == evidence["axis_span_tolerance"] == TOL,
                  "population_rank", event["restart"])
    for key, expected in (("axis_min", x.min(axis=0)), ("axis_max", x.max(axis=0)),
                          ("axis_span", span), ("center", center), ("scale", np.full(16, .2)),
                          ("singular_values", singular)):
        actual = np.asarray(evidence[key], float)
        checks.expect(actual.shape == expected.shape and np.allclose(actual, expected, rtol=1e-12, atol=1e-15),
                      "population_numeric_evidence", key)
    startup = event["restart"] == 0
    checks.expect(evidence["startup"] is startup and evidence["full_space_required"] is startup and
                  event["full_space_required"] is startup and
                  event["archive_count"] == evidence["archive_count"] and
                  event["inherited_count"] == evidence["inherited_count"], "population_origin", event["restart"])
    if startup:
        checks.expect(evidence["archive_count"] == evidence["inherited_count"] == 0 and
                      evidence["archive_origin"] == "none", "cold_start_inheritance", evidence)
    if event["status"] == "completed":
        checks.expect(len(x) == 64 and (not startup or full), "population_startup_gate", event["restart"])
    return x.tolist()


def audit_record(record, raw, *, checks=None, frozen=None):
    """Audit one gracefully finalized worker without trusting any saved endpoint."""
    checks = checks or Checks()
    expect = checks.expect
    policy = record["policy"]
    expect(policy in POLICIES and record["seed"] in SEEDS and record["budget_s"] == 600,
           "record_scope", (policy, record["seed"]))
    expect(record["source"] == "OMP" and record["radius_mm"] == 20 and record["levels"] == LEVELS,
           "record_settings", policy)
    config = record["adaptive_config"]
    expect(config == dict(population=64, projected_seeds=12,
                          initial_population_gate="startup_only_full_affine_dimension"), "adaptive_config", config)
    if frozen is not None:
        expect(record["reference"] == frozen["reference"], "reference_settings", record["reference"])
    expect(raw == record["events"], "raw_event_equality", "JSONL and result event arrays differ")
    active = {"coarse_adaptive": "coarse", "fine_adaptive": "fine", "random_fine": None}[policy]
    registry = {}
    for index, row in enumerate(record["candidates"]):
        identity = candidate_id(row["candidate"], checks)
        expect(row["id"] == identity and row["i"] == index and identity not in registry and
               finite(row["first_seen_s"]) and row["first_seen_s"] >= 0, "candidate_registry", index)
        registry[identity] = row
    previous = record["startup_s"]
    expect(finite(previous) and previous >= 0, "startup_clock", previous)
    if raw:
        expect(raw[0]["stage"] == "cold_inputs", "cold_input_evidence", "first event")
    expect(sum(e["stage"] == "cold_inputs" for e in raw) == int(bool(raw)),
           "cold_input_evidence", "exactly one initial check")
    terminal, baseline, initializations, gate, optimizer = False, None, [], None, None
    pending, coarse_seen, promoted = [], set(), 0
    memo, first_eligible, total, best_seen = {}, {}, 0, False
    for position, event in enumerate(raw):
        stage, status = event["stage"], event["status"]
        expect(stage in ("cold_inputs", "prepare_fine", "baseline", "build_fine", "prepare_coarse",
                         "build_coarse", "initialize", "proposal", "fine", "coarse", "draw", "progress_io"),
               "event_stage", stage)
        start, end = event["start_s"], event["end_s"]
        expect(not terminal and finite(start) and finite(end) and previous <= start < 600 and end >= start,
               "event_clock", position)
        expect(status in ("completed", "late", "failed") and
               (status != "completed" or end <= 600) and (status != "late" or end >= 600),
               "event_status", position)
        previous, terminal = end, status != "completed"
        if position:
            expect(raw[0]["status"] == "completed", "cold_input_evidence", position)
        if "id" in event:
            expect(event["id"] in registry and event.get("candidate") == registry[event["id"]]["candidate"] and
                   registry[event["id"]]["first_seen_s"] <= end, "event_candidate", position)
        if stage == "baseline":
            check_measurement(event["measurement"], checks)
            if status == "completed":
                expect(baseline is None, "baseline_unique", position)
                baseline = event
        if stage == "initialize":
            expect(active is not None and event["level"] == active, "initialization_level", position)
            expect(type(event["restart"]) is int and event["restart"] == len(initializations),
                   "initialization_sequence", event["restart"])
            initializations.append(event)
            vectors = population_evidence(event, checks)
            if status == "completed":
                if event["restart"] == 0:
                    gate = end
                optimizer = dict(vectors=vectors, energies=[None]*len(vectors), index=0, generation=0)
        if stage == "proposal":
            expect(active is not None and event["level"] == active and gate is not None and gate <= start,
                   "proposal_initialization", position)
        if stage not in ("fine", "coarse"):
            continue
        expect(baseline is not None and baseline["measurement"]["eligible"] and baseline["end_s"] <= start,
               "baseline_before_score", position)
        check_measurement(event["measurement"], checks)
        if stage == "coarse":
            expect(policy == "coarse_adaptive", "policy_coarse_stage", position)
            expect(len(pending) < 4, "promotion_before_new_coarse", position)
        if policy == "coarse_adaptive" and stage == "fine":
            expect(event.get("search_stage") == "coarse_promotion" and
                   event.get("proposal_kind") == "external_validation", "policy_fine_stage", position)
        elif active is None:
            expect(event.get("search_stage") == event.get("proposal_kind") == "projected_random",
                   "random_score_kind", position)
            before = raw[position-1]
            expect(before["stage"] == "draw" and before.get("id") == event["id"] and
                   before["status"] == "completed", "random_draw_origin", position)
        if active is not None:
            expect(gate is not None and gate <= start, "score_initialization_gate", position)
        if stage == active:
            before = raw[position-1]
            expect(before["stage"] == "proposal" and before["status"] == "completed" and
                   all(before.get(k) == event.get(k) for k in
                       ("id", "candidate", "search_stage", "proposal_kind", "generation", "index")),
                   "proposal_score_origin", position)
            search_stage = event["search_stage"]
            expect(search_stage in ("differential_evolution", "coordinate_polish"), "adaptive_score_kind", position)
            expected_stage = "coordinate_polish" if best_seen and total % 19 == 0 else "differential_evolution"
            expect(search_stage == expected_stage, "native_polish_schedule", position)
            i = optimizer["index"]
            kind = ("population_initialization" if optimizer["energies"][i] is None else "de_trial")
            if search_stage == "coordinate_polish":
                kind = "coordinate_polish"
            expect(event["proposal_kind"] == kind and event["generation"] == optimizer["generation"] and
                   event["index"] == i, "optimizer_context", position)
            if kind == "population_initialization":
                expect(event["candidate"]["vector"] == optimizer["vectors"][i], "initial_population_request", position)
        if status == "completed":
            expect(type(event.get("cached")) is bool and type(event.get("physics_evaluated")) is bool and
                   finite(event.get("objective")), "completed_score_fields", position)
        key = (stage, event["id"])
        if event.get("cached"):
            prior = memo.get(key)
            expect(prior is not None and event["measurement"] == prior["measurement"] and
                   event.get("objective") == prior.get("objective"), "cache_evidence", position)
        elif status == "completed":
            memo[key] = event
        if status == "completed" and event["measurement"]["eligible"]:
            first_eligible.setdefault(stage, event)
        if stage == "coarse" and status == "completed" and not event.get("cached") and event["measurement"]["eligible"]:
            if event["id"] not in coarse_seen:
                coarse_seen.add(event["id"])
                pending.append(event)
        if event.get("search_stage") == "coarse_promotion":
            expect(policy == "coarse_adaptive" and stage == "fine" and len(pending) == 4,
                   "promotion_quartet", position)
            winner = min(pending, key=lambda e: (e["measurement"]["peak"], e["start_s"]))
            expect(event["id"] == winner["id"] and promoted < len(record["promotions"]), "promotion_choice", position)
            item = record["promotions"][promoted]
            expect(item["ids"] == [e["id"] for e in pending] and item["selected"] == winner["id"] and
                   pending[-1]["end_s"] <= item["at_s"] <= start, "promotion_record", position)
            pending, promoted = [], promoted+1
        if stage == active and status == "completed":
            if event["search_stage"] == "differential_evolution":
                i, energy = optimizer["index"], event["objective"]
                if optimizer["energies"][i] is None or energy <= optimizer["energies"][i]:
                    optimizer["energies"][i] = energy
                    optimizer["vectors"][i] = event["candidate"]["vector"]
                optimizer["index"] = (i+1) % len(optimizer["vectors"])
                optimizer["generation"] += optimizer["index"] == 0
            total += 1
            best_seen |= event["measurement"]["eligible"]
    expect(finite(record["elapsed_s"]) and record["elapsed_s"] >= previous, "final_clock", record["elapsed_s"])
    expect(record["stop"] in ("deadline", "failed", "search_stopped"), "stop_reason", record["stop"])
    if record["stop"] == "deadline":
        expect(record["elapsed_s"] >= 600, "deadline_stop_clock", record["elapsed_s"])
    if baseline is None:
        expect(record["baseline"] is None and record["baseline_ready_s"] is None, "baseline_record", "unavailable")
    else:
        expect(record["baseline"] == baseline["measurement"] and record["baseline_ready_s"] == baseline["end_s"],
               "baseline_record", "value and completion")
    extra = record["promotions"][promoted:]
    if extra:
        expect(policy == "coarse_adaptive" and len(extra) == 1 and len(pending) == 4 and
               record["stop"] in ("deadline", "failed"), "terminal_promotion", extra)
        winner = min(pending, key=lambda e: (e["measurement"]["peak"], e["start_s"]))
        expect(extra[0]["ids"] == [e["id"] for e in pending] and extra[0]["selected"] == winner["id"] and
               pending[-1]["end_s"] <= extra[0]["at_s"] <= record["elapsed_s"], "terminal_promotion_choice", extra)
    expect(policy == "coarse_adaptive" or not record["promotions"], "policy_promotions", policy)
    draws = [e for e in raw if e["stage"] == "draw" and ("id" in e or "rejection" in e)]
    expect(record["completed_draws"] == len(draws) and record["rejected"] ==
           dict(Counter(e["rejection"] for e in draws if e.get("rejection"))), "draw_counts", policy)
    expect(set(record["parity"]) == {"fine", "coarse"}, "parity_levels", policy)
    for level in ("fine", "coarse"):
        parity = record["parity"].get(level)
        if level in first_eligible:
            expect(parity is not None and parity.get("passed") is True and
                   parity.get("candidate") == first_eligible[level]["candidate"], "parity_gate", level)
        if parity is not None:
            expect(parity.get("passed") is True and all(finite(parity.get(k)) and parity[k] >= 0 for k in
                   ("peak_relative_difference", "density_scaled_max_difference", "quadrature_absolute_difference",
                    "cpu_seconds", "gpu_seconds")) and parity["peak_relative_difference"] <= 1e-10,
                   "parity_values", level)
            identity = candidate_id(parity["candidate"], checks)
            expect(any(e["stage"] == level and e.get("id") == identity and
                       e.get("measurement", {}).get("eligible") for e in raw), "parity_candidate", level)
    endpoints = {str(t): checkpoint(raw, t) for t in CHECKPOINTS}
    expect(record["checkpoints"] == endpoints and record["final"] == endpoints["600"], "checkpoint_equality", policy)
    return dict(checkpoints=endpoints, coverage=observed_coverage(raw), startup_event_count=len(initializations[:1]),
                restart_event_count=max(0, len(initializations)-1),
                stop=record["stop"], censoring_time_s=min(600., record["elapsed_s"]),
                first_5pct_right_censored=endpoints["600"]["first_5pct_s"] is None,
                late_score_events=[dict(id=e["id"], stage=e["stage"], status=e["status"], end_s=e["end_s"],
                    measurement=e["measurement"], returned_deposition=returned_deposition(e))
                    for e in raw if e["stage"] in ("coarse", "fine") and e["end_s"] > 600],
                low_observed_rank_is_failure=False)


def omp_endpoints(review):
    """Only fixed-wall bookkeeping; no magnetic or thermal calculation."""
    wall = np.asarray(review["wall"], float)
    limits = {}
    for row in review["optimization_footprints"]:
        segment = row["segment"]
        a, b = wall[segment:segment+2]
        direction = b-a
        alpha = np.clip((np.asarray(row["points"])-a) @ direction / (direction @ direction), 0, 1)
        if np.ptp(alpha) > TOL:
            limits[segment] = (float(min(alpha)), float(max(alpha)))
    limits.update({i: (0., 1.) for i in review["optimization_policy"]["existing_design_segments"]})
    return {str(i): [(wall[i]+u*(wall[i+1]-wall[i])).tolist() for u in limits[i]]
            for i in (58, 59, 60, 68) if i in limits}


def audit_frozen(study, checks):
    frozen = read(study / "frozen.json")
    expect = checks.expect
    expect(frozen["seeds"] == list(SEEDS) and frozen["policies"] == list(POLICIES) and
           frozen["checkpoints_s"] == list(CHECKPOINTS) and frozen["budget_s"] == 600 and
           frozen["population"] == 64 and frozen["projected_seeds"] == 12 and
           frozen["source"] == "OMP" and frozen["radius_mm"] == 20 and frozen["levels"] == LEVELS,
           "frozen_protocol", "registered scope")
    expect(Path(frozen["code_root"]).resolve() == CODE_ROOT and
           frozen["python_version"] == platform.python_version(), "runtime_identity", frozen["code_root"])
    expect(frozen["reference"]["power_W"] == 1e6, "assigned_power", frozen["reference"])
    verify_hash(study / "protocol.ko.md", frozen["study_protocol_sha256"], checks)
    for name, expected in (frozen["source_hashes"] | frozen["original_reference_hashes"]).items():
        verify_hash(safe_child(CODE_ROOT, name, checks), expected, checks)
        verify_hash(safe_child(study / "code-snapshot", name, checks), expected, checks)
    source = Path(frozen["source_study"]).resolve()
    for name, expected in frozen["source_records"].items():
        verify_hash(safe_child(source, name, checks), expected, checks)
    preflight, downloads = read(source / "preflight.json"), read(source / "downloads.json")["files"]
    cases = frozen["cases"]
    expect(len(cases) == 3 and len({case["case_id"] for case in cases}) == 3,
           "case_count", "three distinct frozen equilibria")
    original = read(CODE_ROOT / REVIEW)
    for case in cases:
        cid, root = case["case_id"], Path(case["root"]).resolve()
        expect(cid and all(c.isalnum() or c in "_-" for c in cid) and root.drive.lower() == "d:",
               "case_path", case)
        expect(set(case["hashes"]) == set(INPUTS), "input_inventory", cid)
        for name, expected in case["hashes"].items():
            verify_hash(safe_child(root, name, checks), expected, checks)
        review, manifest = read(root / REVIEW), read(root / MANIFEST)
        expect(review["wall"] == original["wall"] and
               review["optimization_policy"] == original["optimization_policy"] and
               omp_endpoints(review) == omp_endpoints(original) and
               case["hashes"][GEOMETRY] == frozen["original_reference_hashes"][GEOMETRY],
               "unchanged_original_wall", cid)
        expect(manifest["shot_sha256"] == case["hashes"][SHOT] and
               manifest["time_ms"] == case["time_ms"] == 1400, "equilibrium_identity", cid)
        entries = [r for r in preflight["cases"] if r["case_id"] == cid]
        expect(len(entries) == 1, "preflight_case", cid)
        entry = entries[0]
        expected_root = Path(entry["root"])
        expected_root = expected_root.resolve() if expected_root.is_absolute() else (source / expected_root).resolve()
        expect((entry.get("status") == "passed" or entry.get("passed") is True) and expected_root == root and
               {k.replace("\\", "/"): v for k, v in entry["inputs"].items()} == case["hashes"], "preflight_inputs", cid)
        raw = entry["source"]
        expect(sum(item == raw for item in downloads) == 1 and Path(raw["local"]).stat().st_size == raw["bytes"] and
               raw["sha256"] == case["hashes"][SHOT], "download_provenance", cid)
        verify_hash(raw["local"], raw["sha256"], checks)
        forbidden = [root / "campaign_protected_candidates.json", root / "campaign_warm_candidates.json"]
        forbidden += [root / "results/omp-tangent-curvature-20261029" / f"R{mm:03d}" / "manifest.json"
                      for mm in (100, 50, 20)]
        expect(not any(path.exists() for path in forbidden), "cold_input_root", cid)
    return frozen


def audit_study(study):
    study, checks, rows, failures, hashes = Path(study).resolve(), Checks(), [], [], {}
    all_terminal, previous_exit = True, None
    report = dict(schema="nf-anytime-final-independent-audit.v1", scope=SCOPE,
                  created_utc=datetime.now(timezone.utc).isoformat(), study=str(study),
                  auditor_sha256=digest(__file__), endpoint_source="raw events.jsonl only",
                  finalized=False, passed=False, expected_arm_count=27)
    try:
        frozen = audit_frozen(study, checks)
        campaign_hash = digest(study / "frozen.json")
        hashes["frozen.json"] = campaign_hash
        report["campaign_sha256"] = campaign_hash
        baselines = {}
        registered_folders = set()
        for case_index, case in enumerate(frozen["cases"]):
            for seed_index, seed in enumerate(SEEDS):
                shift = (case_index+seed_index) % 3
                for policy in POLICIES[shift:]+POLICIES[:shift]:
                    label = dict(case_id=case["case_id"], seed=seed, policy=policy)
                    folder = study / "live" / case["case_id"] / str(seed) / policy
                    registered_folders.add(folder.relative_to(study).as_posix())
                    row = dict(label, status="not_terminal", certified_checkpoints=False)
                    rows.append(row)
                    if not (folder / "process.json").exists():
                        all_terminal = False
                        continue
                    terminal_observed = False
                    try:
                        process = read(folder / "process.json")
                        launch, end = process["launch_ns"], process["exit_observed_ns"]
                        checks.expect(type(launch) is int and type(end) is int and 0 <= launch <= end and
                                      finite(process["elapsed_s"]) and
                                      abs(process["elapsed_s"]-(end-launch)/1e9) <= 1e-9,
                                      "process_clock", label)
                        checks.expect(process["returncode"] is None or type(process["returncode"]) is int,
                                      "process_returncode", label)
                        terminal_observed = True
                        checks.expect(previous_exit is None or previous_exit <= launch, "serial_order", label)
                        previous_exit = end
                        env = process["environment"]
                        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE"):
                            checks.expect(env[name] == "1", "worker_environment", name)
                        for name in ("CUPY_CACHE_DIR", "CUDA_CACHE_PATH", "TEMP", "TMP", "MPLCONFIGDIR", "HF_HOME"):
                            checks.expect(Path(env[name]).resolve().is_relative_to(folder / "cache"), "private_cache_path", name)
                        for name in ("process.json", "events.jsonl", "result.json", "stdout.log"):
                            path = folder / name
                            if path.exists():
                                hashes[path.relative_to(study).as_posix()] = digest(path)
                        process_failed = process["returncode"] not in (0, 1)
                        raw, partial = (read_events(folder / "events.jsonl", allow_partial_tail=process_failed)
                                        if (folder / "events.jsonl").exists() else ([], False))
                        if process_failed:
                            row.update(status="process_failed", process=process, observed_event_count=len(raw),
                                       partial_event_tail=partial, checkpoints=None,
                                       censoring="process failure; no certified final worker record")
                            continue
                        record = read(folder / "result.json")
                        checks.expect(record["case_id"] == case["case_id"] and record["seed"] == seed and record["policy"] == policy and
                                      record["input_root"] == case["root"] and record["campaign_sha256"] == campaign_hash and
                                      record["launch_ns"] == launch, "worker_provenance", label)
                        checks.expect(record["elapsed_s"] <= process["elapsed_s"] and
                                      (process["returncode"] == 1) == (record["stop"] == "failed"), "worker_exit", label)
                        result = audit_record(record, raw, checks=checks, frozen=frozen)
                        status = "arm_failed" if record["stop"] == "failed" else (
                            "early_stopped" if record["stop"] == "search_stopped" else "completed")
                        row.update(status=status, certified_checkpoints=True, error=record.get("error"),
                                   process=process, **result)
                        if record["baseline"] is not None:
                            old = baselines.setdefault(case["case_id"], record["baseline"])
                            checks.expect(record["baseline"] == old, "baseline_across_all_seeds_methods", label)
                    except (AuditError, OSError, ValueError, KeyError, TypeError, IndexError, np.linalg.LinAlgError) as exc:
                        if not terminal_observed:
                            all_terminal = False
                        row.update(status="audit_failed", certified_checkpoints=False)
                        failures.append(dict(label, check=getattr(exc, "code", type(exc).__name__), detail=str(exc)))
        checks.expect(len(rows) == 27, "arm_count", len(rows))
        observed_folders = {path.parent.relative_to(study).as_posix()
                            for path in (study / "live").glob("*/*/*/process.json")}
        checks.expect(observed_folders.issubset(registered_folders), "unregistered_arm", sorted(observed_folders-registered_folders))
        for name, expected in hashes.items():
            verify_hash(study / name, expected, checks)
        audit_frozen(study, checks)  # Detect evidence changes during this read-only audit.
    except (AuditError, OSError, ValueError, KeyError, TypeError, IndexError, np.linalg.LinAlgError) as exc:
        all_terminal = False
        failures.append(dict(check=getattr(exc, "code", type(exc).__name__), detail=str(exc)))
    report.update(finalized=all_terminal and len(rows) == 27, passed=all_terminal and len(rows) == 27 and not failures,
                  all_arms_terminal=all_terminal and len(rows) == 27, arm_count=len(rows), rows=rows,
                  all_arms_completed_normally=len(rows) == 27 and all(row["status"] == "completed" for row in rows),
                  observed_status_counts=dict(Counter(row["status"] for row in rows)),
                  check_counts=dict(checks.counts), failures=failures, evidence_hashes=hashes,
                  diagnostics=dict(
                      audit_pass="Passed means the audit found no evidence inconsistency after all 27 processes terminated. It does not mean every experiment succeeded; inspect observed_status_counts.",
                      coverage="Unique vectors, centered and scaled by fixed 0.2 m bounds; rank tolerance 1e-12. Low observed rank is not an audit failure.",
                      score_counts="Fresh score requests include geometry rejections. Returned deposition is not a precise count of thermal-solver entry. Eligibility combines geometry and physics checks.",
                      censoring="Missing/process-failed arms have no fabricated checkpoints. Recorded worker failures retain independently reconstructable earlier outcomes.",
                      finite_panel="Exact expected uniform best is conditional on a fixed coarse-eligible outcome table, with k=floor(N/4) requests without replacement and baseline fallback. It assumes no independence of physical predictions and is not an equal-wall-time or global-space claim."))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    study, out = args.study.resolve(), args.out.resolve()
    if study.drive.lower() != "d:" or out.drive.lower() != "d:":
        raise ValueError("study and audit output must stay on D")
    if any(out.is_relative_to(study / name) for name in ("live", "code-snapshot", "reference", "uncertainty")):
        raise ValueError("write the independent report outside sealed evidence directories")
    report = audit_study(study)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({key: report[key] for key in ("finalized", "passed", "arm_count", "observed_status_counts")}))
    raise SystemExit(0 if report["passed"] else 2)


if __name__ == "__main__":
    main()
