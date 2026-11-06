"""Post-trial finite-panel allocation reference; no equal-wall-time claim.

Panel rules and preset seeds are fixed before live search. Membership is selected
by geometry alone and sealed before the additional panel thermal calculations.
The reused seeds/generator may reproduce live or earlier-study candidates; this
is not held-out validation. Every accepted member receives both registered scores,
including failures. The exact random comparator uses the same coarse-eligible
universe as the grouped policy and retains the original-wall fallback.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
from time import perf_counter


CODE_ROOT = Path(__file__).resolve().parents[1]
SEEDS = (20611, 20627, 20639)
POLICIES = ("coarse_adaptive", "fine_adaptive", "random_fine")
CHECKPOINTS = (200, 400, 600)
LEVELS = {"coarse": ("cubic", 320, .0005, 161, .001),
          "fine": ("cubic", 1280, .000125, 641, .00025)}
PANEL_SIZE = 24
MAX_ATTEMPTS = 6000
BATCH = 4
PANEL_SCOPE = ("deterministic rules and preset seeds fixed before live search; geometry-only membership sealed "
               "before additional panel thermal calculations; reused seeds/generator may reproduce live or "
               "earlier-study candidates; no historical thermal results enter the sampler; not held-out validation")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def object_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def finite_positive(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value > 0)


def validate_measurement(state):
    require(isinstance(state, dict) and isinstance(state.get("eligible"), bool),
            "measurement requires explicit boolean eligibility")
    require(isinstance(state.get("errors", []), list) and isinstance(state.get("ledger", {}), dict),
            "malformed measurement errors or ledger")
    if state["eligible"]:
        require(all(finite_positive(state.get(key)) for key in ("peak", "A2")),
                "eligible measurement must have positive finite peak and A2")
    return state


def measurement_record(value):
    """Keep eligibility, failures and observables together, also for rejected geometry."""
    require(isinstance(value, dict), "measurement is not a mapping")
    state = dict(eligible=value.get("eligible"), errors=list(value.get("errors", [])),
                 ledger=value.get("ledger", {}), peak=value.get("peak"), A2=value.get("A2"))
    validate_measurement(state)
    # Failed numerical output is evidence, but JSON must not encode NaN/Infinity.
    for key in ("peak", "A2"):
        if state[key] is not None and not finite_positive(state[key]):
            state["errors"].append("nonpositive_or_nonfinite_" + key)
            state[key] = None
    json.dumps(state, allow_nan=False)
    return state


def validate_candidate(identity, candidate):
    """Match Candidate.record/identity without importing the physical evaluator."""
    require(isinstance(candidate, dict) and set(candidate) == {"source", "radius_mm", "vector"},
            "malformed candidate record")
    require(candidate["source"] == "OMP" and candidate["radius_mm"] == 20,
            "candidate source or radius changed")
    vector = candidate["vector"]
    require(isinstance(vector, list) and len(vector) == 16 and
            all(isinstance(x, (int, float)) and not isinstance(x, bool)
                and math.isfinite(x) and abs(x) <= .1 + 1e-12 for x in vector),
            "candidate requires sixteen finite bounded OMP coefficients")
    require(identity == object_digest(candidate), "candidate identity differs from coefficients")
    return candidate


def expected_uniform_best(values, k, baseline):
    """Exact order-statistic expectation for k draws without replacement.

    None denotes an ineligible fine request. Its outcome is the baseline, as is
    any eligible outcome worse than the baseline. For sorted outcome index i,
    C(n-i-1,k-1)/C(n,k) is the probability it is the first selected index. Ties
    are separate population members, so no special tie approximation is needed.
    """
    require(finite_positive(baseline), "positive finite baseline required")
    require(isinstance(k, int) and not isinstance(k, bool) and 0 <= k <= len(values),
            "integer sample size outside the finite population")
    require(all(v is None or finite_positive(v) for v in values), "invalid finite-population outcome")
    if k == 0:
        return float(baseline)
    outcomes = sorted(baseline if v is None else min(v, baseline) for v in values)
    denominator = math.comb(len(outcomes), k)
    return math.fsum(value * (math.comb(len(outcomes)-i-1, k-1) / denominator)
                     for i, value in enumerate(outcomes[:len(outcomes)-k+1]))


def summarize_panel(rows, baseline):
    """Apply the preregistered eligible-groups-of-four rule to a scored panel."""
    validate_measurement(baseline)
    require(baseline["eligible"], "eligible original-wall fallback required")
    require([r["i"] for r in rows] == list(range(len(rows))), "panel order must be unique and contiguous")
    require(len({r["id"] for r in rows}) == len(rows), "duplicate panel identity")
    for row in rows:
        validate_measurement(row["coarse"])
        validate_measurement(row["fine"])
    universe = [row for row in rows if row["coarse"]["eligible"]]
    complete = len(universe) // BATCH * BATCH
    groups = [universe[i:i+BATCH] for i in range(0, complete, BATCH)]
    promoted = [min(group, key=lambda row: (row["coarse"]["peak"], row["i"])) for group in groups]

    def best(pool):
        chosen = dict(kind="original_wall", id=None, i=None, peak=baseline["peak"], A2=baseline["A2"])
        for row in pool:
            if row["fine"]["eligible"] and row["fine"]["peak"] < chosen["peak"]:
                chosen = dict(kind="new_design", id=row["id"], i=row["i"],
                              peak=row["fine"]["peak"], A2=row["fine"]["A2"])
        return chosen

    selected, reference = best(promoted), best(rows)
    k = len(promoted)  # Failed/ineligible fine requests still consume allocation.
    values = [row["fine"]["peak"] if row["fine"]["eligible"] else None for row in universe]
    expected = expected_uniform_best(values, k, baseline["peak"])
    return dict(coarse_eligible_ids=[r["id"] for r in universe],
                promotion_groups=[dict(ids=[r["id"] for r in group], selected=winner["id"])
                                  for group, winner in zip(groups, promoted)],
                promoted_ids=[r["id"] for r in promoted], fine_requests_k=k,
                unused_eligible_tail_ids=[r["id"] for r in universe[complete:]],
                selected=selected, reference_best=reference,
                regret_absolute=selected["peak"]-reference["peak"],
                regret_relative=selected["peak"]/reference["peak"]-1,
                uniform_expected_peak=expected,
                uniform_expected_regret_absolute=expected-reference["peak"],
                uniform_expected_regret_relative=expected/reference["peak"]-1,
                uniform_population_size=len(universe),
                uniform_subset_count=math.comb(len(universe), k),
                uniform_expectation="exact combinatorial order statistics; no Monte Carlo",
                scope="allocation after a common coarse screen on this locked finite panel only",
                equal_wall_time_comparison=False, global_optimum_claim=False)


class EvidenceStore:
    """Exclusive artifacts with byte-level seals; an interrupted unsealed write fails closed."""
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)

    def path(self, name):
        path = Path(name)
        require(not path.is_absolute() and ".." not in path.parts and path.parts,
                "unsafe evidence path")
        target = (self.root / path).resolve()
        require(target.is_relative_to(self.root) and target != self.root, "unsafe evidence path")
        return target

    @staticmethod
    def seal(path):
        return path.with_name(path.name + ".sha256")

    def exists(self, name):
        path = self.path(name)
        if self.seal(path).exists() and not path.exists():
            raise ValueError("sealed artifact is missing: " + str(name))
        return path.exists()

    def write_bytes(self, name, payload):
        path = self.path(name)
        require(not path.name.endswith(".sha256"), "hash sidecar cannot be an artifact")
        path.parent.mkdir(parents=True, exist_ok=True)
        if self.seal(path).exists():
            raise FileExistsError(path)
        with path.open("xb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        with self.seal(path).open("x", encoding="ascii") as stream:
            stream.write(hashlib.sha256(payload).hexdigest() + "\n")
            stream.flush()
            os.fsync(stream.fileno())

    def write(self, name, value):
        payload = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf8")
        self.write_bytes(name, payload)

    def read_bytes(self, name):
        path = self.path(name)
        payload = path.read_bytes()
        expected = self.seal(path).read_text(encoding="ascii").strip()
        require(hashlib.sha256(payload).hexdigest() == expected, "artifact byte hash differs: " + str(name))
        return payload

    def read(self, name):
        return json.loads(self.read_bytes(name))

    def sha256(self, name):
        return hashlib.sha256(self.read_bytes(name)).hexdigest()

    def artifacts(self, pattern="*"):
        return sorted(p.relative_to(self.root).as_posix() for p in self.root.rglob(pattern)
                      if p.is_file() and not p.name.endswith((".sha256", ".lock")))

    def verify_all(self):
        for name in self.artifacts():
            self.read_bytes(name)
        for seal in self.root.rglob("*.sha256"):
            require(seal.with_name(seal.name[:-7]).is_file(), "orphan artifact hash: " + str(seal))


def validate_frozen(study):
    """Verify frozen live code, its snapshot, inputs and preregistration before reuse."""
    study = Path(study).resolve()
    frozen = read_json(study / "frozen.json")
    require(frozen["seeds"] == list(SEEDS) and frozen["policies"] == list(POLICIES)
            and frozen["checkpoints_s"] == list(CHECKPOINTS) and frozen["budget_s"] == 600,
            "live protocol seeds, arms or checkpoints changed")
    require(frozen["source"] == "OMP" and frozen["radius_mm"] == 20, "live physical scope changed")
    require(frozen["levels"] == {key: list(value) for key, value in LEVELS.items()}, "native levels changed")
    require(Path(frozen["code_root"]).resolve() == CODE_ROOT, "run from the frozen code root")
    require(frozen["python_version"] == platform.python_version(), "Python version changed")
    require(digest(study / "protocol.ko.md") == frozen["study_protocol_sha256"], "study protocol changed")
    for name, expected in (frozen["source_hashes"] | frozen["original_reference_hashes"]).items():
        path = (CODE_ROOT / name).resolve()
        require(path.is_relative_to(CODE_ROOT), "unsafe frozen source path")
        require(digest(path) == expected, "frozen source/input changed: " + name)
        require(digest(study / "code-snapshot" / name) == expected, "frozen snapshot changed: " + name)
    for name, expected in frozen["source_records"].items():
        require(digest(Path(frozen["source_study"]) / name) == expected, "source evidence changed: " + name)
    require(len(frozen["cases"]) == 3 and len({c["case_id"] for c in frozen["cases"]}) == 3,
            "three distinct frozen equilibria required")
    for case in frozen["cases"]:
        require(case["case_id"] and all(ch.isalnum() or ch in "_-" for ch in case["case_id"]), "unsafe case ID")
        require(case["time_ms"] == 1400, "equilibrium time changed")
        for name, expected in case["hashes"].items():
            require(digest(Path(case["root"]) / name) == expected, "staged input changed: " + case["case_id"])
    return frozen


def freeze_post(study, kind, details, sources):
    frozen = validate_frozen(study)
    store = EvidenceStore(Path(study) / kind)
    store.verify_all()
    source_hashes = {Path(p).resolve().relative_to(CODE_ROOT).as_posix(): digest(p) for p in sources}
    protocol = dict(schema="nf-anytime-" + kind + "-protocol.v1",
                    campaign_sha256=digest(Path(study) / "frozen.json"),
                    study_protocol_sha256=frozen["study_protocol_sha256"],
                    source_hashes=source_hashes, details=details,
                    reference=frozen["reference"], backend="gpu", source="OMP", radius_mm=20,
                    scope="post-trial calculation; additional cost outside the live competition")
    if store.exists("protocol.json"):
        old = store.read("protocol.json")
        require({k: v for k, v in old.items() if k != "created_utc"} == protocol,
                "post-trial protocol or code changed; mixed continuation refused")
    else:
        protocol["created_utc"] = datetime.now(timezone.utc).isoformat()
        store.write("protocol.json", protocol)
    for name, expected in source_hashes.items():
        snapshot = "code-snapshot/" + name
        if not store.exists(snapshot):
            store.write_bytes(snapshot, (CODE_ROOT / name).read_bytes())
        require(store.sha256(snapshot) == expected, "post-trial code snapshot differs")
    return frozen, store


def load_audited_arms(study, frozen):
    """Read only finished campaigns and recheck the independent audit's evidence.

    A supervisor-killed arm has no certified checkpoints. Its complete recorded
    events may still identify a last overrun design, but never enter the native
    comparison pool without a validated result.
    """
    from verification.nf_anytime_campaign import reconstruct, validate_decisions, validate_record
    study = Path(study)
    audit_path = study / "audit.json"
    audit = read_json(audit_path)
    require(audit["campaign_sha256"] == digest(study / "frozen.json"), "audit belongs to another campaign")
    expected = {(case["case_id"], seed, policy) for case in frozen["cases"] for seed in SEEDS for policy in POLICIES}
    keys = [(r["case_id"], r["seed"], r["policy"]) for r in audit["rows"]]
    require(len(keys) == len(expected) and set(keys) == expected, "audit does not contain exactly the frozen arms")
    hashes, arms = {"frozen.json": digest(study / "frozen.json"), "audit.json": digest(audit_path)}, []
    for row in audit["rows"]:
        require(row["status"] in ("completed", "arm_failed", "process_failed"),
                "live campaign has not finished or fails audit")
        prefix = Path("live") / row["case_id"] / str(row["seed"]) / row["policy"]
        process = read_json(study / prefix / "process.json")
        require(process == row["process"], "audited process evidence changed")
        hashes[(prefix / "process.json").as_posix()] = digest(study / prefix / "process.json")
        event_path = study / prefix / "events.jsonl"
        raw, partial_tail = [], False
        if event_path.exists():
            hashes[(prefix / "events.jsonl").as_posix()] = digest(event_path)
            lines = event_path.read_text(encoding="utf8").splitlines()
            for index, line in enumerate(lines):
                try:
                    raw.append(json.loads(line))
                except json.JSONDecodeError:
                    require(row["status"] == "process_failed" and index == len(lines)-1,
                            "malformed audited event stream")
                    partial_tail = True
        if row["status"] == "process_failed":
            arms.append(dict(case_id=row["case_id"], seed=row["seed"], policy=row["policy"],
                             audit_status=row["status"], budget_s=600, baseline=None,
                             candidates=[], candidate_registry_complete=False,
                             events=raw, native_pool_eligible=False, partial_event_tail=partial_tail,
                             checkpoints={str(t): dict(available=False) for t in CHECKPOINTS}))
            continue
        result_path = study / prefix / "result.json"
        require(digest(result_path) == row["result_sha256"], "audited result bytes changed")
        record = read_json(result_path)
        validate_record(record, raw)
        validate_decisions(record)
        checkpoints = {str(t): reconstruct(record, t) for t in CHECKPOINTS}
        require(checkpoints == row["checkpoints"], "independently audited checkpoints differ")
        require(record["case_id"] == row["case_id"] and record["seed"] == row["seed"]
                and record["policy"] == row["policy"], "audited arm identity differs")
        hashes[(prefix / "result.json").as_posix()] = digest(result_path)
        arms.append(dict(case_id=row["case_id"], seed=row["seed"], policy=row["policy"],
                         audit_status=row["status"], budget_s=record["budget_s"], baseline=record["baseline"],
                         candidates=record["candidates"], candidate_registry_complete=True,
                         events=raw, checkpoints=checkpoints, native_pool_eligible=True))
    return arms, hashes


def bind_evidence(store, name, value):
    if store.exists(name):
        require(store.read(name) == value, "reused evidence changed: " + name)
    else:
        store.write(name, value)


def verify_external_hashes(study, hashes):
    for name, expected in hashes.items():
        require(digest(Path(study) / name) == expected, "live evidence changed during post-trial calculation: " + name)


def panel_live_overlap(panel, arms):
    """Describe identity overlap after generation; never use outcomes to select members.

    The measured scope is this campaign's live arms at the same equilibrium,
    across all seeds and methods. Earlier study archives are not scanned.
    """
    case_arms = [arm for arm in arms if arm["case_id"] == panel["case_id"]]
    require(bool(case_arms), "live overlap requires audited arms for the panel equilibrium")
    panel_ids = [row["id"] for row in panel["rows"]]
    require(len(set(panel_ids)) == len(panel_ids), "duplicate panel identity in overlap provenance")
    for row in panel["rows"]:
        validate_candidate(row["id"], row["candidate"])
    proposed, fine_observed, by_arm = set(), set(), []
    for arm in case_arms:
        live_ids, fine_ids = set(), set()
        for row in arm.get("candidates", []):
            validate_candidate(row["id"], row["candidate"])
            live_ids.add(row["id"])
        for event in arm["events"]:
            if "id" not in event:
                continue
            validate_candidate(event["id"], event["candidate"])
            live_ids.add(event["id"])
            if (event.get("stage") != "fine" or event.get("status") != "completed"
                    or event.get("cached") or not isinstance(event.get("measurement"), dict)):
                continue
            end = event.get("end_s")
            require(isinstance(end, (int, float)) and math.isfinite(end) and end >= 0,
                    "invalid completion clock in fine-observation overlap")
            if end <= arm["budget_s"]:
                fine_ids.add(event["id"])
        proposed.update(live_ids)
        fine_observed.update(fine_ids)
        registry_complete = bool(arm.get("candidate_registry_complete", "candidates" in arm))
        events_complete = arm.get("audit_status") != "process_failed" and not arm.get("partial_event_tail", False)
        proposed_overlap = [identity for identity in panel_ids if identity in live_ids]
        fine_overlap = [identity for identity in panel_ids if identity in fine_ids]
        by_arm.append(dict(seed=arm["seed"], policy=arm["policy"], status=arm.get("audit_status"),
                           live_proposed_count=len(live_ids), live_on_time_fine_observed_count=len(fine_ids),
                           proposed_intersection_ids=proposed_overlap, proposed_intersection_count=len(proposed_overlap),
                           on_time_fine_intersection_ids=fine_overlap, on_time_fine_intersection_count=len(fine_overlap),
                           complete_live_proposal_coverage=registry_complete,
                           complete_live_fine_observation_coverage=events_complete))
    proposed_overlap = [identity for identity in panel_ids if identity in proposed]
    fine_overlap = [identity for identity in panel_ids if identity in fine_observed]
    return dict(case_id=panel["case_id"], panel_seed=panel["seed"], panel_size=len(panel_ids),
                scope="current campaign's recorded live candidates in the same equilibrium across all seeds/methods",
                live_proposed_ids=sorted(proposed), live_proposed_count=len(proposed),
                live_on_time_fine_observed_ids=sorted(fine_observed), live_on_time_fine_observed_count=len(fine_observed),
                proposed_intersection_ids=proposed_overlap, proposed_intersection_count=len(proposed_overlap),
                on_time_fine_intersection_ids=fine_overlap, on_time_fine_intersection_count=len(fine_overlap),
                by_arm=by_arm,
                complete_live_proposal_coverage=all(row["complete_live_proposal_coverage"] for row in by_arm),
                complete_live_fine_observation_coverage=all(row["complete_live_fine_observation_coverage"] for row in by_arm),
                fine_observation_definition="completed noncached fine measurement by 600 seconds, including ineligible outcomes and completed geometry decisions; late outcomes excluded",
                incomplete_coverage_note="failed processes expose only recorded identities/events, not unlogged proposals or unfinished measurements",
                historical_archive_overlap_checked=False, held_out_validation_claim=False,
                live_outcomes_used_for_generation=False, membership_scope=PANEL_SCOPE)


def generate_panel(geometry, seed, *, draw_fn=None, rng=None):
    import numpy as np
    if draw_fn is None:
        from verification.nf_base_live import draw as draw_fn
    rng = np.random.default_rng(seed) if rng is None else rng
    accepted, attempts, rejected = [], [], Counter()
    start = perf_counter()
    for attempt in range(1, MAX_ATTEMPTS + 1):
        candidate, reason = draw_fn(geometry, rng, accepted)
        if reason:
            rejected[reason] += 1
            attempts.append(dict(attempt=attempt, rejection=reason))
            continue
        require(candidate is not None, "generator returned no candidate and no rejection")
        record = candidate.record()
        validate_candidate(candidate.identity, record)
        coefficients = list(record["vector"]) + [0.] * 8
        accepted.append(dict(i=len(accepted), id=candidate.identity, candidate=record,
                             coefficients=coefficients, coefficients_sha256=object_digest(coefficients),
                             cumulative_draws=attempt))
        attempts.append(dict(attempt=attempt, id=candidate.identity))
        if len(accepted) == PANEL_SIZE:
            break
    return dict(seed=seed, rows=accepted, attempt_records=attempts, attempts=len(attempts),
                rejected=dict(rejected), complete=len(accepted) == PANEL_SIZE,
                generation_s=perf_counter()-start, thermal_results_used_for_selection=False)


def validate_panel(panel, seed):
    require(panel["seed"] == seed and panel["complete"] == (len(panel["rows"]) == PANEL_SIZE), "panel seed or completeness differs")
    require(0 <= len(panel["rows"]) <= PANEL_SIZE and 1 <= panel["attempts"] <= MAX_ATTEMPTS, "invalid panel size or proposal count")
    require(panel["complete"] or panel["attempts"] == MAX_ATTEMPTS, "incomplete panel stopped before its proposal limit")
    require(panel["thermal_results_used_for_selection"] is False, "panel selected using thermal scores")
    require(len(panel["attempt_records"]) == panel["attempts"], "missing proposal evidence")
    require([r["attempt"] for r in panel["attempt_records"]] == list(range(1, panel["attempts"]+1)),
            "panel attempt order differs")
    require(dict(Counter(r["rejection"] for r in panel["attempt_records"] if "rejection" in r)) == panel["rejected"],
            "panel rejection counts differ from attempt evidence")
    require([r["i"] for r in panel["rows"]] == list(range(len(panel["rows"]))), "panel order differs")
    require([r["id"] for r in panel["attempt_records"] if "id" in r] == [r["id"] for r in panel["rows"]],
            "accepted attempt order differs from panel")
    identities = set()
    for row in panel["rows"]:
        validate_candidate(row["id"], row["candidate"])
        require(row["id"] not in identities, "duplicate panel identity")
        identities.add(row["id"])
        require(row["coefficients"] == row["candidate"]["vector"] + [0.] * 8 and
                row["coefficients_sha256"] == object_digest(row["coefficients"]), "panel coefficients changed")
        require(panel["attempt_records"][row["cumulative_draws"]-1] ==
                dict(attempt=row["cumulative_draws"], id=row["id"]), "candidate proposal provenance differs")


def evaluate_score(evaluator, candidate, score_fn=None):
    if score_fn is None:
        from verification.nf_base_budget import score as score_fn
    start = perf_counter()
    try:
        value = score_fn(evaluator, candidate)
        measured = measurement_record(value)
        status, error = "completed", None
    except Exception as exc:
        error = type(exc).__name__ + ": " + str(exc)
        measured = dict(eligible=False, errors=[error], ledger={}, peak=None, A2=None)
        status = "failed"
    return dict(measurement=measured, status=status, error=error, seconds=perf_counter()-start)


def score_panel(store, prefix, panel, evaluators, baseline, *, score_fn=None):
    """Resume individual requests without replacing failures or selecting new members."""
    from d3d_optimization.campaign_geometry import Candidate
    panel_name = prefix + "/panel.json"
    require(store.read(panel_name) == panel, "candidate panel must be sealed before scoring")
    panel_hash, protocol_hash = store.sha256(panel_name), store.sha256("protocol.json")
    rows = [dict(item) for item in panel["rows"]]
    evidence = {}
    for level in ("coarse", "fine"):
        for row in rows:
            name = prefix + "/scores/" + row["id"] + "." + level + ".json"
            identity = dict(id=row["id"], candidate=row["candidate"], level=level,
                            panel_sha256=panel_hash, protocol_sha256=protocol_hash)
            if store.exists(name):
                scored = store.read(name)
                require(all(scored[key] == val for key, val in identity.items()), "reused score provenance differs")
                validate_measurement(scored["measurement"])
            else:
                scored = dict(identity, **evaluate_score(evaluators[level], Candidate.from_record(row["candidate"]), score_fn))
                store.write(name, scored)
            row[level] = scored["measurement"]
            evidence[name] = store.sha256(name)
    return dict(schema="nf-anytime-reference-panel.v1", panel_sha256=panel_hash,
                protocol_sha256=protocol_hash, complete=panel["complete"], baseline=baseline,
                rows=rows, score_hashes=evidence, reference_fine_evaluations=len(rows),
                assessment=summarize_panel(rows, baseline))


def reference_details():
    return dict(seeds=list(SEEDS), panel_size=PANEL_SIZE, max_geometry_attempts=MAX_ATTEMPTS,
                levels={key: list(value) for key, value in LEVELS.items()}, batch=BATCH,
                generation="nf_base_live.draw on fine geometry; membership sealed before additional panel heat calculations",
                membership_scope=PANEL_SCOPE, held_out_validation_claim=False,
                overlap_provenance="record live proposed and on-time fine-observed identity intersections after generation; all seeds/methods in the same equilibrium; never use overlap or heat values for selection",
                promotion="complete groups of four fresh coarse-eligible rows in original order; minimum (peak,index); no tail flush",
                random_comparator="exact uniform without-replacement k from all coarse-eligible members; ineligible fine outcome retains baseline",
                reference="original wall plus every fine-eligible locked panel member", tuning=False)


def configure_environment(study, phase):
    """Keep runtime artifacts on the study drive, outside sealed result evidence."""
    root = Path(study).resolve() / "posttrial-runtime" / phase
    locations = {key: root / leaf for key, leaf in
                 (("CUPY_CACHE_DIR", "cupy"), ("CUDA_CACHE_PATH", "cuda"),
                  ("TEMP", "tmp"), ("TMP", "tmp"), ("MPLCONFIGDIR", "matplotlib"),
                  ("HF_HOME", "huggingface"))}
    for path in locations.values():
        path.mkdir(parents=True, exist_ok=True)
    values = {key: str(path) for key, path in locations.items()}
    values.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  PYTHONDONTWRITEBYTECODE="1")
    os.environ.update(values)
    sys.dont_write_bytecode = True
    return values


def start_invocation(store):
    number = len(store.artifacts("*.started.json"))
    name = "runs/" + str(number).zfill(4)
    variables = ("CUPY_CACHE_DIR", "CUDA_CACHE_PATH", "TEMP", "TMP", "MPLCONFIGDIR", "HF_HOME",
                 "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE")
    store.write(name + ".started.json", dict(created_utc=datetime.now(timezone.utc).isoformat(), pid=os.getpid(),
                                             environment={key: os.environ.get(key) for key in variables}))
    return name


def cost_summary(store):
    preparations = [store.read(p) for p in store.artifacts("*.json") if p.startswith("preparation/")]
    scores = [store.read(p) for p in store.artifacts("*.json") if "/scores/" in p]
    panels = [store.read(p) for p in store.artifacts("panel.json")]
    starts = store.artifacts("*.started.json")
    finishes = store.artifacts("*.finished.json")
    return dict(preparation_seconds=sum(p["seconds"] for p in preparations),
                candidate_scoring_seconds=sum(p["seconds"] for p in scores),
                panel_generation_seconds=sum(p["generation_s"] for p in panels),
                completed_invocation_seconds=sum(store.read(p)["seconds"] for p in finishes),
                invocation_count=len(starts), unfinished_invocation_count=len(starts)-len(finishes),
                note="extra post-trial costs; all recorded preparations (including resume), proposals and scores are retained; unfinished process wall time is unknown")


def run_reference(study):
    configure_environment(study, "reference")
    from dataclasses import asdict
    from d3d_optimization.campaign_state import keep_awake, process_lock
    from d3d_optimization.reference_transport import ReferencePoint
    from verification import nf_base_budget as physics
    from verification.nf_anytime_campaign import versions
    start = perf_counter()
    with process_lock(Path(study) / "campaign.lock"), process_lock(Path(study) / "posttrial.lock"), keep_awake():
        frozen, store = freeze_post(study, "reference", reference_details(), [Path(__file__)])
        require(versions() == frozen["packages"], "runtime package versions changed")
        arms, hashes = load_audited_arms(study, frozen)
        bind_evidence(store, "live-evidence.json", hashes)
        if store.exists("summary.json"):
            return store.read("summary.json")
        invocation = start_invocation(store)
        ref = ReferencePoint(**frozen["reference"])
        require(asdict(ref) == frozen["reference"], "reference settings changed")
        try:
            for case in frozen["cases"]:
                prefixes = [case["case_id"] + "/" + str(seed) for seed in SEEDS]
                if all(store.exists(prefix + "/result.json") for prefix in prefixes):
                    continue
                validate_frozen(study)
                evaluators, baselines = {}, {}
                for level, index in (("fine", 2), ("coarse", 0)):
                    prepared_at = perf_counter()
                    try:
                        original, _ = physics.prep(Path(case["root"]), "OMP", ref, index)
                        baseline = measurement_record(physics.measurement(original.baseline, ref))
                        require(baseline["eligible"], "original-wall baseline is ineligible")
                        evaluators[level] = physics.evaluator(original, 20, "gpu")
                        baselines[level] = baseline
                        preparation = dict(status="completed", baseline=baseline, seconds=perf_counter()-prepared_at)
                    except Exception as exc:
                        preparation = dict(status="failed", error=type(exc).__name__ + ": " + str(exc), seconds=perf_counter()-prepared_at)
                        store.write("preparation/" + case["case_id"] + "/" + invocation.split("/")[-1] + "/" + level + ".json", preparation)
                        raise
                    store.write("preparation/" + case["case_id"] + "/" + invocation.split("/")[-1] + "/" + level + ".json", preparation)
                    del original
                # Seal every seed panel before any additional panel scoring.
                # Fixed seeds can reproduce candidates that the live arms scored.
                for seed, prefix in zip(SEEDS, prefixes):
                    name = prefix + "/panel.json"
                    if not store.exists(name):
                        panel = dict(generate_panel(evaluators["fine"], seed), case_id=case["case_id"],
                                     protocol_sha256=store.sha256("protocol.json"),
                                     membership_scope=PANEL_SCOPE, held_out_validation_claim=False)
                        panel["live_overlap"] = dict(panel_live_overlap(panel, arms),
                                                     live_evidence_sha256=store.sha256("live-evidence.json"))
                        store.write(name, panel)
                    panel = store.read(name)
                    validate_panel(panel, seed)
                    require(panel["case_id"] == case["case_id"] and panel["protocol_sha256"] == store.sha256("protocol.json"), "panel protocol differs")
                    require(panel["live_overlap"] == dict(panel_live_overlap(panel, arms),
                                                           live_evidence_sha256=store.sha256("live-evidence.json")),
                            "panel live-overlap provenance differs")
                for seed, prefix in zip(SEEDS, prefixes):
                    if store.exists(prefix + "/result.json"):
                        continue
                    panel = store.read(prefix + "/panel.json")
                    result = score_panel(store, prefix, panel, evaluators, baselines["fine"])
                    result.update(case_id=case["case_id"], seed=seed, coarse_baseline=baselines["coarse"],
                                  live_overlap=panel["live_overlap"], membership_scope=PANEL_SCOPE,
                                  held_out_validation_claim=False)
                    store.write(prefix + "/result.json", result)
                    print("REFERENCE", case["case_id"], seed, result["assessment"]["regret_relative"], flush=True)
                del evaluators
            validate_frozen(study)
            verify_external_hashes(study, hashes)
            store.write(invocation + ".finished.json", dict(status="completed", seconds=perf_counter()-start))
        except Exception as exc:
            if not store.exists(invocation + ".finished.json"):
                store.write(invocation + ".finished.json", dict(status="failed", seconds=perf_counter()-start,
                                                               error=type(exc).__name__ + ": " + str(exc)))
            raise
        results = []
        for case in frozen["cases"]:
            for seed in SEEDS:
                name = case["case_id"] + "/" + str(seed) + "/result.json"
                result = store.read(name)
                require(result["assessment"] == summarize_panel(result["rows"], result["baseline"]), "stored panel assessment differs")
                results.append(dict(case_id=case["case_id"], seed=seed, complete=result["complete"],
                                    result_sha256=store.sha256(name), live_overlap=result["live_overlap"],
                                    held_out_validation_claim=False, **result["assessment"]))
        summary = dict(schema="nf-anytime-reference-summary.v1", protocol_sha256=store.sha256("protocol.json"),
                       rows=results, costs=cost_summary(store), equal_wall_time_comparison=False,
                       global_optimum_claim=False, held_out_validation_claim=False, membership_scope=PANEL_SCOPE)
        store.write("summary.json", summary)
        store.verify_all()
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--freeze-only", action="store_true", help="seal the post-trial code/protocol without any physical calculation")
    args = parser.parse_args()
    study = args.study.resolve()
    require(study.drive.lower() == "d:", "study must remain on D as requested")
    if args.freeze_only:
        from d3d_optimization.campaign_state import process_lock
        with process_lock(study / "posttrial.lock"):
            freeze_post(study, "reference", reference_details(), [Path(__file__)])
        print("Frozen", study / "reference" / "protocol.json", flush=True)
    else:
        result = run_reference(study)
        print("REFERENCE COMPLETE", len(result["rows"]), "panels", flush=True)


if __name__ == "__main__":
    main()
