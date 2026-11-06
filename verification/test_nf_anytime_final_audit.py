"""Synthetic final-audit checks. Execute only after the timed campaign finishes."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest

from verification import nf_anytime_final_audit as audit


def measurement(peak=10., eligible=True):
    return dict(eligible=eligible, peak=peak, A2=2., errors=[],
                ledger=dict(prescribed=1., deposited=1., unhit=0., unresolved=0., balance_error=0.))


def candidate(vector):
    value = dict(source="OMP", radius_mm=20, vector=list(vector))
    identity = hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()
    return identity, value


def population():
    x = np.random.default_rng(17).uniform(-.03, .03, (64, 16))
    singular = np.linalg.svd((x-x.mean(axis=0))/.2, compute_uv=False)
    return dict(vectors=x.tolist(), count=64, dimension=16, rank=16,
                axis_min=x.min(axis=0).tolist(), axis_max=x.max(axis=0).tolist(),
                axis_span=np.ptp(x, axis=0).tolist(), center=x.mean(axis=0).tolist(),
                scale=[.2]*16, singular_values=singular.tolist(),
                rank_method="centered_and_fixed_bound_width_scaled", rank_tolerance=1e-12,
                axis_span_tolerance=1e-12, full_space_passed=True,
                full_space_required=True, startup=True, archive_count=0,
                inherited_count=0, archive_origin="none")


def record():
    initial = population()
    identity, design = candidate(initial["vectors"][0])
    common = dict(id=identity, candidate=design, search_stage="differential_evolution",
                  proposal_kind="population_initialization", generation=0, index=0)
    events = [
        dict(stage="cold_inputs", start_s=.1, end_s=.2, status="completed"),
        dict(stage="prepare_fine", start_s=.3, end_s=.8, status="completed"),
        dict(stage="baseline", start_s=1., end_s=2., status="completed", measurement=measurement()),
        dict(stage="build_fine", start_s=2.1, end_s=2.2, status="completed"),
        dict(stage="initialize", start_s=3., end_s=4., status="completed", level="fine", restart=0,
             archive_count=0, inherited_count=0, full_space_required=True, initial_population=initial),
        dict(stage="proposal", start_s=5., end_s=6., status="completed", level="fine", **common),
        dict(stage="fine", start_s=6.1, end_s=7., status="completed", measurement=measurement(8.),
             cached=False, physics_evaluated=True, objective=.8, **common),
    ]
    counts = dict(fine=1, coarse=0, eligible_fine=1, eligible_coarse=0,
                  cached=0, proposals=1, rejected=0, late=0, failed=0)
    checkpoints = {str(t): dict(t_s=float(t), available=True,
        incumbent=dict(kind="new_design", id=identity, candidate=design, peak=8., A2=2.),
        first_5pct_s=7., counts=deepcopy(counts)) for t in (200, 400, 600)}
    return dict(policy="fine_adaptive", seed=20611, budget_s=600., startup_s=0., source="OMP", radius_mm=20,
                reference={"power_W": 1e6}, levels=deepcopy(audit.LEVELS),
                adaptive_config=dict(population=64, projected_seeds=12,
                                     initial_population_gate="startup_only_full_affine_dimension"),
                events=events, candidates=[dict(i=0, id=identity, candidate=design, first_seen_s=5.5)],
                promotions=[], completed_draws=0, rejected={},
                baseline=measurement(), baseline_ready_s=2., stop="deadline", error=None,
                elapsed_s=600., checkpoints=checkpoints, final=deepcopy(checkpoints["600"]),
                parity=dict(coarse=None, fine=dict(candidate=design, passed=True,
                    peak_relative_difference=0., density_scaled_max_difference=0.,
                    quadrature_absolute_difference=0., cpu_seconds=1., gpu_seconds=1.)))


def test_completed_valid_arm_has_independent_endpoints_and_low_coverage_is_diagnostic():
    r = record()
    result = audit.audit_record(r, deepcopy(r["events"]))
    assert result["checkpoints"]["600"]["incumbent"]["peak"] == 8.
    assert result["coverage"]["600"]["fine"]["eligible"]["rank"] == 0
    assert result["coverage"]["600"]["fine"]["eligible"]["unique_designs"] == 1


@pytest.mark.parametrize("change,match", [
    ("wrong_level", "initialization_level"), ("duplicate_startup", "initialization_sequence"),
    ("rank", "population"), ("cache_without_origin", "cache_evidence"),
    ("parity", "parity"), ("future_checkpoint", "checkpoint"),
    ("raw_changed", "raw_event_equality"), ("postdeadline_start", "event_clock"),
])
def test_audit_rejects_material_evidence_corruption(change, match):
    r = record()
    raw = deepcopy(r["events"])
    if change == "wrong_level":
        r["events"][4]["level"] = "coarse"
    elif change == "duplicate_startup":
        extra = deepcopy(r["events"][4])
        extra.update(start_s=4.1, end_s=4.2)
        r["events"].insert(5, extra)
    elif change == "rank":
        r["events"][4]["initial_population"]["vectors"] = np.zeros((64, 16)).tolist()
    elif change == "cache_without_origin":
        r["events"][-1]["cached"] = True
    elif change == "parity":
        r["parity"]["fine"] = None
    elif change == "future_checkpoint":
        r["checkpoints"]["200"]["incumbent"]["peak"] = 1.
    elif change == "raw_changed":
        raw[-1]["measurement"]["peak"] = 9.
    else:
        r["events"][-1].update(start_s=601., end_s=610., status="late")
        r["elapsed_s"] = 611.
    if change != "raw_changed":
        raw = deepcopy(r["events"])
    with pytest.raises(audit.AuditError, match=match):
        audit.audit_record(r, raw)


def test_late_fine_value_is_retained_but_cannot_win_even_after_deadline():
    events = record()["events"]
    events[-1].update(end_s=610., status="late", measurement=measurement(1.))
    for cutoff in (200, 400, 600, 1000):
        result = audit.checkpoint(events, cutoff, 600.)
        assert result["incumbent"]["kind"] == "original_wall"
        assert result["first_5pct_s"] is None
    assert audit.checkpoint(events, 1000, 600.)["counts"]["late"] == 1
    assert events[-1]["measurement"]["peak"] == 1.


def test_exact_checkpoint_boundary_is_inclusive():
    events = record()["events"]
    events[-1]["end_s"] = 200.
    assert audit.checkpoint(events, 199.9, 600.)["incumbent"]["kind"] == "original_wall"
    assert audit.checkpoint(events, 200., 600.)["first_5pct_s"] == 200.


def test_random_policy_cannot_reuse_an_adaptive_initialization():
    r = record()
    r["policy"] = "random_fine"
    with pytest.raises(audit.AuditError, match="initialization_level"):
        audit.audit_record(r, deepcopy(r["events"]))


def test_coarse_policy_cannot_relabel_a_direct_fine_score():
    r = record()
    r["policy"] = "coarse_adaptive"
    r["events"][4]["level"] = "coarse"
    r["events"][5]["level"] = "coarse"
    with pytest.raises(audit.AuditError, match="policy_fine_stage"):
        audit.audit_record(r, deepcopy(r["events"]))


def test_real_quartet_order_and_minimum_are_checked_before_accepting_fine_choice():
    r = record()
    r.update(policy="coarse_adaptive", candidates=[], events=r["events"][:5])
    r["events"][4]["level"] = "coarse"
    rows = []
    for index in range(4):
        identity, design = candidate(r["events"][4]["initial_population"]["vectors"][index])
        start = 5.+3*index
        common = dict(id=identity, candidate=design, search_stage="differential_evolution",
                      proposal_kind="population_initialization", generation=0, index=index)
        r["candidates"].append(dict(i=index, id=identity, candidate=design, first_seen_s=start+.2))
        r["events"].append(dict(stage="proposal", level="coarse", status="completed",
                                start_s=start, end_s=start+.4, **common))
        score = dict(stage="coarse", status="completed", start_s=start+.5, end_s=start+1,
                     cached=False, physics_evaluated=True, objective=(8.-index)/10.,
                     measurement=measurement(8.-index), **common)
        rows.append(score)
        r["events"].append(score)
    winner = rows[-1]
    r["events"].append(dict(stage="fine", status="completed", start_s=15.5, end_s=16.5,
                            id=winner["id"], candidate=winner["candidate"], search_stage="coarse_promotion",
                            proposal_kind="external_validation", generation=None, index=None,
                            cached=False, physics_evaluated=True, objective=.45, measurement=measurement(4.5)))
    r["promotions"] = [dict(ids=[row["id"] for row in rows], selected=winner["id"], at_s=15.3)]
    r["parity"]["coarse"] = deepcopy(r["parity"]["fine"])
    r["parity"]["fine"]["candidate"] = winner["candidate"]
    for value in r["checkpoints"].values():
        value.update(first_5pct_s=16.5,
                     incumbent=dict(kind="new_design", id=winner["id"], candidate=winner["candidate"], peak=4.5, A2=2.))
        value["counts"].update(coarse=4, eligible_coarse=4, proposals=4)
    r["final"] = deepcopy(r["checkpoints"]["600"])
    assert audit.audit_record(r, deepcopy(r["events"]))["checkpoints"]["600"]["incumbent"]["peak"] == 4.5
    r["events"][-1].update(id=rows[0]["id"], candidate=rows[0]["candidate"])
    with pytest.raises(audit.AuditError, match="promotion_choice"):
        audit.audit_record(r, deepcopy(r["events"]))


def test_deadline_placeholder_is_not_counted_as_returned_deposition():
    r = record()
    item = r["events"][-1]
    item.update(status="late", end_s=601.,
                measurement=dict(eligible=False, errors=["Deadline: "], ledger={}, peak=None, A2=None))
    item.pop("objective")
    item.pop("physics_evaluated")
    assert not audit.returned_deposition(item)
    assert audit.returned_deposition(record()["events"][-1])


def test_byte_hash_check_rejects_a_changed_artifact(tmp_path):
    path = tmp_path / "source.py"
    path.write_bytes(b"original\n")
    expected = audit.digest(path)
    path.write_bytes(b"changed\n")
    with pytest.raises(audit.AuditError, match="file_hash"):
        audit.verify_hash(path, expected, audit.Checks())


def test_raw_reader_keeps_only_complete_prefix_of_killed_process(tmp_path):
    path = tmp_path / "events.jsonl"
    path.write_text('{"stage":"cold_inputs"}\n{"unfinished":', encoding="utf8")
    events, partial = audit.read_events(path, allow_partial_tail=True)
    assert partial and events == [{"stage": "cold_inputs"}]
    with pytest.raises((ValueError, audit.AuditError)):
        audit.read_events(path, allow_partial_tail=False)


def test_incomplete_campaign_cannot_finalize_or_fabricate_checkpoints(tmp_path, monkeypatch):
    (tmp_path / "frozen.json").write_text("{}", encoding="utf8")
    frozen = {"cases": [dict(case_id=name, root=str(tmp_path / name)) for name in ("a", "b", "c")]}
    monkeypatch.setattr(audit, "audit_frozen", lambda *_: frozen)
    result = audit.audit_study(tmp_path)
    assert result["arm_count"] == 27
    assert not result["finalized"] and not result["passed"]
    assert result["observed_status_counts"] == {"not_terminal": 27}
    assert all(not row["certified_checkpoints"] and "checkpoints" not in row for row in result["rows"])


def test_terminal_process_failures_are_retained_without_inventing_endpoints(tmp_path, monkeypatch):
    (tmp_path / "frozen.json").write_text("{}", encoding="utf8")
    frozen = {"cases": [dict(case_id=name, root=str(tmp_path / name)) for name in ("a", "b", "c")]}
    monkeypatch.setattr(audit, "audit_frozen", lambda *_: frozen)
    sequence = 0
    for case_index, case in enumerate(frozen["cases"]):
        for seed_index, seed in enumerate(audit.SEEDS):
            shift = (case_index+seed_index) % 3
            for policy in audit.POLICIES[shift:]+audit.POLICIES[:shift]:
                folder = tmp_path / "live" / case["case_id"] / str(seed) / policy
                folder.mkdir(parents=True)
                environment = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE")}
                environment.update({key: str(folder / "cache" / key) for key in
                                    ("CUPY_CACHE_DIR", "CUDA_CACHE_PATH", "TEMP", "TMP", "MPLCONFIGDIR", "HF_HOME")})
                process = dict(launch_ns=sequence*2000, exit_observed_ns=sequence*2000+1000,
                               elapsed_s=.000001, returncode=9, error=None, environment=environment)
                (folder / "process.json").write_text(json.dumps(process), encoding="utf8")
                sequence += 1
    result = audit.audit_study(tmp_path)
    assert result["finalized"] and result["passed"]  # Integrity succeeds; the experiments failed.
    assert result["observed_status_counts"] == {"process_failed": 27}
    assert all(row["checkpoints"] is None and not row["certified_checkpoints"] for row in result["rows"])
