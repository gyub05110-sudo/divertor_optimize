"""Protocol immutability and serial order are properties of the experiment."""
import json
from pathlib import Path
import pytest
from copy import deepcopy


def driver():
    from verification import nf_anytime_campaign
    return nf_anytime_campaign


def test_each_method_has_every_serial_position_within_each_equilibrium():
    d = driver()
    for case in range(3):
        orders = [d.arm_order(case, seed) for seed in d.SEEDS]
        for position in range(3):
            assert {order[position] for order in orders} == set(d.POLICIES)
        assert all(set(order) == set(d.POLICIES) for order in orders)


def test_checkpoint_requires_actual_completed_measurement_and_retains_late_data():
    d = driver()
    record = dict(budget_s=600, startup_s=1., baseline={"eligible": True, "peak": 10., "A2": 2.},
                  baseline_ready_s=120., stop="deadline", events=[
        dict(stage="fine", start_s=190., end_s=210., status="completed",
             id="a", candidate={}, measurement={"eligible": True, "peak": 6., "A2": 3.}),
        dict(stage="fine", start_s=590., end_s=610., status="late",
             id="b", candidate={}, measurement={"eligible": True, "peak": 2., "A2": 5.})])
    result = d.reconstruct(record, 200)
    assert result["peak"] == 10. and result["fine_count"] == 0
    assert d.reconstruct(record, 400)["peak"] == 6.
    assert d.reconstruct(record, 600)["peak"] == 6.
    assert record["events"][-1]["measurement"]["peak"] == 2.


def test_no_valid_output_before_baseline_ready():
    d = driver()
    record = dict(budget_s=600, baseline={"eligible": True, "peak": 10., "A2": 2.},
                  baseline_ready_s=201., events=[])
    assert not d.reconstruct(record, 200)["available"]


def test_available_original_is_not_treated_as_search_failure():
    d = driver()
    record = dict(budget_s=600, baseline={"eligible": True, "peak": 10., "A2": 2.},
                  baseline_ready_s=100., events=[])
    got = d.reconstruct(record, 200)
    assert got["available"] and got["kind"] == "original_wall"
    assert got["peak_reduction_pct"] == 0


def test_unchanged_freeze_reuses_protocol_but_changed_code_is_rejected(tmp_path, monkeypatch):
    d = driver()
    code = tmp_path / "code"
    code.mkdir()
    file = code / "runner.py"
    file.write_text("frozen\n", encoding="utf8")
    study = tmp_path / "study"
    study.mkdir()
    (study / "protocol.ko.md").write_text("fixed study\n", encoding="utf8")
    source = tmp_path / "source"
    source.mkdir()
    for name in d.SOURCE_RECORDS:
        (source / name).write_text("{}", encoding="utf8")
    monkeypatch.setattr(d, "CODE_ROOT", code)
    monkeypatch.setattr(d, "INPUTS", ())
    monkeypatch.setattr(d, "source_files", lambda: [file])
    monkeypatch.setattr(d, "load_cases", lambda path: [{"case_id": "one", "root": str(source), "hashes": {}}])
    original = d.freeze(study, source)
    assert original == d.freeze(study, source)
    file.write_text("changed\n", encoding="utf8")
    with pytest.raises(ValueError, match="protocol|frozen|changed"):
        d.freeze(study, source)


def test_bad_peak_or_clock_cannot_be_selected():
    d = driver()
    record = dict(budget_s=600, baseline={"eligible": True, "peak": 10., "A2": 2.},
                  baseline_ready_s=100., events=[dict(stage="fine", start_s=150., end_s=149.,
                    status="completed", measurement={"eligible": True, "peak": 1., "A2": 4.})])
    with pytest.raises(ValueError):
        d.reconstruct(record, 200)
    record["events"][0].update(end_s=151.)
    record["events"][0]["measurement"]["peak"] = float("nan")
    with pytest.raises(ValueError):
        d.reconstruct(record, 200)


def valid_record():
    from d3d_optimization.campaign_geometry import Candidate
    c = Candidate("OMP", 20, (.001,) * 16)
    base = dict(eligible=True, errors=[], ledger={}, peak=10., A2=2.)
    events = [dict(stage="baseline", start_s=1., end_s=100., status="completed", measurement=base),
              dict(stage="fine", start_s=190., end_s=210., status="completed", id=c.identity,
                   candidate=c.record(), measurement=dict(eligible=True, peak=6., A2=3.))]
    return dict(budget_s=600, startup_s=1., baseline=base, baseline_ready_s=100., events=events,
                candidates=[dict(id=c.identity, candidate=c.record())], stop="deadline", elapsed_s=600.)


@pytest.mark.parametrize("mutation", ["postdeadline", "overlap", "bad_status", "wrong_id", "wrong_baseline", "raw_changed"])
def test_event_audit_rejects_causal_and_provenance_corruption(mutation):
    d = driver()
    record = valid_record()
    raw = deepcopy(record["events"])
    if mutation == "postdeadline":
        record["events"][-1].update(start_s=601., end_s=610., status="late")
        record["elapsed_s"] = 611.
    elif mutation == "overlap":
        record["events"][-1]["start_s"] = 99.
    elif mutation == "bad_status":
        record["events"][-1]["status"] = "late"
    elif mutation == "wrong_id":
        record["events"][-1]["id"] = "wrong"
    elif mutation == "wrong_baseline":
        record["baseline_ready_s"] = 99.
    elif mutation == "raw_changed":
        raw[-1]["measurement"]["peak"] = 7.
    if mutation != "raw_changed":
        raw = deepcopy(record["events"])
    with pytest.raises(ValueError):
        d.validate_record(record, raw)


def test_event_audit_allows_retained_terminal_late_measurement():
    d = driver()
    record = valid_record()
    record["events"][-1].update(start_s=590., end_s=610., status="late")
    record["elapsed_s"] = 611.
    d.validate_record(record, deepcopy(record["events"]))
    assert d.reconstruct(record, 600)["peak"] == 10.


def promotion_record(late=False):
    from d3d_optimization.campaign_geometry import Candidate
    from verification.nf_anytime_search import snapshot
    r = valid_record()
    r.update(policy="coarse_adaptive", events=r["events"][:1], candidates=[], promotions=[])
    for i in range(4):
        c = Candidate("OMP", 20, (.001 * (i + 1),) * 16)
        r["candidates"].append(dict(id=c.identity, candidate=c.record()))
        r["events"].append(dict(stage="coarse", start_s=110.+10*i, end_s=115.+10*i,
                                status="completed", id=c.identity, candidate=c.record(),
                                measurement=dict(eligible=True, peak=8.-i, A2=3.)))
    selected = r["candidates"][-1]
    r["promotions"].append(dict(ids=[c["id"] for c in r["candidates"]], selected=selected["id"], at_s=150.))
    r["events"].append(dict(stage="fine", start_s=590. if late else 155., end_s=610. if late else 170.,
                            status="late" if late else "completed", search_stage="coarse_promotion",
                            **selected, measurement=dict(eligible=True, peak=4., A2=3.)))
    r["elapsed_s"] = 611. if late else 600.
    r["checkpoints"] = {str(t): snapshot(r,t) for t in (200,400,600)}
    r["final"] = r["checkpoints"]["600"]
    return r


@pytest.mark.parametrize("late", [False, True])
def test_promotion_and_checkpoint_audit_uses_coarse_evidence_and_completion_time(late):
    d = driver()
    r = promotion_record(late)
    assert d.validate_record(r, deepcopy(r["events"]))
    assert d.validate_decisions(r)
    assert d.reconstruct(r,600)["peak"] == (10. if late else 4.)


def test_promotion_choice_cannot_be_substituted_after_scoring():
    d = driver()
    r = promotion_record()
    r["events"][-1].update(r["candidates"][0])
    with pytest.raises(ValueError, match="coarse choice"):
        d.validate_decisions(r)


def test_saved_checkpoint_is_checked_against_raw_events():
    d = driver()
    r = promotion_record()
    r["checkpoints"]["200"]["incumbent"]["peak"] = 1.
    with pytest.raises(ValueError, match="incumbent"):
        d.validate_decisions(r)


def test_runtime_source_closure_follows_transitive_local_imports(tmp_path, monkeypatch):
    d = driver()
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg/start.py").write_text("from .nested import helper\nimport outside\n", encoding="utf8")
    (tmp_path / "pkg/nested.py").write_text("from last import end\n", encoding="utf8")
    (tmp_path / "outside.py").write_text("import os\n", encoding="utf8")
    (tmp_path / "last.py").write_text("end=1\n", encoding="utf8")
    monkeypatch.setattr(d, "CODE_ROOT", tmp_path)
    actual = d.local_source_closure([tmp_path / "pkg/start.py"])
    assert {p.relative_to(tmp_path).as_posix() for p in actual} == {
        "pkg/start.py", "pkg/nested.py", "outside.py", "last.py"}


def test_coarse_policy_cannot_bypass_promotions_by_relabelling_a_fine_event():
    from verification.nf_anytime_search import snapshot
    d = driver()
    r = valid_record()
    r.update(policy="coarse_adaptive", promotions=[])
    r["events"][-1]["search_stage"] = "differential_evolution"
    r["checkpoints"] = {str(t): snapshot(r,t) for t in (200,400,600)}
    r["final"] = r["checkpoints"]["600"]
    with pytest.raises(ValueError, match="policy|promotion"):
        d.validate_decisions(r)


def test_initialization_gate_is_recomputed_from_vectors():
    import numpy as np
    from verification.nf_anytime_search import _population_evidence
    d = driver()
    vectors = np.vstack([np.zeros((1,16)), .01*np.eye(16), np.zeros((47,16))])
    evidence = _population_evidence(vectors, True, 0, 0)
    r = valid_record()
    r.update(policy="fine_adaptive", adaptive_config={"population":64,"projected_seeds":12})
    r["events"].insert(1, dict(stage="initialize", start_s=101., end_s=102., status="completed",
                                level="fine", restart=0, initial_population=evidence))
    assert d.validate_initialization(r)
    r["events"][1]["initial_population"]["vectors"] = np.zeros((64,16)).tolist()
    with pytest.raises(ValueError, match="population|diversity"):
        d.validate_initialization(r)
