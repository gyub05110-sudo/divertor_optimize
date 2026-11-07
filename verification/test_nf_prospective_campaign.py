"""Regression checks for the prospective split, budget and process evidence."""
from pathlib import Path
import json
import pytest

from verification import nf_prospective_campaign as campaign
from verification.nf_anytime_reference import EvidenceStore
from d3d_optimization.broadening_study import digest


def cases(prefix):
    return [dict(case_id=f"{prefix}{i}", root=f"D:/inputs/{prefix}{i}",
                 hashes={campaign.INPUTS[0]:f"{prefix}{i}hash"}) for i in range(3)]


def test_primary_and_long_pairs_are_fixed_before_results():
    selected = dict(policy="fine_adaptive", population=24, projected_seeds=8)
    jobs = campaign.holdout_jobs(cases("new"), selected)
    assert len(jobs) == 18
    assert len({j["relative_dir"] for j in jobs}) == 18
    assert sum(j["budget_s"] for j in jobs) == 18000
    for i in range(0, len(jobs), 2):
        left, right = jobs[i:i+2]
        assert (left["case_id"], left["seed"], left["budget_s"]) == (right["case_id"], right["seed"], right["budget_s"])
        assert {left["policy"], right["policy"]} == {"coarse_adaptive", "fine_adaptive"}
        assert left["budget_s"] == (1800 if left["seed"] == 20711 else 600)
    assert jobs[0]["policy"] != jobs[2]["policy"] != jobs[4]["policy"]
    assert jobs[0]["policy"] != jobs[6]["policy"]


def test_only_registered_development_choices_can_enter_holdout():
    with pytest.raises(ValueError):
        campaign.holdout_jobs(cases("new"), dict(policy="fine_adaptive", population=64, projected_seeds=12))
    jobs = campaign.development_jobs(cases("old"))
    assert len(jobs) == 6 and {r["seed"] for r in jobs} == {20611}
    assert {(r["population"],r["projected_seeds"]) for r in jobs} == {(24,8),(32,12)}


def test_previously_used_equilibrium_cannot_be_relabelled():
    old = cases("old")
    new = cases("new")
    campaign.check_disjoint(old,new)
    new[1]["hashes"] = old[0]["hashes"]
    with pytest.raises(ValueError):
        campaign.check_disjoint(old,new)
    new = cases("new")
    new[2]["case_id"] = old[0]["case_id"]
    with pytest.raises(ValueError):
        campaign.check_disjoint(old,new)


def test_process_receipt_binds_exact_raw_bytes_even_after_normal_exit(tmp_path):
    out = tmp_path/"arm"
    store = EvidenceStore(out)
    job = {"case_id":"one", "phase":"holdout", "registration_sha256":"aa"*32}
    store.write("job.json",job)
    (out/"events.jsonl").write_bytes(b'{"valid":true}\n')
    (out/"result.json").write_text("{}", encoding="utf-8")
    (out/"stdout.log").write_bytes(b"")
    process = dict(job_sha256=store.sha256("job.json"), registration_sha256=job["registration_sha256"],
                   launch_ns=100,exit_observed_ns=500,elapsed_s=400/1e9,returncode=0,error=None,
                   events_sha256=digest(out/"events.jsonl"),result_sha256=digest(out/"result.json"),
                   stdout_sha256=digest(out/"stdout.log"))
    store.write("process.json",process)
    assert campaign.verify_process(out,job)==process
    (out/"events.jsonl").write_bytes(b'{"valid":false}\n')
    with pytest.raises(ValueError):
        campaign.verify_process(out,job)


def test_parent_observed_timeout_is_retained_and_cannot_fabricate_fine_outcome(tmp_path):
    out=tmp_path/"arm"
    store=EvidenceStore(out)
    job=dict(case_id="one",phase="holdout",registration_sha256="cc"*32,policy="fine_adaptive",
             population=24,projected_seeds=8,budget_s=600,seed=20729,relative_dir="live/one")
    store.write("job.json",job)
    (out/"stdout.log").write_bytes(b"timeout")
    p=dict(job_sha256=store.sha256("job.json"),registration_sha256=job["registration_sha256"],
           launch_ns=10,exit_observed_ns=110,elapsed_s=1e-7,returncode=None,error="watchdog",
           events_sha256=None,result_sha256=None,stdout_sha256=digest(out/"stdout.log"))
    store.write("process.json",p)
    assert campaign.verify_process(out,job)["returncode"] is None
    row=campaign.assess(out,job,{}, {})
    assert row["de_requests"] == 0 and row["normal"] is False and row["parent_elapsed_s"] == p["elapsed_s"]
    (out/"result.json").write_text("{}",encoding="utf8")
    with pytest.raises(ValueError):
        campaign.verify_process(out,job)


def test_child_failure_before_physics_keeps_parent_clock_and_verified_audit(monkeypatch, tmp_path):
    from verification import nf_prospective_worker as worker
    job=dict(case_id="one",phase="development",registration_sha256="a"*64,policy="fine_adaptive",
             population=24,projected_seeds=8,budget_s=600,seed=20611,relative_dir="development/one",input_root="D:/inputs/one")
    out=tmp_path/"arm"
    store=EvidenceStore(out)
    store.write("job.json",job)
    result=dict(case_id="one",input_root=job["input_root"],launch_ns=100,job_sha256=store.sha256("job.json"),
                registration_sha256=job["registration_sha256"],levels={},elapsed_s=600.1,stop="deadline")
    (out/"result.json").write_text(json.dumps(result),encoding="utf8")
    (out/"events.jsonl").write_bytes(b"")
    (out/"stdout.log").write_bytes(b"")
    p=dict(job_sha256=store.sha256("job.json"),registration_sha256=job["registration_sha256"],launch_ns=100,
           exit_observed_ns=601_000_000_100,elapsed_s=601.,returncode=1,error=None,events_sha256=digest(out/"events.jsonl"),
           result_sha256=digest(out/"result.json"),stdout_sha256=digest(out/"stdout.log"))
    store.write("process.json",p)
    monkeypatch.setattr(worker,"audit_arm",lambda *a,**kw: {"normal":False,"elapsed_s":600.1})
    row=campaign.assess(out,job,{"heating":1},{})
    assert row["status"]=="arm_failed" and row["elapsed_s"]==600.1 and row["parent_elapsed_s"]==601.

