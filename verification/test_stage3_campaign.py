import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import pytest

from verification import stage3_campaign as campaign


def cases(tmp_path):
    rows = []
    for i in range(12):
        path = tmp_path / f"case-{i}.parquet"
        data = f"synthetic equilibrium {i}".encode()
        path.write_bytes(data)
        rows.append(dict(id=f"case-{i}", path=str(path), sha256=hashlib.sha256(data).hexdigest(),
                         status="ready", rank=i+1, time_ms=1400.))
    return rows


def test_design_balances_leading_arm_and_preserves_all_slots(tmp_path):
    rows = cases(tmp_path)
    rows[3]["status"] = "time_missing"
    schedule = campaign.build_tasks(rows)
    assert len(schedule) == len({t["relative_output"] for t in schedule}) == 108
    assert all(t["budget_seconds"] == 600 for t in schedule)
    for row in rows:
        group = [t for t in schedule if t["case_id"] == row["id"]]
        assert len(group) == 9
        assert {t["registration_status"] for t in group} == {row["status"]}
        first = [next(t for t in group if t["seed"] == seed)["arm"] for seed in campaign.SEEDS]
        assert Counter(first) == {arm: 1 for arm in campaign.ARMS}
        for seed in campaign.SEEDS:
            assert {t["arm"] for t in group if t["seed"] == seed} == set(campaign.ARMS)
    rows[3]["id"] = "../escape"
    with pytest.raises(ValueError):
        campaign.build_tasks(rows)


def test_registration_exclusion_never_invokes_heat_and_is_resumable(tmp_path):
    rows = cases(tmp_path)
    rows[0].update(status="no_exact_frame", reason="1400 ms missing")
    task = campaign.build_tasks(rows)[0]
    seal = dict(registration_sha256="12"*32)
    assert task["case_id"] == rows[0]["id"]
    report = campaign.excluded_receipt(task, seal)
    dest = tmp_path/task["relative_output"]
    dest.mkdir(parents=True)
    (dest/"receipt.json").write_text(json.dumps(report))
    assert campaign.verified_receipt(dest, task, seal)["status"] == "registration_excluded"
    assert report["elapsed_s"] is None
    assert report["detail"] == "1400 ms missing"
    edited = dict(task, arm="different")
    with pytest.raises(ValueError):
        campaign.verified_receipt(dest, edited, seal)


def test_seal_checks_source_and_registration_content_before_run(tmp_path):
    meta = tmp_path/"inputs"/"registration.json"
    meta.parent.mkdir()
    meta.write_text(json.dumps({"cases": cases(tmp_path)}))
    src = tmp_path/"code.py"
    src.write_text("x=1")
    plan = campaign.make_plan(tmp_path, root=tmp_path, files=[src], revision="synthetic")
    campaign.check_plan(tmp_path, plan, root=tmp_path)
    src.write_text("x=2")
    with pytest.raises(ValueError, match="code.py"):
        campaign.check_plan(tmp_path, plan, root=tmp_path)
    src.write_text("x=1")
    meta.write_text(meta.read_text()+" ")
    with pytest.raises(ValueError, match="registration"):
        campaign.check_plan(tmp_path, plan, root=tmp_path)


def test_real_invalid_parquet_is_a_retained_cli_outcome_and_resume_skips(tmp_path, monkeypatch):
    task = campaign.build_tasks(cases(tmp_path))[0]
    plan = dict(tasks=[task], registration_sha256="5a"*32)
    monkeypatch.setattr(campaign, "load_seal", lambda study: plan)
    receipt = campaign.worker(tmp_path, 0)
    dest = tmp_path/task["relative_output"]
    assert receipt["status"] == "completed" and receipt["returncode"] == 2
    assert json.loads((dest/"result.json").read_text())["status"] == "input_error"
    monkeypatch.setattr(campaign.subprocess, "run", lambda *_a, **_kw: pytest.fail("repeated completed job"))
    assert campaign.worker(tmp_path, 0) == receipt


def test_seal_repairs_interruption_between_plan_and_checksum(tmp_path, monkeypatch):
    src = tmp_path/"inputs"/"registration.json"
    src.parent.mkdir()
    src.write_text(json.dumps({"cases": cases(tmp_path)}))
    old = campaign.make_plan(tmp_path, root=tmp_path, files=[], revision="test")
    (tmp_path/"expected.json").write_text(json.dumps(old))
    assert campaign.seal(tmp_path) == old
    assert campaign.load_seal(tmp_path) == old


def test_clean_worker_import_does_not_preload_scientific_libraries():
    import subprocess
    code = ("import sys; import verification.stage3_campaign; "
            "assert not set(sys.modules).intersection({'numpy','scipy','pyarrow','cupy'})")
    subprocess.run([sys.executable, "-B", "-c", code], cwd=campaign.ROOT, check=True)


def test_recover_original_result_after_wrapper_crash_without_retiming(tmp_path, monkeypatch):
    task = campaign.build_tasks(cases(tmp_path))[0]
    plan = dict(tasks=[task], registration_sha256="11"*32)
    monkeypatch.setattr(campaign, "load_seal", lambda _: plan)
    original_save = campaign.public_save
    def crash(path, value):
        if path.name == "completion.json":
            raise OSError("interrupted after CLI result")
        original_save(path, value)
    monkeypatch.setattr(campaign, "public_save", crash)
    with pytest.raises(OSError, match="interrupted"):
        campaign.worker(tmp_path, 0)
    monkeypatch.setattr(campaign, "public_save", original_save)
    monkeypatch.setattr(campaign.runpy, "run_path", lambda *_a, **_kw: pytest.fail("repeated thermal run"))
    receipt = campaign.worker(tmp_path, 0)
    assert receipt["status"] == "completed"
    assert receipt["timing_source"] == "recovered_cli_result"
    assert receipt["returncode"] is None and receipt["elapsed_s"] is None
    assert receipt["source_attempt"].endswith("/001")
    assert campaign.worker(tmp_path, 0) == receipt


def test_independent_processes_serialize_one_fixture_and_record_parent_clock(tmp_path):
    import subprocess
    from concurrent.futures import ThreadPoolExecutor
    from time import sleep, monotonic

    reg = tmp_path/"inputs"/"registration.json"
    reg.parent.mkdir()
    reg.write_text(json.dumps({"cases": cases(tmp_path)}))
    plan = campaign.make_plan(tmp_path, root=tmp_path, files=[], revision="fixture")
    campaign.public_save(tmp_path/"expected.json", plan)
    campaign.seal(tmp_path)
    task = plan["tasks"][0]
    cmd = [sys.executable,"-B","-u","-m","verification.stage3_campaign","worker","--study",str(tmp_path),"--index","0"]
    with ThreadPoolExecutor(max_workers=1) as pool:
        with campaign.resource_lock(tmp_path/"worker.lock"):
            other = subprocess.Popen(cmd, cwd=campaign.ROOT, env=campaign.environment(tmp_path))
            observed = pool.submit(campaign.observe_worker, tmp_path, task, cmd)
            deadline = monotonic()+10
            while not list((tmp_path/"observations").rglob("*.started.json")) and monotonic() < deadline:
                sleep(.05)
            assert not (tmp_path/"attempts").exists() and other.poll() is None
        assert other.wait(30) == 0
        observed.result(timeout=30)
    folders = list((tmp_path/"attempts").rglob("command.json"))
    assert len(folders) == 1
    receipt = campaign.verified_receipt(tmp_path/task["relative_output"], task, plan)
    assert receipt["returncode"] == 2
    parent = campaign.read(next((tmp_path/"observations").rglob("*.finished.json")))
    assert parent["elapsed_s"] > 0 and parent["returncode"] == 0


def test_receipt_detects_event_tampering_and_binds_task_identity(tmp_path, monkeypatch):
    task = campaign.build_tasks(cases(tmp_path))[0]
    plan = dict(tasks=[task], registration_sha256="34"*32)
    monkeypatch.setattr(campaign, "load_seal", lambda _: plan)
    out = tmp_path/task["relative_output"]
    receipt = campaign.worker(tmp_path, 0)
    assert campaign.verified_receipt(out, task, plan) == receipt
    assert receipt["elapsed_s"] > 0
    (out/"events.jsonl").write_text("changed")
    with pytest.raises(ValueError, match="events"):
        campaign.verified_receipt(out, task, plan)
