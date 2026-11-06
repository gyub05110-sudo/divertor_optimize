"""Regression checks for the additive campaign, with no thermal computation."""
import importlib
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter_ns

import pytest


def driver():
    return importlib.import_module("verification.nf_unused_live")


def staged(tmp_path):
    d = driver()
    study = tmp_path / "study"
    root = study / "cases" / "new_case" / "input-root"
    for name in (d.REVIEW, d.GEOMETRY):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((d.CODE_ROOT / name).read_bytes())
    shot = root / d.SHOT
    shot.parent.mkdir(parents=True, exist_ok=True)
    shot.write_bytes(b"synthetic shot, never loaded by the tests")
    d.save(root / d.MANIFEST, dict(shot_sha256=d.digest(shot), time_ms=1400))
    (study / "protocol.ko.md").write_text("Fixed synthetic test protocol\n", encoding="utf8")
    d.save(study / "selection.json", dict(seeds=list(d.SEEDS), budget_s=200,
           protocol_sha256=d.digest(study / "protocol.ko.md")))
    source = dict(path="synthetic-public-source", local=str(shot),
                  sha256=d.digest(shot), bytes=shot.stat().st_size)
    d.save(study / "downloads.json", dict(files=[source]))
    inputs = {name.replace("/", "\\"): d.digest(root / name) for name in d.INPUTS}
    d.save(study / "preflight.json", dict(cases=[dict(case_id="new_case", status="passed",
           root=str(root), time_ms=1400, source=source, inputs=inputs),
           dict(case_id="excluded", status="outside_domain", reason="inner strike outside segment 69")],
           protocol_sha256=d.digest(study / "protocol.ko.md"),
           selection_sha256=d.digest(study / "selection.json"),
           downloads_sha256=d.digest(study / "downloads.json"),
           driver_sha256=d.digest(d.CODE_ROOT / "verification/nf_unused_inputs.py"),
           geometry_driver_sha256=d.digest(d.CODE_ROOT / "verification/nf_equilibrium_times.py")))
    return d, study, root


def test_original_wall_and_manifest_hash_are_checked(tmp_path):
    d, study, root = staged(tmp_path)
    cases = d.load_cases(study)
    assert len(cases) == 1 and cases[0]["case_id"] == "new_case"
    record = d.input_record(root)
    assert set(record["hashes"]) == set(d.INPUTS)
    changed = d.read(root / d.REVIEW)
    changed["wall"][0][0] += .001
    d.save(root / d.REVIEW, changed)
    with pytest.raises(ValueError, match="supplied material wall"):
        d.input_record(root)
    (root / d.REVIEW).write_bytes((d.CODE_ROOT / d.REVIEW).read_bytes())
    (root / d.SHOT).write_bytes(b"changed source")
    with pytest.raises(ValueError, match="shot hash"):
        d.input_record(root)


def test_worker_passes_input_root_seed_and_original_launch_clock(tmp_path, monkeypatch):
    d, study, root = staged(tmp_path)
    captured = {}
    monkeypatch.setattr(d.live, "SEED", d.live.SEED)
    monkeypatch.setattr(d.live.frozen, "ROOT", d.live.frozen.ROOT)
    monkeypatch.setattr(d.live, "SOURCE", d.live.SOURCE)

    def run(policy, start, *, budget_s):
        captured.update(policy=policy, start=start, seed=d.live.SEED,
                        root=d.live.frozen.ROOT, budget=budget_s)
        return dict(stop="failed", error="synthetic baseline failure")

    monkeypatch.setattr(d.live, "run_arm", run)
    out = study / "worker"
    out.mkdir()
    assert d.run_worker(out, root, 20627, "screen", 123456) == 1
    assert captured == dict(policy="screen", start=123456, seed=20627, root=root,
                            budget=200)
    record = d.read(out / "screen.json")
    assert record["launch_ns"] == 123456 and record["seed"] == 20627
    assert record["input_root"] == str(root.resolve()) and record["pid"] > 0


def test_worker_startup_exhausting_budget_never_prepares(tmp_path, monkeypatch):
    d, study, root = staged(tmp_path)
    monkeypatch.setattr(d.live, "SEED", d.live.SEED)
    monkeypatch.setattr(d.live.frozen, "ROOT", d.live.frozen.ROOT)
    monkeypatch.setattr(d.live, "SOURCE", d.live.SOURCE)
    monkeypatch.setattr(d.live.frozen, "prep", lambda *args: pytest.fail("late preparation started"))
    out = study / "deadline-worker"
    out.mkdir()
    assert d.run_worker(out, root, 20611, "direct", perf_counter_ns() - 201_000_000_000) == 0
    record = d.read(out / "direct.json")
    assert record["stop"] == "deadline" and record["startup_s"] >= 201
    assert record["events"] == [] and record["baseline"] is None
    assert d.audit_arm(record)["valid"]


def test_supervisor_runs_second_fresh_child_after_first_failure(tmp_path, monkeypatch):
    d, study, root = staged(tmp_path)
    bootstrap = """
import runpy, sys
from types import SimpleNamespace
from verification import nf_base_budget as f
def prep(root, source, ref, level):
    if sys.argv[sys.argv.index('--arm')+1] == 'direct':
        raise ValueError('synthetic supplied baseline failure')
    assert str(root) == sys.argv[sys.argv.index('--input-root')+1]
    return SimpleNamespace(baseline={'eligible':True,'peak':10.,'A2':2.}), 0.
f.prep=prep
f.evaluator=lambda *args: object()
f.measurement=lambda row, ref: row
f.propose=lambda *args: None
runpy.run_module('verification.nf_unused_live', run_name='__main__')
"""
    real_run = subprocess.run
    caches = []
    order = []

    def dispatch(args, **kwargs):
        if args[0] == "nvidia-smi":
            return subprocess.CompletedProcess(args, 0, stdout="0, synthetic GPU, driver\n")
        assert args[1:3] == ["-m", "verification.nf_unused_live"]
        assert Path(kwargs["cwd"]) == d.CODE_ROOT
        order.append(args[args.index("--arm") + 1])
        env = kwargs["env"]
        caches.append(env["CUPY_CACHE_DIR"])
        assert env["CUPY_CACHE_DIR"] != env["CUDA_CACHE_PATH"]
        assert env["PYTHONDONTWRITEBYTECODE"] == "1"
        assert Path(env["TEMP"]).is_relative_to(study)
        return real_run([args[0], "-c", bootstrap, *args[3:]], **kwargs)

    monkeypatch.setattr(subprocess, "run", dispatch)
    case = d.load_cases(study)[0]
    d.prepare_campaign(study)
    out = d.run_pair(study, case, 20611, 0)
    assert order == ["direct", "screen"] and len(set(caches)) == 2
    direct, screen = (d.read(out / (arm + ".json")) for arm in order)
    assert direct["stop"] == "failed" and screen["stop"] == "draw_limit"
    assert direct["pid"] != screen["pid"]
    assert screen["completed_draws"] == 600
    processes = [d.read(out / (arm + "_process.json")) for arm in order]
    assert [p["returncode"] for p in processes] == [1, 0]
    assert processes[0]["exit_observed_ns"] <= processes[1]["launch_ns"]
    assert not list(out.glob("*-cache-*"))
    pair = d.read(out / "comparison.json")
    assert pair["status"] == "arm_failed" and pair["shared_prefix"]["length"] == 0
    with pytest.raises(FileExistsError):
        d.run_pair(study, case, 20611, 0)
    summary = d.audit_campaign(study)
    assert len(summary["pairs"]) == 3
    assert summary["pairs"][0]["status"] == "arm_failed"
    assert summary["pairs"][1]["status"] == "not_run"
    changed = d.read(out / "screen_process.json")
    changed["returncode"] = 1
    d.save(out / "screen_process.json", changed)
    assert d.audit_pair(out)["status"] == "audit_failed"
    (out / "screen.json").write_text('{"incomplete":', encoding="utf8")
    malformed = d.audit_pair(out)
    assert malformed["status"] == "audit_failed" and "screen" in malformed["errors"]


def test_order_is_balanced_by_case_seed_index():
    d = driver()
    assert d.arm_order(0, 20611) == ("direct", "screen")
    assert d.arm_order(0, 20627) == ("screen", "direct")
    assert d.arm_order(1, 20611) == ("screen", "direct")


def test_audit_rejects_late_value_in_incumbent(tmp_path, monkeypatch):
    d, study, root = staged(tmp_path)
    monkeypatch.setattr(d.live, "SEED", d.live.SEED)
    monkeypatch.setattr(d.live.frozen, "ROOT", d.live.frozen.ROOT)
    out = study / "audit-worker"
    out.mkdir()
    d.run_worker(out, root, 20611, "direct", perf_counter_ns() - 201_000_000_000)
    record = d.read(out / "direct.json")
    record["incumbent"] = dict(kind="original_wall", peak=10., A2=2.)
    with pytest.raises(ValueError, match="incumbent"):
        d.audit_arm(record)


def test_changed_campaign_input_is_not_silently_resumed(tmp_path):
    d, study, root = staged(tmp_path)
    d.prepare_campaign(study)
    d.prepare_campaign(study)
    (study / "protocol.ko.md").write_text("Changed rules\n", encoding="utf8")
    with pytest.raises(ValueError, match="protocol"):
        d.prepare_campaign(study)


def test_audit_reconstructs_archived_events_and_excludes_late_measurements():
    d = driver()
    folder = d.CODE_ROOT / "results/d3d-base-budget-20261105/live"
    for policy in ("direct", "screen"):
        record = d.read(folder / f"{policy}.json")
        assert d.audit_arm(record)["valid"]
    direct = d.read(folder / "direct.json")
    direct["events"][-1]["measurement"] = {"eligible": True, "peak": .01, "A2": 3.}
    with pytest.raises(ValueError, match="late measurement"):
        d.audit_arm(direct)


def test_consistent_replacement_shot_and_manifest_still_break_preflight(tmp_path):
    d, study, root = staged(tmp_path)
    (root / d.SHOT).write_bytes(b"replacement source")
    manifest = d.read(root / d.MANIFEST)
    manifest["shot_sha256"] = d.digest(root / d.SHOT)
    d.save(root / d.MANIFEST, manifest)
    with pytest.raises(ValueError, match="preflight input hashes"):
        d.load_cases(study)


def test_frozen_preflight_metadata_is_verified(tmp_path):
    d, study, root = staged(tmp_path)
    preflight = d.read(study / "preflight.json")
    preflight["selection_sha256"] = "wrong"
    d.save(study / "preflight.json", preflight)
    with pytest.raises(ValueError, match="preflight.*selection"):
        d.prepare_campaign(study)
