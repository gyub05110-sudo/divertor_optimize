"""Exercise the real policy loop without performing any thermal evaluations."""
from collections import Counter
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

from verification import nf_base_budget as frozen
from verification import nf_base_live as live


class Clock:
    ns = 0

    def add(self, seconds):
        self.ns += round(seconds * 1e9)

    def __call__(self):
        return self.ns


class Geometry:
    metadata = {"source": "OMP"}
    rmin = .020

    def constraints(self, co):
        return {"passed": True}

    def precise_lengths(self, co):
        return {"floor": 1.15}, 0.


@pytest.fixture
def experiment(monkeypatch):
    clock = Clock()
    submitted = []
    draws = []

    def prep(root, source, ref, level):
        clock.add({0: 2, 2: 3}[level])
        return SimpleNamespace(level=level, baseline={"peak": 10., "A2": 2., "eligible": True}), 0.

    def evaluator(original, radius, backend):
        assert radius == 20 and backend == "gpu"
        clock.add(1)
        geom = Geometry()
        geom.level = original.level
        return geom

    def propose(geom, rng):
        assert geom.level == 2
        clock.add(1)
        draws.append(1)
        return np.full(16, .004 * len(draws))

    def score(geom, candidate):
        i = round(candidate.vector[0] / .004) - 1
        submitted.append((geom.level, i))
        if geom.level == 2:
            clock.add(5)
            peak = [11., 10.5, 9., 8.][i % 4]
        else:
            clock.add(2)
            peak = [1., .5, .5, .1][i % 4]
        return dict(eligible=geom.level == 2 or i % 4 != 3, peak=peak, A2=3.)

    monkeypatch.setattr(frozen, "prep", prep)
    monkeypatch.setattr(frozen, "evaluator", evaluator)
    monkeypatch.setattr(frozen, "measurement", lambda row, ref: row)
    monkeypatch.setattr(frozen, "propose", propose)
    monkeypatch.setattr(frozen, "geometry_check", lambda e, co: True)
    monkeypatch.setattr(frozen, "score", score)
    return SimpleNamespace(run=lambda policy, budget: live.run_arm(
        policy, 0, clock_ns=clock, budget_s=budget), clock=clock, submitted=submitted)


def test_on_time_fine_scores_and_original_fallback_are_separate(experiment):
    result = experiment.run("direct", 16)
    assert experiment.submitted == [(2, 0), (2, 1)]
    assert result["fine_ids"] == [0, 1]
    assert result["baseline"]["peak"] == result["incumbent"]["peak"] == 10.
    assert result["incumbent"]["kind"] == "original_wall"
    assert result["best_eligible_new"]["i"] == 1
    assert result["first_5pct_s"] is None
    assert result["elapsed_s"] == 16


def test_late_fine_value_is_discarded_and_overshoot_is_visible(experiment):
    result = experiment.run("direct", 15.99)
    assert result["fine_ids"] == [0]
    assert result["best_eligible_new"]["peak"] == 11.
    assert result["events"][-1]["status"] == "late_discarded"
    assert "measurement" not in result["events"][-1]
    assert result["events"][-1]["end_s"] == 16


def test_startup_and_late_cold_prep_are_charged(experiment):
    experiment.clock.add(2)
    result = experiment.run("direct", 4.5)
    assert result["events"][0]["start_s"] == 2
    assert result["events"][0]["end_s"] == 5
    assert result["baseline"] is None and result["incumbent"] is None
    assert experiment.submitted == []


def test_late_projection_never_admits_a_candidate(experiment):
    result = experiment.run("direct", 4.5)
    assert not result["generated"] and not experiment.submitted
    assert result["events"][-1]["stage"] == "draw"
    assert result["events"][-1]["status"] == "late_discarded"


def test_full_batch_promotes_eligible_minimum_with_generation_order_tie(experiment):
    result = experiment.run("screen", 24)
    assert experiment.submitted == [(0, 0), (0, 1), (0, 2), (0, 3), (2, 1)]
    assert result["coarse_ids"] == [0, 1, 2, 3]
    assert result["fine_ids"] == [1]


def test_deadline_in_partial_batch_never_promotes_early(experiment):
    result = experiment.run("screen", 17)
    assert experiment.submitted == [(0, 0), (0, 1), (0, 2)]
    assert result["coarse_ids"] == [0, 1, 2]
    assert len(result["generated"]) == 4 and result["fine_ids"] == []


def test_first_improvement_clock_and_candidate_cap(experiment):
    result = experiment.run("direct", 200)
    assert len(result["generated"]) == len(result["fine_ids"]) == 12
    assert result["first_5pct_s"] == 22
    assert result["incumbent"]["kind"] == "new_design"
    assert result["best_eligible_new"]["i"] == 3


def test_stream_matches_frozen_panel_including_rng_and_all_rejections(experiment, monkeypatch):
    class Filtered(Geometry):
        def constraints(self, co):
            return {"passed": co[0] != .004}

        def precise_lengths(self, co):
            return {"floor": 1.21 if co[0] == .008 else 1.2}, 0.

    def propose(geom, rng):
        n = rng.integers(25)
        return np.full(16, .004*n) if n else None

    monkeypatch.setattr(frozen, "propose", propose)
    monkeypatch.setattr(frozen, "geometry_check", lambda e, co: co[0] != .012)
    geom = Filtered()
    seed = 7
    expected = frozen.panel(geom, seed)
    rng = np.random.default_rng(seed)
    accepted = []
    rejected = Counter()
    for attempt in range(1, 601):
        candidate, reason = live.draw(geom, rng, accepted)
        if reason:
            rejected[reason] += 1
        else:
            accepted.append(live.candidate_row(candidate, len(accepted), attempt))
            if len(accepted) == 12:
                break
    assert {r["id"] for r in accepted} == {r["id"] for r in expected["rows"]}
    assert [(r["id"], r["cumulative_draws"]) for r in accepted] == [
        (r["id"], r["cumulative_draws"]) for r in expected["rows"]]
    assert dict(rejected) == expected["rejected"]
    assert set(rejected) == {"geometry_limits", "projection", "whole_wall", "near_duplicate"}


def test_final_short_group_only_when_all_draws_exhausted(experiment, monkeypatch):
    count = 0

    def only_one(geom, rng):
        nonlocal count
        count += 1
        return np.full(16, .004) if count == 1 else None

    monkeypatch.setattr(frozen, "propose", only_one)
    result = experiment.run("screen", 200)
    assert count == 600 and result["rejected"] == {"projection": 599}
    assert experiment.submitted == [(0, 0), (2, 0)]


def test_driver_compares_only_shared_prefix_and_detects_drift():
    a = [{"id": "A", "cumulative_draws": 2}, {"id": "B", "cumulative_draws": 4}]
    assert live.prefix_audit(a, a[:1]) == {"length": 1, "equal": True}
    assert not live.prefix_audit(a, [dict(a[0], cumulative_draws=3)])["equal"]
    assert not live.prefix_audit(a, [dict(a[0], id="C")])["equal"]


def test_supervisor_launches_two_fresh_workers_and_audits_results(tmp_path, monkeypatch):
    # Use the real child entry point and file interchange with cheap synthetic
    # setup and an exhausted stream. No physical evaluator can run in this test.
    bootstrap = """
import runpy
from types import SimpleNamespace
from verification import nf_base_budget as f
f.prep = lambda *args: (SimpleNamespace(baseline={"eligible":True,"peak":10.,"A2":2.}), 0.)
f.evaluator = lambda *args: object()
f.measurement = lambda row, ref: row
f.propose = lambda *args: None
runpy.run_module("verification.nf_base_live", run_name="__main__")
"""
    launch = subprocess.run
    caches = []
    order = []

    def dispatch(args, **kwargs):
        if args[0] == "nvidia-smi":
            assert "uuid" not in args[1]
            return subprocess.CompletedProcess(args, 0, stdout="0, test GPU, test driver")
        assert args[1:3] == ["-m", "verification.nf_base_live"]
        order.append(args[args.index("--arm") + 1])
        cache = kwargs["env"]["CUPY_CACHE_DIR"]
        assert cache != kwargs["env"]["CUDA_CACHE_PATH"]
        caches.append(cache)
        return launch([args[0], "-c", bootstrap, *args[3:]], **kwargs)

    real_digest = frozen.digest
    monkeypatch.setattr(frozen, "digest", lambda path: real_digest(path) if str(path).endswith(".py") else "test-digest")
    monkeypatch.setattr(subprocess, "run", dispatch)
    out = tmp_path / "pair"
    live.run_pair(out)
    assert order == ["direct", "screen"] and caches[0] != caches[1]
    assert not list(out.glob("*-cache-*"))
    direct, screen = (frozen.load_json(out/(p+".json")) for p in order)
    assert direct["pid"] != screen["pid"]
    assert direct["completed_draws"] == screen["completed_draws"] == 600
    assert direct["stop"] == screen["stop"] == "draw_limit"
    assert direct["startup_s"] > 0 and screen["startup_s"] > 0
    first, second = (frozen.load_json(out/(p+"_process.json")) for p in order)
    assert first["exit_observed_ns"] <= second["launch_ns"]
    assert first["elapsed_s"] >= direct["elapsed_s"] and second["elapsed_s"] >= screen["elapsed_s"]
    assert frozen.load_json(out/"comparison.json")["shared_prefix"] == {"length": 0, "equal": True}
    protocol = frozen.load_json(out/"protocol.json")
    assert protocol["gpu_inventory"] == ["0, test GPU, test driver"]
    assert "host" not in protocol and "python" not in protocol
    assert "DIVERTOR_GPU_RUNTIME" not in protocol["environment"]
    with pytest.raises(FileExistsError):
        live.run_pair(out)
