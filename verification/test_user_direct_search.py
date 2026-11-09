"""Direct-only search uses native DE and the supplied-input admission clock."""
import importlib.util
import json
from math import inf
from types import SimpleNamespace

import numpy as np
import pytest

import optimize_divertor as cli
from d3d_optimization import campaign_search
from verification.test_user_search import Clock, result


def direct():
    name = "d3d_optimization.user_direct_search"
    assert importlib.util.find_spec(name), "missing direct-fine comparator"
    return __import__(name, fromlist=["search_direct"])


def geometry():
    return SimpleNamespace(metadata={"source": "OMP"}, bound=.1, rmin=.02, gap_min=.04,
        baseline_peak=10., lengths=dict(omp_entry=.18, shelf=.1, omp_connector=.11, floor=.55))


def test_population_reproduces_frozen_cold_32_12_and_spans_every_axis(tmp_path, monkeypatch):
    module = direct()
    fine = geometry()
    project = lambda _, target: target
    monkeypatch.setattr(module, "project_floor", project)
    monkeypatch.setattr(campaign_search, "project_floor", project)
    seed = 20711 + 20 + 10000
    prior = campaign_search.seed_population(fine, tmp_path, dict(population=32, projected_seeds=12,
                search_deadline_epoch=inf), seed)
    actual = module._initialize(fine, seed, [], [], inf)
    assert actual == prior
    assert np.linalg.matrix_rank(np.asarray(actual["population"])-np.mean(actual["population"], axis=0)) == 16


def simulate(monkeypatch, *, budget, fine_score, baseline=None, project=None):
    module = direct()
    clock, fine, events, candidates = Clock(5.), geometry(), [], []
    monkeypatch.setattr(module, "perf_counter", clock)
    monkeypatch.setattr(module, "project_floor", project or (lambda _, target: target))
    def evaluate(e, candidate):
        assert e is fine
        candidates.append(candidate)
        clock.advance(1.)
        raw = fine_score(len(candidates)-1)
        if isinstance(raw, Exception):
            raise raw
        return raw
    def observe(event):
        json.dumps(event, allow_nan=False)
        events.append(event)
    outcome = module.search_direct(fine, radius_mm=20, seed=20711, started=100.,
        deadline=100.+budget, baseline=baseline, observe=observe, evaluate=evaluate)
    return outcome, events, candidates


def test_direct_generation_progresses_beyond_full_initial_population_with_native_accept(monkeypatch):
    outcome, events, candidates = simulate(monkeypatch, budget=50.,
        fine_score=lambda i: result(5.-i*.01, deposition={}))
    init, = [e for e in events if e["stage"] == "initialize"]
    assert init["rank"] == 16 and init["population"] == 32 and init["projected_seeds"] == 12
    assert not any(e["stage"] == "coarse" for e in events)
    fine = [e for e in events if e["stage"] == "fine"]
    assert outcome["counts"]["coarse"] == 0 and len(fine) == outcome["counts"]["fine"]
    trials = [e for e in fine if e["proposal_kind"] == "de_trial"]
    assert trials and all(e["generation"] > 0 and e["physics_evaluated"] for e in trials)
    assert all(e["coarse_peak_q"] is None for e in fine)
    assert outcome["winner"] == candidates[-1]
    assert outcome["fine_finalists"][0]["id"] == outcome["winner"].identity
    assert outcome["first_feasible_elapsed_s"] == 6.


def test_proposal_vectors_replay_archived_searchcase_step(tmp_path, monkeypatch):
    """Replay measured energies through the original optimizer, including interleaved polish."""
    from d3d_optimization.campaign_geometry import Candidate
    fine = geometry()
    _, events, _ = simulate(monkeypatch, budget=90.,
        fine_score=lambda i: result(5. + (i*7 % 11)/10., deposition={}))
    proposals = [e for e in events if e["stage"] == "proposal"]
    scores = {e["id"]: e["objective"] for e in events if e["stage"] == "fine"}
    monkeypatch.setattr(campaign_search, "project_floor", lambda _, target: target)
    opt = campaign_search.seed_population(fine, tmp_path, dict(population=32, projected_seeds=12,
        search_deadline_epoch=inf), 20711+20+10000)
    class Replay:
        folder, source, e = tmp_path, "OMP", fine
        state = dict(optimizer=opt, best=None, best_infeasible=None, archive=[], restart=0,
                     physics_evaluations=0, feasible_evaluations=0, total_proposals=0)
        def checkpoint(self, force=False):
            pass
        def score(self, x, stage):
            i = self.state["total_proposals"]
            c = Candidate("OMP", 20, tuple(map(float, x)))
            assert c.identity == proposals[i]["id"]
            assert (stage == "coordinate_polish") == (proposals[i]["proposal_kind"] == "coordinate_polish")
            return dict(identity=c.identity, candidate=c.record(), objective=scores[c.identity], feasible=True,
                        cached=proposals[i]["status"]=="rejected", stage=stage, evaluation_seconds=0.,
                        physics_evaluated=True)
    replay = Replay()
    for _ in proposals:
        campaign_search.SearchCase.step(replay)
    assert replay.state["optimizer"]["generation"] >= 2


def test_original_reference_ineligible_late_fine_cannot_win(monkeypatch):
    base = result(.001, reference_only=True)
    out, events, candidates = simulate(monkeypatch, budget=8.25, baseline=base,
        fine_score=lambda i: [result(6.), result(.01, feasible=False, geometry_valid=True, penalty=1.),
                              result(7.), result(.001)][i])
    assert out["winner"] == candidates[0] and out["selected"]["peak_q"] == 6.
    assert out["counts"]["eligible_fine"] == 2 and out["counts"]["late"] == 1
    assert out["inflight_overrun_s"] == .75
    assert events[-1]["stage"] == "fine" and events[-1]["status"] == "late"
    assert not events[-1]["incumbent_updated"]


def test_full_direct_cli_skips_coarse_preparation_and_uses_same_fine_reference(tmp_path, monkeypatch):
    module = direct()
    from d3d_optimization import campaign_geometry, user_input, user_physics, user_output, user_precheck
    path, out = tmp_path / "shot.parquet", tmp_path / "out"
    path.touch()
    user = SimpleNamespace(provenance={}, review={k: {} for k in ("axis_RZ_m", "xpoint_RZ_m", "legs")})
    base = result(5., reference_only=True, physics_measurable=True, deposition={})
    original = SimpleNamespace(baseline=base)
    fine = geometry()
    levels = []
    def prepare(_user, reference, level):
        levels.append(level)
        assert level == cli.FINE
        return original
    def search(owner, **kwargs):
        assert owner is fine and kwargs["baseline"] is base and kwargs["seed"] == 20711
        return dict(selected=None, winner=None, status="no_design", termination="deadline")
    monkeypatch.setattr(user_input, "load_user_input", lambda *a: user)
    monkeypatch.setattr(user_precheck, "fixed_wall_precheck", lambda _: {"excluded": False})
    monkeypatch.setattr(user_physics, "shelf_limits", lambda *a: {"admissible": True})
    monkeypatch.setattr(user_physics, "prepare_user", prepare)
    monkeypatch.setattr(campaign_geometry, "evaluator", lambda *a: fine)
    monkeypatch.setattr(user_output, "export_measurement", lambda *a: {"peak_MW_m2": 5.})
    monkeypatch.setattr(module, "search_direct", search)
    assert cli.main(["--equilibrium", str(path), "--output", str(out), "--time-ms", "1400",
                     "--method", "direct", "--seed", "20711"]) == 0
    report = json.loads((out/"result.json").read_text())
    assert report["status"] == "no_design" and report["arguments"]["method"] == "direct"
    assert report["fixed_wall"] == {"excluded": False} and report["surface_policy"] is None
    assert report["levels"] == {"fine": list(cli.FINE)} and levels == [cli.FINE]
    assert "prepare_coarse_search" not in [e["stage"] for e in report["timings"]["preparation"]]


def test_cli_screened_default_and_direct_full_space(tmp_path):
    args = ["--equilibrium", "a.parquet", "--output", str(tmp_path), "--time-ms", "1400"]
    assert getattr(cli.arguments(args), "method", None) == "screened"
    with pytest.raises(SystemExit):
        cli.arguments(args + ["--method", "direct", "--family", "floor-only"])
