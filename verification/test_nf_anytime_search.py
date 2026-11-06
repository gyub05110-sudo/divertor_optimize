"""Deadline and evidence tests; no physical evaluator or GPU may run here."""
import importlib
import importlib.util
import json
from collections import OrderedDict
from types import SimpleNamespace

import numpy as np
import pytest

from d3d_optimization.campaign_geometry import Candidate


def runner():
    assert importlib.util.find_spec("verification.nf_anytime_search") is not None, (
        "the additive anytime runner has not been implemented"
    )
    return importlib.import_module("verification.nf_anytime_search")


class Clock:
    def __init__(self):
        self.ns = 0

    def add(self, seconds):
        self.ns += round(seconds * 1e9)

    def __call__(self):
        return self.ns


def measured(peak=10., eligible=True):
    return dict(eligible=eligible, errors=[] if eligible else ["invalid"],
                ledger={"deposited": 1.}, peak=peak, A2=2.)


def candidate(i):
    return Candidate("OMP", 20, tuple(np.full(16, .001 * (i + 1))))


def event(i, end, peak, stage="fine", status="completed", eligible=True, cached=False):
    c = candidate(i)
    return dict(stage=stage, start_s=end-1, end_s=end, status=status,
                id=c.identity, candidate=c.record(), cached=cached,
                measurement=measured(peak, eligible))


def record(events=(), ready=3.):
    return dict(budget_s=600., baseline=measured(), baseline_ready_s=ready,
                events=list(events), candidates=[])


def test_snapshot_exact_boundary_and_late_result_is_never_admitted():
    live = runner()
    r = record([event(0, 200., 9.), event(1, 600.01, 1., status="late")])
    assert live.snapshot(r, 199.99)["incumbent"]["kind"] == "original_wall"
    at = live.snapshot(r, 200.)
    assert at["incumbent"]["id"] == candidate(0).identity
    assert at["first_5pct_s"] == 200.
    assert live.snapshot(r, 1000.)["incumbent"]["peak"] == 9.


def test_snapshot_requires_baseline_setup_before_requested_time():
    live = runner()
    assert not live.snapshot(record(ready=200.01), 200.)["available"]
    assert not live.snapshot(record(ready=600.01), 1000.)["available"]


def test_failed_baseline_measurement_is_not_a_fallback():
    live = runner()
    r = record(ready=None)
    r["baseline"] = None
    r["events"] = [dict(stage="baseline", start_s=1., end_s=2., status="failed",
                        measurement=measured())]
    assert not live.snapshot(r, 200.)["available"]


def test_coarse_only_and_worse_fine_scores_cannot_displace_baseline():
    live = runner()
    result = live.snapshot(record([event(0, 4., .1, "coarse"), event(1, 5., 11.)]), 200.)
    assert result["incumbent"]["kind"] == "original_wall"
    assert result["first_5pct_s"] is None
    assert result["counts"]["coarse"] == result["counts"]["fine"] == 1


def test_snapshot_recomputes_shape_identity_and_ignores_optimizer_state():
    live = runner()
    row = event(0, 4., 8.)
    row["id"] = "stale-untrusted-id"
    r = record([row])
    r["incumbent"] = dict(peak=.001, id="not-evidence")
    assert live.snapshot(r, 200.)["incumbent"]["id"] == candidate(0).identity


@pytest.fixture
def experiment(monkeypatch, tmp_path):
    live = runner()
    clock = Clock()
    options = dict(score_s={0: 2., 2: 5.}, prep_s={0: 2., 2: 3.},
                   peak=lambda level, i: 9.-.1*i, repeat=False, fail=None,
                   housekeeping_s=0., eligible=lambda level, i: True,
                   draw_rejection=None, restart_after_first=False)
    calls, feedback, initializations = [], [], []

    class Geometry:
        metadata = {"source": "OMP"}
        rmin = .02
        baseline_peak = 10.
        lengths = dict(omp_entry=1., shelf=1., omp_connector=1., floor=1.)

        def __init__(self, level, reference):
            self.level, self.reference = level, reference

        def describe(self, raw):
            return dict(feasible=raw["eligible"], errors=raw["errors"],
                        peak_MW_m2=raw["peak"], ledger=raw["ledger"])

    def prep(root, source, ref, level):
        clock.add(options["prep_s"][level])
        return SimpleNamespace(level=level, reference=ref, baseline=measured(),
                               metadata={"source": source}), 0.

    def evaluator(original, radius, backend):
        clock.add(1.)
        return Geometry(original.level, original.reference)

    def scoring(geom, c):
        i = round(c.vector[0] / .001)-1
        calls.append((geom.level, i))
        clock.add(options["score_s"][geom.level])
        if options["fail"] == (geom.level, i):
            raise ValueError("synthetic scorer failure")
        return measured(options["peak"](geom.level, i), options["eligible"](geom.level, i))

    def init(self, original, mm, folder, root, config):
        self.original, self.mm, self.folder, self.root, self.config = original, mm, folder, root, config
        self.source = "OMP"
        self.e = evaluator(original, mm, "gpu")
        self.cpu = Geometry(original.level, original.reference)
        self.memo = OrderedDict()
        self.state = dict(optimizer=None, total_proposals=0, restart=0, archive=[],
                          best=None, best_infeasible=None, last_improvement=0,
                          active_seconds=0., parity=None)

    def initialize(self, inherited=()):
        initializations.append((self.config.copy(), tuple(inherited)))
        clock.add(2.)
        self.state["optimizer"] = dict(generation=0)

    def score(self, x, stage):
        c = Candidate(self.source, self.mm, tuple(map(float, x)))
        if c.identity in self.memo:
            clock.add(.1)
            return dict(self.memo[c.identity], cached=True)
        raw = scoring(self.e, c)
        row = self.e.describe(raw)
        row.update(candidate=c.record(), identity=c.identity, cached=False,
                   physics_evaluated=True, objective=raw["peak"]/10.,
                   stage=stage, evaluation_seconds=options["score_s"][self.e.level])
        self.memo[c.identity] = row
        return row

    def step(self):
        clock.add(1.)
        n = self.state["total_proposals"]
        i = n//2 if options["repeat"] else n
        row = self.score(candidate(i).vector, "differential_evolution")
        feedback.append(row["identity"])
        self.state["total_proposals"] += 1
        if options["restart_after_first"] and n == 0:
            self.state["optimizer"]["generation"] = 100
        clock.add(options["housekeeping_s"])

    draws = []

    def draw(geom, rng, accepted):
        clock.add(1.)
        i = len(draws)
        draws.append(i)
        reason = options["draw_rejection"]
        return (None, reason) if reason else (candidate(i), None)

    monkeypatch.setattr(live, "prep", prep)
    monkeypatch.setattr(live, "evaluator", evaluator)
    monkeypatch.setattr(live, "measurement", lambda raw, ref: raw)
    monkeypatch.setattr(live, "score", scoring)
    monkeypatch.setattr(live, "draw", draw)
    monkeypatch.setattr(live.SearchCase, "__init__", init)
    monkeypatch.setattr(live.SearchCase, "initialize", initialize)
    monkeypatch.setattr(live.SearchCase, "score", score)
    monkeypatch.setattr(live.SearchCase, "step", step)
    monkeypatch.setattr(live.SearchCase, "checkpoint", lambda *a, **k: None)
    monkeypatch.setattr(live, "save", lambda *a, **k: None)

    def run(policy, budget=40.):
        return live.run_arm(policy, tmp_path/"input", 20611, tmp_path/policy,
                            0, budget_s=budget, clock_ns=clock)

    return SimpleNamespace(run=run, clock=clock, calls=calls, options=options,
                           feedback=feedback, initializations=initializations, live=live,
                           out=tmp_path)


def test_late_score_retained_with_full_measurement_and_no_feedback(experiment):
    r = experiment.run("fine_adaptive", 11.)
    scores = [e for e in r["events"] if e["stage"] == "fine"]
    assert scores[-1]["status"] == "late"
    assert set(scores[-1]["measurement"]) >= {"eligible", "errors", "ledger", "peak", "A2"}
    assert scores[-1]["candidate"] == candidate(0).record()
    assert not experiment.feedback
    assert r["stop"] == "deadline"
    assert experiment.live.snapshot(r, 11.)["incumbent"]["kind"] == "original_wall"


def test_score_at_deadline_is_admitted(experiment):
    r = experiment.run("fine_adaptive", 12.)
    assert experiment.live.snapshot(r, 12.)["incumbent"]["peak"] == 9.
    assert r["events"][-1]["status"] == "completed"


def test_timely_score_survives_late_step_housekeeping(experiment):
    experiment.options["housekeeping_s"] = 20.
    r = experiment.run("fine_adaptive", 15.)
    assert experiment.calls == [(2, 0)]
    assert experiment.live.snapshot(r, 15.)["incumbent"]["peak"] == 9.
    assert r["elapsed_s"] > 15.


def test_cold_start_and_late_setup_count_against_budget(experiment):
    experiment.clock.add(2.)
    r = experiment.run("random_fine", 4.)
    assert r["startup_s"] == 2.
    assert r["events"][0]["start_s"] == 2.
    assert not experiment.live.snapshot(r, 4.)["available"]
    assert not experiment.calls


def test_startup_exhaustion_writes_empty_raw_stream(experiment):
    experiment.clock.add(20.)
    r = experiment.run("fine_adaptive", 10.)
    assert r["stop"] == "deadline" and r["events"] == []
    path = experiment.out/"fine_adaptive/events.jsonl"
    assert path.exists() and path.read_text(encoding="utf-8") == ""


def test_historical_warm_shapes_are_rejected_without_reading_them(experiment):
    historical = experiment.out/"input/campaign_warm_candidates.json"
    historical.parent.mkdir()
    historical.write_text("preserve this historical input", encoding="utf-8")
    r = experiment.run("fine_adaptive", 40.)
    assert r["stop"] == "failed" and "historical shape" in r["error"]
    assert historical.read_text(encoding="utf-8") == "preserve this historical input"
    assert not experiment.calls and not experiment.initializations


def test_coarse_promotes_after_four_unique_eligible_results(experiment):
    experiment.options["repeat"] = True
    r = experiment.run("coarse_adaptive", 40.)
    assert experiment.calls[:5] == [(0, 0), (0, 1), (0, 2), (0, 3), (2, 3)]
    fine = [e for e in r["events"] if e["stage"] == "fine"]
    assert fine[0]["id"] == candidate(3).identity
    coarse_before = [e for e in r["events"] if e["stage"] == "coarse" and e["end_s"] <= fine[0]["start_s"]]
    assert len(coarse_before) == 7
    assert sum(not e["cached"] for e in coarse_before) == 4


def test_ineligible_coarse_scores_do_not_fill_group(experiment):
    experiment.options["eligible"] = lambda level, i: level == 2 or i != 1
    r = experiment.run("coarse_adaptive", 23.)
    assert not [e for e in r["events"] if e["stage"] == "fine"]


def test_cutoff_in_partial_batch_never_flushes_a_promotion(experiment):
    r = experiment.run("coarse_adaptive", 17.)
    assert all(level == 0 for level, _ in experiment.calls)
    assert not [e for e in r["events"] if e["stage"] == "fine"]


def test_failed_scorer_preserves_prior_evidence_and_explicit_error(experiment):
    experiment.options["fail"] = (2, 1)
    r = experiment.run("fine_adaptive", 40.)
    assert r["stop"] == "failed" and "synthetic scorer failure" in r["error"]
    assert r["events"][-1]["status"] == "failed"
    assert len(r["candidates"]) == 2
    assert experiment.live.snapshot(r, 40.)["incumbent"]["peak"] == 9.


@pytest.mark.parametrize("policy", ["fine_adaptive", "random_fine"])
def test_search_is_not_capped_at_twelve_candidates(experiment, policy):
    r = experiment.run(policy, 90.)
    assert len(r["candidates"]) > 12
    assert r["stop"] == "deadline"


def test_random_rejections_are_preserved_and_charged(experiment):
    experiment.options["draw_rejection"] = "whole_wall"
    r = experiment.run("random_fine", 7.)
    assert r["rejected"]["whole_wall"] == 3
    assert not r["candidates"] and not experiment.calls


def test_adaptive_initialization_preserves_twelve_projected_population(experiment):
    experiment.run("fine_adaptive", 8.)
    config, inherited = experiment.initializations[0]
    assert config["population"] == config["projected_seeds"] == 12
    assert config["seed"] == 20611
    assert inherited == ()


@pytest.mark.parametrize("failure", [False, True])
def test_finalized_raw_events_match_returned_evidence(experiment, failure):
    if failure:
        experiment.options["fail"] = (2, 0)
    r = experiment.run("fine_adaptive", 11.)
    path = experiment.out/"fine_adaptive/events.jsonl"
    assert path.exists(), "every completed, late or failed event must be appended"
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    assert rows == r["events"]
    assert all(set(row) >= {"stage", "start_s", "end_s", "status"} for row in rows)
    assert set(r["parity"]) == {"fine", "coarse"}


def test_native_restart_policy_is_kept(experiment):
    experiment.options["restart_after_first"] = True
    experiment.run("fine_adaptive", 20.)
    assert len(experiment.initializations) == 2


def test_deadline_during_proposal_keeps_shape_without_starting_score(experiment):
    r = experiment.run("fine_adaptive", 6.5)
    assert not experiment.calls
    assert r["events"][-1]["stage"] == "proposal"
    assert r["events"][-1]["status"] == "late"
    assert len(r["candidates"]) == 1


@pytest.mark.parametrize("budget", [6., 9.])
def test_real_searchcase_uses_fine_feedback_and_full_parity_without_physics(monkeypatch, tmp_path, budget):
    live = runner()
    from d3d_optimization import campaign_search as search
    clock, calls, feedback = Clock(), [], []
    ref = live.ReferencePoint()
    raw = dict(feasible=True, errors=[], peak_q=9., geometry_valid=True,
               elapsed_seconds=2., deposition=dict(ledger={"deposited": 1.},
               points=np.array([[1., 2.]]), q=np.array([9.]), cell_power=np.array([1.]),
               areas=np.array([1.]), quadrature_integrated_power=1.))

    class Geometry:
        metadata = {"source": "OMP"}
        rmin, baseline_peak = .02, 10.
        lengths = dict(omp_entry=1., shelf=1., omp_connector=1., floor=1.)
        reference = ref

        def constraints(self, co):
            return {"passed": True}

        def describe(self, result):
            return dict(feasible=result["feasible"], errors=result["errors"])

    def prep(*args):
        clock.add(3.)
        return SimpleNamespace(metadata={"source": "OMP"}, baseline=raw), 0.

    def evaluator(original, mm, backend):
        clock.add(1.)
        e = Geometry()
        e.backend = backend
        return e

    def evaluate(e, c):
        e.constraints(c.coefficients)
        calls.append((e.backend, c.identity))
        clock.add(2. if e.backend == "gpu" else 1.)
        return raw

    original_accept = search.accept

    def accept(opt, x, objective):
        feedback.append(objective)
        return original_accept(opt, x, objective)

    monkeypatch.setattr(live, "prep", prep)
    monkeypatch.setattr(search, "evaluator", evaluator)
    monkeypatch.setattr(search, "evaluate_candidate", evaluate)
    monkeypatch.setattr(search, "project_floor", lambda e, x: np.asarray(x))
    monkeypatch.setattr(search, "accept", accept)
    r = live.run_arm("fine_adaptive", tmp_path/"input", 20611, tmp_path/"run", 0,
                     budget_s=budget, clock_ns=clock)
    assert r["stop"] == "deadline", r.get("error")
    assert [backend for backend, _ in calls] == (["gpu"] if budget == 6. else ["gpu", "cpu", "gpu"])
    assert feedback == ([] if budget == 6. else [.9])
    if budget == 9.:
        assert r["parity"]["fine"]["passed"]
        assert r["parity"]["fine"]["density_scaled_max_difference"] == 0.
    else:
        assert r["parity"]["fine"] is None
    assert r["events"][-1]["measurement"]["A2"] == 1.


def test_projected_seeds_obey_injected_monotonic_deadline(experiment, monkeypatch):
    from d3d_optimization import campaign_search as search
    projections = []

    def project(e, target):
        projections.append(target.tolist())
        experiment.clock.add(2.)
        return np.zeros(4)

    def initialize(self, inherited=()):
        self.state["optimizer"] = search.seed_population(self.e, self.root, self.config, 20611)

    monkeypatch.setattr(search, "project_floor", project)
    monkeypatch.setattr(experiment.live.SearchCase, "initialize", initialize)
    r = experiment.run("fine_adaptive", 8.)
    assert len(projections) == 2
    assert not experiment.calls and r["stop"] == "deadline"


def test_event_write_crossing_deadline_keeps_completed_score(experiment, monkeypatch):
    original = experiment.live._Arm.append_event

    def append(self, row):
        original(self, row)
        if row["stage"] == "fine":
            experiment.clock.add(20.)

    monkeypatch.setattr(experiment.live._Arm, "append_event", append)
    r = experiment.run("fine_adaptive", 15.)
    assert experiment.calls == [(2, 0)]
    assert experiment.live.snapshot(r, 15.)["incumbent"]["peak"] == 9.
    assert r["events"][-1]["end_s"] == 12.
    assert r["elapsed_s"] == 32.
