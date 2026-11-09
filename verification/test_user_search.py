"""Deadline and objective tests using a monotonic process clock and native result shapes."""
import importlib
import importlib.util
import json
from types import SimpleNamespace

import numpy as np
import pytest


class Clock:
    def __init__(self, elapsed):
        self.now = 100. + elapsed

    def __call__(self):
        return self.now

    def advance(self, seconds):
        self.now += seconds


def result(peak, feasible=True, **extra):
    return dict(peak_q=peak, feasible=feasible, errors=[] if feasible else ["exposure"],
                profile=np.array([peak]), **extra)


@pytest.fixture
def run(monkeypatch):
    assert importlib.util.find_spec("d3d_optimization.user_search") is not None, "user search is missing"
    module = importlib.import_module("d3d_optimization.user_search")

    def simulate(*, budget=45., setup=5., coarse=None, fine=None, family="full", seed=3,
                 baseline=None, project=None):
        clock = Clock(setup)
        events, targets = [], []
        calls = dict(coarse=[], fine=[])
        coarse = coarse or (lambda i, c: (1., result(10.-i*.001)))
        fine = fine or (lambda i, c: (2., result(8.-i*.001)))

        def geometry(level):
            return SimpleNamespace(level=level, metadata={"source": "OMP"}, bound=.1, rmin=.02,
                reference=SimpleNamespace(power_W=2e6),
                lengths=dict(omp_entry=.18, shelf=.1, omp_connector=.11, floor=.55))

        def evaluate(e, candidate):
            bucket = calls[e.level]
            index = len(bucket)
            bucket.append(candidate)
            elapsed, raw = (coarse if e.level == "coarse" else fine)(index, candidate)
            clock.advance(elapsed)
            if isinstance(raw, Exception):
                raise raw
            return raw

        def project_floor(e, target):
            targets.append(np.array(target))
            elapsed, floor = project(len(targets)-1, target) if project else (0., target)
            clock.advance(elapsed)
            return floor

        def observe(event):
            json.dumps(event, allow_nan=False)
            events.append(event)

        monkeypatch.setattr(module, "project_floor", project_floor)
        monkeypatch.setattr(module, "perf_counter", clock)
        output = module.search_peak(geometry("coarse"), geometry("fine"), radius_mm=20, seed=seed,
            started=100., deadline=100.+budget, baseline=baseline, observe=observe, family=family,
            evaluate=evaluate)
        return SimpleNamespace(output=output, events=events, calls=calls, targets=targets)
    return simulate


def test_winner_is_lowest_exact_eligible_on_time_fine_peak_even_if_reference_lower(run):
    lower = np.nextafter(7., 0.)
    levels = [result(7., A2=99.), result(lower, A2=.01), result(lower, A2=500.),
              result(.01, feasible=False), result(float("nan"))]
    sim = run(baseline=result(1., feasible=False, reference_only=True, physics_measurable=True),
              fine=lambda i, c: (2., levels[i] if i < len(levels) else result(9.)))
    out = sim.output
    assert out["status"] == "qualified"
    assert out["winner"] == sim.calls["fine"][1]
    assert out["selected"] is levels[1]
    admitted = [e for e in sim.events if e["stage"] == "fine" and e["within_budget"] and e["eligible"]]
    assert len(admitted) == out["counts"]["eligible_fine"]
    assert out["first_feasible_elapsed_s"] == admitted[0]["end_s"]
    assert out["fine_finalists"][0]["id"] == out["winner"].identity
    assert {x.identity for x in sim.calls["fine"]} <= {x.identity for x in sim.calls["coarse"]}
    assert len({x.identity for x in sim.calls["fine"]}) == len(sim.calls["fine"])


@pytest.mark.parametrize("fine_duration, status, overrun", [(4., "qualified", 0.), (4.25, "no_design", .25)])
def test_first_fine_eligible_at_deadline_and_not_an_instant_later(run, fine_duration, status, overrun):
    sim = run(budget=10., fine=lambda i, c: (fine_duration, result(1.)))
    assert len(sim.calls["fine"]) == len(sim.calls["coarse"]) == 1
    assert sim.output["status"] == status
    assert sim.output["counts"]["late"] == int(overrun > 0)
    assert sim.output["inflight_overrun_s"] == overrun
    assert sim.output["first_feasible_elapsed_s"] == (10. if not overrun else None)
    event = sim.events[-1]
    assert event["stage"] == "fine" and event["start_s"] == 6.
    assert event["end_s"] == 10.+overrun
    assert event["incumbent_updated"] == (overrun == 0)


@pytest.mark.parametrize("eligible, reference_only, within_budget, expected", [
    (True, False, True, "fallback"), (False, True, True, "no_design"),
    (True, True, True, "no_design"), (True, False, False, "no_design")])
def test_expired_setup_starts_no_work_and_preserves_reference_eligibility(
        run, eligible, reference_only, within_budget, expected):
    base = result(3., eligible, reference_only=reference_only, within_budget=within_budget)
    sim = run(budget=4., setup=5., baseline=base)
    assert sim.output["status"] == expected and sim.output["winner"] is None
    assert sim.output["selected"] is (base if expected == "fallback" else None)
    assert sim.output["inflight_overrun_s"] == 0.
    assert not sim.targets and not sim.events


def test_candidate_failure_is_audited_and_eligible_fallback_survives_worse_fine_scores(run):
    base = result(7.)
    sim = run(baseline=base, fine=lambda i, c: (2., ValueError("synthetic integration failure")
                                               if i == 0 else result(7. + i)))
    assert sim.output["status"] == "fallback" and sim.output["selected"] is base
    assert sim.output["counts"]["failed"] == 1 and sim.output["first_feasible_elapsed_s"] is not None
    failure, = [e for e in sim.events if e["status"] == "failed"]
    assert failure["stage"] == "fine" and "synthetic integration failure" in str(failure["errors"])


def test_late_projection_is_audited_and_never_starts_coarse(run):
    sim = run(budget=6., project=lambda i, x: (1.5, x))
    assert len(sim.events) == 1
    assert sim.events[0]["stage"] == "proposal" and not sim.events[0]["within_budget"]
    assert sim.output["inflight_overrun_s"] == .5
    assert not sim.calls["coarse"]


def test_cold_family_origin_and_matched_seed_with_real_proposals(run):
    full, floor = (run(family=f) for f in ("full", "floor-only"))
    np.testing.assert_equal(full.targets[0], np.zeros(4))
    np.testing.assert_equal(full.targets[:3], floor.targets[:3])
    assert full.calls["coarse"][0] == floor.calls["coarse"][0]
    full_vectors = np.array([c.vector for c in full.calls["coarse"]])
    assert np.all(np.ptp(full_vectors, axis=0) > 0)
    assert all(c.vector[:12] == (0.,)*12 for c in floor.calls["coarse"])
    assert np.all(np.ptp([c.vector[12:] for c in floor.calls["coarse"]], axis=0) > 0)
    assert any(e["proposal_kind"].startswith("local") for e in full.events if e["stage"] == "coarse")


@pytest.mark.parametrize("family, expected", [
    ("full", "1469b4ea8dd25bccd71cbb8102051acae77283c52242dab07b5492a105025a82"),
    ("floor-only", "9a88ce582c79d842245f892c4fb3a1bf884307246b8536c4e3b466ac0f50e8ac"),
])
def test_historical_seeded_proposal_stream(run, family, expected):
    from hashlib import sha256

    sim = run(family=family, budget=90.)
    identities = [e["id"] for e in sim.events if e["stage"] == "proposal"]
    assert len(identities) == 53
    assert sha256("\n".join(identities).encode()).hexdigest() == expected


@pytest.mark.parametrize("family,face", [
    ("floor-entry", "omp_entry"), ("floor-shelf", "shelf"), ("floor-connector", "omp_connector"),
])
def test_focused_proposals_move_only_named_face_and_projected_floor(run, family, face):
    from d3d_optimization.independent_sources import split_design_wall
    from d3d_optimization.user_input import REFERENCE

    # Candidate.bind appends eight zero IMP coefficients. Use the actual
    # engineering order consumed by TangentEvaluator.curves for the rest.
    _, surfaces, _ = split_design_wall(json.loads(REFERENCE.read_text()))
    names = [s["name"] for s in surfaces]
    sim = run(family=family, budget=90., project=lambda i, target: (0., target * .9))
    proposals = [e for e in sim.events if e["stage"] == "proposal"]
    assert {e["proposal_kind"] for e in proposals} == {
        "projected_origin", "projected_restart", "local_projected"}
    vectors = np.array([e["candidate"]["vector"] for e in proposals])
    for j, name in enumerate(names[:4]):
        block = vectors[:, 4*j:4*j+4]
        if name in (face, "floor"):
            assert np.all(np.ptp(block, axis=0) > 0.)
        else:
            np.testing.assert_array_equal(block, 0.)
    floor_start = 4 * names.index("floor")
    np.testing.assert_array_equal(vectors[:, floor_start:floor_start+4], .9*np.array(sim.targets))


def test_next_fine_confirmation_uses_best_pending_coarse_peak(run):
    peaks = [20., 15., 13., 14., 12., 10., 8.]
    sim = run(coarse=lambda i,c: (1., result(peaks[i % len(peaks)], A2=1000.+i)))
    pending = {}
    for event in sim.events:
        if event["stage"] == "coarse" and event["within_budget"] and event["eligible"]:
            pending[event["id"]] = event["peak_q"]
        if event["stage"] == "fine":
            assert pending[event["id"]] == min(pending.values())
            del pending[event["id"]]


def test_rejected_projection_and_coarse_ineligible_do_not_reach_fine(run):
    sim = run(project=lambda i, x: (1., None if i < 2 else x),
              coarse=lambda i,c: (1., result(.01, feasible=False)))
    assert sim.output["status"] == "no_design"
    assert sim.output["counts"]["projection_rejected"] == 2 and not sim.calls["fine"]
    assert any("exposure" in e["errors"] for e in sim.events if e["stage"] == "coarse")


def test_backend_failure_ends_search_with_incumbent_and_auditable_reason(run):
    # A candidate-specific integration ValueError can be skipped. A broken
    # backend must not consume the rest of the campaign silently retrying.
    sim = run(fine=lambda i,c: (2., result(5.) if i == 0 else RuntimeError("backend unavailable")))
    assert sim.output["status"] == "qualified" and sim.output["winner"] == sim.calls["fine"][0]
    assert sim.output["termination"] == "failed"
    assert sim.output["counts"]["failed"] == 1
    assert sim.events[-1]["status"] == "failed" and "backend unavailable" in str(sim.events[-1]["errors"])
