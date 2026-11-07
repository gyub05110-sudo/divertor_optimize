"""Synthetic clocks and ledgers for the registered prospective arm worker."""
import copy
from dataclasses import asdict
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
from time import perf_counter_ns
from types import SimpleNamespace

import numpy as np
import pytest

from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.reference_transport import ReferencePoint
from verification import nf_anytime_search as historical


def test_prospective_worker_has_an_independent_additive_entry_point():
    assert importlib.util.find_spec("verification.nf_prospective_worker"), "configured worker is missing"


@pytest.fixture
def worker():
    return importlib.import_module("verification.nf_prospective_worker")


def measured(peak=10., area=2.):
    if peak is None:
        return dict(eligible=False, peak=None, A2=None, errors=["continuous_shape_constraint_failed"], ledger={})
    return dict(eligible=True, peak=peak, A2=area, errors=[],
                ledger=dict(prescribed=1., deposited=1., unhit=0., balance_error=0.))


def arm(population=24, budget=600, policy="fine_adaptive"):
    projected = 8 if population == 24 else 12
    active = "coarse" if policy == "coarse_adaptive" else "fine"
    cuts = [200, 400, 600] + ([1800] if budget == 1800 else [])
    registration = dict(policy=policy, population=population, projected_seeds=projected, budget_s=budget,
                        checkpoints_s=cuts)
    raw, registry = [], {}

    def event(stage, start, end, status="completed", **data):
        item = dict(stage=stage, start_s=start, end_s=end, status=status, **data)
        raw.append(item)
        return item

    event("cold_inputs", 1., 2.)
    event("prepare_fine", 2., 3.)
    base = measured()
    event("baseline", 3., 4., measurement=base)
    event("build_fine", 4., 5.)
    if active == "coarse":
        event("prepare_coarse", 5., 6.)
        event("build_coarse", 6., 7.)
    rng = np.random.default_rng(1107)
    vectors = rng.uniform(-.09, .09, (population, 16))
    evidence = historical._population_evidence(vectors, True, 0, 0)
    event("initialize", 8., 9., level=active, restart=0, archive_count=0, inherited_count=0,
          full_space_required=True, initial_population=evidence)
    index, generation, total = 0, 0, 0

    def score(vector, state, end, *, kind, cached=False, status="completed", start=None):
        nonlocal index, generation, total
        candidate = Candidate("OMP", 20, tuple(float(v) for v in vector))
        stage = "coordinate_polish" if kind == "coordinate_polish" else "differential_evolution"
        start = end-1. if start is None else start
        common = dict(id=candidate.identity, candidate=candidate.record(), search_stage=stage,
                      proposal_kind=kind, index=index, generation=generation)
        registry.setdefault(candidate.identity, dict(i=len(registry), id=candidate.identity,
                            candidate=candidate.record(), first_seen_s=start))
        event("proposal", start, start+.01, level=active, **common)
        event(active, start+.02, end, status=status, cached=cached,
              physics_evaluated=state["peak"] is not None, measurement=state,
              objective=state["peak"]/10. if state["eligible"] else 11., **common)
        if status == "completed":
            if kind != "coordinate_polish":
                index = (index+1) % population
                generation += index == 0
            total += 1

    for i, vector in enumerate(vectors):
        if total and total % 19 == 0:
            score(rng.uniform(-.09, .09, 16), measured(None), 11.+2*total, kind="coordinate_polish")
        score(vector, base if i == 0 else measured(None), 11.+2*total,
              kind="population_initialization")
    while total % 19 == 0:
        score(rng.uniform(-.09, .09, 16), measured(None), 11.+2*total, kind="coordinate_polish")
    score(rng.uniform(-.09, .09, 16), measured(None), 190., kind="de_trial")
    score(vectors[0], base, 200., kind="de_trial", cached=True)
    score(rng.uniform(-.09, .09, 16), measured(8., 3.), 400., kind="de_trial")
    if budget == 1800:
        score(rng.uniform(-.09, .09, 16), measured(7., 3.5), 900., kind="de_trial")
    score(rng.uniform(-.09, .09, 16), measured(.5, 9.), budget+2., start=budget-2.,
          kind="de_trial", status="late")
    first = next(e for e in raw if e["stage"] == active and e["measurement"]["eligible"])
    parity = dict(fine=None, coarse=None)
    parity[active] = dict(candidate=first["candidate"], peak_relative_difference=0.,
                          density_scaled_max_difference=0., quadrature_absolute_difference=0.,
                          cpu_seconds=.01, gpu_seconds=.01, passed=True)
    record = dict(schema="nf-anytime-search.v1", policy=policy, seed=20611, budget_s=float(budget),
                  startup_s=.5, stop="deadline", elapsed_s=float(budget+3), error=None, source="OMP",
                  radius_mm=20, levels=dict(coarse=list(historical.LEVELS[0]), fine=list(historical.LEVELS[2])),
                  reference=asdict(ReferencePoint()), baseline=base, baseline_ready_s=4., events=raw,
                  candidates=list(registry.values()), promotions=[], rejected={}, completed_draws=0,
                  adaptive_config=dict(population=population, projected_seeds=projected,
                                       initial_population_gate="startup_only_full_affine_dimension"),
                  prospective_config=registration, parity=parity)
    record["checkpoints"] = {str(t): historical.snapshot(record, t) for t in cuts}
    record["final"] = historical.snapshot(record, budget)
    return record


def audit(worker, record, raw=None):
    return worker.audit_arm(record, copy.deepcopy(record["events"]) if raw is None else raw,
                            policy=record["policy"], seed=20611, budget_s=record["budget_s"],
                            population=record["adaptive_config"]["population"],
                            projected_seeds=record["adaptive_config"]["projected_seeds"])


@pytest.mark.parametrize("budget", [600, 1800])
def test_only_fresh_eligible_returned_de_scores_enter_causal_counts(worker, budget):
    result = audit(worker, arm(budget=budget))
    assert result["normal"] and result["startup_full_rank"]
    early, exact = (result["checkpoints"][str(t)] for t in (200, 400))
    assert (early["de_requests"], early["de_proposals"], early["eligible_thermal_de_trials"]) == (2, 2, 0)
    assert (early["peak"], early["A2"]) == (10., 2.)
    assert (exact["eligible_thermal_de_trials"], exact["peak"], exact["peak_ratio"]) == (1, 8., .8)
    assert result["eligible_thermal_de_trials"] == (1 if budget == 600 else 2)
    assert result["de_requests"] == result["completed_de_requests"]+1
    assert result["late_scores"][-1]["measurement"]["peak"] == .5
    assert result["peak"] == (8. if budget == 600 else 7.)
    assert result["A2"] == (3. if budget == 600 else 3.5)


def test_coarse_de_scores_are_counted_but_only_fine_returns_select_a_wall(worker):
    result = audit(worker, arm(population=64, policy="coarse_adaptive"))
    assert result["normal"] and result["peak_ratio"] == 1. and result["A2"] == 2.
    assert result["selected"]["kind"] == "original_wall" and result["eligible_thermal_de_trials"] == 1


def test_terminal_proposal_without_a_score_is_only_a_proposal(worker):
    record = arm()
    record["events"].pop()
    result = audit(worker, record)
    assert result["de_proposals"] == 4 and result["de_requests"] == 3
    assert result["eligible_thermal_de_trials"] == 1


@pytest.mark.parametrize("mutation", ["moments", "rank", "warm", "configuration", "parity", "late", "checkpoint", "initial_as_de", "thermal_flag", "event_stream"])
def test_inconsistent_independent_evidence_rejected(worker, mutation):
    r = arm()
    raw = r["events"]
    init = next(e for e in raw if e["stage"] == "initialize")
    if mutation == "moments": init["initial_population"]["center"][0] += .01
    if mutation == "rank":
        x = np.asarray(init["initial_population"]["vectors"])
        x[:, 15] = x[:, 14]
        init["initial_population"] = historical._population_evidence(x, True, 0, 0)
    if mutation == "warm": init["archive_count"] = 1
    if mutation == "configuration": r["prospective_config"]["projected_seeds"] = 12
    if mutation == "parity": r["parity"]["fine"]["passed"] = False
    if mutation == "late": raw[-1]["status"] = "completed"
    if mutation == "checkpoint": r["checkpoints"]["400"]["incumbent"]["peak"] = .5
    if mutation == "initial_as_de":
        i = next(i for i, e in enumerate(raw) if e.get("proposal_kind") == "population_initialization")
        raw[i]["proposal_kind"] = raw[i+1]["proposal_kind"] = "de_trial"
    if mutation == "thermal_flag": next(e for e in raw if e["stage"] == "fine")["physics_evaluated"] = False
    if mutation == "event_stream": raw = raw[:-1]
    with pytest.raises(ValueError): audit(worker, r, raw)


def test_failure_retains_confirmed_endpoint_and_disqualifies_arm(worker):
    r = arm()
    end = r["events"][-1]
    end.update(end_s=599.9, status="failed", error="synthetic failure")
    r.update(elapsed_s=599.95, stop="failed", error="synthetic failure")
    r["checkpoints"] = {str(t): historical.snapshot(r, t) for t in (200, 400, 600)}
    r["final"] = historical.snapshot(r, 600)
    result = audit(worker, r)
    assert not result["normal"] and result["reasons"]
    assert result["peak"] == 8. and result["eligible_thermal_de_trials"] == 1


@pytest.mark.parametrize("policy,population,projected,budget", [("coarse_adaptive",64,12,600),("fine_adaptive",24,8,600),("fine_adaptive",32,12,1800)])
def test_isolated_runner_uses_parent_clock_and_exhausted_budget(worker, tmp_path, policy, population, projected, budget):
    started = perf_counter_ns()-(budget+1)*10**9
    record = worker.run_configured(policy, population, projected, tmp_path/"input", 20611, tmp_path/"arm", started, budget)
    assert record["adaptive_config"]["population"] == population
    assert record["adaptive_config"]["projected_seeds"] == projected
    assert record["stop"] == "deadline" and record["events"] == []
    assert record["startup_s"] >= budget and (tmp_path/"arm/events.jsonl").read_bytes() == b""
    assert list(record["checkpoints"]) == (["200","400","600"] if budget == 600 else ["200","400","600","1800"])
    assert not audit(worker, record)["normal"]
    assert (historical.POPULATION, historical.PROJECTED_SEEDS, historical.CHECKPOINTS) == (64, 12, (200,400,600))


def test_worker_refuses_unregistered_config_and_historical_seed_file(worker, tmp_path):
    with pytest.raises(ValueError): worker.run_configured("fine_adaptive", 64, 12, tmp_path, 20611, tmp_path/"a", perf_counter_ns(), 600)
    with pytest.raises(ValueError): worker.run_configured("fine_adaptive", 24, 8, tmp_path, 20611, tmp_path/"a", perf_counter_ns(), 601)
    (tmp_path/"campaign_warm_candidates.json").write_text("historical")
    r = worker.run_configured("fine_adaptive", 24, 8, tmp_path, 20611, tmp_path/"a", perf_counter_ns(), 600)
    assert r["stop"] == "failed" and "historical" in r["error"]
    assert not audit(worker, r)["normal"]


@pytest.mark.parametrize("policy,size,n_seeds", [("coarse_adaptive",64,12),("fine_adaptive",24,8),("fine_adaptive",32,12)])
@pytest.mark.parametrize("rank_deficient", [False, True])
def test_worker_uses_native_cold_initializer_with_effective_gate(worker, tmp_path, monkeypatch,
                                                               policy, size, n_seeds, rank_deficient):
    from d3d_optimization import campaign_search

    started = perf_counter_ns()
    clock = [started]
    projections = []
    original_loader = worker._engine

    def loader(config):
        engine = original_loader(config)
        native_run = engine.run_arm
        engine.run_arm = lambda *a, **kw: native_run(*a, **kw, clock_ns=lambda: clock[0])
        engine.prep = lambda root, source, ref, level: (SimpleNamespace(
            metadata={"source": source}, reference=ref, baseline=measured()), 0.)
        engine.measurement = lambda value, ref: value
        return engine

    def evaluator(original, radius, backend):
        return SimpleNamespace(metadata=original.metadata, rmin=radius/1000, reference=original.reference,
                               lengths=dict(omp_entry=1., shelf=1., omp_connector=1., floor=1.))

    def project(geom, target):
        projections.append(target.copy())
        return target

    def single_initialization(self, seconds, inherited=()):
        self.initialize(inherited)
        clock[0] = started+600*10**9

    monkeypatch.setattr(worker, "_engine", loader)
    monkeypatch.setattr(campaign_search, "evaluator", evaluator)
    monkeypatch.setattr(campaign_search, "project_floor", project)
    monkeypatch.setattr(campaign_search.SearchCase, "run_slice", single_initialization)
    if rank_deficient:
        native_population = campaign_search.new_population
        monkeypatch.setattr(campaign_search, "new_population",
                            lambda vectors, seed: native_population(np.zeros((len(vectors),16)), seed))
    r = worker.run_configured(policy, size, n_seeds, tmp_path/"input", 20611, tmp_path/"run", started, 600)
    raw = [json.loads(line) for line in (tmp_path/"run/events.jsonl").read_bytes().splitlines()]
    result = audit(worker, r, raw)
    assert len(projections) == n_seeds
    init = next(e for e in raw if e["stage"] == "initialize")
    assert init["initial_population"]["count"] == size
    assert init["initial_population"]["archive_origin"] == "none" and init["inherited_count"] == 0
    if rank_deficient:
        assert r["stop"] == "failed" and not result["normal"] and init["initial_population"]["rank"] == 0
        assert result["startup_population"]["rank"] == 0 and result["startup_population"]["status"] == "failed"
        assert f"{size} members" in init["error"] and f"{size} members" in r["error"]
    else:
        assert r["stop"] == "deadline" and result["normal"] and init["initial_population"]["rank"] == 16


def test_cli_authenticates_job_before_work_and_persists_failed_evidence(worker, tmp_path, monkeypatch):
    out, root = tmp_path/"arm", tmp_path/"case"
    out.mkdir(); root.mkdir()
    (root/"campaign_warm_candidates.json").write_text("historical")
    job = dict(case_id="case", input_root=str(root.resolve()), seed=20611, policy="fine_adaptive", population=24,
               projected_seeds=8, budget_s=600, phase="development", registration_sha256="a"*64)
    payload = json.dumps(job).encode()
    (out/"job.json").write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    started = perf_counter_ns()
    args = ["--policy","fine_adaptive","--population","24","--projected-seeds","8","--input-root",str(root),
            "--seed","20611","--out",str(out),"--started-ns",str(started),"--budget-s","600.0"]
    monkeypatch.setenv("NF_PROSPECTIVE_JOB_SHA256", "0"*64)
    with pytest.raises(ValueError): worker.main(args)
    assert not (out/"events.jsonl").exists()
    monkeypatch.setenv("NF_PROSPECTIVE_JOB_SHA256", digest)
    with pytest.raises(ValueError): worker.main(["999" if a == "20611" else a for a in args])
    assert not (out/"events.jsonl").exists()
    assert worker.main(args) != 0
    result = json.loads((out/"result.json").read_bytes())
    assert result["job_sha256"] == digest and result["case_id"] == "case" and result["launch_ns"] == started
    assert result["registration_sha256"] == job["registration_sha256"] and result["pid"] > 0
    assert Path(result["input_root"]) == root and result["events"] == [json.loads((out/"events.jsonl").read_bytes())]


CASES = ("deb842cf1e", "fad1624516", "4d14a35fcb")


def observations(peaks24=(8.,8.,8.), peaks32=(9.,9.,9.), de24=(1,1,0), de32=(1,1,1)):
    out = []
    for pop, proj, peaks, counts in [(24,8,peaks24,de24), (32,12,peaks32,de32)]:
        for case, peak, n in zip(CASES, peaks, counts):
            out.append(dict(policy="fine_adaptive", population=pop, projected_seeds=proj, seed=20611, budget_s=600,
                            case_id=case, normal=True, startup_full_rank=True, baseline_peak=10., peak=peak, A2=2.,
                            peak_ratio=peak/10., de_requests=12, eligible_thermal_de_trials=n))
    return out


def test_selection_requires_thermal_de_in_two_inputs_and_all_three_normal_endpoints(worker):
    result = worker.choose_fine(observations(de24=(0,1,0)))
    assert result["qualified"] and result["selected"]["population"] == 32
    obs = observations(de24=(1,1,0), de32=(0,0,0))
    obs[1]["normal"] = False
    failure = worker.choose_fine(obs)
    assert not failure["qualified"] and failure["selected"] is None and failure["reasons"]
    assert all(len(c["by_case"]) == 3 and c["reasons"] for c in failure["configurations"])
    assert not worker.choose_fine(observations()[:-1])["qualified"]


def test_selection_is_paired_geometric_mean_then_de_cases_then_population(worker):
    # Arithmetic means rank these configurations the opposite way.
    obs = observations(peaks24=(1.,9.,9.), peaks32=(5.,5.,5.))
    assert worker.choose_fine(obs)["selected"]["population"] == 24
    tied = observations(peaks32=(8.,8.,8.))
    assert worker.choose_fine(tied)["selected"]["population"] == 32
    tied[-1]["eligible_thermal_de_trials"] = 0
    assert worker.choose_fine(tied)["selected"]["population"] == 24
    with pytest.raises(ValueError): worker.choose_fine(tied + [tied[0]])
