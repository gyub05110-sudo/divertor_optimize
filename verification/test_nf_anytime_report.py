"""Synthetic report checks; execute only after the live timing campaign stops."""
import hashlib
from itertools import combinations
import json

import numpy as np
import pytest

from verification import nf_anytime_report as report


def measurement(peak=10., area=2., eligible=True):
    return dict(eligible=eligible, peak=peak, A2=area, errors=[], ledger={})


def shape(i=0):
    candidate = dict(source="OMP", radius_mm=20, vector=[float(i)/1000.]*16)
    identity = hashlib.sha256(json.dumps(candidate, sort_keys=True).encode()).hexdigest()
    return dict(id=identity, candidate=candidate)


def baseline(end=5.):
    return dict(stage="baseline", start_s=end-1., end_s=end, status="completed",
                measurement=measurement())


def score(i, end, peak, *, stage="fine", status="completed", cached=False,
          physics=True, eligible=True, area=3., kind="population_initialization"):
    return dict(stage=stage, start_s=end-1., end_s=end, status=status, **shape(i),
                cached=cached, physics_evaluated=physics, proposal_kind=kind,
                measurement=measurement(peak, area, eligible))


def arm(events, policy="fine_adaptive", **kwargs):
    return dict(case_id="case_a", seed=20611, policy=policy, budget_s=600.,
                events=events, elapsed_s=600., process=dict(elapsed_s=601.),
                audit_status="completed", native_pool_eligible=True,
                stop="deadline", parity={}, **kwargs)


def test_curve_starts_only_when_baseline_is_available_and_keeps_exact_cutoff():
    a = arm([baseline(40.), score(1, 200., 9.), score(2, 601., 1., status="late")])
    assert not report.reconstruct(a, 39.99)["available"]
    assert report.reconstruct(a, 200.)["peak"] == 9.
    curve = report.anytime_curve(a)
    assert curve[0]["t_s"] == 40.
    assert curve[-1]["t_s"] == 600. and curve[-1]["peak"] == 9.
    assert report.reconstruct(a, 900.)["peak"] == 9.


def test_coarse_ineligible_cached_and_worse_results_do_not_win():
    a = arm([baseline(), score(1, 6., .1, stage="coarse"), score(2, 7., .2, eligible=False),
             score(3, 8., .3, cached=True), score(4, 9., 11.)])
    assert report.reconstruct(a, 200.)["kind"] == "original_wall"


def test_first_five_percent_keeps_time_or_explicit_right_censoring():
    successful = arm([baseline(), score(1, 17., 9.5)])
    observed = report.first_improvement(successful)
    assert observed["event_observed"] and observed["time_s"] == 17.
    failed = arm([baseline(), score(1, 10., 9.8)])
    failed.update(elapsed_s=37., stop="failed", audit_status="arm_failed")
    censored = report.first_improvement(failed)
    assert censored["time_s"] is None and not censored["event_observed"]
    assert censored["censor_time_s"] == 37.
    assert censored["reason"] == "arm_failed"
    failed.update(stop="search_stopped", audit_status="completed")
    assert report.first_improvement(failed)["reason"] == "search_stopped_before_budget"
    assert report.anytime_curve(failed)[-1]["t_s"] == 37.


def test_uncertified_process_failure_cannot_supply_a_curve_or_incumbent():
    a = arm([baseline(), score(1, 20., 1.)])
    a.update(audit_status="process_failed", native_pool_eligible=False)
    assert report.anytime_curve(a) == []
    assert not report.reconstruct(a, 600.)["available"]
    assert report.first_improvement(a)["time_s"] is None


def test_late_fine_diagnostic_preserves_strict_and_tolerance_comparisons_without_admission():
    a = arm([baseline(), score(1, 600., 8.), score(2, 601., 8.-1e-11, status="late"),
             score(3, 602., 7., status="late"), score(4, 603., .1, stage="coarse", status="late")])
    result = report.late_result_diagnostic(a)
    assert result["late_coarse_count"] == 1 and result["late_fine_count"] == 2
    assert result["genuine_eligible_late_fine_count"] == 2
    assert result["strict_improvement_count"] == 2
    assert result["improvement_beyond_replay_tolerance_count"] == 1
    assert result["fine_results"][0]["strict_improvement"]
    assert result["fine_results"][0]["replay_tolerance_ordering"] == "within_replay_tolerance"
    assert result["fine_results"][1]["relative_improvement"] == pytest.approx(.125)
    assert result["incumbent_at_cutoff"]["peak"] == report.reconstruct(a, 600.)["peak"] == 8.
    assert report.first_improvement(a)["time_s"] == 600.
    assert report.anytime_curve(a)[-1]["peak"] == 8.
    assert {r["id"] for r in result["fine_results"]} == {shape(2)["id"], shape(3)["id"]}


def test_late_diagnostic_keeps_missing_ineligible_cached_and_partial_results_uncompared():
    partial = score(2, 600., 1., status="late", physics=False)
    missing = score(3, 602., None, status="late", physics=False, eligible=False)
    missing.pop("measurement")
    a = arm([baseline(), score(1, 600., 9.), partial, missing,
             score(4, 603., 1., status="late", eligible=False),
             score(5, 604., 1., status="late", cached=True)])
    result = report.late_result_diagnostic(a)
    assert result["late_fine_count"] == 4
    assert result["genuine_eligible_late_fine_count"] == result["compared_fine_count"] == 0
    assert all(r["strict_improvement"] is None and r["relative_improvement"] is None
               for r in result["fine_results"])
    assert report.reconstruct(a, 600.)["peak"] == 9.
    assert report.first_improvement(a)["time_s"] == 600.


def test_late_fine_without_available_baseline_has_no_fabricated_comparison():
    a = arm([score(1, 601., 1., status="late")])
    result = report.late_result_diagnostic(a)
    assert result["genuine_eligible_late_fine_count"] == 1
    assert result["compared_fine_count"] == result["strict_improvement_count"] == 0
    assert result["fine_results"][0]["comparison_reason"] == "incumbent_unavailable"
    assert result["fine_results"][0]["relative_improvement"] is None
    assert report.first_improvement(a)["time_s"] is None


def test_fine_adaptive_phase_disclosure_uses_completed_checkpoint_counts():
    arms = []
    for seed in report.SEEDS:
        a = arm([baseline(), score(1, 20., 9.), score(2, 30., 8., kind="coordinate_polish"),
                 score(3, 601., 7., status="late", kind="de_trial")])
        a["seed"] = seed
        arms.append(a)
    summary = report.summarize_live(arms)
    phase = summary["fine_adaptive_phase_at_600"]
    assert phase["arm_count"] == phase["zero_completed_de_trial_arm_count"] == 3
    assert phase["completed_proposal_kinds"] == {"population_initialization": 3, "coordinate_polish": 3}
    assert phase["only_initialization_and_polish"]
    text = report.fine_adaptive_phase_text(phase)
    assert "3/3" in text and "DE trial 완료 0건" in text
    assert "최적 튜닝" in text and "점근적" in text
    arms[0]["events"].append(score(4, 599., 7., kind="de_trial"))
    changed = report.summarize_live(arms)["fine_adaptive_phase_at_600"]
    assert changed["zero_completed_de_trial_arm_count"] == 2
    assert not changed["only_initialization_and_polish"]
    assert "DE trial 완료 0건" not in report.fine_adaptive_phase_text(changed)


def test_requests_physics_cache_and_geometry_rejections_have_distinct_counts():
    rejected = score(3, 40., None, physics=False, eligible=False, kind="de_trial")
    rejected["measurement"]["errors"] = ["continuous_shape_constraint_failed"]
    a = arm([baseline(), score(1, 20., 9.), score(1, 30., 9., cached=True), rejected,
             score(4, 201., 8., kind="coordinate_polish")])
    counts = report.checkpoint_metrics(a, 200.)["counts"]["fine"]
    assert counts["requests_started"] == 4
    assert counts["requests_completed"] == 3
    assert counts["in_flight"] == 1
    assert counts["deposition_result_evaluations"] == 1
    assert counts["cached_requests"] == counts["geometry_rejections"] == 1
    assert counts["proposal_kinds"] == {"population_initialization": 2, "de_trial": 1}


def test_active_coverage_uses_actual_vectors_and_excludes_fine_promotion():
    coarse = score(1, 20., 8., stage="coarse")
    proposed = dict(stage="proposal", level="coarse", start_s=18., end_s=19.,
                    status="completed", **shape(1))
    fine = score(2, 30., 7., kind="external_validation")
    a = arm([baseline(), proposed, coarse, fine], "coarse_adaptive")
    metrics = report.checkpoint_metrics(a, 200.)
    assert metrics["active_level"] == "coarse"
    assert metrics["coverage"]["proposed"]["unique_count"] == 1
    assert metrics["coverage"]["deposition"]["unique_count"] == 1
    assert metrics["coverage"]["eligible"]["unique_count"] == 1
    assert metrics["coverage"]["deposition"]["rank"] == 0


def test_diversity_uses_centered_fixed_bounds_and_does_not_count_initial_archive():
    x = np.vstack([np.zeros(16), np.eye(16)*.01])
    actual = report.vector_coverage(x.tolist())
    assert actual["rank"] == 16 and actual["unique_count"] == 17
    assert actual["axis_span"] == [.01]*16 and actual["scale"] == [.2]*16
    flat = x[:, 12:]
    floor_only = np.column_stack([np.zeros((len(x), 12)), flat])
    assert report.vector_coverage(floor_only.tolist())["rank"] == 4
    a = arm([baseline(), dict(stage="initialize", start_s=5., end_s=6., status="completed",
                             initial_population=dict(vectors=x.tolist(), rank=16))])
    assert report.checkpoint_metrics(a, 200.)["coverage"]["deposition"]["rank"] == 0


def test_evaluated_coverage_excludes_late_and_cached_scores():
    a = arm([baseline(), score(1, 20., 9.), score(2, 30., 8., cached=True),
             score(3, 601., 1., status="late")])
    pool = report.checkpoint_metrics(a, 600.)["coverage"]["eligible"]
    assert pool["unique_count"] == 1


def test_equal_k_expectation_matches_exhaustive_subsets_and_retains_failures():
    values, baseline_peak, k = [7., None, 12., 8., 6.], 10., 2
    expected = sum(min([baseline_peak] + [values[i] for i in ids if values[i] is not None])
                   for ids in combinations(range(len(values)), k))/10.
    assert report.exact_equal_k_expectation(values, k, baseline_peak) == pytest.approx(expected)
    assert report.exact_equal_k_expectation(values, 0, baseline_peak) == baseline_peak


def test_finite_panel_comparator_uses_common_coarse_universe_and_counts_failed_requests():
    rows = []
    for i in range(6):
        rows.append(dict(i=i, **shape(i+1), coarse=measurement(1.+i, eligible=i != 5),
                         fine=measurement(8.+i, eligible=i != 0)))
    rows[5]["fine"] = measurement(1.)  # Best design excluded by coarse screen.
    result = report.panel_metrics(rows, measurement())
    assert result["fine_requests_k"] == 1
    assert result["uniform_population_size"] == 5
    assert result["selected_peak"] == 10.  # Failed promoted request retains wall.
    assert result["reference_peak"] == 1.
    assert result["uniform_expected_peak"] == pytest.approx((10.+9.+10.+10.+10.)/5.)
    assert not result["equal_wall_time_comparison"]


def test_paired_settings_distinguish_ordering_consistency_from_range_overlap():
    first = {s: measurement(v) for s, v in zip(report.SETTINGS, [8., 9., 10., 11.])}
    second = {s: measurement(v) for s, v in zip(report.SETTINGS, [9., 10., 11., 12.])}
    pair = report.setting_pair(first, second)
    assert pair["paired_ordering"] == "first_lower_all_settings"
    assert pair["range_relation"] == "overlap_or_touch"
    assert not pair["fine_ordering_resolved_by_ranges"]


def test_paired_settings_report_reversal_tie_and_missing_eligibility():
    first = {s: measurement(9.) for s in report.SETTINGS}
    second = {s: measurement(9.) for s in report.SETTINGS}
    assert report.setting_pair(first, second)["paired_ordering"] == "exact_tie_all_settings"
    second[report.SETTINGS[0]] = measurement(8.)
    second[report.SETTINGS[1]] = measurement(10.)
    assert report.setting_pair(first, second)["paired_ordering"] == "setting_order_reversal"
    second[report.SETTINGS[2]] = measurement(None, eligible=False)
    assert report.setting_pair(first, second)["paired_ordering"] == "incomplete_or_ineligible"


def test_tiny_strict_native_difference_is_a_nominal_win_with_separate_tolerance_diagnostic():
    rows = []
    for seed in (20611, 20627):
        for policy, peak in zip(report.POLICIES, (8., 8.+1e-11, 9.)):
            a = arm([baseline(), score(1, 20., peak)], policy)
            a["seed"] = seed
            rows.append(a)
    pairs = report.paired_comparisons(rows, (200,))
    assert len(pairs) == 6
    assert all(p["case_id"] == "case_a" for p in pairs)
    chosen = [p for p in pairs if p["first_policy"] == report.POLICIES[0]
              and p["second_policy"] == report.POLICIES[1]]
    assert [p["ordering"] for p in chosen] == ["first_lower", "first_lower"]
    assert [p["replay_tolerance_ordering"] for p in chosen] == ["within_replay_tolerance"]*2
    tiny = arm([baseline(), score(1, 20., 10.-1e-11, area=2.+1e-11)])
    metrics = report.checkpoint_metrics(tiny, 200.)
    assert metrics["q_down"] and metrics["q_down_A2_up"]
    assert not metrics["q_down_beyond_replay_tolerance"]


def test_exact_setting_range_separation_is_preserved_before_tolerance_diagnostic():
    first = {s: measurement(9.) for s in report.SETTINGS}
    second = {s: measurement(9.+1e-11) for s in report.SETTINGS}
    result = report.setting_pair(first, second)
    assert result["range_relation"] == "first_below_second"
    assert result["paired_ordering"] == "first_lower_all_settings"
    assert result["replay_tolerance_range_relation"] == "overlap_or_touch"


def test_native_replay_coverage_separates_unobserved_and_eligibility_only_candidates():
    def replay(previous, observations, observables, eligibility_only):
        return dict(passed=True, all_available_observations_matched=True,
                    previously_scored=previous, compared_count=observations,
                    observable_compared_count=observables, eligibility_only_count=eligibility_only)

    genuine = replay(True, 2, 4, 0)
    geometry_only = replay(True, 1, 0, 1)
    unobserved = replay(False, 0, 0, 0)
    candidates = [dict(id=str(i), native_reproducibility=row)
                  for i, row in enumerate((genuine, geometry_only, unobserved))]
    coverage = report.native_replay_coverage(candidates)
    assert coverage["locked_candidate_count"] == 3
    assert coverage["previously_scored_candidate_count"] == 2
    assert coverage["candidates_with_observable_comparisons"] == 1
    assert coverage["eligibility_only_candidate_count"] == 1
    assert coverage["unobserved_candidate_count"] == 1
    assert coverage["compared_observation_count"] == 3
    assert coverage["compared_observable_count"] == 4
    assert coverage["eligibility_only_observation_count"] == 1
    assert report.observable_native_replay_matched(genuine)
    assert not report.observable_native_replay_matched(geometry_only)
    assert not report.observable_native_replay_matched(unobserved)
    text = report.native_replay_intro(dict(all_available_native_observations_matched=True,
        cases=[dict(case_id="case_a", native_replay_coverage=coverage)]))
    assert "대조 가능한 이전 native 관측" in text
    assert "기존 native 관측이 없는 후보" in text
    assert "적격성만" in text
    assert "고정 후보 전체의 열결과 일치" not in text


def test_numerical_lead_keeps_range_separation_distinct_from_paired_setting_consistency():
    pairs = []
    for other in ("fine_adaptive", "random_fine"):
        for relation, ordering in (("first_below_second", "first_lower_all_settings"),
                                   ("overlap_or_touch", "first_lower_all_settings"),
                                   ("overlap_or_touch", "exact_tie_all_settings")):
            pairs.append(dict(cutoff_s=600, first_policy="coarse_adaptive", second_policy=other,
                              range_relation=relation, paired_ordering=ordering, available=True,
                              native_replay_verified_for_pair=True))
    text = report.numerical_lead(dict(available=True, cases=[dict(method_checkpoint_pairs=pairs)]))
    assert text.count("범위 분리 1/3") == 2
    assert text.count("같은 설정마다 저비용 피크가 낮음 2/3") == 2
    assert "미세 우열" in text and "신뢰구간" in text
    assert report.numerical_lead(dict(available=False)) == ""


def rejection_fixture(tmp_path, *, corrupt=None):
    from verification import nf_anytime_rejections as audit

    def seal(name, value):
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, sort_keys=True), encoding="utf-8")
        checksum = report.digest(path)
        path.with_name(path.name+".sha256").write_text(checksum, encoding="ascii")
        return checksum

    frozen = dict(cases=[dict(case_id="case_a")], reference={})
    seal("frozen.json", frozen)
    coarse = score(1, 20., 11., stage="coarse", eligible=False)
    coarse["measurement"]["errors"] = ["receiving_angle_constraint_failed"]
    live = arm([baseline(), coarse], "coarse_adaptive")
    live["checkpoints"] = {"600": report.reconstruct(live, 600.)}
    prefix = "live/case_a/20611/coarse_adaptive/"
    live_hashes = {"frozen.json": report.digest(tmp_path/"frozen.json"),
                   prefix+"result.json": seal(prefix+"result.json", live),
                   prefix+"events.jsonl": seal(prefix+"events.jsonl", live["events"])}
    kind = "screen-rejections/"
    sources = {}
    for path in audit.source_paths():
        name = path.relative_to(report.CODE_ROOT).as_posix()
        target = tmp_path/kind/"code-snapshot"/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
        sources[name] = report.digest(path)
        target.with_name(target.name+".sha256").write_text(sources[name], encoding="ascii")
    protocol = dict(campaign_sha256=live_hashes["frozen.json"], source_hashes=sources,
                    details=audit.protocol_details(), reference={}, source="OMP", radius_mm=20, backend="gpu")
    protocol_sha = seal(kind+"protocol.json", protocol)
    seal(kind+"live-evidence.json", live_hashes)
    selected = dict(protocol_sha256=protocol_sha, live_evidence_hashes=live_hashes,
                    **audit.select_rejections([live], live_hashes))
    if corrupt == "selection":
        selected["cases"][0]["candidates"][1]["qualifying_observations"][0]["incumbent_600"]["peak"] = 99.
    selection_sha = seal(kind+"selection.json", selected)
    chosen = selected["cases"][0]
    scores, score_hashes = {}, {}
    for row in chosen["candidates"]:
        record = dict(id=row["id"], candidate=row["candidate"], case_id="case_a", setting="cubic1280",
                      setting_tuple=["cubic", 1280, .000125, 641, .00025],
                      protocol_sha256=protocol_sha, selection_sha256=selection_sha,
                      status="completed", error=None, seconds=.01,
                      measurement=measurement(10. if row["id"] == report.BASELINE_ID else 12.))
        if corrupt == "score_protocol" and row["id"] != report.BASELINE_ID:
            record["protocol_sha256"] = "wrong"
        if corrupt == "score_setting" and row["id"] != report.BASELINE_ID:
            record["setting_tuple"][1] = 640
        scores[row["id"]] = record
        name = "case_a/scores/"+row["id"]+".cubic1280.json"
        score_hashes[name] = seal(kind+name, record)
    result = dict(protocol_sha256=protocol_sha, selection_sha256=selection_sha,
                  score_hashes=score_hashes, **audit.assess_rejections(chosen, scores))
    result_sha = seal(kind+"case_a/result.json", result)
    counts = {key: result[key] for key in ("false_negative_count", "fine_ineligible_count", "fine_evaluation_failed_count")}
    summary = dict(protocol_sha256=protocol_sha, selection_sha256=selection_sha,
                   selected_candidate_count=1, qualifying_observation_count=1, complete_live_coverage=True,
                   post_hoc=True, primary_winners_changed=False, scope=audit.SCOPE, costs={"seconds": .02},
                   cases=[dict(case_id="case_a", result_sha256=result_sha, **counts)], **counts)
    if corrupt == "summary_count":
        summary["false_negative_count"] = 0
    seal(kind+"summary.json", summary)
    return frozen, [live], live_hashes, chosen, scores


def test_rejection_metrics_distinguish_false_negative_from_gain_ineligibility_and_failure(tmp_path):
    _, _, _, selected, scores = rejection_fixture(tmp_path)
    identity = selected["candidates"][1]["id"]
    result = report.rejection_case_metrics(selected, scores)
    row = result["rows"][0]
    assert result["false_negative_count"] == 1 and row["fine_eligible"]
    assert not row["vs_baseline"]["eligible_peak_improvement"]  # Eligible even though peak worsens.
    assert row["vs_live_coarse_incumbents"][0]["reference"]["peak"] == 10.
    scores[identity]["measurement"] = measurement(8., eligible=False)
    scores[identity]["measurement"]["errors"] = ["still_ineligible"]
    result = report.rejection_case_metrics(selected, scores)
    assert result["fine_ineligible_count"] == 1 and result["false_negative_count"] == 0
    assert result["rows"][0]["vs_baseline"]["strictly_lower_observed_peak"]
    assert not result["rows"][0]["vs_baseline"]["eligible_peak_improvement"]
    scores[identity].update(status="failed", error="synthetic solver failure")
    failed = report.rejection_case_metrics(selected, scores)
    assert failed["fine_evaluation_failed_count"] == 1
    assert failed["rows"][0]["coarse_eligibility_false_negative"] is None
    assert failed["rows"][0]["fine"]["error"] == "synthetic solver failure"


def test_rejection_loader_reconstructs_locked_selection_and_keeps_separate_costs(tmp_path):
    frozen, arms, hashes, _, _ = rejection_fixture(tmp_path)
    loaded = report.load_rejections(tmp_path, frozen, arms, dict(hashes))
    assert loaded["available"] and loaded["false_negative_count"] == 1
    assert loaded["selected_candidate_count"] == loaded["qualifying_observation_count"] == 1
    assert not loaded["primary_winners_changed"] and loaded["costs"] == {"seconds": .02}
    text = report.rejections_text(loaded)
    assert "전역 false-negative 비율" in text and "cubic1280" in text
    assert "receiving_angle_constraint_failed" in text
    assert "알 수 없음" in text
    score_path = next((tmp_path/"screen-rejections/case_a/scores").glob("*.json"))
    score_path.write_bytes(score_path.read_bytes()+b" ")
    with pytest.raises(ValueError, match="byte seal"):
        report.load_rejections(tmp_path, frozen, arms, dict(hashes))


@pytest.mark.parametrize("corrupt", ["selection", "score_protocol", "score_setting", "summary_count"])
def test_rejection_loader_rejects_sealed_but_inconsistent_provenance_and_counts(tmp_path, corrupt):
    frozen, arms, hashes, _, _ = rejection_fixture(tmp_path, corrupt=corrupt)
    with pytest.raises(ValueError, match="rejection"):
        report.load_rejections(tmp_path, frozen, arms, dict(hashes))


def test_rejection_loader_allows_a_missing_postphase_without_inventing_zero_results(tmp_path):
    result = report.load_rejections(tmp_path, dict(cases=[]), [], {})
    assert not result["available"] and "false_negative_count" not in result


def test_duration_breakdown_retains_charged_gaps_and_does_not_double_count_parity():
    a = arm([baseline(5.), score(1, 12., 9.)])
    a.update(startup_s=1., elapsed_s=20., parity=dict(fine=dict(cpu_seconds=.4, gpu_seconds=.5)),
             process=dict(elapsed_s=21.))
    costs = report.stage_costs(a)
    assert costs["stage_seconds"] == {"baseline": 1., "fine": 1.}
    assert costs["clock_gaps_seconds"] == 17.
    assert costs["reported_cpu_parity_seconds"] == .4
    assert costs["admitted_scoring_event_seconds"] == 1.
    assert costs["process_finalization_and_observation_seconds"] == 1.
    assert costs["parity_is_already_inside_score_events"]


def test_timeout_parent_duration_remains_visible_without_inventing_worker_time():
    arms = []
    for case in ("a", "b", "c"):
        for seed in report.SEEDS:
            for policy in report.POLICIES:
                row = arm([baseline()], policy)
                row.update(case_id=case, seed=seed)
                arms.append(row)
    timed_out = arms[0]
    timed_out.pop("elapsed_s")
    timed_out.update(events=[], stop=None, audit_status="process_failed", native_pool_eligible=False,
                     process=dict(elapsed_s=1200.25, returncode=None, error="timeout"))
    costs = report.stage_costs(timed_out)
    assert costs["elapsed_s"] is costs["overrun_seconds"] is None
    assert costs["process_elapsed_s"] == 1200.25
    assert costs["process_overrun_seconds"] == 600.25
    summary = dict(live=report.summarize_live(arms), frozen_case_ids=["a", "b", "c"],
                   finite_panel=dict(available=False, reason="not run"),
                   numerical_settings=dict(available=False, reason="not run"),
                   study="synthetic", generated_utc="synthetic", source_sha256="synthetic")
    text = report.korean_report(summary)
    assert "부모 관측 종료" in text and "부모 관측 초과" in text
    cost_line = next(line for line in text.splitlines()
                     if line.startswith("| 저비용 적응 탐색 |") and "1200.25" in line)
    assert "600.25" in cost_line
    assert "미확인 1/9" in cost_line


def test_read_only_sealed_reader_rejects_changed_bytes(tmp_path):
    folder = tmp_path/"reference"
    folder.mkdir()
    path = folder/"summary.json"
    path.write_text('{"value":1}', encoding="utf-8")
    path.with_name("summary.json.sha256").write_text(report.digest(path), encoding="ascii")
    before = set(tmp_path.rglob("*"))
    assert report._sealed(tmp_path, "reference", "summary.json", {}) == {"value": 1}
    assert set(tmp_path.rglob("*")) == before
    path.write_text('{"value":2}', encoding="utf-8")
    with pytest.raises(ValueError, match="byte seal"):
        report._sealed(tmp_path, "reference", "summary.json", {})


def test_panel_overlap_reports_reused_live_candidates_without_holdout_claim():
    overlap = dict(live_proposed_ids=["a", "b"], live_proposed_count=2,
                   live_on_time_fine_observed_ids=["b"], live_on_time_fine_observed_count=1,
                   proposed_intersection_ids=["b"], proposed_intersection_count=1,
                   on_time_fine_intersection_ids=["b"], on_time_fine_intersection_count=1,
                   held_out_validation_claim=False, historical_archive_overlap_checked=False,
                   live_outcomes_used_for_generation=False)
    compact = report.panel_overlap(overlap, ["b", "c"])
    assert compact["proposed_intersection_count"] == compact["on_time_fine_intersection_count"] == 1
    assert not compact["held_out_validation_claim"]
    assert "live_proposed_ids" not in compact
    overlap["on_time_fine_intersection_count"] = 0
    with pytest.raises(ValueError, match="overlap"):
        report.panel_overlap(overlap, ["b", "c"])


def reference_panel(case, seed, entries, baseline_peak=10.):
    return dict(case_id=case, seed=seed, baseline=measurement(baseline_peak),
                panel_sha256=f"panel-{case}-{seed}", result_sha256=f"result-{case}-{seed}",
                rows=[dict(i=i, **shape(identity), fine=measured)
                      for i, (identity, measured) in enumerate(entries)])


def test_panel_shape_counts_distinguish_measurement_rows_from_unique_geometries():
    panels = [reference_panel(case, seed, [(i, measurement()) for i in identities])
              for case, seed, identities in (("a", 1, (1, 2)), ("a", 2, (1, 3)),
                                             ("b", 1, (1, 2)), ("b", 2, (1, 4)),
                                             ("c", 1, (1, 5)), ("c", 2, (1, 6)))]
    counts = report.panel_diagnostics(panels)["shape_counts"]
    assert counts["fine_measurement_rows"] == 12
    assert counts["candidate_case_count"] == 9
    assert counts["global_unique_shape_count"] == 6
    assert counts["common_to_all_shape_ids"] == [shape(1)["id"]]
    assert counts["common_to_all_shape_count"] == 1
    assert counts["shared_across_equilibria_shape_count"] == 2
    assert all(row["unique_shape_count"] == 3 and row["duplicate_measurement_rows"] == 1
               for row in counts["per_case"])


def test_posthoc_transfer_preserves_worsening_missing_and_ineligible_destinations():
    panels = [reference_panel("a", 1, [(1, measurement(8., 3.)), (2, measurement(9.))]),
              reference_panel("a", 2, [(1, measurement(8., 3.))]),
              reference_panel("b", 1, [(1, measurement(12., 1.)), (3, measurement(7.))]),
              reference_panel("c", 1, [(1, measurement(6., 4., eligible=False))]),
              reference_panel("d", 1, [(4, measurement(11.))])]
    transfer = report.panel_diagnostics(panels)["posthoc_transfer"]
    choices = {row["case_id"]: row for row in transfer["source_choices"]}
    assert choices["a"]["id"] == shape(1)["id"]
    assert choices["a"]["peak"] == 8. and choices["a"]["A2"] == 3.
    assert len(choices["a"]["observations"]) == 2
    assert choices["a"]["selected_observation"]["panel_sha256"] == "panel-a-1"
    rows = {row["destination_case_id"]: row for row in transfer["rows"]
            if row["source_case_id"] == "a"}
    assert rows["b"]["status"] == "eligible"
    assert rows["b"]["peak_reduction_pct"] == pytest.approx(-20.)
    assert rows["b"]["A2_change_pct"] == pytest.approx(-50.)
    assert rows["c"]["status"] == "ineligible" and rows["c"]["peak"] == 6.
    assert rows["c"]["peak_reduction_pct"] is None
    assert rows["d"]["status"] == "unavailable" and rows["d"]["measurement"] is None
    assert all(row["status"] == "baseline" and row["peak_reduction_pct"] == 0.
               for row in transfer["rows"] if row["source_case_id"] == "d")
    assert transfer["new_physics_evaluations"] == 0 and transfer["posthoc"]
    assert not transfer["held_out_validation_claim"] and not transfer["global_optimum_claim"]
    assert not transfer["live_checkpoint_selection"]
    text = report.panel_diagnostics_text(report.panel_diagnostics(panels))
    assert "사후 입력 간 재사용" in text and "holdout" in text
    assert "미관측" in text and "부적격" in text


def test_posthoc_transfer_retains_wall_on_exact_tie_and_checks_duplicate_observations():
    first = reference_panel("a", 1, [(1, measurement(10.))])
    repeat = reference_panel("a", 2, [(1, measurement(10.))])
    result = report.panel_diagnostics([first, repeat])
    assert result["posthoc_transfer"]["source_choices"][0]["id"] == report.BASELINE_ID
    repeat["rows"][0]["fine"] = measurement(8.)
    with pytest.raises(ValueError, match="duplicate panel observation"):
        report.panel_diagnostics([first, repeat])
    repeat["rows"][0]["fine"] = measurement(10., eligible=False)
    with pytest.raises(ValueError, match="duplicate panel observation"):
        report.panel_diagnostics([first, repeat])


def test_all_twenty_seven_terminal_audited_arms_are_required():
    frozen = dict(cases=[dict(case_id=c) for c in ("a", "b", "c")],
                  seeds=list(report.SEEDS), policies=list(report.POLICIES))
    rows = [dict(case_id=c["case_id"], seed=s, policy=p, status="completed")
            for c in frozen["cases"] for s in report.SEEDS for p in report.POLICIES]
    report.require_terminal_audit(frozen, dict(rows=rows))
    with pytest.raises(ValueError, match="27"):
        report.require_terminal_audit(frozen, dict(rows=rows[:-1]))
    rows[0]["status"] = "not_run"
    with pytest.raises(ValueError, match="terminal"):
        report.require_terminal_audit(frozen, dict(rows=rows))


def test_group_summary_counts_three_equilibria_not_nine_independent_discharges():
    arms = []
    for case in ("a", "b", "c"):
        for seed in report.SEEDS:
            for policy in report.POLICIES:
                row = arm([baseline(), score(1, 20., 8.)], policy)
                row.update(case_id=case, seed=seed)
                arms.append(row)
    summary = report.summarize_live(arms)
    assert summary["equilibrium_count"] == 3
    assert summary["seed_repeats_per_equilibrium"] == 3
    assert summary["independent_discharge_count_claim"] is False
    assert len(summary["paired_comparisons"]) == 81
