"""Independent Stage 3 end-to-end clock, selection and denominator checks."""
from copy import deepcopy
import json
import math

import pytest

from d3d_optimization.campaign_geometry import Candidate
from optimize_divertor import COARSE, FINE
from verification import stage3_report as report


def measurement(q, eligible):
    return dict(peak_MW_m2=q, physics_measurable=True, candidate_eligible=eligible, A2_m2=1.0)


def journal(peaks=(10.2, 9.5, 9., 1.), direct=False):
    events = []
    for i, (q, end) in enumerate(zip(peaks, (195., 220., 598., 600.5))):
        wall = Candidate("OMP", 20, (0.,)*12 + (.01*i,)*4)
        row = dict(index=i, candidate=wall.record(), id=wall.identity, duration_s=.25,
                   within_budget=True, status="completed")
        start = end - 1 if i < 3 else 599.1
        events.append(dict(row, stage="proposal", start_s=start, end_s=start+.25))
        if not direct:
            events.append(dict(row, stage="coarse", eligible=True, peak_q=q/3, start_s=start+.25, end_s=start+.50))
        begin = start + (.25 if direct else .5)
        events.append(dict(row, stage="fine", eligible=True, peak_q=q, start_s=begin, end_s=end,
                           duration_s=end-begin, within_budget=end<=600, status="completed" if end<=600 else "late",
                           incumbent_updated=i in (0,1,2) and q < min(peaks[:i], default=math.inf)))
    best = min((e for e in events if e["stage"]=="fine" and e["within_budget"]), key=lambda e:e["peak_q"])
    data = dict(status="qualified", arguments=dict(family="full" if direct else "floor-only", radius_mm=20, budget_seconds=600,
                seed=20711, backend="cpu", time_ms=1400., method="direct" if direct else "screened"),
                reference_SI=dict(power_W=2e6), levels=dict(coarse=COARSE, fine=FINE),
                timings=dict(baseline_ready_s=130.,search_finished_s=601.,elapsed_before_result_write_s=604.,
                             deadline_overrun_s=1.),
                fixed_shelf=dict(admissible=True),
                baseline=dict(measurement(20.,False), within_budget=True, reference_only=True),
                selected=measurement(best["peak_q"]*2,True), candidate=dict(id=best["id"],**best["candidate"]),
                search=dict(counts=dict(proposals=4, coarse=0 if direct else 4, fine=4, eligible_fine=3,
                                        eligible_coarse=0 if direct else 4, failed=0, late=1),
                            termination="deadline", first_feasible_elapsed_s=195.))
    return data,events


def test_end_to_end_checkpoints_first_hit_and_strict_on_time_selection():
    data, events = journal()
    row = report.audit_run(data, events)
    snap = row["checkpoints"]
    assert snap["200"]["reduction_pct"] == pytest.approx(-2.)
    assert snap["400"]["reduction_pct"] == pytest.approx(5.)
    assert snap["600"]["q_native"] == 9.
    assert row["first_5pct"]["observed_s"] == 220.
    assert row["first_5pct"]["censor_s"] is None
    assert row["late_evaluations"] == 1


def test_no_threshold_and_unstarted_registration_have_different_exposure():
    data, events = journal((11.,10.8,10.5,.9))
    out = report.audit_run(data,events)
    assert out["first_5pct"] == dict(observed_s=None, right_censored=True, censor_s=600., reason="deadline_without_5pct")
    assert out["checkpoints"]["600"]["reduction_pct"] < 0
    absent = report.no_run("registration_excluded", "unsupported grid")
    assert absent["first_5pct"]["censor_s"] is None
    assert not absent["first_5pct"]["right_censored"]


@pytest.mark.parametrize("change", ["deadline", "geometry_id", "winner", "export", "fine_count", "coarse_score", "clock"])
def test_mismatched_native_evidence_is_rejected(change):
    data, events = journal()
    e = next(e for e in events if e["stage"]=="fine")
    if change=="deadline": events[-1]["within_budget"] = True
    if change=="geometry_id": e["id"] = "wrong"
    if change=="winner": data["candidate"] = dict(id=e["id"], **e["candidate"])
    if change=="export": data["selected"]["peak_MW_m2"] = 7.
    if change=="fine_count": data["search"]["counts"]["fine"] = 2
    if change=="coarse_score": e["stage"] = "coarse"
    if change=="clock": e["end_s"] = e["start_s"]-1
    with pytest.raises(ValueError):
        report.audit_run(data,events)


def test_equal_fine_scores_keep_first_completion_not_proposal_or_id_order():
    data, events = journal((9.,9.,9.,1.))
    out = report.audit_run(data,events)
    assert out["checkpoints"]["600"]["id"] == events[2]["id"]


def test_pairing_within_seed_and_shot_even_with_failed_comparators():
    data,events = journal()
    ready = report.audit_run(data,events)
    failed = report.no_run("registration_excluded","unsupported")
    rows = []
    for shot in ("a","b"):
        for seed in (20711,20729):
            for arm in report.ARMS:
                metrics=deepcopy(ready if shot=="a" else failed)
                if shot=="a" and seed==20711 and arm=="full-screened":
                    metrics["checkpoints"]["600"]["q_native"] = math.nextafter(9.,0.)
                    metrics["nonfloor_max_displacement_m"] = .02
                    metrics["has_nonfloor_coefficients"] = True
                if shot=="a" and seed==20729 and arm=="floor-screened":
                    metrics = deepcopy(failed)
                rows.append(dict(case_id=shot, seed=seed, arm=arm, **metrics))
    pairs=report.paired(rows)
    selected=[p for p in pairs if p["arm"]=="full-screened" and p["comparator"]=="floor-screened" and p["cutoff_s"]==600]
    assert len(selected)==4
    assert sum(p["signed_outcome"]=="arm_lower" for p in selected)==1
    assert sum(p["signed_outcome"]=="unavailable" for p in selected)==3
    assert sum(p["replay_trigger"] for p in selected)==1
    by_shot=report.per_shot(rows)
    a=next(v for v in by_shot if v["case_id"]=="a" and v["arm"]=="floor-screened")
    assert a["scheduled_seeds"]==2 and a["qualified_seeds"]==1 and a["reached_5pct_seeds"]==1



def test_registered_pair_statistics_have_distinct_full_and_speed_denominators():
    data,events=journal()
    row=report.audit_run(data,events)
    rows=[dict(deepcopy(row),case_id="a",seed=20711,arm=arm) for arm in report.ARMS]
    for r,q in zip(rows,(10.,8.,9.)):
        r["checkpoints"]["600"]["q_native"]=q
    pairs=[p for p in report.paired(rows) if p["cutoff_s"]==600]
    shape=next(p for p in pairs if p["arm"]=="full-screened" and p["comparator"]=="floor-screened")
    speed=next(p for p in pairs if p["arm"]=="full-screened" and p["comparator"]=="full-direct")
    assert shape["registered_gap_pct"]==pytest.approx(25.)
    assert speed["registered_gap_pct"]==pytest.approx(100/9)
    assert shape["replay_trigger"] and not shape["nonfloor_changed"]
    assert shape["practical_outcome"]=="arm_lower_over_1pct"


def test_denominators_use_files_and_seed_counts_including_structural_exclusions():
    data,events=journal()
    good=report.audit_run(data,events)
    fail=report.no_run("registration_excluded","unsupported")
    rows=[dict(case_id=c,seed=i,arm=a,**deepcopy(good if c=="a" and i in (1,2) else fail))
          for c in ("a","b") for i in (1,2,3) for a in report.ARMS]
    grouped=report.per_shot(rows)
    stats=report.by_arm(grouped,rows)
    assert stats["floor-screened"]["scheduled_shots"]==2
    assert stats["floor-screened"]["scheduled_shot_seed_pairs"]==6
    assert stats["floor-screened"]["shots_improving_any_seed"]==1
    assert stats["floor-screened"]["shots_improving_all_seeds"]==0
    assert stats["floor-screened"]["end_to_end_improved_pct"]==pytest.approx(100*2/6)


def test_native_export_checks_use_saved_boundary_and_reconstructed_geometry(tmp_path):
    from verification.test_user_shape_results import fixture
    path=tmp_path/"baseline"
    path.mkdir()
    base=fixture(path)
    data=dict(reference_SI=dict(power_W=2e6),baseline=base,selected=None)
    row={}
    report.check_exports(tmp_path,data,row)
    assert row["nonfloor_max_displacement_m"] is None
    base["peak_MW_m2"]+=.1
    with pytest.raises(ValueError,match="peak"):
        report.check_exports(tmp_path,data,row)


def test_incomplete_report_uses_separate_destination_and_never_completes(tmp_path):
    cases=[]
    from verification.stage3_campaign import make_plan
    for i in range(12):
        path=tmp_path/f"{i}.parquet"
        path.write_text(str(i))
        cases.append(dict(id=str(i),rank=i+1,time_ms=1400.,status="ready",path=str(path),sha256=report.file_hash(path)))
    study=tmp_path/"study"
    (study/"inputs").mkdir(parents=True)
    reg=study/"inputs"/"registration.json"
    reg.write_text(json.dumps(dict(cases=cases)))
    plan=make_plan(study,files=[],revision="fixture")
    (study/"expected.json").write_text(json.dumps(plan))
    (study/"expected.sha256").write_text(report.file_hash(study/"expected.json"))
    with pytest.raises(ValueError,match="receipt"):
        report.run(study)
    out=report.run(study,allow_incomplete=True)
    assert not out["complete"] and len(out["runs"])==108
    assert not (study/"stage3_report.json").exists()
    assert (study/"interim"/"stage3_report.json").exists()



def test_recovered_terminal_cli_does_not_acquire_process_time_or_returncode():
    data,events=journal()
    row=report.audit_run(data,events)
    report.check_process(dict(timing_source="recovered_cli_result",returncode=None,elapsed_s=None),row)
    with pytest.raises(ValueError):
        report.check_process(dict(returncode=None,elapsed_s=None),row)


def test_unmeasurable_but_exported_failed_baseline_does_not_count_as_measured(tmp_path):
    from verification.test_user_shape_results import fixture
    path=tmp_path/"baseline"
    path.mkdir()
    raw=fixture(path)
    raw.update(physics_measurable=False,errors=["unresolved_heat_peak"],peak_MW_m2=None)
    report.check_exports(tmp_path,dict(reference_SI=dict(power_W=2e6),baseline=raw),{})
    row=report.no_run("calculation_error","unmeasurable")
    row.update(case_id="a",baseline_peak_MW_m2=2.,baseline_measurable=False)
    summary=report.case_flow({"a":dict(id="a",rank=1,sha256="x",status="ready")},[row])[0]
    assert not summary["measured_baseline_any"]


def test_checkpoint_states_do_not_leak_future_threshold_or_an_unready_reference():
    data,events=journal()
    first=[e for e in events if e.get("index")==0 and e["stage"]=="fine"][0]
    first.update(start_s=201.,end_s=201.25,duration_s=.25)
    data["search"]["first_feasible_elapsed_s"]=201.25
    snaps=report.audit_run(data,events)["checkpoints"]
    assert snaps["200"]["availability"]=="reference_only"
    assert snaps["200"]["first_5pct"]["observed_s"] is None
    assert snaps["200"]["first_5pct"]["censor_s"]==200.
    assert snaps["400"]["first_5pct"]["observed_s"]==220.
    data["timings"]["baseline_ready_s"]=200.5
    snaps=report.audit_run(data,events)["checkpoints"]
    assert snaps["200"]["availability"]=="reference_pending"


def test_nonadmissible_measurable_wall_is_separate_from_observed_miss():
    data,_=journal()
    data.update(status="constraint_excluded",candidate=None,selected=None,fixed_shelf={"admissible":False})
    del data["search"]
    data["timings"].update(search_finished_s=150.,elapsed_before_result_write_s=151.,deadline_overrun_s=0.)
    out=report.audit_run(data,[])
    assert out["baseline_measurable"]
    assert out["first_5pct"]["observed_s"] is None and out["first_5pct"]["censor_s"] is None
    assert out["checkpoints"]["600"]["availability"]=="constraint_excluded"
    assert report.no_run("registration_excluded","missing")["checkpoints"]["600"]["availability"]=="registration_excluded"
