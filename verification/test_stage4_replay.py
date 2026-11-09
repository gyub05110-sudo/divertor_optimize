"""Synthetic checks of Stage 4 fixed-wall replay planning, thresholds and restart integrity."""
from importlib import import_module, util
import json

import pytest

from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.user_output import file_hash
from verification import stage4_report as report


def replay():
    name = "verification.stage4_replay"
    assert util.find_spec(name), "Stage 4 post-budget replay is required"
    return import_module(name)


def fixture_report(study, scores, predicted="floor-only"):
    rows = []
    for arm, q in zip(report.design.ALL, scores, strict=True):
        folder = f"runs/test/30713/{arm}"
        candidate = Candidate("OMP",20,(0.,)*12 + (round(q/10000,6),)*4) if q else None
        record = dict(id=candidate.identity, **candidate.record()) if candidate else None
        row = dict(case_id="test", arm=arm, seed=30713, relative_output=folder,
                   status="qualified" if q else "no_design", candidate=record,
                   predicted_family=predicted, predicted_floor_sufficient=predicted=="floor-only",
                   abstained=predicted=="full", baseline_peak_MW_m2=100.,
                   checkpoints={"600":dict(q_native=q,peak_MW_m2=q,id=record["id"] if record else None,
                   reduction_pct=None if q is None else 100*(100-q)/100)})
        data=dict(arguments=dict(family=arm,seed=30713,radius_mm=20,backend="cpu"), candidate=record,
                  inputs=dict(time_ms=1400.), reference_SI=dict(power_W=1e6), software=dict(code_sha256={}),
                  baseline=dict(peak_MW_m2=100., A2_m2=2.), selected=dict(peak_MW_m2=q,A2_m2=3.) if q else None)
        path = study / folder / "result.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(data))
        row["result_sha256"] = file_hash(path)
        rows.append(row)
    return dict(runs=rows, pairs=[report.product_pair(rows[0],rows[1],600)], predictions=[report.prediction_case(rows)])


def test_first_seed_contrasts_have_deduped_candidates_and_ordered_best_ties(tmp_path):
    source = fixture_report(tmp_path,[96.,92.,96.,92.,92.,94.])
    source["runs"].append(dict(case_id="test", seed=30731, arm="full", status="no_design"))
    plan = replay().plan_case(tmp_path,"test",source,"audit")
    assert len(plan["candidates"]) == 2
    product, forced = plan["pairs"]
    assert product["trigger"] and forced["trigger"]
    assert (product["left"],product["right"]) == ("hotspot","full")
    assert (forced["left"],forced["right"]) == ("floor-only","floor-entry")
    assert plan["levels"]["cubic640"][1] == 640
    path=tmp_path/source["runs"][0]["relative_output"]/"result.json"
    path.write_text("{}")
    with pytest.raises(ValueError,match="digest"):
        replay().plan_case(tmp_path,"test",source,"audit")


@pytest.mark.parametrize("scores,predicted,expected", [
    ([95.1,94.9,94.9,95.,95.,95.],"floor-only",[True,False]),
    ([100.,100.,101.,100.,100.,100.],"floor-only",[False,False]),
    ([100.,100.,99.,None,95.,99.],"floor-only",[False,False]),
    ([None,91.,95.,92.,93.,94.],"full",[False,False]),
    ([None,None,None,None,None,None],None,[False,False]),
])
def test_predeclared_thresholds_and_retained_unavailable_controls(tmp_path,scores,predicted,expected):
    plan=replay().plan_case(tmp_path,"test",fixture_report(tmp_path,scores,predicted),"audit")
    assert [p["trigger"] for p in plan["pairs"]] == expected
    if not any(expected):
        folder=tmp_path/"replays/test"
        result=replay().case_replay(folder,plan)
        assert result["complete"] and result["trigger_count"]==0 and result["results"]=={}
        assert len(result["comparisons"])==2 and replay().case_replay(folder,plan)==result


def test_grid_reversal_threshold_crossing_failed_constraints_and_native_area_check(tmp_path):
    mod = replay()
    plan = mod.plan_case(tmp_path,"test",fixture_report(tmp_path,[96.,92.,96.,92.,92.,94.]),"audit")
    pair=plan["pairs"][0]
    def grid(l,r,a2=3.,ok=True):
        score=lambda q,area:dict(q_native=q,A2_m2=area,physics_measurable=True,candidate_eligible=True)
        return dict(status="measured",fixed_shelf=dict(admissible=ok),baseline=score(100.,2.),
                    candidates={pair["native"]["hotspot"]["id"]:score(l,a2), pair["native"]["full"]["id"]:score(r,3.)})
    grids=dict(cubic1280=grid(96.,92.),quintic1280=grid(93.,96.),cubic640=grid(96.,92.,ok=False))
    out=mod.adjudicate(pair,grids)
    assert out["native_reproduction_passed"]
    assert out["by_level"]["cubic1280"]["sign"]=="right_lower"
    assert out["by_level"]["quintic1280"]["sign"]=="left_lower"
    assert out["by_level"]["quintic1280"]["left_5pct"] and not out["by_level"]["quintic1280"]["right_5pct"]
    assert out["by_level"]["cubic640"]["gap_pct"] is None
    grids["cubic1280"]=grid(96.,92.,a2=2.99)
    assert not mod.adjudicate(pair,grids)["native_reproduction_passed"]


def test_case_resume_validates_baseline_and_grid_receipts_without_remeasuring(tmp_path,monkeypatch):
    mod=replay()
    plan=mod.plan_case(tmp_path,"test",fixture_report(tmp_path,[96.,92.,96.,92.,92.,94.]),"audit")
    root=tmp_path/"replays/test"
    calls=[]
    monkeypatch.setattr(mod,"load_replay_user",lambda cfg: None)
    def run_level(folder,name,plan,user):
        calls.append(name)
        folder.mkdir()
        data= dict(status="calculation_error", error="synthetic nonconvergence")
        (folder/"setup.json").write_text(json.dumps(data))
        return dict(data,baseline=None,candidates={})
    monkeypatch.setattr(mod.fixed,"run_level",run_level)
    a=mod.case_replay(root,plan)
    assert a==mod.case_replay(root,plan) and len(calls)==3 and a["trigger_count"]==2
    (root/"cubic640/setup.json").write_text("{}")
    with pytest.raises(ValueError,match="digest"):
        mod.case_replay(root,plan)


def test_no_heldout_replay_before_complete_report_and_export_matches_summary(tmp_path,monkeypatch):
    mod=replay()
    source=fixture_report(tmp_path,[None]*6,None)
    audit_file=tmp_path/"stage4_report.json"
    (tmp_path/"expected.json").write_text("{}")
    (tmp_path/"inputs").mkdir()
    (tmp_path/"inputs/registration.json").write_text(json.dumps(dict(cases=[dict(id="test")])))
    monkeypatch.setattr(mod.design,"load_seal",lambda _:dict(registration_sha256="registered"))
    monkeypatch.setattr(mod.report,"audit_jobs",lambda *args, **kwargs:source["runs"])
    source.update(complete=False, plan_sha256=file_hash(tmp_path/"expected.json"),registration_sha256="registered",
                  report_driver_sha256=file_hash(report.__file__))
    audit_file.write_text(json.dumps(source))
    with pytest.raises(ValueError,match="complete audited"):
        mod.run(tmp_path)
    assert not (tmp_path/"replays").exists()
    source["complete"]=True
    audit_file.write_text(json.dumps(source))
    mod.run(tmp_path)
    assert (tmp_path/"stage4_replay.json").read_bytes()==(tmp_path/"replays/summary.json").read_bytes()
    cases=json.loads((tmp_path/"stage4_replay.json").read_text())["cases"]
    assert len(cases)==1 and cases[0]["trigger_count"]==0 and cases[0]["first_seed_statuses"]["hotspot"]=="no_design"


def test_exact_one_percent_closed_boundary_never_triggers_despite_float_gap(tmp_path):
    plan = replay().plan_case(tmp_path,"test",fixture_report(tmp_path,[2.02,2.,2.02,2.,2.01,2.01]),"audit")
    assert [p["trigger"] for p in plan["pairs"]] == [False,False]
    assert plan["pairs"][1]["timed_gap_pct"] > 1.  # Percent display can round above one.


def test_exact_closed_replay_boundaries_preserve_at_least_five_and_within_one(tmp_path):
    mod = replay()
    plan = mod.plan_case(tmp_path, "test", fixture_report(tmp_path,[96.,92.,96.,92.,92.,94.]),"audit")
    pair = plan["pairs"][0]
    score=lambda q,area:dict(q_native=q,A2_m2=area,physics_measurable=True,candidate_eligible=True)
    def level(left,right,base):
        return dict(status="measured",fixed_shelf=dict(admissible=True),baseline=score(base,2.),
                    candidates={pair["native"]["hotspot"]["id"]:score(left,3.),pair["native"]["full"]["id"]:score(right,3.)})
    levels = mod.adjudicate(pair, dict(cubic1280=level(96.,92.,100.),cubic640=level(2.02,2.,19.),
                                      quintic1280=level(19.*.95,19.,19.)))["by_level"]
    assert levels["cubic640"]["within_1pct"] is True
    assert levels["quintic1280"]["left_5pct"] is True


def test_partial_result_with_a_verified_timely_winner_can_trigger(tmp_path):
    source = fixture_report(tmp_path, [96.,92.,96.,92.,92.,94.])
    source["runs"][0]["status"] = "partial_result"
    source["runs"][2]["status"] = "partial_result"
    source["pairs"] = [report.product_pair(source["runs"][0], source["runs"][1], 600)]
    source["predictions"] = [report.prediction_case(source["runs"])]
    plan = replay().plan_case(tmp_path, "test", source, "audit")
    assert all(p["trigger"] for p in plan["pairs"]) and len(plan["candidates"]) == 2
