"""Fixed-wall Stage 3 replay inventory, eligibility, and persistence decisions."""
import json

import pytest

from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.user_output import file_hash
from verification import stage3_replay as replay


def emit(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj))


def test_freeze_deduplicates_fixed_walls_and_retains_failed_comparator(tmp_path):
    shot="one"
    floor=Candidate("OMP",20,(0.,)*12 + (.01,)*4)
    full=Candidate("OMP",20,(.002,)*12 + (.01,)*4)
    base=dict(inputs=dict(time_ms=1400.,equilibrium_sha256="sha"),reference_SI={"power_W":1e6},
              arguments=dict(radius_mm=20,backend="cpu"),software=dict(code_sha256={}),
              baseline=dict(peak_MW_m2=10.,A2_m2=.5))
    rows=[]
    pairs=[]
    for seed in (1,2):
        for arm,c,q in (("floor-screened",floor,9.),("full-screened",full,8.),("full-direct",None,None)):
            folder=f"runs/{shot}/{seed}/{arm}"
            path=tmp_path/folder/"result.json"
            record=dict(id=c.identity,**c.record()) if c else None
            data=dict(base,candidate=record,selected=dict(peak_MW_m2=q,A2_m2=.6) if c else None)
            emit(path,data)
            rows.append(dict(case_id=shot,seed=seed,arm=arm,relative_output=folder,result_sha256=file_hash(path),
                             status="qualified" if c else "no_design",candidate=record,
                             checkpoints={"600":dict(q_native=q,id=c.identity if c else None)}))
        for arm in ("full-screened","full-direct"):
            pairs.append(dict(case_id=shot,seed=seed,arm=arm,comparator="floor-screened",cutoff_s=600,
                              replay_trigger=arm=="full-screened",signed_outcome="arm_lower" if arm=="full-screened" else "unavailable",
                              nonfloor_changed=True if arm=="full-screened" else None))
    protocol=replay.plan_case(tmp_path,shot,rows,pairs,"frozen")
    assert len(protocol["pairs"])==4
    assert len(protocol["candidates"])==2
    arm=protocol["pairs"][0]
    assert arm["full_id"]==full.identity and arm["ablation_id"]==floor.identity and arm["floor_id"]==floor.identity
    assert protocol["levels"]["cubic1280"][1]==1280
    assert set(protocol["levels"])=={"cubic1280","quintic1280","cubic640"}
    (tmp_path/rows[0]["relative_output"]/"result.json").write_text("{}")
    with pytest.raises(ValueError,match="digest"):
        replay.plan_case(tmp_path,shot,rows,pairs,"frozen")


def test_replay_reports_reversal_failed_ablation_and_native_reproduction_separately():
    roles=dict(replay_trigger=True,full_id="f",floor_id="c",ablation_id="a",
               full_native_q=8.,floor_native_q=9.,full_native_A2_m2=1.,floor_native_A2_m2=1.,
               baseline_native_q=10.,baseline_native_A2_m2=.5)
    def score(q,eligible=True):
        return dict(physics_measurable=True,candidate_eligible=eligible,peak_MW_m2=q,A2_m2=1.,q_native=q)
    def grid(f,c,a,eligible=True):
        return dict(status="measured",fixed_shelf=dict(admissible=True),baseline=dict(score(10.,False),A2_m2=.5),
                    candidates=dict(f=score(f),c=score(c),a=score(a,eligible)))
    result=replay.adjudicate(roles,dict(cubic1280=grid(8.,9.,8.5),cubic640=grid(9.2,9.,8.),
                                        quintic1280=grid(8.,9.,7.9,False)))
    assert result["native_reproduction_passed"]
    assert result["full_lower_than_floor_all_levels"] is False
    assert result["nonfloor_benefit_all_levels"] is None
    assert result["by_level"]["cubic1280"]["ablation_minus_full_q"]==.5
    assert result["by_level"]["quintic1280"]["ablation_minus_full_q"] is None
    assert not replay.adjudicate(roles,dict(cubic1280=grid(7.,9.,8.5)))["native_reproduction_passed"]


def test_floor_only_full_winner_has_zero_direct_effect_when_replay_is_admissible():
    roles=dict(replay_trigger=True,full_id="f",floor_id="c",ablation_id="f",full_native_q=8.,floor_native_q=9.,
               full_native_A2_m2=1.,floor_native_A2_m2=1.,baseline_native_q=10.,baseline_native_A2_m2=1.)
    m=lambda q: dict(physics_measurable=True,candidate_eligible=True,q_native=q,peak_MW_m2=q,A2_m2=1.)
    level=dict(status="measured",fixed_shelf=dict(admissible=True),baseline=m(10.),candidates=dict(f=m(8.),c=m(9.)))
    assessed=replay.adjudicate(roles,{n:level for n in replay.LEVELS})
    assert assessed["full_lower_than_floor_all_levels"] and not assessed["nonfloor_benefit_all_levels"]
    assert assessed["by_level"]["cubic1280"]["ablation_minus_full_q"]==0.


def test_replay_cache_preserves_rejected_geometries_and_skips_second_prepare(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from dataclasses import asdict
    from d3d_optimization.reference_transport import ReferencePoint

    wall=Candidate("OMP",20,(0.,)*16)
    protocol=dict(common=dict(reference_SI=asdict(ReferencePoint()),radius_mm=20,backend="cpu"),
                  levels=dict(cubic1280=list(replay.LEVELS["cubic1280"])), candidates={wall.identity:wall.record()})
    count={"prepare":0,"evaluate":0}
    def fake_prepare(*args):
        count["prepare"]+=1
        return SimpleNamespace(baseline=dict(physics_measurable=False,errors=["failed_baseline"]))
    def fake_evaluate(*args):
        count["evaluate"]+=1
        return dict(errors=["angle_constraint"],geometry_screen=dict(passed=False))
    monkeypatch.setattr(replay,"shelf_limits",lambda *args: dict(admissible=False))
    monkeypatch.setattr(replay,"prepare_user",fake_prepare)
    monkeypatch.setattr(replay,"evaluate_user_candidate",fake_evaluate)
    monkeypatch.setattr(replay,"evaluator",lambda *args: object())
    a=replay.run_level(tmp_path/"grid","cubic1280",protocol,None)
    b=replay.run_level(tmp_path/"grid","cubic1280",protocol,None)
    assert a==b and count=={"prepare":1,"evaluate":1}
    assert a["candidates"][wall.identity]["errors"]==["angle_constraint"]
    assert a["candidates"][wall.identity]["candidate_eligible"] is False


def test_effect_decomposition_reports_same_floor_and_changed_floor_percentages():
    roles=dict(replay_trigger=True,full_id="f",floor_id="c",ablation_id="a",full_native_q=8.,floor_native_q=9.,
               full_native_A2_m2=1.,floor_native_A2_m2=1.,baseline_native_q=10.,baseline_native_A2_m2=1.)
    m=lambda q: dict(physics_measurable=True,candidate_eligible=True,q_native=q,peak_MW_m2=q,A2_m2=1.)
    level=dict(status="measured",fixed_shelf=dict(admissible=True),baseline=m(10.),candidates=dict(f=m(8.),a=m(8.75),c=m(9.)))
    compare=replay.adjudicate(roles,{n:level for n in replay.LEVELS})["by_level"]["cubic1280"]
    assert compare["nonfloor_reduction_vs_ablated_pct"]==pytest.approx(100*.75/8.75)
    assert compare["changed_floor_reduction_vs_floor_winner_pct"]==pytest.approx(100*.25/9)
    assert compare["family_gap_pct"]==pytest.approx(12.5)
