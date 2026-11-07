"""Baseline diagnostics use retained native results and never admit a design."""
from copy import deepcopy
import importlib
import importlib.util
import json
from types import SimpleNamespace

import numpy as np
import pytest

from d3d_optimization.frozen_mapping import HeatMap


def test_additive_diagnostic_entry_point_exists():
    assert importlib.util.find_spec("verification.nf_prospective_baseline_diagnostic"), "missing diagnostic"


@pytest.fixture
def diag():
    return importlib.import_module("verification.nf_prospective_baseline_diagnostic")


WALL = np.array([[1.,0.],[2.,0.],[2.,1.],[1.,0.]])
WIDTHS = (1,2,3,5,10)
MESSAGE = "ValueError: continuous-boundary baseline is infeasible"


def direct_maps(missing=(0., 0., .0005, .002, .02), angle=2.):
    surfaces = [dict(name="shelf",segment_index=0,original_segment=0,angle_floor_deg=1.5),
                dict(name="floor",segment_index=1,original_segment=1,angle_floor_deg=None)]
    maps, sources = [], []
    for width, loss in zip(WIDTHS, missing):
        weights=np.array([.4,.6-loss,loss])
        sources.append(SimpleNamespace(width=width/1000,weights=weights,x=np.array([.001,.002,.003])))
        maps.append(HeatMap(points=np.array([[1.9,0.],[2.,.2],[np.nan,np.nan]]), segment=np.array([0,1,-1]),
                            q=np.array([1.,2.,0.]), angles=np.array([angle,4.,np.nan]),
                            segment_power=np.array([.4,.6-loss,0.]), areas=np.ones(3),
                            mapped_power=1.-loss, unmapped_power=loss, balance_error=0.,
                            hit_power=np.array([.4,.6-loss,0.])))
    parts=dict(shelf=np.array([0]),floor=np.array([1]))
    angles={}
    for s in surfaces:
        floor=s["angle_floor_deg"]
        a=maps[0].angle_application(floor or 0.,parts[s["name"]])
        a["angle_floor_deg"]=floor
        if floor is None: a.update(status="diagnostic_only_no_angle_floor",angle_constraint_passed=None,applied_hit_count=0)
        angles[s["name"]]=a
    result=dict(geometry_valid=True,feasible=max(missing)<=.001 and angle>=1.5,penalty=0.,
                wall=WALL.copy(),original_segment=np.arange(3),surface_segments=parts,angles=angles,
                missing=max(missing),power_balance_error=0.,maps=maps,
                surface_power={s["name"]:[float(m.segment_power[parts[s["name"]]].sum()) for m in maps] for s in surfaces})
    return SimpleNamespace(sources=sources,surfaces=surfaces,dimension=8),result


def test_cross_width_guard_and_prescribed_direct_map_are_independent(diag):
    geom,data = direct_maps()
    s = diag.direct_report(geom, data, WALL)
    assert s["widths_mm"] == list(WIDTHS) and s["failed_predicates"] == ["maximum_unmapped_power"]
    assert s["prescribed_direct"]["direct_predicates_passed"] and s["prescribed_direct"]["width_mm"] == 3
    assert s["prescribed_direct"]["unmapped_power"] == .0005
    assert s["thresholds"]["maximum_unmapped_power"] == .001
    assert s["power_units"] == "fraction of assigned source power"
    assert [row["width_mm"] for row in s["by_width"] if not row["predicates"]["unmapped_power"]] == [5,10]
    assert s["prescribed_direct"]["surfaces"]["shelf"]["minimum_hit"]["nearest_corner"]["index"] == 1
    assert all(row["predicates"]["illuminated"] for row in s["by_width"])
    assert "thermal" not in s and np.isnan(data["maps"][0].points[-1,0])


def test_inactive_angle_illumination_and_balance_follow_native_scope(diag):
    geom,data=direct_maps(missing=(0.,)*5,angle=.4)
    s=diag.direct_report(geom,data,WALL)
    assert s["failed_predicates"] == ["powered_surface_angles"]
    assert s["prescribed_direct"]["surfaces"]["shelf"]["minimum_hit"]["donor"] == 0
    assert s["prescribed_direct"]["surfaces"]["floor"]["angle_floor_deg"] is None
    for m in data["maps"]:
        m.segment[:2] = 2
        m.segment_power[:] = [0.,0.,1.]
        m.balance_error=.00011
    data["power_balance_error"]=.00011
    data["surface_power"]={n:[0.]*5 for n in data["surface_power"]}
    data["angles"]["shelf"].update(status="target_unilluminated",angle_constraint_passed=None,minimum_angle_deg=None)
    s=diag.direct_report(geom,data,WALL)
    assert set(s["failed_predicates"]) == {"any_design_surface_illuminated","maximum_balance_error"}
    assert s["prescribed_direct"]["surfaces"]["shelf"]["minimum_hit"] is None


def test_native_constructor_exception_and_unreached_transport_are_preserved(diag,monkeypatch):
    from d3d_optimization.continuous_boundary import BoundaryAwareEvaluator, BoundaryGuard
    from d3d_optimization.multisurface import MultiSurfaceEvaluator
    geom,data=direct_maps()
    base=object.__new__(MultiSurfaceEvaluator)
    base.__dict__.update(vars(geom))
    monkeypatch.setattr(MultiSurfaceEvaluator,"evaluate",lambda *args:dict(data))
    monkeypatch.setattr(BoundaryGuard,"__init__",lambda *args:None)
    monkeypatch.setattr(BoundaryGuard,"witness",lambda *args:dict(resolved=True,angle_passed=True,angle_deg=3.,
                        surface="shelf",split_segment=0,point=[1.9,0.],q=[2.]*5,zero_power_limit=True))
    observer=diag.NativeObserver(WALL,.003)
    with observer.observing(), pytest.raises(ValueError,match="continuous-boundary baseline is infeasible"):
        BoundaryAwareEvaluator(base,None,None,"OMP",None,.0005)
    result=observer.report()
    assert result["continuous"][-1]["native_feasible"] is False
    assert result["continuous"][-1]["witness"]["location"]["original_segment"] == 0
    assert result["continuous"][-1]["prescribed_direct_with_witness_passed"] is True
    assert result["thermal_paths"] == [] and result["original_transport"] == [] and not result["observer_errors"]
    basis=diag.failure_basis(result,{"type":"ValueError","message":"continuous-boundary baseline is infeasible",
        "frames":[{"file":"d3d_optimization/continuous_boundary.py","function":"__init__","line":129}]})
    assert basis["stage"]=="continuous_boundary_constructor" and basis["failed_predicates"]==["maximum_unmapped_power"]
    assert basis["prescribed_direct_with_witness_passed"] and not basis["thermal_boundary_measurement_observed"]
    assert not hasattr(base,"boundary_guard")


def test_unresolved_witness_uses_recorded_intersection_and_never_awards_3mm(diag):
    geom,data=direct_maps(missing=(0.,)*5)
    data.update(feasible=False,boundary_witness=dict(resolved=False,reason="first_hit_flux_root_not_bracketed"))
    report=diag.continuous_report(geom,data,WALL,.003,
        {"wall_segment":0,"point_RZ_m":[1.99999,0.],"path_distance_m":.2})
    assert report["failed_predicates"] == ["separatrix_witness_resolved"]
    assert report["prescribed_direct_with_witness_passed"] is False
    assert report["witness"]["native_first_intersection"]["location"]["nearest_corner"]["distance_m"] == pytest.approx(1e-5)
    assert report["witness"]["native_first_intersection"]["location"]["surface"] == "shelf"


def test_caught_deposition_failure_keeps_known_power_and_location(diag):
    def synthetic_deposition():
        deposited, unhit, prescribed, balance, maxclosure=.976,.022,1.,.002,.01
        raise FloatingPointError("periodic kernel failed power closure")
    try: synthetic_deposition()
    except FloatingPointError as exc: facts=diag.exception_evidence(exc)
    assert facts["type"] == "FloatingPointError"
    assert facts["frames"][-1]["function"] == "synthetic_deposition"
    assert facts["frames"][-1]["numeric_locals"] == dict(deposited=.976,unhit=.022,prescribed=1.,balance=.002,maxclosure=.01)


def test_final_thermal_observables_survive_an_ineligible_wall(diag):
    dep=dict(ledger=dict(prescribed=1.,deposited=.998,unhit=.002,balance_error=0.),q=np.array([2.,5.]),
             cell_power=np.array([.6,.398]),areas=np.array([.3,.4]),points=np.array([[1.5,0.],[2.,.2]]),
             sample_parent=np.array([0,1]),diagnostics={"max_donor_closure_error":0.},quadrature_integrated_power=.998)
    data=dict(geometry_valid=True,feasible=False,errors=["unhit_power"],precise_length_ratios={"shelf":1.},
              receiving_angles={"shelf":{"passed":True,"minimum_angle_deg":2.}},deposition=dep,peak_q=5.,
              wall=WALL.copy(),original_segment=np.arange(3))
    out=diag.transport_report(data,WALL,1.e6,True)
    assert out["failed_predicates"] == ["deposition_unhit_power"]
    assert out["measurement"]["peak_MW_m2"] == 5. and out["measurement"]["A2_m2"]>0
    assert out["measurement"]["eligible"] is False and out["ledger"]["unhit"]==.002


def test_thermal_receiving_angle_failure_localizes_the_actual_powered_sample(diag):
    dep=dict(ledger=dict(prescribed=1.,deposited=1.,unhit=0.,balance_error=0.),q=np.array([3.,5.,2.]),
             cell_power=np.array([.25,0.,.75]),areas=np.ones(3),points=np.array([[1.9,0.],[1.95,0.],[1.8,0.]]),
             sample_cell=np.array([0,1,2]),sample_parent=np.array([0,0,0]),diagnostics={},quadrature_integrated_power=1.)
    data=dict(geometry_valid=True,feasible=False,errors=["receiving_angle_constraint_failed"],precise_length_ratios={"shelf":1.},
              receiving_angles={"shelf":{"passed":False,"minimum_angle_deg":.8}},deposition=dep,peak_q=5.,
              surface_segments={"shelf":np.array([0])},wall=WALL.copy(),original_segment=np.arange(3))
    out=diag.transport_report(data,WALL,1e6,True,incidence=np.array([.8,.01,2.]))
    assert out["failed_predicates"]==["receiving_angles"]
    m=out["powered_receiving_angle_locations"]["shelf"]
    assert m["sample"]==0 and m["surface"]=="shelf" and m["incidence_deg"]==.8
    assert m["nearest_corner"]["index"]==1


def records():
    jobs=[]; rows=[]; raw={}
    for i,(case,policy,stage) in enumerate([("casea","coarse_adaptive","prepare_fine"),
             ("casea","fine_adaptive","prepare_fine"),("caseb","coarse_adaptive","prepare_coarse"),
             ("caseb","fine_adaptive",None)]):
        folder=f"live/{case}/{i}/{policy}"
        job=dict(case_id=case,input_root=case,seed=i,policy=policy,phase="holdout",budget_s=600,relative_dir=folder,
                 population=64 if policy=="coarse_adaptive" else 24,projected_seeds=12 if policy=="coarse_adaptive" else 8)
        jobs.append(job)
        rows.append(dict(**job,normal=stage is None,status="completed" if stage is None else "arm_failed",
                         launch_ns=i+1,job_sha256="a"*64,result_sha256="b"*64,events_sha256="c"*64,
                         process_sha256="d"*64,error=None if stage is None else MESSAGE))
        raw[folder]=[] if stage is None else [dict(stage=stage,status="failed",error=MESSAGE,start_s=1.,end_s=2.)]
    frozen=dict(phase="holdout",jobs=jobs,reference={"lambda_q_m":.003},source="OMP",radius_mm=20,
                levels={"fine":["cubic",1280,.000125,641,.00025],"coarse":["cubic",320,.0005,161,.001]},
                cases=[dict(case_id=case,root=case,hashes={"shot":"x"},time_ms=1400.) for case in ("casea","caseb")])
    return frozen,dict(phase="holdout",rows=rows,planned_arm_count=4,completed_arm_count=1),raw


def test_final_audit_registration_deduplicates_only_same_failed_stage(diag):
    frozen,audit,raw=records()
    lock=diag.select_failures(frozen,audit,raw)
    assert [(v["case_id"],v["stage"],len(v["failure_evidence"])) for v in lock["targets"]] == [
        ("casea","prepare_fine",2),("caseb","prepare_coarse",1)]
    assert not any(v["case_id"]=="casea" and v["level"]=="coarse" for v in lock["targets"])
    assert lock["registered_arms"]==4 and lock["failed_arms"]==3 and lock["completed_arms"]==1
    assert lock["targets"][0]["failure_evidence"][0]["event_index"]==0
    assert lock["targets"][0]["failure_evidence"][0]["events_sha256"]=="c"*64
    with pytest.raises(ValueError):
        diag.select_failures(frozen,dict(audit,rows=audit["rows"][:-1]),raw)


def test_successful_study_has_no_replay_and_other_failures_remain_visible(diag):
    frozen,audit,raw=records()
    for row in audit["rows"]: row.update(normal=True,status="completed",error=None)
    audit["completed_arm_count"]=4
    assert diag.select_failures(frozen,audit,{})["targets"]==[]
    row=audit["rows"][0]
    row.update(normal=False,status="arm_failed",error="GPU score failed")
    raw[row["relative_dir"]]=[dict(stage="coarse",status="failed",error=row["error"],start_s=40.,end_s=41.)]
    audit["completed_arm_count"]=3
    result=diag.select_failures(frozen,audit,raw)
    assert result["targets"]==[] and result["other_failed_arms"][0]["stage"]=="coarse"


def test_replay_requires_preexisting_sealed_registration_before_preparing(diag,tmp_path,monkeypatch):
    from verification.nf_anytime_reference import EvidenceStore
    store=EvidenceStore(tmp_path)
    with pytest.raises((ValueError,FileNotFoundError)):
        diag.replay_study(tmp_path)
    assert not (tmp_path/"audit.json").exists()
    with pytest.raises((ValueError,FileNotFoundError)):
        diag.register_study(tmp_path)
    assert not (tmp_path/"baseline-diagnostic/registration.json").exists()


def test_exact_native_prepare_arguments_and_earlier_stop_cannot_match_baseline_failure(diag,tmp_path):
    from verification.nf_anytime_reference import digest,object_digest
    from d3d_optimization.reference_transport import ReferencePoint
    from dataclasses import asdict
    wall_file=tmp_path/diag.REVIEW
    wall_file.parent.mkdir()
    wall_file.write_text(json.dumps({"wall":WALL.tolist()}))
    target=dict(case_id="example",stage="baseline",level="fine",fidelity=["cubic",1280,.000125,641,.00025],
                source="OMP",input_root=str(tmp_path),time_ms=1400.,reference=asdict(ReferencePoint()),
                input_sha256={diag.REVIEW:digest(wall_file)},original_wall_coordinates_sha256=object_digest(WALL.tolist()),
                failure_evidence=[dict(event_id="abc#1",event_error=MESSAGE)])
    calls=[]
    def failing(root,source,ref,index):
        calls.append((root,source,ref,index))
        raise ValueError("continuous-boundary baseline is infeasible")
    outcome=diag.run_target(target,prepare=failing)
    assert calls==[(tmp_path,"OMP",ReferencePoint(),2)]
    assert not outcome["all_primary_native_errors_reproduced"] and outcome["matching_primary_failure_ids"]==[]
    assert not outcome["full_thermal_deposition_returned"] and outcome["process_cpu_seconds"]>=0
    outcome=diag.run_target(dict(target,stage="prepare_fine"),prepare=failing)
    assert outcome["all_primary_native_errors_reproduced"] and outcome["native_error"]==MESSAGE


def test_sealed_all_pass_replay_makes_zero_native_calls_and_is_idempotent(diag,tmp_path,monkeypatch):
    from verification.nf_anytime_reference import digest
    plan=dict(diagnostic_source_sha256=digest(diag.__file__),registered_arms=18,completed_arms=18,failed_arms=0,
              targets=[],other_failed_arms=[])
    monkeypatch.setattr(diag,"_registered_payload",lambda _:deepcopy(plan))
    monkeypatch.setattr(diag,"run_target",lambda _:pytest.fail("all-pass study scheduled physics"))
    diag.register_study(tmp_path)
    first=diag.replay_study(tmp_path)
    assert first["original_wall_replay_count"]==0 and first["primary_denominator_unchanged"]
    assert diag.replay_study(tmp_path)==first


def test_registration_binds_every_primary_receipt_and_rejects_bytes_changed_afterward(diag,tmp_path,monkeypatch):
    from verification.nf_anytime_reference import EvidenceStore,digest
    from verification import nf_prospective_campaign as campaign
    from d3d_optimization.reference_transport import ReferencePoint
    from dataclasses import asdict
    frozen,audit,raw=records()
    live=EvidenceStore(tmp_path)
    live.write_bytes("protocol.ko.md",b"original protocol")
    for c in frozen["cases"]:
        c["root"]=str(tmp_path/"inputs"/c["case_id"])
        review=diag.Path(c["root"])/diag.REVIEW
        review.parent.mkdir(parents=True)
        review.write_text(json.dumps({"wall":WALL.tolist()}))
        c["hashes"]={diag.REVIEW:digest(review)}
    roots={c["case_id"]:c["root"] for c in frozen["cases"]}
    for j in frozen["jobs"]: j["input_root"]=roots[j["case_id"]]
    frozen.update(code_root=str(diag.ROOT),source_hashes={},original_reference_hashes={},packages={},python_version="3",
                  protocol_sha256=live.sha256("protocol.ko.md"),reference=asdict(ReferencePoint()))
    live.write("frozen.json",frozen)
    registration=live.sha256("frozen.json")
    for j,row in zip(frozen["jobs"],audit["rows"]):
        receipt=EvidenceStore(tmp_path/j["relative_dir"])
        receipt.write("job.json",dict(j,registration_sha256=registration))
        arm=receipt.root
        (arm/"events.jsonl").write_text("".join(json.dumps(e)+"\n" for e in raw[j["relative_dir"]]))
        rec=dict(events=raw[j["relative_dir"]],error=row["error"],launch_ns=row["launch_ns"],job_sha256=receipt.sha256("job.json"))
        (arm/"result.json").write_text(json.dumps(rec))
        (arm/"stdout.log").write_bytes(b"")
        p=dict(launch_ns=row["launch_ns"],exit_observed_ns=row["launch_ns"]+10**9,elapsed_s=1.,returncode=int(not row["normal"]),
               registration_sha256=registration,error=None,job_sha256=receipt.sha256("job.json"),events_sha256=digest(arm/"events.jsonl"),
               result_sha256=digest(arm/"result.json"),stdout_sha256=digest(arm/"stdout.log"))
        receipt.write("process.json",p)
        row.update({k:p[k] for k in ("job_sha256","events_sha256","result_sha256")},process_sha256=receipt.sha256("process.json"))
    live.write("audit.json",dict(audit,registration_sha256=registration))
    monkeypatch.setattr(campaign,"verify_registration",lambda *_:None)
    plan=diag.register_study(tmp_path)
    assert len(plan["targets"])==2 and plan["registered_arms"]==4
    success=frozen["jobs"][-1]["relative_dir"]+"/result.json"
    assert plan["primary_sha256"][success]==digest(tmp_path/success)
    (tmp_path/success).write_bytes(b"changed successful result")
    with pytest.raises(ValueError):
        diag.replay_study(tmp_path)
    assert not (tmp_path/"baseline-diagnostic/results").exists()
