"""Post-budget fixed-wall replay for strict Stage 3 full-family wins, grouped by discharge."""
import argparse
from math import isclose
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

from optimize_divertor import serial, software
from d3d_optimization.campaign_geometry import Candidate, evaluator
from d3d_optimization.campaign_state import process_lock, public_save
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.user_input import load_user_input
from d3d_optimization.user_output import file_hash
from d3d_optimization.user_physics import evaluate_user_candidate, prepare_user, shelf_limits
from verification.stage3_report import positive, read, require, REQUIRED_FILES
from verification.user_shape_replay import INPUT_KEYS, LEVELS, _measurement, _replay_input_paths


def plan_case(study,case,rows,pairs,search_report_sha256):
    eligible_runs={(r["seed"],r["arm"]):r for r in rows if r["case_id"]==case}
    inventory=[p for p in pairs if p["case_id"]==case and p["cutoff_s"]==600 and p["comparator"]=="floor-screened"]
    planned=[]
    candidates={}
    common=None
    for pair in inventory:
        seed=pair["seed"]
        full=eligible_runs[seed,pair["arm"]]
        floor=eligible_runs[seed,pair["comparator"]]
        q_full=full["checkpoints"]["600"]["q_native"]
        q_floor=floor["checkpoints"]["600"]["q_native"]
        trigger=(full["status"]=="qualified" and floor["status"]=="qualified" and q_full<q_floor)
        require(trigger==pair["replay_trigger"], "reported replay scope disagrees with strict native peak selection")
        row=dict(seed=seed,arm=full["arm"],comparator=floor["arm"],
                 contrast_scope=pair.get("contrast_scope"),signed_outcome=pair["signed_outcome"],replay_trigger=trigger,
                 full_status=full["status"],floor_status=floor["status"],full_native_q=q_full,floor_native_q=q_floor,
                 nonfloor_changed=pair["nonfloor_changed"])
        if trigger:
            for role,run in (("floor",floor),("full",full)):
                path=study/run["relative_output"]/"result.json"
                require(file_hash(path)==run["result_sha256"], "search result digest changed")
                result=read(path)
                cand=Candidate.from_record(result["candidate"])
                require(cand.identity==run["candidate"]["id"]==run["checkpoints"]["600"]["id"] and
                        run["candidate"]==result["candidate"],"winner identity changed")
                require(role!="floor" or not any(cand.vector[:12]),"nonfloor geometry in floor comparator")
                cfg=dict(inputs=result["inputs"],reference_SI=result["reference_SI"],radius_mm=result["arguments"]["radius_mm"],
                         backend=result["arguments"]["backend"],search_code_sha256=result["software"]["code_sha256"])
                require(common is None or common==cfg, "different fixed inputs or numerical software within discharge")
                common=cfg
                scale=cfg["reference_SI"]["power_W"]/1e6
                baseline=result["baseline"]
                q_base=baseline["peak_MW_m2"]/scale
                require("baseline_native_q" not in row or row["baseline_native_q"]==q_base, "paired baseline differs")
                row.update(baseline_native_q=q_base,baseline_native_A2_m2=baseline["A2_m2"])
                candidates[cand.identity]=cand.record()
                row[role+"_id"]=cand.identity
                row[role+"_native_A2_m2"]=result["selected"]["A2_m2"]
                row[role+"_result_sha256"]=run["result_sha256"]
                if role=="full":
                    ablated=Candidate("OMP",cand.radius_mm,(0.,)*12 + cand.vector[12:])
                    candidates[ablated.identity]=ablated.record()
                    row["ablation_id"]=ablated.identity
                    row["zero_nonfloor_full"]=not any(cand.vector[:12])
        planned.append(row)
    return dict(schema="nf-stage3-fixed-replay.v1",case_id=case,search_report_sha256=search_report_sha256,
                scope="Fixed original wall and two within-budget winners plus full winner with 12 nonfloor coefficients zeroed",
                replay_driver_sha256=file_hash(__file__),common=common,pairs=planned,candidates=candidates,
                levels={k:list(v) for k,v in LEVELS.items()})


def native_matches(row,q,area,reference=False):
    ok=row.get("physics_measurable") if reference else row.get("candidate_eligible") and row.get("physics_measurable")
    return bool(ok and positive(row.get("q_native")) and positive(row.get("A2_m2")) and
                isclose(row["q_native"],q,rel_tol=1e-9) and isclose(row["A2_m2"],area,rel_tol=1e-9))


def adjudicate(pair,grids):
    if not pair["replay_trigger"]:
        return dict(adjudicated=False,reason=pair["signed_outcome"])
    levels={}
    reproduction={}
    for name in LEVELS:
        grid=grids.get(name) or {}
        constraints=grid.get("fixed_shelf",{}).get("admissible") is True
        base=grid.get("baseline") or {}
        records=grid.get("candidates",{})
        full=records.get(pair["full_id"],{})
        floor=records.get(pair["floor_id"],{})
        ablated=records.get(pair["ablation_id"],{})
        reference=base.get("q_native") if base.get("physics_measurable") else None
        def eligible_q(row):
            q=row.get("q_native")
            return q if row.get("physics_measurable") and row.get("candidate_eligible") and constraints and positive(q) else None
        x,y,z=map(eligible_q,(full,floor,ablated))
        levels[name]=dict(status=grid.get("status","missing"),fixed_shelf_admissible=grid.get("fixed_shelf",{}).get("admissible"),
                          reference_q=reference,full_q=x,floor_q=y,ablation_q=z,
                          full_reduction_vs_reference_pct=100*(reference-x)/reference if reference and x else None,
                          family_gap_pct=100*(y-x)/x if x and y else None,
                          floor_minus_full_q=y-x if x and y else None,
                          full_lower_than_floor=(x<y) if x and y else None,
                          ablation_minus_full_q=z-x if x and z else None,
                          nonfloor_reduction_vs_ablated_pct=100*(z-x)/z if x and z else None,
                          nonfloor_lowers_peak=(x<z) if x and z else None,
                          changed_floor_reduction_vs_floor_winner_pct=100*(y-z)/y if z and y else None,
                          ablation_minus_floor_q=z-y if z and y else None)
        if name=="cubic1280":
            reproduction=dict(baseline=native_matches(base,pair["baseline_native_q"],pair["baseline_native_A2_m2"],True),
                              full=constraints and native_matches(full,pair["full_native_q"],pair["full_native_A2_m2"]),
                              floor=constraints and native_matches(floor,pair["floor_native_q"],pair["floor_native_A2_m2"]))
    def persists(key):
        vals=[v[key] for v in levels.values()]
        return all(vals) if all(v is not None for v in vals) and all(reproduction.values()) else None
    return dict(adjudicated=True,native_reproduction_relative_tolerance=1e-9,native_reproduction=reproduction,
                native_reproduction_passed=all(reproduction.values()),by_level=levels,
                full_lower_than_floor_all_levels=persists("full_lower_than_floor"),
                nonfloor_benefit_all_levels=persists("nonfloor_lowers_peak"))


def measure(directory,owner,raw,elapsed):
    """Publish one complete fixed wall and its metadata together. An interrupted write can be retried."""
    if "deposition" in raw:
        row=_measurement(directory,owner,raw)
        row["q_native"]=raw.get("peak_q")
    else:
        directory.mkdir()
        row=dict(physics_measurable=False,candidate_eligible=False,q_native=None,errors=raw.get("errors",[]))
        for key in ("geometry_screen","shape_constraints","precise_length_ratios"):
            if key in raw:
                row[key]=raw[key]
    row["evaluation_s"]=elapsed
    public_save(directory/"measurement.json",serial(row))
    return row


def load_measurement(directory):
    if not directory.exists():
        return None
    record=read(directory/"measurement.json")
    hashes=record.get("sha256",{})
    require(not hashes or (set(hashes)==REQUIRED_FILES and directory.name==record["directory"]),
            "invalid replay measurement inventory")
    for name,digest in hashes.items():
        require(file_hash(directory/name)==digest, "replay measurement changed at "+str(directory/name))
    return record


def run_level(folder,name,protocol,user):
    folder.mkdir(exist_ok=True)
    cfg=protocol["common"]
    resolution=tuple(protocol["levels"][name])
    reference=ReferencePoint(**cfg["reference_SI"])
    roles=["baseline",*sorted(protocol["candidates"])]
    scores={role:load_measurement(folder/role) for role in roles}
    missing=[k for k,v in scores.items() if v is None]
    setup_file=folder/"setup.json"
    if missing or not setup_file.exists():
        t=perf_counter()
        setup=dict(resolution=list(resolution))
        try:
            setup["fixed_shelf"]=shelf_limits(user,resolution[0])
            original=prepare_user(user,reference,resolution)
            solver=evaluator(original,cfg["radius_mm"],cfg["backend"])
        except (ValueError, ArithmeticError) as exc:
            setup.update(status="calculation_error",error=f"{type(exc).__name__}: {exc}",prepare_s=perf_counter()-t)
            public_save(setup_file,serial(setup))
            return dict(setup,baseline=None,candidates={})
        setup.update(status="measured",prepare_s=perf_counter()-t)
        public_save(setup_file,serial(setup))
        for role in missing:
            t=perf_counter()
            if role=="baseline":
                owner,raw=original,original.baseline
                elapsed=None # Already included in prepare_s.
            else:
                owner=solver
                try:
                    raw=evaluate_user_candidate(solver,Candidate.from_record(protocol["candidates"][role]))
                except (ValueError,ArithmeticError) as exc:
                    raw=dict(feasible=False,errors=[f"{type(exc).__name__}: {exc}"])
                elapsed=perf_counter()-t
            with TemporaryDirectory(dir=folder,prefix=".partial-") as scratch:
                staging=Path(scratch)/role
                measure(staging,owner,raw,elapsed)
                staging.rename(folder/role)
            scores[role]=load_measurement(folder/role)
            del raw
    setup=read(setup_file)
    return dict(setup,baseline=scores["baseline"],candidates={k:scores[k] for k in roles[1:]})


def case_replay(folder,protocol):
    folder.mkdir(parents=True,exist_ok=True)
    with process_lock(folder/"replay.lock"):
        protocol_file=folder/"protocol.json"
        if protocol_file.exists():
            require(protocol==read(protocol_file), "fixed replay protocol changed")
        else:
            public_save(protocol_file,protocol)
        report_file=folder/"report.json"
        saved=read(report_file) if report_file.exists() else {}
        if saved:
            require(saved["protocol_sha256"]==file_hash(protocol_file),"replay report protocol digest mismatch")
            for relative,digest in saved.get("measurement_receipts_sha256",{}).items():
                require(file_hash(folder/relative)==digest,"replay measurement receipt changed: "+relative)
        common=protocol["common"]
        grids={}
        if common:
            source=common["inputs"]
            paths=_replay_input_paths(source)
            user=load_user_input(paths["equilibrium"],source["time_ms"],paths["geometry"],paths["design_policy"])
            require(all(source[k]==user.provenance[k] for k in INPUT_KEYS), "fixed physical input changed")
            current=software()
            require(current["code_sha256"]==common["search_code_sha256"], "solver code differs from timed search")
            for name in protocol["levels"]:
                earlier=saved.get("results",{}).get(name)
                if earlier and earlier.get("status")=="calculation_error":
                    grids[name]=earlier
                else:
                    grids[name]=run_level(folder/name,name,protocol,user)
                manifests={p.relative_to(folder).as_posix():file_hash(p) for p in folder.glob("*/*/measurement.json")}
                interim=dict(complete=False,protocol_sha256=file_hash(protocol_file),results=grids,
                             measurement_receipts_sha256=manifests)
                public_save(report_file,serial(interim))
        result=dict(complete=True,case_id=protocol["case_id"],protocol_sha256=file_hash(protocol_file),
                    scope=protocol["scope"],trigger_count=sum(p["replay_trigger"] for p in protocol["pairs"]),
                    fixed_unique_geometry_count=len(protocol["candidates"]),results=grids,
                    comparisons=[dict(p,**adjudicate(p,grids)) for p in protocol["pairs"]],
                    measurement_receipts_sha256={p.relative_to(folder).as_posix():file_hash(p)
                                                 for p in folder.glob("*/*/measurement.json")})
        public_save(report_file,serial(result))
        return result


def run(study,case=None):
    from verification.stage3_campaign import load_seal
    study=Path(study).resolve()
    plan=load_seal(study)
    report_file=study/"stage3_report.json"
    source=read(report_file)
    require(source["complete"] and source["expected_sha256"]==file_hash(study/"expected.json") and
            source["registration_sha256"]==plan["registration_sha256"], "complete audited timed campaign required")
    report_digest=file_hash(report_file)
    cases=[v["id"] for v in read(study/"inputs"/"registration.json")["cases"]]
    require(case is None or case in cases,"unknown registered discharge")
    dest=study/"replays"
    for case_id in [case] if case else cases:
        protocol=plan_case(study,case_id,source["runs"],source["paired"],report_digest)
        result=case_replay(dest/case_id,protocol)
        print(case_id,result["trigger_count"],"strict wins;",result["fixed_unique_geometry_count"],"unique geometries",flush=True)
    if case is None:
        cases_done=[]
        for case_id in cases:
            path=dest/case_id/"report.json"
            row=read(path)
            require(row["complete"],"incomplete discharge replay")
            cases_done.append(dict(case_id=case_id,report_sha256=file_hash(path),trigger_count=row["trigger_count"],
                                   fixed_unique_geometry_count=row["fixed_unique_geometry_count"],comparisons=row["comparisons"]))
        public_save(dest/"summary.json",dict(complete=True,search_report_sha256=report_digest,cases=cases_done))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",required=True,type=Path)
    p.add_argument("--case",help="single registered discharge; omit to replay all and write summary.json")
    run(**vars(p.parse_args(argv)))


if __name__=="__main__":
    main()
