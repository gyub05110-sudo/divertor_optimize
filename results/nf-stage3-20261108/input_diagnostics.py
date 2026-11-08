"""Explain terminal Stage 3 registration and loader failures without preparing a heat calculation."""
import argparse
import json
from pathlib import Path
import sys
import traceback

import numpy as np
from shapely.geometry import Point, Polygon

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from d3d_optimization.campaign_state import public_save
from d3d_optimization.user_input import load_user_input
from d3d_optimization.user_output import file_hash
from optimize_divertor import serial
from verification.stage3_campaign import read, require, verified_receipt


def point_envelope(points, field):
    points = np.asarray(points).reshape(-1, 2)
    low = np.array((field.r[0], field.z[0]))
    high = np.array((field.r[-1], field.z[-1]))
    excess = np.maximum.reduce([low - points, points - high, np.zeros_like(points)])
    outside = np.flatnonzero((excess > 0).any(axis=1))
    return dict(count=len(points), min_RZ_m=points.min(axis=0), max_RZ_m=points.max(axis=0),
                outside_count=len(outside), max_excess_RZ_m=excess.max(axis=0),
                outside_examples=[dict(index=int(i), RZ_m=points[i]) for i in outside[:5]])


def diagnose(case):
    """Run only the normal input loader. Read the actual rejected query and converged root from traceback locals."""
    path=Path(case["path"])
    require(file_hash(path)==case["sha256"], "registered source bytes changed")
    try:
        load_user_input(path, case["time_ms"])
    except (ValueError, OSError) as failure:
        chain=list(traceback.walk_tb(failure.__traceback__))
        project=[(f,n) for f,n in chain if Path(f.f_code.co_filename).resolve().is_relative_to(ROOT/"d3d_optimization")]
        frames={f.f_code.co_name:f.f_locals for f,n in project}
        initial=frames["load_user_input"]
        source=[dict(file=Path(f.f_code.co_filename).name, function=f.f_code.co_name, line=n) for f,n in project]
        out=dict(reproduced_reason=str(failure), failure_type=type(failure).__name__, call_path=source,
                 stage=source[-1]["function"], classification="input_validation")
        field=initial.get("field")
        if field is None:
            return out
        lcfs,axis=initial["lcfs"],initial["axis"]
        polygon=Polygon(lcfs)
        out.update(grid=dict(min_RZ_m=[field.r[0],field.z[0]], max_RZ_m=[field.r[-1],field.z[-1]],
                             points_RZ=[len(field.r),len(field.z)]),
                   magnetic_axis_RZ_m=axis, LCFS=point_envelope(lcfs,field), wall=point_envelope(initial["wall"],field),
                   polygon_valid=polygon.is_valid, LCFS_contains_axis=polygon.contains(Point(axis)))
        if "separatrix" in frames:
            local=frames["separatrix"]
            measurement=dict(LCFS_lowest_RZ_m=local["guess"],distance_limit_m=.15,gradient_limit=1e-7)
            if "solution" in local:
                solution=local["solution"]
                displacement=float(np.linalg.norm(solution.x-local["guess"]))
                residual=float(np.linalg.norm(solution.fun))
                measurement.update(root_RZ_m=solution.x,root_solver_success=bool(solution.success),
                                   root_solver_message=str(solution.message),root_evaluations=int(solution.nfev),
                                   root_distance_from_guess_m=displacement,root_gradient_norm=residual,
                                   distance_passed=displacement<=.15,gradient_passed=residual<=1e-7)
            if "hessian" in local:
                measurement["hessian_determinant"]=float(np.linalg.det(local["hessian"]))
            out["lower_X_search"]=measurement
        if source[-1]["function"]=="checked":
            query=frames["checked"]["p"]
            out["rejected_field_query"]=point_envelope(query,field)
            if out["LCFS"]["outside_count"]:
                label="LCFS_vertices_outside_grid"
            elif out["wall"]["outside_count"]:
                label="fixed_material_vertices_outside_grid"
            elif any(s["function"]=="<lambda>" and s["file"]=="case.py" for s in source) and "solution" not in frames.get("separatrix",{}):
                label="lower_X_root_trial_outside_grid"
            elif "_strike_legs" in frames:
                label="separatrix_leg_query_outside_grid"
            else:
                label="field_query_outside_grid"
            out["classification"]=label
        elif "lower X-point cannot be validated" in str(failure):
            out["classification"]="lower_X_root_rejected"
        elif "stationary point is not a magnetic saddle" in str(failure):
            out["classification"]="lower_stationary_point_not_saddle"
        elif "separatrix" in str(failure):
            out["classification"]="separatrix_reconstruction_failed"
        return out
    return dict(classification="not_reproduced",reproduced_reason=None)


def run(study):
    study=Path(study).resolve()
    expected=study/"expected.json"
    plan_sha=file_hash(expected)
    require((study/"expected.sha256").read_text().strip()==plan_sha,"expected plan digest mismatch")
    plan=read(expected)
    registration=study/"inputs"/"registration.json"
    require(file_hash(registration)==plan["registration_sha256"],"input registration changed")
    for name in ("d3d_optimization/user_input.py","d3d_optimization/case.py","d3d_optimization/frozen_field.py",
                 "d3d_optimization/frozen_mapping.py","d3d_optimization/continuous_boundary.py",
                 "d3d_optimization/independent_sources.py","baseline_review/baseline_review.json",
                 "d3d_optimization/data/DIIID_geom.json","baseline_review/optimization_region_policy.json"):
        require(file_hash(ROOT/name)==plan["files"][name], "sealed loader source or material geometry changed")
    cases={c["id"]:c for c in read(registration)["cases"]}
    destination=study/"input-diagnostics.json"
    prior=read(destination) if destination.exists() else {}
    if prior:
        require(prior["expected_sha256"]==plan_sha and prior["diagnostic_sha256"]==file_hash(__file__),
                "previous diagnostic used a different plan or diagnostic script")
    previous={c["case_id"]:c for c in prior.get("cases",[])}
    evidence={}
    terminal=0
    for task in plan["tasks"]:
        folder=study/task["relative_output"]
        receipt_file=folder/"receipt.json"
        if not receipt_file.exists():
            continue
        terminal+=1
        receipt=read(receipt_file)
        excluded=receipt["status"]=="registration_excluded"
        if not excluded and (receipt["status"]!="completed" or receipt["returncode"] not in (2,None)):
            continue
        if not excluded:
            require(file_hash(folder/"result.json")==receipt["result_sha256"], "terminal result digest mismatch")
            result=read(folder/"result.json")
            if result["status"]!="input_error":
                continue
        verified_receipt(folder,task,plan)
        reason=receipt["detail"] if excluded else result["reason"]
        item=evidence.setdefault(task["case_id"],dict(excluded=excluded,reasons=set(),receipts=[]))
        require(excluded==item["excluded"],"registration and input_error mixed in same discharge")
        item["reasons"].add(reason)
        item["receipts"].append(dict(seed=task["seed"],arm=task["arm"],receipt_sha256=file_hash(receipt_file)))
    output=[]
    for identity,records in evidence.items():
        case=cases[identity]
        reason=sorted(records["reasons"])
        if identity in previous:
            row=previous[identity]
            require(reason==row["observed_reasons"] and row["input_sha256"]==case["sha256"],"structural selection changed")
        else:
            details=(dict(classification="registration_"+case["status"],exact_frame_count=case.get("exact_frame_count"),
                          reproduced_reason=None) if records["excluded"] else diagnose(case))
            row=dict(case_id=identity,rank=case["rank"],input_sha256=case["sha256"],time_ms=case["time_ms"],
                     registration_status=case["status"],observed_reasons=reason,diagnosis=details,
                     loader_reason_matches=details["reproduced_reason"] in reason if not records["excluded"] else None)
        row["source_receipts"]=records["receipts"]
        output.append(row)
    report=dict(schema="stage3-input-domain-diagnostic.v1",expected_sha256=plan_sha,
                diagnostic_sha256=file_hash(__file__),terminal_slots=terminal,scheduled_slots=len(plan["tasks"]),
                all_slots_terminal=terminal==len(plan["tasks"]),cases=output,
                scope="Post-run CPU diagnosis of completed input_error and registration-exclusion records only. "
                      "One input replay per distinct input-error file. Registration exclusions use sealed metadata only. "
                      "The checks establish support by this loader and fixed lower-divertor model; physical infeasibility "
                      "of alternative geometries or equilibrium topologies is untested.")
    public_save(destination,serial(report))
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path)
    report=run(parser.parse_args().study)
    print(f"{report['terminal_slots']}/{report['scheduled_slots']} terminal slots, {len(report['cases'])} structural inputs")
    for c in report["cases"]:
        print(c["case_id"],c["diagnosis"]["classification"],c["loader_reason_matches"])

