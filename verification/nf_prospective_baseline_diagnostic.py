"""Post-trial CPU diagnosis of original-wall failures in the sealed holdout audit.

Register case/failed native stage pairs before replay. Observers retain native
returns and exceptions from the original preparation call. They perform no new
tracing, source fits, geometry selection, optimization or alternate-width solve.
"""
import argparse
from contextlib import contextmanager, ExitStack
from datetime import datetime, timezone
from functools import wraps
import hashlib
import math
import os
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time
from unittest.mock import patch

import numpy as np

from verification import nf_anytime_final_audit as audited
from verification.nf_anytime_reference import EvidenceStore, bind_evidence, object_digest, require


ROOT = Path(__file__).resolve().parents[1]
REVIEW = "baseline_review/baseline_review.json"
STAGES = {"prepare_fine":"fine", "baseline":"fine", "prepare_coarse":"coarse"}
DEFAULT_STUDY = Path("D:/DivertorOptimization/studies/nf-prospective-search-20261107")


def _json(value):
    """Native inactive entries may contain NaN. JSON null marks an undefined observable."""
    if isinstance(value, np.ndarray):
        return _json(value.tolist())
    if isinstance(value, np.generic):
        return _json(value.item())
    if isinstance(value, dict):
        return {str(k):_json(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [_json(v) for v in value]
    return None if isinstance(value,float) and not math.isfinite(value) else value


def _location(point, k, result, wall):
    """Use the receiving chord midpoint to disambiguate the two faces at a corner."""
    point, wall = np.asarray(point,float), np.asarray(wall,float)
    if point.shape != (2,) or not np.isfinite(point).all() or k is None or k < 0:
        return None
    native_wall = np.asarray(result["wall"])
    midpoint = (native_wall[k]+native_wall[k+1])/2
    edges = np.diff(wall,axis=0)
    sq = np.einsum("ij,ij->i",edges,edges)
    projection = np.divide(np.sum((midpoint-wall[:-1])*edges,axis=1),sq,
                           where=sq>0,out=np.zeros(len(edges)))
    distances = np.linalg.norm(wall[:-1]+np.clip(projection,0.,1.)[:,None]*edges-midpoint,axis=1)
    original = int(np.argmin(distances))
    corners = np.linalg.norm(wall[[original,original+1]]-point,axis=1)
    corner = original+int(np.argmin(corners))
    if corner==len(wall)-1 and np.array_equal(wall[0],wall[-1]):
        corner=0
    surface=next((name for name,ids in result.get("surface_segments",{}).items() if k in ids),
                 "fixed" if "surface_segments" in result else None)
    return dict(point_RZ_m=point.tolist(),wall_segment=int(k),split_segment=int(result["original_segment"][k]),
                original_segment=original,surface=surface,chord_midpoint_matching_residual_m=float(distances[original]),
                nearest_corner=dict(index=corner,point_RZ_m=wall[corner].tolist(),distance_m=float(min(corners))))


def _geometric_failures(geom, result):
    from shapely.geometry import LineString, Polygon
    lines = {name:LineString(t) for name,t in result.get("targets",{}).items()}
    failures = []
    for name,t in result.get("targets",{}).items():
        line = lines[name]
        checks = dict(material_length=float(np.linalg.norm(np.diff(t,axis=0),axis=1).sum()) <= 1.2*geom.lengths[name],
                      simple=line.is_simple,fixed_clear=not LineString(t[1:-1]).intersects(geom.fixed_geometry),
                      lcfs_clear=not line.intersects(Polygon(geom.lcfs)),
                      vessel_clear=geom.vessel_solids is None or not line.intersects(geom.vessel_solids))
        failures.extend(dict(surface=name,predicate=key) for key,ok in checks.items() if not ok)
    targets = list(result.get("targets",{}).items())
    for j,(name,t) in enumerate(targets):
        for other,previous in targets[:j]:
            if LineString(t[1:-1]).intersects(LineString(previous[1:-1])):
                failures.append(dict(surface=name,other=other,predicate="design_interior_intersection"))
    return failures


def direct_report(geom, result, wall, prescribed_width=.003):
    """Read the per-width maps and reconstruct the exact MultiSurface predicates."""
    widths = [float(src.width*1000) for src in geom.sources]
    out = dict(widths_mm=widths,native_feasible=bool(result["feasible"]),geometry_valid=bool(result["geometry_valid"]),
               power_units="fraction of assigned source power",thresholds=dict(maximum_unmapped_power=.001,
               maximum_balance_error=.0001),by_width=[],prescribed_direct=None)
    if not result["geometry_valid"]:
        return dict(out,failed_predicates=["geometry_valid"],geometry_failures=_geometric_failures(geom,result))
    for src,m,width in zip(geom.sources,result["maps"],widths,strict=True):
        surface_data = {}
        for face in geom.surfaces:
            name,floor = face["name"],face["angle_floor_deg"]
            ids = result["surface_segments"][name]
            angle = m.angle_application(floor or 0.,ids)
            donors = np.flatnonzero(np.isin(m.segment,ids)&(np.asarray(m.hit_power)>0))
            hit = None
            if len(donors):
                invalid = donors[~np.isfinite(m.angles[donors])]
                i = int(invalid[0] if len(invalid) else donors[np.argmin(m.angles[donors])])
                hit = _location(m.points[i],m.segment[i],result,wall)
                if hit is not None:
                    hit.update(donor=i,source_x_m=float(src.x[i]),incidence_deg=_json(m.angles[i]))
            surface_data[name] = dict(power=float(np.asarray(m.segment_power)[ids].sum()), angle_floor_deg=floor,
                powered_hit_count=len(donors),minimum_angle_deg=angle["minimum_angle_deg"],minimum_hit=hit,
                angle_status=angle["status"] if floor is not None else "diagnostic_only_no_angle_floor",
                angle_passed=angle["angle_constraint_passed"] if floor is not None else None)
        checks = dict(unmapped_power=m.unmapped_power<=.001,balance=m.balance_error<=.0001,
                      illuminated=any(v["power"]>0 for v in surface_data.values()),
                      powered_angles=all(v["angle_passed"] is not False for v in surface_data.values()))
        missed = np.flatnonzero(np.asarray(m.segment)<0)
        row = dict(width_mm=width,prescribed_power=float(np.sum(src.weights)),mapped_power=float(m.mapped_power),
                   unmapped_power=float(m.unmapped_power),balance_error=float(m.balance_error),
                   unmapped_donor_count=len(missed),unmapped_source_x_range_m=(
                   [float(src.x[missed].min()),float(src.x[missed].max())] if len(missed) else None),
                   surfaces=surface_data,predicates=checks,direct_predicates_passed=all(checks.values()))
        out["by_width"].append(row)
        if math.isclose(src.width,prescribed_width,rel_tol=0.,abs_tol=1e-14):
            out["prescribed_direct"]=row
    native_angles = result["angles"]
    predicates = dict(maximum_unmapped_power=result["missing"]<=.001,
                      maximum_balance_error=result["power_balance_error"]<=.0001,
                      any_design_surface_illuminated=any(any(p>0 for p in ps) for ps in result["surface_power"].values()),
                      powered_surface_angles=not any(v["status"] in ("target_angle_violation","invalid_target_angle")
                                                     for v in native_angles.values()))
    out.update(native_maximum_unmapped_power=result["missing"],native_maximum_balance_error=result["power_balance_error"],
               native_union_angles=native_angles,predicates=predicates,
               failed_predicates=[k for k,ok in predicates.items() if not ok])
    return _json(out)


def continuous_report(geom, result, wall, prescribed_width, intersection):
    direct = direct_report(geom,result,wall,prescribed_width)
    failed = list(direct["failed_predicates"])
    native = result.get("boundary_witness")
    witness = dict(native) if native is not None else None
    if witness is not None:
        first_location = None
        if intersection is not None and intersection.get("point_RZ_m") is not None:
            first_location = _location(intersection["point_RZ_m"],intersection["wall_segment"],result,wall)
        witness["native_first_intersection"] = (dict(intersection,location=first_location) if intersection is not None else None)
        if not native["resolved"]:
            failed.append("separatrix_witness_resolved")
        else:
            if not native["angle_passed"]:
                failed.append("separatrix_witness_angle")
            if intersection is not None and intersection["wall_segment"]>=0:
                k=intersection["wall_segment"]
            else:
                ids=np.flatnonzero(np.asarray(result["original_segment"])==native["split_segment"])
                k=int(min(ids,key=lambda i:np.linalg.norm(np.asarray(result["wall"])[i:i+2].mean(axis=0)-native["point"]))) if len(ids) else None
            witness["location"]=_location(native["point"],k,result,wall)
    both = None
    if direct["prescribed_direct"] is not None and witness is not None:
        both=bool(direct["prescribed_direct"]["direct_predicates_passed"] and witness["resolved"] and witness["angle_passed"])
    return _json(dict(widths_mm=direct["widths_mm"],native_feasible=bool(result["feasible"]),failed_predicates=failed,
                      witness=witness,prescribed_direct=direct["prescribed_direct"],
                      prescribed_direct_with_witness_passed=both))


def transport_report(result, wall, power_W, direct_feasible, incidence=None):
    """Record the original thermal score only when deposition actually returned."""
    ratios=result.get("precise_length_ratios",{})
    dep=result.get("deposition")
    gates=dict(material_length=(max(ratios.values()) <= 1.2+1e-12) if ratios else None,
               continuous_direct=direct_feasible,geometry_valid=result.get("geometry_valid"))
    ledger=dep.get("ledger") if dep is not None else None
    if result.get("receiving_angles") is not None:
        gates["receiving_angles"]=all(v["passed"] for v in result["receiving_angles"].values())
    if ledger is not None:
        gates.update(deposition_unhit_power=ledger["unhit"]<=.001,deposition_balance=ledger["balance_error"]<=2e-12)
    outcome=None
    locations={}
    if dep is not None:
        cell_power=np.asarray(dep["cell_power"])
        denominator=float(np.sum(cell_power**2/np.asarray(dep["areas"])))
        imax=int(np.argmax(dep["q"]))
        outcome=dict(eligible=bool(result["feasible"]),peak_MW_m2=float(result["peak_q"])*power_W/1e6,
                     A2_m2=float(cell_power.sum()**2/denominator) if denominator>0 else None,
                     peak_location=_location(dep["points"][imax],dep["sample_parent"][imax],result,wall))
        if incidence is not None:
            powered=cell_power[dep["sample_cell"]]>0
            for name in result.get("receiving_angles",{}):
                ids=np.flatnonzero(powered & np.isin(dep["sample_parent"],result["surface_segments"][name]))
                if not len(ids):
                    locations[name]=None
                    continue
                invalid=ids[~np.isfinite(incidence[ids])]
                k=int(invalid[0] if len(invalid) else ids[np.argmin(incidence[ids])])
                point=_location(dep["points"][k],dep["sample_parent"][k],result,wall)
                locations[name]=dict(point,sample=k,incidence_deg=incidence[k])
    return _json(dict(native_feasible=bool(result["feasible"]),errors=result.get("errors",[]),gates=gates,
                      failed_predicates=[k for k,ok in gates.items() if ok is not None and not ok],
                      precise_length_ratios=ratios,receiving_angles=result.get("receiving_angles"),ledger=ledger,
                      thresholds=dict(deposition_unhit_power=.001,deposition_balance=2e-12,material_length=1.2+1e-12),
                      powered_receiving_angle_locations=locations,
                      measurement=outcome,diagnostics=dep.get("diagnostics") if dep is not None else None,
                      quadrature_integrated_power=dep.get("quadrature_integrated_power") if dep is not None else None))


def exception_evidence(exc):
    tb=exc.__traceback__
    frames=[]
    fields=("deposited","unhit","prescribed","balance","maxclosure","grad_sep","i","j","k","distance","delta")
    while tb is not None:
        frame=tb.tb_frame
        values={key:_json(frame.f_locals[key]) for key in fields if key in frame.f_locals and
                isinstance(frame.f_locals[key],(int,float,np.number)) and not isinstance(frame.f_locals[key],bool)}
        path=Path(frame.f_code.co_filename)
        frames.append(dict(file=path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path),
                           function=frame.f_code.co_name,line=tb.tb_lineno,numeric_locals=values))
        tb=tb.tb_next
    return dict(type=type(exc).__name__,message=str(exc),frames=frames)


def failure_basis(observations,native_exception):
    """Attribute the native stop using recorded layers, including caught inner exceptions."""
    final=observations["original_transport"]
    guards=observations["continuous"]
    initial=next((r for r in guards if len(r["widths_mm"])>1),None)
    maps=next((r for r in observations["direct_maps"] if len(r["widths_mm"])>1),None)
    site=native_exception["frames"][-1] if native_exception and native_exception["frames"] else None
    stage="native_call_completed" if native_exception is None else "native_exception"
    predicates=[]
    inner=observations["native_errors"]
    if site and Path(site["file"]).name=="continuous_boundary.py" and site["function"]=="__init__" and initial:
        stage,predicates="continuous_boundary_constructor",initial["failed_predicates"]
    elif final and not final[-1]["native_feasible"]:
        stage,predicates="original_boundary_transport",final[-1]["failed_predicates"]
        for item in inner:
            if item["call"] in ("ThermalPaths.evaluate","deposit_boundary"):
                stage=item["call"]
                break
    elif native_exception and observations["direct_maps"] and not observations["direct_maps"][-1]["geometry_valid"]:
        stage,predicates="original_geometry",observations["direct_maps"][-1]["failed_predicates"]
    return dict(stage=stage,native_raise_site=site,failed_predicates=predicates,
                guard_widths_mm=initial["widths_mm"] if initial else None,
                guard_direct_width_failures=({key:[r["width_mm"] for r in maps["by_width"] if not r["predicates"][key]]
                    for key in ("unmapped_power","balance","illuminated","powered_angles")} if maps else None),
                prescribed_direct_with_witness_passed=initial["prescribed_direct_with_witness_passed"] if initial else None,
                thermal_boundary_measurement_observed=any(r["measurement"] is not None for r in final))


class NativeObserver:
    """Wrappers return the very same native object or re-raise the same exception."""
    def __init__(self,wall,prescribed_width,power_W=1.e6):
        self.wall=np.asarray(wall,float)
        self.width=prescribed_width
        self.power_W=power_W
        self.data={k:[] for k in ("launch_grids","direct_maps","continuous","thermal_paths","broadening",
                                  "original_transport","native_errors","observer_errors")}
        self.intersection=None
        self.incidence_points=self.incidence_values=None
        self.at_transport_entry=0

    def report(self):
        return self.data

    @contextmanager
    def observing(self):
        from d3d_optimization import boundary_evaluator as boundary, continuous_boundary as continuous
        from d3d_optimization import frozen_field, independent_sources, multisurface, reference_transport

        def observe(label,callback,*args):
            try:
                callback(*args)
            except Exception as exc:
                self.data["observer_errors"].append(dict(call=label,error=exception_evidence(exc)))

        with ExitStack() as stack:
            def install(owner,name,callback,pre=None):
                native=getattr(owner,name)
                label=native.__qualname__
                @wraps(native)
                def wrapper(*args,**kwargs):
                    if pre is not None:
                        observe(label,pre,args,kwargs)
                    try:
                        value=native(*args,**kwargs)
                    except Exception as exc:
                        observe(label,lambda:self.data["native_errors"].append(dict(call=label,error=exception_evidence(exc))))
                        raise
                    if callback is not None:
                        observe(label,callback,args,value)
                    return value
                stack.enter_context(patch.object(owner,name,wrapper))

            def launch(args,kwargs):
                starts=np.ascontiguousarray(args[1],dtype="<f8")
                self.data["launch_grids"].append(dict(shape=list(starts.shape),dtype="little-endian float64",step_m=kwargs.get("step"),
                    max_length_m=kwargs.get("max_length"),first_RZ_m=starts[0].tolist(),last_RZ_m=starts[-1].tolist(),
                    sha256=hashlib.sha256(starts.tobytes()).hexdigest()))

            def first_hit(args,hits):
                k=int(hits.segment[0])
                self.intersection=dict(wall_segment=k,point_RZ_m=hits.points[0].tolist() if k>=0 else None,
                                       path_distance_m=float(hits.distance[0]) if k>=0 else None)

            def direct(args,value):
                geom,c=args
                require(np.all(np.asarray(c)==0),"diagnostic received a non-original shape")
                row=direct_report(geom,value,self.wall,self.width)
                row["zero_coefficient_dimension"]=len(c)
                self.data["direct_maps"].append(row)

            def guard(args,value):
                geom,c=args
                require(np.all(np.asarray(c)==0),"diagnostic received a non-original shape")
                self.data["continuous"].append(continuous_report(geom.base,value,self.wall,self.width,self.intersection))

            def thermal(args,value):
                obj=args[0]
                self.data["thermal_paths"].append(_json(reference_transport.thermal_summary(value,obj.source.weights,obj.reference)))

            def deposition(args,value):
                cell_power=np.asarray(value["cell_power"])
                denom=float(np.sum(cell_power**2/np.asarray(value["areas"])))
                observed={k:value[k] for k in ("ledger","diagnostics","quadrature_integrated_power")}
                observed.update(peak_MW_m2=float(np.max(value["q"]))*self.power_W/1e6,
                                A2_m2=float(cell_power.sum()**2/denom) if denom>0 else None,
                                receiving_angle_gate_at_this_point="pending")
                self.data["broadening"].append(_json(observed))

            def incidence(args,value):
                self.incidence_points,self.incidence_values=args[1],value

            def original(args,value):
                geom,c=args
                require(np.all(np.asarray(c)==0),"diagnostic received a non-original shape")
                prior=self.data["continuous"]
                direct_ok=prior[-1]["native_feasible"] if len(prior)>self.at_transport_entry else None
                dep=value.get("deposition")
                samples=self.incidence_values if dep is not None and self.incidence_points is dep["points"] else None
                self.data["original_transport"].append(transport_report(value,self.wall,geom.reference.power_W,direct_ok,samples))

            install(independent_sources,"trace_paths",None,launch)
            install(multisurface.MultiSurfaceEvaluator,"evaluate",direct)
            install(continuous,"first_hits",first_hit)
            install(continuous.BoundaryAwareEvaluator,"evaluate",guard,
                    lambda args,kwargs:setattr(self,"intersection",None))
            install(reference_transport.ThermalPaths,"evaluate",thermal)
            install(boundary,"deposit_boundary",deposition)
            install(frozen_field.FrozenField,"incidence",incidence)
            install(boundary.BoundaryEvaluator,"evaluate",original,
                    lambda args,kwargs:setattr(self,"at_transport_entry",len(self.data["continuous"])))
            yield


def select_failures(frozen,audit,events):
    """Lock distinct failed original-wall calls and keep every registered launch in the denominator."""
    jobs,rows=frozen["jobs"],audit["rows"]
    require(frozen["phase"]==audit["phase"]=="holdout" and len(jobs)==len(rows)==audit["planned_arm_count"],
            "final audited held-out denominator differs")
    require(sum(row["normal"] is True for row in rows)==audit["completed_arm_count"],"audited completion count differs")
    cases={case["case_id"]:case for case in frozen["cases"]}
    targeted,other={},[]
    for job,row in zip(jobs,rows,strict=True):
        require(all(row[k]==job[k] for k in ("case_id","seed","policy","relative_dir","population","projected_seeds","budget_s")),
                "failed arm identity differs from its registered job")
        case=cases[job["case_id"]]
        require(job["phase"]=="holdout" and job["input_root"]==case["root"],"registered launch input root differs")
        if row["normal"]:
            continue
        folder=row["relative_dir"]
        failed=[(i,e) for i,e in enumerate(events.get(folder,())) if e.get("status")=="failed"]
        require(len(failed)<=1,"more than one terminal failure in an arm")
        i,event=failed[0] if failed else (None,{})
        stage=event.get("stage")
        reference={k:row.get(k) for k in ("case_id","seed","policy","population","projected_seeds","budget_s",
                   "relative_dir","status","error","launch_ns","process_sha256","job_sha256","result_sha256","events_sha256")}
        reference.update(event_index=i,stage=stage,event_error=event.get("error"),start_s=event.get("start_s"),
                         end_s=event.get("end_s"),event_id=f"{row['events_sha256']}#{i}" if i is not None else None)
        if stage not in STAGES or row["status"]!="arm_failed":
            other.append(reference)
            continue
        require(row["error"]==event["error"] and row["result_sha256"] and row["events_sha256"],
                "original-wall failure lacks its matching finalized native error")
        level=STAGES[stage]
        require(level!="coarse" or job["policy"]=="coarse_adaptive","coarse failure attached to fine-only policy")
        target=targeted.setdefault((case["case_id"],stage),dict(case_id=case["case_id"],stage=stage,level=level,
            fidelity=frozen["levels"][level],source=frozen["source"],radius_mm=frozen["radius_mm"],reference=frozen["reference"],
            input_root=case["root"],input_sha256=case["hashes"],time_ms=case["time_ms"],failure_evidence=[]))
        target["failure_evidence"].append(reference)
    return dict(registered_arms=len(jobs),completed_arms=audit["completed_arm_count"],
                failed_arms=len(jobs)-audit["completed_arm_count"],targets=list(targeted.values()),other_failed_arms=other)


def _registered_payload(study):
    from verification.nf_anytime_campaign import versions
    from verification.nf_prospective_campaign import verify_process, verify_registration

    live=EvidenceStore(study)
    audit=live.read("audit.json")
    frozen=live.read("frozen.json")
    frozen_hash=live.sha256("frozen.json")
    require(audit["registration_sha256"]==frozen_hash and frozen["source"]=="OMP" and
            frozen["radius_mm"]==20 and frozen["levels"]==audited.LEVELS, "native holdout settings or audit link differ")
    require(Path(frozen["code_root"]).resolve()==ROOT,"diagnostic must use the frozen code root")
    verify_registration(live,frozen)
    inputs={case["case_id"]:case for case in frozen["cases"]}
    streams={}
    primary={"audit.json":live.sha256("audit.json"),"frozen.json":frozen_hash,"protocol.ko.md":live.sha256("protocol.ko.md")}
    for job,row in zip(frozen["jobs"],audit["rows"],strict=True):
        folder=live.path(job["relative_dir"])
        proc=verify_process(folder,dict(job,registration_sha256=frozen_hash))
        for field,name in (("job_sha256","job.json"),("events_sha256","events.jsonl"),
                           ("result_sha256","result.json"),("process_sha256","process.json")):
            path=folder/name
            digest=audited.digest(path) if field=="process_sha256" else proc[field]
            require(digest==row.get(field),"final audit hash differs for "+str(path))
            if digest is not None:
                primary[path.relative_to(live.root).as_posix()]=digest
        primary[(folder/"stdout.log").relative_to(live.root).as_posix()]=proc["stdout_sha256"]
        require(row["launch_ns"]==proc["launch_ns"],"failure launch differs from process receipt")
        if row["normal"]:
            continue
        if row["status"]=="arm_failed":
            rec=audited.read(folder/"result.json")
            raw,_=audited.read_events(folder/"events.jsonl")
            require(raw==rec["events"] and rec["error"]==row["error"] and
                    rec["job_sha256"]==proc["job_sha256"] and rec["launch_ns"]==proc["launch_ns"],
                    "failure raw stream differs from bound finalized result")
            streams[job["relative_dir"]]=raw
    selection=select_failures(frozen,audit,streams)
    for target in selection["targets"]:
        case=inputs[target["case_id"]]
        wall=audited.read(Path(case["root"])/REVIEW)["wall"]
        target["original_wall_coordinates_sha256"]=object_digest(wall)
    return dict(schema="nf-prospective-baseline-diagnostic-registration.v1",study_root=str(live.root),
                phase="posthoc_after_final_holdout_audit",primary_sha256=primary,protocol_sha256=frozen["protocol_sha256"],
                native_source_sha256=frozen["source_hashes"],native_reference_input_sha256=frozen["original_reference_hashes"],
                diagnostic_source_sha256=audited.digest(Path(__file__)),python_version=platform.python_version(),
                registered_python_version=frozen["python_version"],registered_packages=frozen["packages"],runtime_packages=versions(),
                original_wall_only=True,solver="native original-wall prep on CPU",geometry_search=False,
                launch_grid_evidence="same native prep, input bytes, operating point and numerical settings; primary streams did not store launch-grid arrays",
                diagnostic_scope="direct baseline predicates by width; a complete deposited thermal measurement requires reaching the original boundary transport",
                changed_primary_scores=False,substitute_or_add_arms=False,**selection)


def register_study(study):
    plan=_registered_payload(study)
    store=EvidenceStore(Path(study)/"baseline-diagnostic")
    code=Path(__file__).read_bytes()
    if not store.exists("source/nf_prospective_baseline_diagnostic.py"):
        store.write_bytes("source/nf_prospective_baseline_diagnostic.py",code)
    require(store.sha256("source/nf_prospective_baseline_diagnostic.py")==plan["diagnostic_source_sha256"],
            "diagnostic source snapshot differs")
    bind_evidence(store,"registration.json",plan)
    return plan


def run_target(target,prepare=None):
    from d3d_optimization.reference_transport import ReferencePoint
    from verification.nf_search_benchmark import measurement, prep
    from verification.nf_anytime_search import _measurement
    from verification.validate_boundary_transport import LEVELS

    index=2 if target["level"]=="fine" else 0
    require(target["stage"] in STAGES and target["level"]==STAGES[target["stage"]] and
            list(LEVELS[index])==target["fidelity"],"native diagnostic fidelity differs from failed call")
    root=Path(target["input_root"])
    wall=audited.read(root/REVIEW)["wall"]
    require(object_digest(wall)==target["original_wall_coordinates_sha256"],"native fixed-wall coordinates changed")
    reference=ReferencePoint(**target["reference"])
    observer=NativeObserver(wall,reference.lambda_q_m,reference.power_W)
    start_utc=datetime.now(timezone.utc).isoformat()
    start,cpu=perf_counter(),process_time()
    reported_error=native_exception=normalized=None
    call="verification.nf_search_benchmark.prep"
    with observer.observing():
        try:
            original,_=(prepare or prep)(root,target["source"],reference,index)
            if target["stage"]=="baseline":
                call="verification.nf_anytime_search._measurement(verification.nf_search_benchmark.measurement)"
                normalized=_measurement(measurement(original.baseline,reference))
        except Exception as exc:
            native_exception=exception_evidence(exc)
            reported_error=f"{type(exc).__name__}: {exc}"
    cpu_s,wall_s=process_time()-cpu,perf_counter()-start
    observed=observer.report()
    matching=[r["event_id"] for r in target["failure_evidence"] if r["event_error"]==reported_error]
    if target["stage"]=="baseline" and call=="verification.nf_search_benchmark.prep":
        matching=[]
    for filename,expected in target["input_sha256"].items():
        require(audited.digest(root/filename)==expected,"staged input changed during diagnosis: "+filename)
    return dict(schema="nf-prospective-baseline-diagnostic-result.v1",case_id=target["case_id"],stage=target["stage"],
                native_entry_point=call,level=target["level"],fidelity=target["fidelity"],source=target["source"],
                time_ms=target["time_ms"],input_root=str(root),reference=target["reference"],
                input_sha256=target["input_sha256"],original_wall_coordinates_sha256=target["original_wall_coordinates_sha256"],
                pid=os.getpid(),interpreter=sys.executable,threads={k:os.environ.get(k) for k in
                ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS")},start_utc=start_utc,
                finish_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=wall_s,process_cpu_seconds=cpu_s,
                native_call_completed=reported_error is None,native_error=reported_error,exception=native_exception,
                primary_failure_count=len(target["failure_evidence"]),matching_primary_failure_ids=matching,
                all_primary_native_errors_reproduced=len(matching)==len(target["failure_evidence"]),
                failure_basis=failure_basis(observed,native_exception),
                normalized_baseline_measurement=normalized,native_observations=observed,optimizer_started=False,
                full_thermal_deposition_returned=bool(observed["broadening"]),
                native_wall_thermal_measurement_returned=any(r["measurement"] is not None for r in observed["original_transport"]),
                input_hashes_rechecked_after=True)


def replay_study(study):
    from d3d_optimization.campaign_state import process_lock

    study=Path(study).resolve()
    path=study/"baseline-diagnostic"
    require((path/"registration.json").is_file(),"register from final sealed holdout audit before replay")
    store=EvidenceStore(path)
    with process_lock(path/"replay.lock"):
        plan=store.read("registration.json")
        require(_registered_payload(study)==plan and store.sha256("source/nf_prospective_baseline_diagnostic.py")==
                plan["diagnostic_source_sha256"],"frozen trial, diagnostic registration or source differs")
        registration_hash=store.sha256("registration.json")
        results=[]
        for target in plan["targets"]:
            name=f"results/{target['case_id']}/{target['stage']}.json"
            if not store.exists(name):
                result=run_target(target)
                result.update(diagnostic_registration_sha256=registration_hash,
                              native_source_registry_sha256=object_digest(plan["native_source_sha256"]),
                              primary_sha256=plan["primary_sha256"],failure_evidence=target["failure_evidence"])
                store.write(name,result)
            result=store.read(name)
            require(result["diagnostic_registration_sha256"]==registration_hash and
                    result["failure_evidence"]==target["failure_evidence"] and
                    result["original_wall_coordinates_sha256"]==target["original_wall_coordinates_sha256"],
                    "retained diagnostic outcome differs from its registration")
            results.append(dict(case_id=target["case_id"],stage=target["stage"],file=name,sha256=store.sha256(name),
                reproduced=result["all_primary_native_errors_reproduced"],native_error=result["native_error"],
                failure_basis=result["failure_basis"],
                primary_failures=result["primary_failure_count"],observer_errors=result["native_observations"]["observer_errors"],
                cpu_seconds=result["process_cpu_seconds"],wall_seconds=result["wall_seconds"]))
        require(_registered_payload(study)==plan,"primary evidence or native code changed during diagnosis")
        summary=dict(schema="nf-prospective-baseline-diagnostic-summary.v1",diagnostic_registration_sha256=registration_hash,
                     registered_arms=plan["registered_arms"],completed_arms=plan["completed_arms"],failed_arms=plan["failed_arms"],
                     other_failed_arms=plan["other_failed_arms"],original_wall_replay_count=len(results),results=results,
                     cpu_seconds=sum(r["cpu_seconds"] for r in results),wall_seconds=sum(r["wall_seconds"] for r in results),
                     primary_source_inputs_and_evidence_unchanged=True,primary_denominator_unchanged=True,
                     all_registered_diagnoses_reproduced=all(r["reproduced"] and not r["observer_errors"] for r in results))
        bind_evidence(store,"summary.json",summary)
        return summary


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,default=DEFAULT_STUDY)
    parser.add_argument("--action",choices=("register","replay"),required=True)
    args=parser.parse_args(argv)
    if args.action=="register":
        plan=register_study(args.study)
        print("REGISTERED",len(plan["targets"]),"original-wall replay targets from",plan["registered_arms"],"held-out arms")
        return 0
    result=replay_study(args.study)
    print("DIAGNOSED",result["original_wall_replay_count"],"case/stage pairs in",args.study/"baseline-diagnostic")
    return int(not result["all_registered_diagnoses_reproduced"])


if __name__=="__main__":
    raise SystemExit(main())
