"""Independent saved-output audit of the 1,200/1,600 ms test. No heat evaluator or GPU calls."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from d3d_optimization.case import load_case, separatrix
from d3d_optimization.continuous_boundary import continuous_design
from d3d_optimization.independent_sources import split_design_wall
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.tangent_geometry import Curve
from verification.validate_boundary_transport import LEVELS

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("D:/DivertorOptimization")
SUBDIR = "results/nf-finalist-robustness-20261106"
SHOT = Path("input/D3D/d3d/d3d_shot_1a3155370c.parquet")
REVIEW = Path("baseline_review/baseline_review.json")
GEOM = Path("input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json")
PANEL = Path("results/d3d-base-budget-20261105/paired/OMP/5161/panel.json")
ACTIVE = ["omp_entry", "shelf", "omp_connector", "floor"]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def separatrix_control(field, xp, psi, wall):
    # Independently locate two psi=psi_X roots just below the X-point. Follow the
    # descending field using adaptive RK45 with a signed-distance material event.
    polygon=Polygon(wall)
    z=xp[1]-.002
    radii=np.linspace(xp[0]-.045, xp[0]+.045,241)
    values=field.values(np.c_[radii,np.full(len(radii),z)])-psi
    pairs=np.where(values[1:]*values[:-1]<0)[0]
    assert len(pairs)==2
    hits={}
    for side,i in zip(("inner","outer"),pairs):
        radius=brentq(lambda r:field.values([[r,z]])[0]-psi,radii[i],radii[i+1],xtol=5e-15)
        start=np.array([radius,z])
        assert polygon.contains(Point(start)) and polygon.boundary.distance(Point(xp))>np.linalg.norm(start-xp)
        sign=-np.sign(field.gradient([start])[0,0])
        def rhs(s,p):
            g=field.gradient([p])[0]
            return sign*np.array([-g[1],g[0]])/np.linalg.norm(g)
        def material(s,p):
            q=Point(p)
            return polygon.boundary.distance(q)*(1 if polygon.covers(q) else -1)
        material.terminal=True
        material.direction=-1
        sol=solve_ivp(rhs,(0,1),start,rtol=1e-10,atol=1e-12,events=material,max_step=.001)
        assert sol.success and len(sol.y_events[0])==1
        hit=sol.y_events[0][0]
        d=np.array([LineString([a,b]).distance(Point(hit)) for a,b in zip(wall,wall[1:])])
        seg=int(np.argmin(d)); assert d[seg]<1e-11
        hits[side]=dict(RZ_m=hit.tolist(),segment=seg,label_error_Wb_rad=float(abs(field.values([hit])[0]-psi)))
    return hits


def audit_shapes(review,field,psi,lcfs,rows,vessel):
    checked,_=continuous_design(review,field,psi)
    split,surfaces,parent=split_design_wall(checked)
    indices={s["name"]:s["segment_index"] for s in surfaces}
    assert [s["name"] for s in surfaces[:4]]==ACTIVE
    assert [s["original_segment"] for s in surfaces[:4]]==[58,59,60,68]
    assert all(parent[i] not in (53,54,55) for i in indices.values())
    assert surfaces[1]["angle_floor_deg"]==1.5 and all(s["angle_floor_deg"] is None for s in surfaces if s["name"]!="shelf")
    ends={s["name"]:s["endpoints_RZ_m"] for s in surfaces[:4]}
    floor_i=indices["floor"]
    direction=split[floor_i+2]-split[floor_i+1]
    original=np.array(review["wall"])
    original_dir=original[70]-original[69]
    np.testing.assert_allclose(direction/np.linalg.norm(direction),original_dir/np.linalg.norm(original_dir),rtol=0,atol=2e-14)
    geometry={}; profiles=[]; all_angles=[]; min_radii=[]; length_max=[]; lip_gaps=[]
    lip=np.array([1.37256,-1.3247])
    min_gap=Curve(*ends["floor"],np.zeros(4)).minimum_distance(lip)[0]
    u=np.linspace(0,1,1001)
    for row in rows:
        coeff=np.array(row["candidate"]["vector"]).reshape(4,4)
        curves=[Curve(*ends[n],c,direction if n=="floor" else None) for n,c in zip(ACTIVE,coeff)]
        samples=[curve.point(u) for curve in curves]
        profiles.append(samples)
        floor=curves[-1]
        der=floor.tangent(1); angle=np.rad2deg(np.arctan2(abs(der[0]*direction[1]-der[1]*direction[0]),der@direction))
        ratio=max(curve.length_ratio() for curve in curves)
        gap=floor.minimum_distance(lip)[0]
        rad=1/max(c.max_curvature() for c in curves)
        # Inspect whole physical wall; candidate curves replace only the four selected intervals.
        replacements={indices[n]:s for n,s in zip(ACTIVE,samples)}
        outline=np.vstack([split[j:j+1] if j not in replacements else
                           replacements[j][:-1] for j in range(len(split)-1)]+[split[-1:]])
        line=LineString(outline)
        passes=angle<1e-8 and gap>=min_gap-1e-10 and rad>=.02-1e-10 and ratio<=1.2+1e-12
        passes &= line.is_simple and Polygon(outline).is_valid and not line.intersects(Polygon(lcfs)) and not line.intersects(vessel)
        geometry[row["i"]]=bool(passes)
        all_angles.append(angle); min_radii.append(rad); length_max.append(ratio); lip_gaps.append(gap)
    return ends, np.array(profiles), geometry, dict(max_G1_error_deg=max(all_angles),min_radius_m=min(min_radii),
                max_length_ratio=max(length_max),min_lip_gap_m=min(lip_gaps),required_min_gap_m=min_gap,
                geometry_ok_indices=[i for i,x in geometry.items() if x],
                note="1,001-vertex/surface sampled complete-wall check; production used its adaptive canonical screen")


def decision(folder,ids,geometry):
    r=read(folder/"result.json");sl=read(folder/"shortlist.json")
    assert r["status"]=="measured" and not r["reason"] and r["preflight_sha256"]==sha(folder/"preflight.json")
    assert r["coarse"]==sl["data"] and (folder/"shortlist.json").stat().st_mtime <= (folder/"result.json").stat().st_mtime
    summary={}
    for level in ("coarse","fine"):
        v=r[level]; rows=v["rows"]; base=v["baseline"]; lerr=[]
        assert len(rows)==12 and [x["i"] for x in rows]==list(range(12)) and [x["id"] for x in rows]==ids
        assert base["eligible"] and not base["errors"]
        for x in [base]+rows:
            if x["eligible"]:
                assert not x["errors"] and x["A2"]>0 and x["peak"]>0
                if "i" in x:
                    assert geometry[x["i"]]
            if x["ledger"]:
                g=x["ledger"]
                err=max(abs(g["prescribed"]-1),abs(g["deposited"]+g["unhit"]+g["unresolved"]-1),g["balance_error"])
                assert err<2e-12
                lerr.append(err)
        summary[level]=dict(eligible=[x["i"] for x in rows if x["eligible"]],
                            exclusions={str(x["i"]):x["errors"] for x in rows if not x["eligible"]},
                            max_power_closure_error=max(lerr),baseline_peak_MW_m2=base["peak"],baseline_A2_m2=base["A2"])
    coarse=r["coarse"]["rows"]; pick=[]
    for i in range(3):
        group=[(x["peak"],x["i"]) for x in coarse[4*i:4*i+4] if x["eligible"]]
        if group: pick.append(min(group)[1])
    assert pick==sl["ids"]==r["shortlisted"]
    base=r["fine"]["baseline"]; options=[dict(base,i=-1)]+[x for x in r["fine"]["rows"] if x["eligible"]]
    selected=min((x for x in options if x["i"] in [-1]+pick),key=lambda x:(x["peak"],x["i"]))
    oracle=min(options,key=lambda x:(x["peak"],x["i"]))
    regret=selected["peak"]/oracle["peak"]-1
    np.testing.assert_allclose([selected["peak"],oracle["peak"],base["peak"],regret],
                               [r["selected_peak"],r["best"],r["baseline_peak"],r["selected_excess_fraction"]],atol=1e-14,rtol=0)
    def score(x):
        return dict(i=x["i"],peak_MW_m2=x["peak"],A2_m2=x["A2"],peak_reduction_pct=100*(1-x["peak"]/base["peak"]))
    summary.update(shortlisted=pick,selected=score(selected),oracle=score(oracle),regret_pct=100*regret,
                   result_sha256=sha(folder/"result.json"),shortlist_sha256=sha(folder/"shortlist.json"))
    return summary


def audit(study,out):
    review=read(ROOT/REVIEW); geom=read(ROOT/GEOM); panel=read(ROOT/PANEL)
    rows=panel["rows"]; ids=[x["id"] for x in rows]
    assert len(set(ids))==12
    for i,row in enumerate(rows):
        c=row["candidate"];assert row["i"]==i and c["radius_mm"]==20 and c["source"]=="OMP" and len(c["vector"])==16
        assert row["id"]==hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()
    assert sha(ROOT/PANEL)==sha(DATA/"studies/d3d-base-budget-20261105/OMP/5161/panel.json")
    assert review["time_ms"]==1400 and review["shot_sha256"]==sha(ROOT/SHOT) and review["baseline_geometry_sha256"]==sha(ROOT/GEOM)
    wall=np.array(review["wall"]);np.testing.assert_array_equal(wall[:-1],geom["limiter"]);np.testing.assert_array_equal(wall[0],wall[-1])
    vv=unary_union([Polygon(v["pts"]) for v in geom["vv"]])
    times=pq.read_table(ROOT/SHOT,columns=["efit_times"]).to_pylist()[0]["efit_times"]
    report=dict(archived_shot_sha256=sha(ROOT/SHOT),public_geometry_sha256=sha(ROOT/GEOM),panel_sha256=sha(ROOT/PANEL),
                script_sha256=sha(Path(__file__)),protocol_sha256=sha(ROOT/SUBDIR/"time-protocol.ko.md"),
                source="single OMP",reference=asdict(ReferencePoint()),production_levels=[LEVELS[i] for i in (0,2)],
                active_original_segments=[58,59,60,68],fixed_OMP_exclusions=[53,54,55],time_slices={},
                limitations=["Two neighboring time slices from the same discharge; one 12-candidate seed panel.",
                             "Scores and ledger checks use emitted summaries; no independent heat recalculation or speed comparison.",
                             "Both shortlisted and full-panel choices use the stored cubic-1280 values. No separate blind evaluator rerun."])
    # Archived (1400 ms) control and an explicitly post-hoc fixed-winner transfer.
    archive_dir=(ROOT/PANEL).parent
    archived=[read(archive_dir/f"pair_{i:02d}.json") for i in range(12)]
    assert [r["id"] for r in archived]==ids
    assert read(archive_dir.parent.parent/"protocol.json")["reference"]==asdict(ReferencePoint())
    archive_base=read(archive_dir.parent/"setup.json")["baseline"]
    archive_pool=[(archive_base["peak"],-1)]+[(r["fine"]["peak"],r["i"]) for r in archived if r["fine"]["eligible"]]
    ref_peak,ref_best=min(archive_pool)
    assert ref_best==read(archive_dir/"result.json")["best_new"] and ref_best>=0
    report["archived_1400_control"]=dict(baseline_peak_MW_m2=archive_base["peak"],baseline_A2_m2=archive_base["A2"],
        best_i=ref_best,best_peak_MW_m2=ref_peak,best_A2_m2=archived[ref_best]["fine"]["A2"])
    for t in (1400,1200,1600):
        assert list(times).count(t)==1
        folder=study/str(t); staged=folder/"input-root"
        base=(ROOT if t==1400 else staged)
        raw=read(base/REVIEW); material=np.array(raw["wall"])
        np.testing.assert_array_equal(wall,material)
        assert raw["optimization_footprints"]==review["optimization_footprints"]
        assert raw["optimization_policy"]==review["optimization_policy"]
        assert read(base/GEOM.with_name("manifest.json"))["time_ms"]==raw["time_ms"]==t
        assert sha(base/SHOT)==sha(ROOT/SHOT) and sha(base/GEOM)==sha(ROOT/GEOM)
        field,lcfs,axis=load_case(base/SHOT,t); sep,xp,psi=separatrix(field,lcfs,axis)
        hits=separatrix_control(field,xp,psi,wall)
        ends,profiles,geook,shape=audit_shapes(raw,field,psi,lcfs,rows,vv)
        if t==1400:
            reference_points=profiles; ref_xp=xp; ref_hits=hits; ref_sep=sep; ref_psi=field.psi-psi
            continue
        info=read(folder/"preflight.json")
        assert info["time_ms"]==t and info["passed"] and info["reason"] is None
        staged_files={str(p.relative_to(staged).as_posix()) for p in staged.rglob("*") if p.is_file()}
        assert staged_files=={x.as_posix() for x in (SHOT,REVIEW,GEOM,GEOM.with_name("manifest.json"))}
        assert raw.keys()=={"wall","optimization_footprints","optimization_policy","sources","time_ms","legs"} and raw["sources"]=={"OMP":{},"IMP":{}}
        assert raw["legs"]==info["legs"] and ends==info["active_OMP_endpoints_RZ_m"]
        assert info["review_sha256"]==sha(ROOT/REVIEW) and info["driver_sha256"]==sha(ROOT/"verification/nf_equilibrium_times.py")
        assert info["shot_sha256"]==sha(ROOT/SHOT) and info["geometry_sha256"]==sha(ROOT/GEOM) and info["panel_sha256"]==sha(ROOT/PANEL)
        np.testing.assert_allclose(xp,info["X_RZ_m"],atol=1e-12,rtol=0)
        np.testing.assert_allclose(axis,info["axis_RZ_m"],atol=1e-12,rtol=0)
        np.testing.assert_allclose(sep,info["omp_sep_R_m"],atol=1e-12,rtol=0)
        for k,h in hits.items():
            assert h["segment"]==info["legs"][k]["segment"]==ref_hits[k]["segment"]
            h["vs_preregistered_trace_m"]=float(np.linalg.norm(np.subtract(h["RZ_m"],info["legs"][k]["strike_RZ_m"])))
            assert h["vs_preregistered_trace_m"]<1e-6
            h["displacement_from_1400_m"]=float(np.linalg.norm(np.subtract(h["RZ_m"],ref_hits[k]["RZ_m"])))
        shift=float(np.max(np.linalg.norm(profiles-reference_points,axis=-1)));assert shift<5e-14
        r=dict(exact_parquet_index=list(times).index(t),X_RZ_m=xp.tolist(),X_shift_from_1400_m=float(np.linalg.norm(xp-ref_xp)),
               LCFS_OMP_R_m=sep,LCFS_OMP_R_shift_from_1400_m=float(sep-ref_sep),separatrix_psi_Wb_rad=psi,
               grid_rms_delta_psi_minus_sepx_vs_1400_Wb_rad=float(np.sqrt(np.mean(((field.psi-psi)-ref_psi)**2))),
               independent_first_hits=hits,max_active_candidate_RZ_change_from_1400_m=shift,
               material_wall_identical_to_original=True,shape_constraints=shape,staged_inputs={s:sha(staged/s) for s in sorted(staged_files)})
        r["decision"]=decision(folder,ids,geook)
        frozen=read(folder/"result.json")["fine"]["rows"][ref_best]
        assert frozen["id"]==ids[ref_best]
        r["posthoc_1400_best_at_new_time"]=dict(i=ref_best,eligible=frozen["eligible"],
            peak_MW_m2=frozen.get("peak"),A2_m2=frozen.get("A2"))
        if frozen["eligible"]:
            r["posthoc_1400_best_at_new_time"].update(
                reduction_pct=100*(1-frozen["peak"]/r["decision"]["fine"]["baseline_peak_MW_m2"]),
                excess_vs_own_time_best_pct=100*(frozen["peak"]/r["decision"]["oracle"]["peak_MW_m2"]-1))
        report["time_slices"][str(t)]=r
        print(t,r["decision"]["selected"],r["decision"]["regret_pct"],flush=True)
    out.write_text(json.dumps(report,indent=2,allow_nan=False),encoding="utf-8")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",type=Path,default=DATA/"studies/nf-finalist-robustness-20261106/times")
    p.add_argument("--output",type=Path,default=ROOT/SUBDIR/"time-audit.json")
    a=p.parse_args();audit(a.study,a.output)
