"""Independent reconstruction of the scientific comparison from saved tables."""
import csv
import gzip
import importlib.util
import json

import numpy as np
import pytest

from d3d_optimization.user_output import file_hash
from verification.test_user_shape_replay import pair


def api():
    assert importlib.util.find_spec("verification.user_shape_results"), "post-run verifier is missing"
    from verification import user_shape_results
    return user_shape_results


def table(path, headers, rows):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers.split(","))
        writer.writerows(rows)


def fixture(folder, wall=None, reg=None):
    if wall is None:
        wall = np.array([[1., 0.], [2., 0.], [2., 1.], [1., 0.]])
        reg = ["fixed", "floor", "fixed"]
    n=len(reg)
    idx=reg.index("floor")
    assigned, frac = 2e6, np.zeros(n)
    frac[[0,idx,-1]]=[.2,.5,.2]
    mid=(wall[:-1]+wall[1:])/2
    ds = np.linalg.norm(np.diff(wall,axis=0),axis=1)
    s = np.r_[0., np.cumsum(ds)]
    areas = np.pi * (wall[:-1,0]+wall[1:,0])*ds
    mean = frac*assigned/(1e6*areas)
    table(folder/"geometry.csv", "vertex,R_m,Z_m,outgoing_region",
          [(i,*p,(reg+[""])[i]) for i,p in enumerate(wall)])
    table(folder/"profile.csv.gz", "wall_segment,region,s_start_m,s_end_m,R_start_m,Z_start_m,R_end_m,Z_end_m,area_m2,deposited_W,cell_mean_MW_m2",
          [(i,reg[i],s[i],s[i+1],*wall[i],*wall[i+1],areas[i],frac[i]*assigned,mean[i]) for i in range(n)])
    table(folder/"donors.csv.gz", "assigned_W,first_region,first_R_m,first_Z_m,source_x_m",
          [[1.8e6,"floor",*mid[idx],.001],[.2e6,"unhit",np.nan,np.nan,.002]])
    peak_samples=np.zeros(n)
    peak_samples[idx]=.7
    np.savez(folder/"boundary.npz", assigned_power_W=assigned, wall=wall, points=mid, q_MW_m2=peak_samples,
             receiving_wall=wall, cell_parent=np.arange(n), cell_power=frac, area_edges=np.r_[0.,areas.cumsum()],
             hit_segment=[idx,-1], hit_points=[mid[idx],[np.nan,np.nan]], source_x_m=[.001,.002])
    return dict(directory=folder.name, sha256={p.name:file_hash(p) for p in folder.iterdir()},
                peak_MW_m2=.7, peak_RZ_m=mid[idx].tolist(), peak_region="floor", A2_m2=frac.sum()**2/(frac**2/areas).sum(),
                ledger=dict(prescribed=1.,deposited=.9,unhit=.1,balance_error=0.),
                surfaces={k:dict(deposited_W=sum(frac[np.array(reg)==k])*assigned,
                                 first_hit_W=1.8e6 if k=="floor" else 0.) for k in set(reg)})


def test_reconstructs_peak_effective_area_closure_and_regions_independently(tmp_path):
    manifest=fixture(tmp_path)
    result=api().audit_measurement(tmp_path,manifest,2e6)
    assert result.metrics["peak_MW_m2"] == .7
    assert result.metrics["unhit_W"] == pytest.approx(.2e6)
    assert result.metrics["regions"]["floor"]["first_hit_W"] == 1.8e6
    assert result.metrics["regions"]["fixed"]["deposited_W"] == pytest.approx(.8e6)


@pytest.mark.parametrize("tamper,reason", [("hash","checksum"),("A2_m2","A2_m2"),("power","deposited_W")])
def test_bad_measurement_evidence_rejected(tmp_path,tamper,reason):
    manifest=fixture(tmp_path)
    if tamper=="hash":
        manifest["sha256"]["boundary.npz"]="wrong"
    elif tamper=="A2_m2":
        manifest["A2_m2"] *= 1.02
    else:
        manifest["surfaces"]["floor"]["deposited_W"] += 10.
    with pytest.raises(ValueError,match=reason):
        api().audit_measurement(tmp_path,manifest,2e6)


def test_material_alignment_uses_fraction_along_matching_reference_segment():
    # The second deformed parent is 2 m long, but its baseline counterpart is 1 m.
    original=np.array([[1.,0.],[2.,0.],[2.,1.]])
    moved=np.array([[1.,0.],[2.,0.],[2.,2.]])
    mapped=api().material_coordinate(original,moved,np.array([1,1]), np.array([[2.,.5],[2.,1.5]]))
    np.testing.assert_allclose(mapped, [1.25,1.75])
    with pytest.raises(ValueError,match="parent"):
        api().material_coordinate(original,moved,np.array([1]),np.array([[2.2,.4]]))


def test_ineligible_time_flag_cannot_hide_an_on_time_event(tmp_path):
    import json
    from d3d_optimization.campaign_geometry import Candidate
    c = Candidate("OMP",20,(0.,)*16)
    run = dict(arguments={"budget_seconds":8}, candidate={"id":c.identity}, search={"counts":{"eligible_fine":1}})
    record=dict(stage="fine",status="completed",within_budget=False,start_s=1.,end_s=5.,
                eligible=True,peak_q=3.,id=c.identity,candidate=c.record())
    (tmp_path/"events.jsonl").write_text(json.dumps(record)+"\n")
    fn = getattr(api(),"fine_receipt",None)
    assert fn is not None
    with pytest.raises(ValueError,match="deadline"):
        fn(tmp_path,run)


def test_coefficients_must_match_exported_material_wall():
    from types import SimpleNamespace
    from d3d_optimization.campaign_geometry import Candidate
    from d3d_optimization.tangent_geometry import Curve
    labels=["omp_entry","shelf","omp_connector","floor","fixed"]
    baseline=np.vstack([np.array([[1.+i+s,0.] for s in (0.,.25,.5,.75)]) for i in range(5)]+[[[6.,0.]]])
    regions=np.repeat(labels,4)
    ref=SimpleNamespace(wall=baseline,regions=regions)
    c=Candidate("OMP",20,(.001,)*16)
    geom=baseline.copy()
    for i in range(4):
        geom[4*i:4*i+5]=Curve(baseline[i*4],baseline[i*4+4],c.vector[i*4:i*4+4]).point(np.linspace(0,1,5))
    check=getattr(api(),"verify_geometry",None)
    assert check is not None
    check(ref,SimpleNamespace(wall=geom,regions=regions),c)
    geom[6,1] += .0002
    with pytest.raises(ValueError,match="wall.*coefficients"):
        check(ref,SimpleNamespace(wall=geom,regions=regions),c)

def test_attribution_uses_nonfloor_finalist_when_full_winner_is_floor_only():
    roles={key: {"id":value} for key,value in {
        "full_winner":"winner","full_winner_floor_ablation":"winner",
        "full_nonfloor_finalist":"distinct","full_nonfloor_finalist_floor_ablation":"ablation"}.items()}
    fn=getattr(api(),"attribution_roles",None)
    assert fn is not None
    assert fn(roles)==("full_nonfloor_finalist","full_nonfloor_finalist_floor_ablation")
    roles["full_nonfloor_finalist"]["id"]=None
    assert fn(roles) is None

def test_protocol_statistic_uses_full_denominator_and_all_levels():
    metrics=getattr(api(),"pair_metrics",None)
    decision=getattr(api(),"protocol_decision",None)
    assert metrics is not None and decision is not None
    def row(q):
        return dict(candidate_eligible=True, peak_MW_m2=q)
    fine=metrics(150.,row(101.005),row(100.))
    assert fine["signed_floor_minus_full_pct"]==pytest.approx(1.005)
    assert fine["floor_reduction_vs_reference_pct"]==pytest.approx(100*(1-101.005/150))
    assert fine["full_reduction_vs_reference_pct"]==pytest.approx(100/3)
    levels={k:fine for k in ("cubic640","cubic1280","quintic1280")}
    assert decision(levels)["status"]=="full_lower_over_1pct"
    levels["cubic640"]=metrics(150.,row(99.9),row(100.))
    assert decision(levels)["status"]=="mixed"
    levels["cubic640"]=metrics(99.,row(99.9),row(100.))
    assert decision(levels)["status"]=="not_adjudicated"
    assert decision({"cubic1280":fine})["status"]=="not_adjudicated"

def test_refined_peak_material_region_must_match_its_position(tmp_path):
    manifest=fixture(tmp_path)
    manifest["peak_region"]="fixed"
    with pytest.raises(ValueError,match="peak region"):
        api().audit_measurement(tmp_path,manifest,2e6)


def test_pair_analysis_consumes_only_saved_measurements(pair,tmp_path):
    from d3d_optimization.campaign_geometry import Candidate
    from d3d_optimization.tangent_geometry import Curve
    from verification.user_shape_replay import prepare_protocol, _reproduction
    api_module=api()
    names=["fixed", *api_module.ACTIVE, "fixed"]
    fixed=np.vstack([[[1+i+u,0.] for u in (0.,.25,.5,.75)] for i in range(6)]
                    + [[[7.,0.],[7.,1.],[1.,1.],[1.,0.]]])
    regions=np.repeat(names,4).tolist()+["fixed"]*3
    def evidence(folder,c=None):
        folder.mkdir()
        wall=fixed.copy()
        if c:
            for j in range(4):
                a,b=(j+1)*4,(j+2)*4
                wall[a:b+1]=Curve(fixed[a],fixed[b],c.vector[j*4:j*4+4]).point(np.linspace(0.,1.,5))
        return dict(fixture(folder,wall,regions),candidate_eligible=c is not None,physics_measurable=True,
                    reference_only=c is None)
    for folder in pair:
        report=json.loads((folder/"result.json").read_text())
        c=Candidate.from_record(report["candidate"])
        report["reference_SI"]["power_W"]=2e6
        report["baseline"]=evidence(folder/"baseline")
        report["selected"]=evidence(folder/"selected",c)
        report["search"]["counts"]={"eligible_fine":3}
        (folder/"result.json").write_text(json.dumps(report))
        events=[json.loads(line) for line in (folder/"events.jsonl").read_text().splitlines()]
        for e in events:
            if e.get("peak_q"):
                e["peak_q"]*=.35/4.
        (folder/"events.jsonl").write_text("".join(json.dumps(e)+"\n" for e in events))
    replay=tmp_path/"replay"
    replay.mkdir()
    protocol=prepare_protocol(*pair,["cubic1280"],"cpu")
    (replay/"protocol.json").write_text(json.dumps(protocol))
    grid={"resolution":list(protocol["levels"]["cubic1280"]),"candidates":{}}
    target=replay/"cubic1280"
    target.mkdir()
    grid["baseline"]=evidence(target/"baseline")
    for identity,record in protocol["candidates"].items():
        grid["candidates"][identity]=dict(evidence(target/identity,Candidate.from_record(record)),id=identity)
    grid["reproduction"]=_reproduction(protocol,grid)
    results=dict(complete=True,protocol_sha256=file_hash(replay/"protocol.json"),results={"cubic1280":grid})
    (replay/"report.json").write_text(json.dumps(results))
    out=tmp_path/"analysis"
    summary=api_module.run(*pair,replay,out)
    assert summary["preregistered_decision"]["status"]=="not_adjudicated"
    assert summary["search_pair"]["signed_floor_minus_full_pct"]==0.
    assert summary["attribution_plot_roles"]["with_nonfloor"]=="full_nonfloor_finalist"
    assert all(file_hash(out/name)==checksum for name,checksum in summary["analysis_sha256"].items())
    assert (out/"profiles.png").stat().st_size>4096
    assert (out/"regional-power.csv").is_file()
    grid["reproduction"]["passed"]=False
    (replay/"report.json").write_text(json.dumps(results))
    with pytest.raises(ValueError,match="canonical reproduction"):
        api_module.run(*pair,replay,tmp_path/"tampered-analysis")
