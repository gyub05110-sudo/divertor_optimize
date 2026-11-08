"""Verify matched search and replay exports without running a physics evaluation."""
import argparse
import csv
from dataclasses import dataclass
import gzip
import json
from pathlib import Path

import numpy as np

from d3d_optimization.user_output import file_hash
from verification.user_shape_replay import prepare_protocol, _reproduction


FILES = {"geometry.csv", "profile.csv.gz", "donors.csv.gz", "boundary.npz"}
ACTIVE = ("omp_entry", "shelf", "omp_connector", "floor")
LABELS = {"baseline":"Supplied wall", "floor":"Floor only", "full":"Four surfaces"}
STYLES = {"baseline":("#4f5560", ":"), "floor":("#10689a", "--"), "full":("#a24713", "-")}


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def match(actual, expected, label, *, atol=1e-12, rtol=1e-11, equal_nan=False):
    a, b = np.asarray(actual), np.asarray(expected)
    require(a.shape == b.shape and np.allclose(a, b, atol=atol, rtol=rtol, equal_nan=equal_nan),
            f"reconstruction differs for {label}")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_columns(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        header, *rows = csv.reader(stream)
    return dict(zip(header, zip(*rows, strict=True), strict=True))


def coord(cols, r, z):
    return np.array([cols[r], cols[z]], dtype=float).T


def arclength(wall):
    return np.r_[0., np.linalg.norm(np.diff(wall, axis=0), axis=1).cumsum()]


def material_coordinate(reference_wall, wall, parent, points):
    """Same material segment and local longitudinal fraction on the baseline."""
    require(reference_wall.shape == wall.shape and np.all((parent >= 0) & (parent < len(wall)-1)),
            "incompatible material parent indexing")
    v = np.diff(wall, axis=0)[parent]
    delta = points - wall[parent]
    frac = np.einsum("ij,ij->i", delta, v) / np.einsum("ij,ij->i", v, v)
    require(np.isfinite(frac).all() and np.all((frac >= -1e-8) & (frac <= 1+1e-8))
            and np.all(np.linalg.norm(delta-frac[:,None]*v, axis=1) <= 1e-9),
            "coordinate is outside its material parent")
    arc = arclength(reference_wall)
    return arc[parent] + np.clip(frac, 0., 1.) * np.diff(arc)[parent]


@dataclass
class Measurement:
    metrics: dict
    wall: np.ndarray
    regions: np.ndarray
    parent: np.ndarray
    midpoints: np.ndarray
    mean_q: np.ndarray


def audit_measurement(folder, record, power):
    folder = Path(folder)
    hashes = record["sha256"]
    require(set(hashes) == FILES and record["directory"] == folder.name, "incomplete measurement file inventory")
    for name, expected in hashes.items():
        require(file_hash(folder/name) == expected, f"checksum mismatch: {folder/name}")

    geom = read_columns(folder/"geometry.csv")
    wall = coord(geom, "R_m", "Z_m")
    region = np.array(geom["outgoing_region"][:-1])
    require(np.array_equal(np.asarray(geom["vertex"], int), np.arange(len(wall)))
            and np.array_equal(wall[0],wall[-1]) and geom["outgoing_region"][-1] == "", "invalid material vertices")
    cols = read_columns(folder/"profile.csv.gz")
    start = coord(cols, "R_start_m", "Z_start_m")
    end = coord(cols, "R_end_m", "Z_end_m")
    parent = np.asarray(cols["wall_segment"], int)
    midpoint = (start+end)/2
    material_coordinate(wall, wall, parent, midpoint)
    require(np.array_equal(region[parent], cols["region"]), "receiving region differs from material parent")
    require(np.array_equal(start[1:], end[:-1]) and np.array_equal(start[0], end[-1]), "receiving wall not closed")
    lengths = np.linalg.norm(end-start, axis=1)
    area = np.pi*(start[:,0]+end[:,0])*lengths
    watts = np.asarray(cols["deposited_W"], float)
    require(np.all(area>0) and np.all(watts>=0), "invalid receiving area or power")
    mean_q = watts/area/1e6
    match(area, np.asarray(cols["area_m2"],float), "toroidal area")
    match(mean_q, np.asarray(cols["cell_mean_MW_m2"],float), "cell mean heat flux")
    arc = np.r_[0.,lengths.cumsum()]
    match(arc[:-1],np.asarray(cols["s_start_m"],float),"profile arc start")
    match(arc[1:],np.asarray(cols["s_end_m"],float),"profile arc end")

    donors = read_columns(folder/"donors.csv.gz")
    first = np.asarray(donors["first_region"])
    assigned = np.asarray(donors["assigned_W"], float)
    require(np.isfinite(assigned).all() and np.all(assigned>=0), "invalid assigned donor power")
    with np.load(folder/"boundary.npz", allow_pickle=False) as raw:
        match(raw["assigned_power_W"], power, "assigned_power_W")
        match(raw["wall"], wall, "material wall")
        match(raw["receiving_wall"], np.vstack([start, end[-1]]), "receiving wall")
        require(np.array_equal(raw["cell_parent"], parent), "NPZ material parent differs")
        match(raw["cell_power"]*power, watts, "cell watts")
        match(np.diff(raw["area_edges"]), area, "area edges", atol=3e-14)
        hit = raw["hit_segment"]
        require(np.all((hit >= -1) & (hit < len(region))), "invalid first-hit segment")
        hit_names = np.where(hit >= 0, region[np.maximum(hit,0)], "unhit")
        require(np.array_equal(hit_names,first), "first-hit region mismatch")
        match(raw["source_x_m"], np.asarray(donors["source_x_m"],float), "source offset")
        match(raw["hit_points"],coord(donors,"first_R_m","first_Z_m"),"first-hit coordinate",equal_nan=True)
        q = raw["q_MW_m2"]
        require(np.isfinite(q).all() and np.all(q>=0), "invalid refined whole-wall heat flux")
        peak_idx = int(np.argmax(q))
        peak = float(q[peak_idx])
        peak_point=raw["points"][peak_idx]
        match(peak_point,record["peak_RZ_m"],"refined peak position")
        delta=np.diff(wall,axis=0)
        alpha=np.clip(np.sum((peak_point-wall[:-1])*delta,axis=1)/np.sum(delta**2,axis=1),0.,1.)
        distance=np.linalg.norm(peak_point-(wall[:-1]+alpha[:,None]*delta),axis=1)
        require(record["peak_region"] in region[distance<1e-9],"refined peak region and material position disagree")

    area2 = float(watts.sum()**2 / np.sum(watts**2/area))
    match(peak,record["peak_MW_m2"],"peak_MW_m2")
    match(area2,record["A2_m2"],"A2_m2")
    unhit = float(assigned[first=="unhit"].sum())
    prescribed, deposited = float(assigned.sum()), float(watts.sum())
    residual = abs(deposited+unhit-prescribed)
    match(prescribed, power, "integrated assigned power")
    require(residual <= 2e-12*power, "first-hit and deposition power do not close")
    for key,value in dict(prescribed=prescribed, deposited=deposited, unhit=unhit, balance_error=residual).items():
        match(value/power,record["ledger"][key],"ledger "+key)
    require(set(region) == set(record["surfaces"]), "regional inventory differs")
    regions = {}
    for name in sorted(set(region.tolist())):
        received = float(watts[region[parent] == name].sum())
        hit_power = float(assigned[first == name].sum())
        match(received,record["surfaces"][name]["deposited_W"],name+" deposited_W")
        match(hit_power,record["surfaces"][name]["first_hit_W"],name+" first_hit_W")
        regional = dict(first_hit_W=hit_power, deposited_W=received, net_spreading_W=received-hit_power)
        on_region=(first == name) & (assigned > 0)
        for column,key in (("physical_width_m","hit_weighted_physical_width_m"),
                           ("parallel_length_m","hit_weighted_parallel_length_m")):
            value=None
            if column in donors and on_region.any():
                values=np.asarray(donors[column],float)[on_region]
                require(np.isfinite(values).all(), "nonfinite powered donor "+column)
                value=float(np.average(values,weights=assigned[on_region]))
            if key in record["surfaces"][name]:
                stored=record["surfaces"][name][key]
                if stored is None or value is None:
                    require(stored is None and value is None, name+" missing "+key)
                else:
                    match(value,stored,name+" "+key)
            regional[key]=value
        regions[name] = regional
    metrics = dict(peak_MW_m2=peak, A2_m2=area2, assigned_W=prescribed, deposited_W=deposited,
                   unhit_W=unhit, closure_residual_W=residual, peak_RZ_m=record["peak_RZ_m"],
                   peak_region=record["peak_region"], regions=regions)
    return Measurement(metrics,wall,region,parent,midpoint,mean_q)


def verify_geometry(reference, actual, candidate):
    """Rebuild only spline coordinates and bind the returned arrays to the hashed vector."""
    from d3d_optimization.tangent_geometry import Curve
    require(np.array_equal(reference.regions,actual.regions), "material region sequence changed")
    expected = reference.wall.copy()
    for i,name in enumerate(ACTIVE):
        ids = np.flatnonzero(reference.regions==name)
        require(len(ids)>0 and np.all(np.diff(ids)==1), "unresolved material portion "+name)
        a,b = int(ids[0]),int(ids[-1]+1)
        direction = expected[b+1]-expected[b] if name=="floor" else None
        curve = Curve(reference.wall[a],reference.wall[b],candidate.vector[4*i:4*i+4],direction)
        expected[a:b+1] = curve.point(np.linspace(0.,1.,b-a+1))
    match(actual.wall,expected,"exported wall from coefficients",atol=1e-11,rtol=0)


def fine_receipt(folder, report):
    events=[json.loads(line) for line in (folder/"events.jsonl").read_text(encoding="utf8").splitlines() if line]
    fine=[e for e in events if e["stage"]=="fine"]
    budget=report["arguments"]["budget_seconds"]
    for row in fine:
        require(row["within_budget"] == (row["end_s"]<=budget), "fine event deadline flag conflicts with its clock")
    admitted=[e for e in fine if e["status"]=="completed" and e["within_budget"] and e.get("eligible")]
    require(admitted and len(admitted)==report["search"]["counts"]["eligible_fine"], "admitted fine count differs")
    best=min(admitted,key=lambda e:e["peak_q"])
    require(report["candidate"]["id"]==best["id"], "returned candidate differs from on-time minimum")
    return dict(id=best["id"],unrounded_peak_q=best["peak_q"],completed_s=best["end_s"],
                eligible_fine=len(admitted),fine_attempts=len(fine),budget_s=budget)


def _csv(path, rows):
    with path.open("x",newline="",encoding="utf8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def attribution_roles(roles):
    for whole, ablated in (("full_winner","full_winner_floor_ablation"),
                           ("full_nonfloor_finalist","full_nonfloor_finalist_floor_ablation")):
        a=roles.get(whole,{})
        b=roles.get(ablated,{})
        if a.get("id") and b.get("id") and a["id"] != b["id"]:
            return whole, ablated
    return None


def pair_metrics(reference_peak, floor, full):
    q_floor, q_full = floor.get("peak_MW_m2"), full.get("peak_MW_m2")
    ratio = q_floor/q_full if q_floor is not None and q_full is not None else None
    return dict(reference_peak_MW_m2=reference_peak, floor_peak_MW_m2=q_floor, full_peak_MW_m2=q_full,
                floor_eligible=floor["candidate_eligible"], full_eligible=full["candidate_eligible"],
                floor_reduction_vs_reference_pct=100*(1-q_floor/reference_peak) if q_floor is not None else None,
                full_reduction_vs_reference_pct=100*(1-q_full/reference_peak) if q_full is not None else None,
                floor_over_full_ratio=ratio, signed_floor_minus_full_pct=100*(ratio-1) if ratio is not None else None)


def protocol_decision(levels):
    names=("cubic1280","quintic1280","cubic640")
    ok=all(k in levels and levels[k]["floor_eligible"] and levels[k]["full_eligible"]
           and levels[k]["floor_reduction_vs_reference_pct"]>0
           and levels[k]["full_reduction_vs_reference_pct"]>0 for k in names)
    status="not_adjudicated"
    if ok:
        ratios=[levels[k]["floor_over_full_ratio"] for k in names]
        if all(.99<=r<=1.01 for r in ratios):
            status="within_1pct"
        elif all(r<.99 for r in ratios):
            status="floor_lower_over_1pct"
        elif all(r>1.01 for r in ratios):
            status="full_lower_over_1pct"
        else:
            status="mixed"
    return dict(status=status, required_levels=list(names), operational_margin_pct=1.,
                requires_both_eligible_and_improved_on_each_level=True,
                statistic="100*(q_floor/q_full-1); positive means lower full-family peak")


def plots(out, series, prefix="", labels=LABELS):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ref=series["baseline"]
    base_arc=arclength(ref.wall)
    changed=np.isin(ref.regions,ACTIVE)
    use=np.r_[changed,False] | np.r_[False,changed]
    region_names={"omp_entry":"OMP entry", "shelf":"Shelf", "omp_connector":"Connector", "floor":"Floor"}
    fig,axes=plt.subplots(2,1,figsize=(9,6.4),height_ratios=(2.1,1),layout="constrained")
    for label,profile in series.items():
        color,line=STYLES[label]
        axes[0].plot(*profile.wall.T,color=color,ls=line,lw=1.7,label=labels[label],zorder={"full":2,"floor":3,"baseline":4}[label])
        if label!="baseline":
            for name in ACTIVE:
                indices=np.flatnonzero(ref.regions==name)
                nodes=np.arange(indices[0],indices[-1]+2)
                chord=ref.wall[nodes[-1]]-ref.wall[nodes[0]]
                normal=np.array([-chord[1],chord[0]])/np.linalg.norm(chord)
                disp=(profile.wall[nodes]-ref.wall[nodes])@normal*1000
                axes[1].plot(base_arc[nodes],disp,color=color,ls=line,lw=1.5,zorder={"full":2,"floor":3}[label])
    pts=ref.wall[use]
    margin=.042
    axes[0].set(xlim=(pts[:,0].min()-margin,pts[:,0].max()+margin),
                ylim=(min(p.wall[use,1].min() for p in series.values())-margin,
                      max(p.wall[use,1].max() for p in series.values())+margin),
                xlabel="R (m)",ylabel="Z (m)",aspect="equal",title="Material wall in the lower divertor")
    axes[0].legend(loc="best",ncol=3,fontsize=9)
    axes[1].axhline(0,color="#888888",lw=.6)
    axes[1].set(xlabel="Arclength on supplied wall (m)",ylabel="Normal displacement (mm)",
                title="Displacement on adjustable portions")
    fig.savefig(out/(prefix+"geometry.png"),dpi=180)
    plt.close(fig)

    fig=plt.figure(figsize=(11,6.3),layout="constrained")
    grid=fig.add_gridspec(2,4,height_ratios=[1.2,1])
    whole=fig.add_subplot(grid[0,:])
    local=[fig.add_subplot(grid[1,i]) for i in range(4)]
    for label,profile in series.items():
        s=material_coordinate(ref.wall,profile.wall,profile.parent,profile.midpoints)
        color,line=STYLES[label]
        whole.plot(s,profile.mean_q,color=color,ls=line,lw=1.35,label=labels[label],zorder={"full":2,"floor":3,"baseline":4}[label])
        for ax,name in zip(local,ACTIVE):
            keep=profile.regions[profile.parent]==name
            orig=np.flatnonzero(ref.regions==name)
            x=(s[keep]-base_arc[orig[0]])*1000
            ax.plot(x,profile.mean_q[keep],color=color,ls=line,lw=1.35,zorder={"full":2,"floor":3,"baseline":4}[label])
    whole.set(xlabel="Arclength on supplied wall (m)",ylabel="Cell mean (MW/m²)",
              title="Whole-wall surface-normal heat flux",xlim=(0,base_arc[-1]),ylim=(0,None))
    whole.legend(loc="best",ncol=3)
    for i,(name,ax) in enumerate(zip(ACTIVE,local)):
        ax.set(title=region_names[name],xlabel="Material distance (mm)",ylim=(0,None))
        ax.ticklabel_format(axis="y",style="sci",scilimits=(-2,2))
        if i==0:
            ax.set_ylabel("Cell mean (MW/m²)")
    fig.savefig(out/(prefix+"profiles.png"),dpi=180)
    plt.close(fig)


def run(full, floor, replay, output):
    from d3d_optimization.campaign_geometry import Candidate
    full,floor,replay,output=map(Path,(full,floor,replay,output))
    require(not output.exists(),"analysis output must be a new directory")
    source=read_json(replay/"protocol.json")
    replay_doc=read_json(replay/"report.json")
    require(replay_doc["complete"] and replay_doc["protocol_sha256"]==file_hash(replay/"protocol.json"),
            "incomplete replay or protocol checksum mismatch")
    # Reuse the manifest, candidate-ID, event and matched-setting checks that
    # sealed this pair before replay. This call performs no physics evaluation.
    current=prepare_protocol(full,floor,list(source["levels"]),source["backend"])
    for key in ("inputs","matched_search","reference_SI","search_code_sha256","roles","candidates","levels","search_levels"):
        require(json.loads(json.dumps(current[key]))==source[key], "replay protocol no longer matches "+key)
    for arm in ("full","floor"):
        for key in ("result_sha256","events_sha256"):
            require(current["sources"][arm][key]==source["sources"][arm][key], "search source checksum mismatch")
    require(set(source["levels"])==set(replay_doc["results"]),"replay level inventory differs")
    power=source["reference_SI"]["power_W"]
    series={}
    receipts={}
    regional=[]
    replay_rows=[]
    replay_gaps={}
    contrast=attribution_roles(source["roles"])
    control_series={}
    def append_regions(level,role,identity,metrics):
        for region,values in metrics["regions"].items():
            regional.append(dict(level=level,role=role,id=identity,region=region,**values))
    for arm in ("floor","full"):
        result=Path(current["sources"][arm]["result"])
        report=read_json(result)
        receipts[arm]=fine_receipt(result.parent,report)
        base=audit_measurement(result.parent/"baseline",report["baseline"],power)
        if "baseline" not in series:
            series["baseline"]=base
            append_regions("search_fine","baseline","",base.metrics)
        else:
            match(base.wall,series["baseline"].wall,"common reference geometry",rtol=0,atol=0)
            match(base.mean_q,series["baseline"].mean_q,"common reference heat profile")
        selected=audit_measurement(result.parent/"selected",report["selected"],power)
        verify_geometry(base,selected,Candidate.from_record(report["candidate"]))
        series[arm]=selected
        append_regions("search_fine",arm,report["candidate"]["id"],selected.metrics)

    for level,grid in replay_doc["results"].items():
        require(grid["resolution"]==source["levels"][level],"replay resolution differs")
        require(set(grid["candidates"])==set(source["candidates"]), "replay candidate inventory differs")
        require(grid["baseline"].get("physics_measurable"),"replayed reference lacks measurable heat")
        baseline=audit_measurement(replay/level/"baseline",grid["baseline"],power)
        append_regions(level,"baseline","",baseline.metrics)
        if contrast and level=="cubic1280":
            control_series["baseline"]=baseline
        replay_rows.append(dict(level=level,role="baseline",id="",eligible="",peak_MW_m2=baseline.metrics["peak_MW_m2"],
                                A2_m2=baseline.metrics["A2_m2"],errors=""))
        for identity,record in grid["candidates"].items():
            require(identity==record["id"], "replay candidate identity mismatch")
            measured=None
            if "sha256" in record:
                measured=audit_measurement(replay/level/identity,record,power)
                verify_geometry(baseline,measured,Candidate.from_record(source["candidates"][identity]))
            else:
                require(not record["candidate_eligible"] and "peak_MW_m2" not in record,"replay heat artifacts missing")
            for role,item in source["roles"].items():
                if item["id"]!=identity:
                    continue
                replay_rows.append(dict(level=level,role=role,id=identity,eligible=record["candidate_eligible"],
                                   peak_MW_m2=record.get("peak_MW_m2"),A2_m2=record.get("A2_m2"),
                                   errors="|".join(record.get("errors",[]))))
                if measured:
                    append_regions(level,role,identity,measured.metrics)
                    if contrast and level=="cubic1280" and role in contrast:
                        control_series["full" if role==contrast[0] else "floor"]=measured
        if level=="cubic1280":
            require(_reproduction(source,grid)==grid.get("reproduction"),"canonical reproduction report conflicts with scores")
        floor_score=grid["candidates"][source["roles"]["floor_winner"]["id"]]
        full_score=grid["candidates"][source["roles"]["full_winner"]["id"]]
        base_peak=baseline.metrics["peak_MW_m2"]
        replay_gaps[level]=pair_metrics(base_peak,floor_score,full_score)
        for row in replay_rows:
            if row["level"]==level:
                peak=row["peak_MW_m2"]
                row["reduction_vs_reference_pct"]=100*(1-peak/base_peak) if peak is not None else None
                row["signed_vs_full_winner_pct"]=(100*(peak/full_score["peak_MW_m2"]-1)
                    if peak is not None and full_score.get("peak_MW_m2") is not None else None)
    search_pair=pair_metrics(series["baseline"].metrics["peak_MW_m2"],
                             dict(series["floor"].metrics,candidate_eligible=True),
                             dict(series["full"].metrics,candidate_eligible=True))
    diff=np.linalg.norm(series["floor"].wall-series["full"].wall,axis=1)
    ref=series["baseline"]
    region_diff={}
    for name in ACTIVE:
        ids=np.flatnonzero(ref.regions==name)
        region_diff[name]=float(1000*np.max(diff[ids[0]:ids[-1]+2]))
    decision=protocol_decision(replay_gaps)
    reproduction=replay_doc["results"].get("cubic1280",{}).get("reproduction")
    if reproduction is not None and not reproduction["passed"]:
        decision.update(status="not_adjudicated",reason="canonical fine reproduction failed")
    summary=dict(schema="user-shape-comparison.v1",source_search=current["sources"],source_replay=str(replay.resolve()),
                 protocol_sha256=file_hash(replay/"protocol.json"),replay_report_sha256=file_hash(replay/"report.json"),
                 verifier_sha256=file_hash(__file__),inputs=source["inputs"],matched_search=source["matched_search"],
                 reference_SI=source["reference_SI"],fine_selection=receipts,
                 search_metrics={k:v.metrics for k,v in series.items()},search_pair=search_pair,
                 preregistered_decision=decision,
                 replay_levels=replay_gaps,canonical_reproduction=reproduction,
                 full_vs_floor_max_point_separation_mm=region_diff,roles=source["roles"],
                 profile_coordinate="Baseline parent vertex arclength plus local fraction on the corresponding material chord. Curves show receiving-cell means; reported peaks use refined point samples.")
    output.mkdir(parents=True,exist_ok=False)
    _csv(output/"regional-power.csv",regional)
    _csv(output/"replay.csv",replay_rows)
    plots(output,series)
    if len(control_series)==3:
        plots(output,control_series,"attribution-",dict(baseline="Supplied wall",floor="Matched floor ablation",
              full="With nonfloor displacement"))
        summary["attribution_plot_roles"]={"with_nonfloor":contrast[0],"matched_floor_ablation":contrast[1],
                                           "level":"cubic1280"}
    summary["analysis_sha256"]={p.name:file_hash(p) for p in sorted(output.iterdir())}
    (output/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf8")
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("full","floor","replay","output"):
        parser.add_argument("--"+arg,type=Path,required=True)
    run(**vars(parser.parse_args()))


if __name__=="__main__":
    main()
