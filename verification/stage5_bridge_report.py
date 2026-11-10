"""Recheck native launch positions, material vertices, field probes and strike receipts and tabulate ingress."""
import argparse
import csv
import hashlib
from pathlib import Path

import h5py
import numpy as np

from verification import stage5_campaign as campaign, stage5_native as native
from verification.stage3_campaign import read, require, sha, public_save


def case_report(study, case):
    cid=case["case_id"]
    inp=study/"inputs"/cid
    out=study/"native/cases"/cid
    ready=read(out/"ready.json")
    require(sha(inp/"case.json")==ready["case_sha256"] and read(inp/"case.json")==case and
            sha(out/"mesh/mmesh.nc")==ready["mesh_sha256"],"input or mesh differs")
    d=np.genfromtxt(inp/"OMP.csv",names=True,delimiter=",")
    source=np.loadtxt(inp/"OMP.source",skiprows=1)
    field=read(out/"field_audit.json")
    require(field==ready["field_audit"] and sha(inp/"OMP.source")==ready["source_sha256"],
            "changed field audit or CDF")
    require(len(d)==len(source)==1280 and (np.diff(source[:,2])>=0).all() and source[-1,2]==1 and
            np.array_equal(source[:,:2],np.column_stack([d["R_m"],d["Z_m"]])),"native launch distribution differs")
    assigned=float(d["power_W"].sum())
    cdf_delta=float(np.max(np.abs(np.cumsum(d["power_W"])/assigned-source[:,2])))
    require(abs(assigned-1e6)<=1e-5 and cdf_delta<=1e-12,"unequal inlet power")
    roles=case["roles"]
    wall_rows=[]
    for ident, vr in ready["walls"].items():
        ix=[r for r in roles if r["wall_id"]==ident]
        require(ix and sha(out/(ident+".pfc.nc"))==vr["pfc_sha256"],"PFC not bound to registered wall")
        with h5py.File(out/(ident+".pfc.nc"),"r") as file:
            coords=np.ascontiguousarray(file["wall/rz"][:],dtype="<f8")
        require(hashlib.sha256(coords.tobytes()).hexdigest()==ident,"native PFC material vertex bytes differ")
        stats=native.strike_audit(study,cid,ident)
        require(stats==vr["strikes"] and stats["pass"],"lower-leg check failed")
        rays=np.loadtxt(out/(ident+".strikes.dat"))
        direct=np.loadtxt(inp/"walls"/(ident+".hits"))
        col=int(stats["toroidal_direction"]>0)
        hits=rays[:,3+2*col:5+2*col]
        lower=np.isfinite(hits).all(axis=1) & (hits[:,0]>0) & (hits[:,1]<d["Z_m"])
        tail=lower & (np.linalg.norm(hits-direct[:,:2],axis=1)>.005)
        wall_rows.append(dict(wall_id=ident,count_over_5mm=int(tail.sum()),
                              donor_W_over_5mm=float(d["power_W"][tail].sum()),
                              roles=[(r["role"],r["seed"]) for r in ix],
                              vertices=len(coords),native_pfc_sha256=vr["pfc_sha256"],**stats))
    require(set(ready["walls"])=={r["wall_id"] for r in roles},"incomplete wall inventory")
    t=case["topology"]
    return dict(case_id=cid,representative=cid in ("3264fdd71f","70481becc0"),baseline_hotspot=case["baseline_peak_region"],
                field_probes=field["count"],max_abs_error_B_T=field["max_abs_error_B_T"],
                max_abs_error_psiN=field["max_abs_error_psiN"],signed_F_T_m=t["F_T_m"],
                axis_RZ_m=t["axis_RZ_m"],xpoint_RZ_m=t["xpoint_RZ_m"],lcfs_sha256=t["lcfs_sha256"],
                equilibrium_sha256=case["raw_sha256"],mesh_sha256=ready["mesh_sha256"],ready_sha256=sha(out/"ready.json"),
                source_sha256=ready["source_sha256"],source_W=assigned,donors=len(d),cdf_max_difference=cdf_delta,
                distinct_walls=len(wall_rows),role_count=len(roles),
                minimum_captured_lower_W=min(r["captured_lower_W"] for r in wall_rows),
                maximum_unmapped_count=max(r["unmapped"] for r in wall_rows),
                maximum_p95_hit_error_mm=1000*max(r["weighted_p95_hit_error_m"] for r in wall_rows),
                maximum_hit_error_mm=1000*max(r["max_hit_error_m"] for r in wall_rows),
                maximum_mean_abs_relative_length_error_pct=100*max(r["weighted_mean_abs_relative_length_error"] for r in wall_rows),
                walls=wall_rows)


def run(study):
    p=study/"inputs/panel.json"
    panel=read(p)
    registered=read(study/"wall-registration.json")
    campaign.verify_panel(registered["roles"],panel)
    cases=[case_report(study,c) for c in panel["cases"]]
    out=dict(schema="nf-stage5-input-bridge.v1",panel_sha256=sha(p),registration_sha256=sha(study/"wall-registration.json"),
             representative_cases=["3264fdd71f","70481becc0"],input_cdf_correction="Round-off in terminal donor CDF only; exact Stage 4 donor watts preserved.",
             historical_smoke=read(study/"native/historical_smoke/audit.json"),cases=cases)
    public_save(study/"bridge.json",out)
    cols=("case_id","representative","baseline_hotspot","signed_F_T_m","field_probes","max_abs_error_psiN","source_W",
          "donors","distinct_walls","role_count","cdf_max_difference","minimum_captured_lower_W","maximum_unmapped_count",
          "maximum_p95_hit_error_mm","maximum_hit_error_mm","maximum_mean_abs_relative_length_error_pct","ready_sha256")
    file=study/"bridge-cases.csv"
    tmp=file.with_suffix(".pending")
    with tmp.open("w",newline="",encoding="utf-8") as f:
        wr=csv.DictWriter(f,fieldnames=cols,extrasaction="ignore",lineterminator="\n")
        wr.writeheader()
        wr.writerows(cases)
    tmp.replace(file)
    print(f"Verified {len(cases)} equilibrium inputs and {sum(c['distinct_walls'] for c in cases)} material walls")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,default=campaign.STUDY)
    run(parser.parse_args().study.resolve())
