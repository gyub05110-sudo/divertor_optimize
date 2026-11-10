"""Exercise the sealed numerical analyzer with an independently relocated miniature dataset."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from verification import stage5_metrics as metrics
from verification.test_stage5_metrics import WALL, REGIONS, native_fixture


def put(path, value):
    path.parent.mkdir(exist_ok=True,parents=True)
    path.write_text(json.dumps(value),encoding="utf-8")
    return metrics.sha(path)


def test_relocated_geometry_and_absent_mesh_replay_identical_sealed_cases(tmp_path):
    from verification import stage5_portable_metrics as replay
    src=tmp_path/"before"
    restored=tmp_path/"elsewhere"
    written=tmp_path/"replayed"
    cid="000000000a"
    inp=src/"inputs"/cid
    base=src/"native/cases"/cid
    w=hashlib.sha256(WALL.astype("<f8").tobytes()).hexdigest()
    inp.mkdir(parents=True)
    regdir="runs/"+cid+"/30713/hotspot/selected"
    measured=src/"stage4"/regdir
    measured.mkdir(parents=True)
    area=np.pi*(WALL[:-1,0]+WALL[1:,0])*np.linalg.norm(np.diff(WALL,axis=0),axis=1)
    fraction=np.full(4,.25)
    np.savez(measured/"boundary.npz",wall=WALL,assigned_power_W=1e6,receiving_wall=WALL,cell_parent=np.arange(4),
             cell_power=fraction,area_edges=np.r_[0.,area.cumsum()],q_MW_m2=[2.,7.,3.,1.],sample_parent=np.arange(4),
             points=.5*(WALL[1:]+WALL[:-1]))
    (measured/"geometry.csv").write_text("vertex,R_m,Z_m,outgoing_region\n"+
        "".join(f"{i},{r},{z},{name}\n" for i,((r,z),name) in enumerate(zip(WALL,[*REGIONS,""]))))
    file_hashes={f.name:metrics.sha(f) for f in measured.iterdir()}
    (inp/"walls").mkdir()
    wall_rel=f"walls/{w}.wall"
    np.savetxt(inp/wall_rel,WALL)
    (inp/"OMP.source").write_text("synthetic")
    files={n:metrics.sha(inp/n) for n in (wall_rel,"OMP.source")}
    roles=[]; registration=[]
    for name in ("baseline",*metrics.FORCED,"hotspot","full"):
        roles.append(dict(role=name,seed=30713,wall_id=w,stage4_directory=str(measured),files_sha256=file_hashes,
                          wall_file=wall_rel,stage4_peak_MW_m2=7.,stage4_peak_region="shelf",
                          stage4_A2_m2=1/np.sum(fraction**2/area),candidate_id=None))
        registration.append(dict(case_id=cid,role=name,optimizer_seed=None if name=="baseline" else 30713,wall_id=w,
                                 measurement_sha256=file_hashes,stage4_directory=regdir))
    case=dict(case_id=cid,roles=roles,reference_SI={"power_W":1e6},predicted_family="floor-only",predicted_region="floor",
              baseline_peak_region="shelf",files_sha256=files)
    case_sha=put(inp/"case.json",case)
    (base/"mesh").mkdir(parents=True)
    mesh=base/"mesh/mmesh.nc"
    mesh.write_bytes(b"omitted native mesh")
    mesh_sha=metrics.sha(mesh)
    receipt=put(base/"mesh.json",dict(mesh_sha256=mesh_sha,mesh_bytes=mesh.stat().st_size))
    ready=put(base/"ready.json",dict(mesh_sha256=mesh_sha,mesh_receipt_sha256=receipt,case_sha256=case_sha))
    (base/(w+".pfc.nc")).write_bytes(b"native PFC")
    for seed in metrics.SEEDS:
        raw=src/"native/results"/cid/w/f"{seed}.nc"
        raw.parent.mkdir(exist_ok=True,parents=True)
        native_fixture(raw)
        _,_,ledger=metrics.read_native(raw,1e6)
        put(raw.with_suffix(".json"),dict(case_id=cid,wall_id=w,seed=seed,histories=4096,ranks=4,raw_sha256=metrics.sha(raw),
            ready_sha256=ready,mesh_sha256=mesh_sha,source_sha256=files["OMP.source"],wall_sha256=files[wall_rel],
            pfc_sha256=metrics.sha(base/(w+".pfc.nc")),**ledger))
    weights,freq=metrics.bootstrap_counts(6)
    original=metrics.analyze_case(src,case,case_sha,weights,freq)
    panel_sha=put(src/"inputs/panel.json",dict(cases=[case],case_record_sha256={cid:case_sha}))
    reg_sha=put(src/"wall-registration.json",dict(roles=registration))
    engine="verification/stage5_metrics.py"
    dest=src/"source-code"/engine; dest.parent.mkdir(parents=True)
    shutil.copyfile(metrics.__file__,dest)
    plan=dict(code={engine:metrics.sha(dest)},registration_sha256=reg_sha,files={
              "inputs/panel.json":panel_sha,f"native/cases/{cid}/ready.json":ready})
    seal=put(src/"execution.json",plan)
    (src/"execution.sha256").write_text(seal+"\n")
    put(src/"metrics.json",dict(cases=[original],analyzer_sha256=metrics.sha(dest),registration_sha256=reg_sha,
                                panel_sha256=panel_sha))
    src.rename(restored)
    (restored/"native/cases"/cid/"mesh/mmesh.nc").unlink()
    assert not measured.exists()
    audit=replay.replay(restored,written,registration_sha=reg_sha,seal_sha=seal)
    assert audit["case_count"]==1 and audit["case_wall_count"]==1 and audit["seed_jobs"]==6
    assert audit["exact_cases"] and audit["exact_profile_bytes"]==1 and audit["mesh_digests_from_receipts"]==[cid]
    original_profile=restored/"profiles"/cid/f"{w}.csv.gz"
    assert (written/"profiles"/cid/f"{w}.csv.gz").read_bytes()==original_profile.read_bytes()
    with pytest.raises(ValueError,match="exists"):
        replay.replay(restored,written,registration_sha=reg_sha,seal_sha=seal)
