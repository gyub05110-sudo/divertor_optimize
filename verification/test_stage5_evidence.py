"""Portability, job identity, native digest and incomplete-publication regression tests."""
from importlib import import_module
import json
import shutil

import pytest


def evidence():
    return import_module("verification.stage5_evidence")


def jobs():
    return [[f"{i:010x}",f"{j:064x}",s] for i,count in enumerate((7,7,8,6,6,6,6))
            for j in range(count) for s in (51031,51047,51059,51071,51089,51109)]


def test_native_inventory_rejects_equal_size_duplicate_and_wrong_seed():
    pack=evidence()
    rows=jobs()
    assert len(pack.job_inventory(rows))==276
    for bad in ([*rows[:-1],rows[0]], [*rows[:-1],[*rows[-1][:2],42]], rows[1:]):
        with pytest.raises(ValueError,match="job"):
            pack.job_inventory(bad)


def test_receipts_bind_raw_bytes_not_old_absolute_paths(tmp_path):
    pack=evidence()
    case,wall,seed=jobs()[0]
    job=(case,wall,seed)
    raw=f"native/results/{case}/{wall}/{seed}"
    attempt=f"dispatch-attempts/{case}/{wall}/{seed}/001"
    files={}
    def write(name,data):
        path=tmp_path/name
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(data if isinstance(data,bytes) else json.dumps(data).encode())
        files[name]=path
        return pack.codec.digest(path)
    h=write(raw+".nc",b"native heat")
    log_hash=write(raw+".log",b"native stdout")
    native_attempt=raw.rsplit("/",1)[0]+f"/attempts/{seed}/001"
    runtime={"ranks":4,"library_mappings":["/work/install/lib/libFIREFLY.so"]}
    runtime_hash=write(native_attempt+"/native.runtime.json",runtime)
    bind=dict(raw_sha256=h,case_id=case,wall_id=wall,seed=seed,source_sha256="source",
              wall_sha256="wall",pfc_sha256="pfc",mesh_sha256="mesh",ready_sha256="ready",
              histories=4096,ranks=4,log_sha256=log_hash,runtime=runtime,runtime_sha256=runtime_hash,
              attempt="001",raw_file=f"{seed}.nc")
    write(native_attempt+"/run.log",b"native stdout")
    write(native_attempt+"/request.json",{"seed":seed})
    write(native_attempt+"/publish.json",bind)
    receipt=write(raw+".json",bind)
    launch=write(attempt+"/launch.json",{"job":list(job),"execution_sha256":"sealed"})
    end=write(attempt+"/completion.json",{"returncode":0})
    write(attempt+"/stdout.log",b"log")
    write(f"receipts/{case}/{wall}/{seed}.json",dict(job=list(job),execution_sha256="sealed",
          returncode=0,raw_path=raw+".nc",raw_sha256=h,native_receipt_path=raw+".json",
          native_receipt_sha256=receipt,attempt=attempt,launch_sha256=launch,completion_sha256=end))
    copy=tmp_path/"bundle"; copy.mkdir()
    pack.codec.write_payload(copy,files,{})
    for p in files.values():
        p.unlink()
    with pack.codec.open_payload(copy) as reader:
        rec=pack.check_job(reader,job,"sealed",{k:bind[k] for k in
                   ("source_sha256","wall_sha256","pfc_sha256","mesh_sha256","ready_sha256")})
        assert rec==dict(seed=seed,raw_sha256=h,receipt_sha256=receipt)
        before=reader.files[raw+".log"]["sha256"]
        reader.files[raw+".log"]["sha256"]="wrong"
        with pytest.raises(ValueError,match="evidence"):
            pack.check_job(reader,job,"sealed",{})
        reader.files[raw+".log"]["sha256"]=before
        del reader.files[raw+".nc"]
        with pytest.raises(ValueError,match="evidence"):
            pack.check_job(reader,job,"sealed",{})


def test_refuses_to_package_without_metrics(tmp_path):
    pack=evidence()
    study=tmp_path/"incomplete";study.mkdir()
    with pytest.raises(ValueError,match="metrics"):
        pack.build(study,tmp_path/"bundle","full")
    assert not (tmp_path/"bundle").exists()


@pytest.mark.parametrize("mode",["full","compact"])
def test_portable_study_verifies_after_original_tree_is_removed(tmp_path,monkeypatch,mode):
    """Build all 276 synthetic receipts, 46 wall records and 51 source aliases."""
    pack=evidence()
    src=tmp_path/"original"; src.mkdir()
    catalog={}
    def put(name,value=b"sample"):
        data=value if isinstance(value,bytes) else json.dumps(value,sort_keys=True).encode()
        h=pack.codec.hashed(data)
        file=src/h
        if not file.exists():
            file.write_bytes(data)
        catalog[name]=file
        return h
    tasks=jobs()
    cids=sorted({t[0] for t in tasks})
    roles=[]
    cases=[]
    metrics_cases=[]
    fixed={}
    ready_for={}
    protocol=put("publication/protocol.ko.md")
    for name in ("firefly-forecast.csv","firefly-products.csv","firefly-comparisons.pdf","firefly-comparisons.png",
                 "REPORT.ko.md","summarize.py","test_summarize.py"):
        put("publication/"+name,b"publication")
    baseline_report=put("stage4/stage4_report.json",{})
    driver=put("analysis-code/verification/stage5_metrics.py")
    outputs={}
    for ci,cid in enumerate(cids):
        wall_ids=sorted({w for c,w,s in tasks if c==cid})
        names=[(30713,r) for r in ("baseline","hotspot","full","floor-only","floor-entry","floor-shelf","floor-connector")]
        if ci==2:
            names.extend([(30757,"hotspot"),(30757,"full")])
        inp=f"inputs/{cid}/"
        base=f"native/cases/{cid}/"
        inputs={"OMP.csv":put(inp+"OMP.csv"),"OMP.source":put(inp+"OMP.source")}
        mappings=[]
        assets={}
        geometries={}
        raw_source=put(f"stage4/inputs/raw/d3d_shot_{cid}.parquet")
        for i,(seed,role) in enumerate(names):
            w=wall_ids[i%len(wall_ids)]
            path=f"runs/{cid}/{seed}/{role}"
            record={"boundary.npz":put("stage4/"+path+"/selected/boundary.npz")}
            roles.append(dict(case_id=cid,role=role,optimizer_seed=None if role=="baseline" else seed,
                              wall_id=w,stage4_directory=path+"/selected",measurement_sha256=record,
                              result_sha256=put("stage4/"+path+"/result.json",{}),events_sha256=put("stage4/"+path+"/events.jsonl"),
                              receipt_sha256=put("stage4/"+path+"/receipt.json",{})))
            mappings.append(dict(role=role,seed=seed,wall_id=w,files_sha256=record))
        mesh_sha=put(base+"mesh/mmesh.nc",b"omitted mesh")
        field={"max_abs_error_psiN":0}
        put(base+"field_audit.json",field)
        put(base+"field_native.txt")
        mesh=dict(mesh_sha256=mesh_sha,mesh_bytes=len(b"omitted mesh"),settings=dict(
            mesh_control_sha256=put(base+"mesh/mmesh.ctr"),casing_sha256=put(base+"model/wall.txt")))
        for w in wall_ids:
            inputs["walls/"+w+".wall"]=put(inp+"walls/"+w+".wall")
            inputs["walls/"+w+".hits"]=put(inp+"walls/"+w+".hits")
            strikes=dict(strikes_sha256=put(base+w+".strikes.dat"),reference_sha256=inputs["walls/"+w+".hits"])
            strikes["pass"]=True
            put(base+w+".strikes.json",strikes)
            assets[w]=dict(pfc_sha256=put(base+w+".pfc.nc"),open_pfc_sha256=put(base+w+".open.nc"),
                           wall_sha256=inputs["walls/"+w+".wall"],strikes=strikes)
            profile=f"profiles/{cid}/{w}.csv.gz"
            geometries[w]=dict(seeds=[],profile=profile,profile_sha256=put(profile))
            outputs[cid,w]=geometries[w]
        case=dict(case_id=cid,roles=mappings,files_sha256=inputs,raw_sha256=raw_source)
        cases.append(case)
        case_hash=put(inp+"case.json",case)
        ready=dict(case_id=cid,case_sha256=case_hash,source_sha256=inputs["OMP.source"],source_csv_sha256=inputs["OMP.csv"],
                   material_inputs=inputs,native_libraries={},mesh_receipt_sha256=put(base+"mesh.json",mesh),
                   mesh_sha256=mesh_sha,field_audit=field,walls=assets)
        ready_hash=put(base+"ready.json",ready)
        fixed[base+"ready.json"]=ready_hash
        ready_for[cid]=ready,ready_hash
        metrics_cases.append(dict(case_id=cid,geometries=geometries))
    reg=dict(jobs=tasks,roles=roles,seeds=list(pack.SEEDS),histories_per_seed=4096,mpi_ranks=4,protocol_sha256=protocol,
             stage4_report_sha256=baseline_report)
    reg_sha=put("wall-registration.json",reg);put("publication/wall-registration.json",reg)
    monkeypatch.setattr(pack,"REGISTRATION",reg_sha)
    panel=dict(stage4_source_sha256={"stage4_report.json":baseline_report},cases=cases,
               case_record_sha256={c["case_id"]:pack.codec.digest(catalog[f"inputs/{c['case_id']}/case.json"]) for c in cases})
    panel_hash=put("inputs/panel.json",panel)
    seal=put("execution.json",dict(jobs=tasks,registration_sha256=reg_sha,protocol_sha256=protocol,code={},files=fixed))
    put("execution.sha256",(seal+"\n").encode())
    for cid,w,seed in tasks:
        folder=f"native/results/{cid}/{w}/"
        stem=folder+str(seed)
        sub=folder+f"attempts/{seed}/001/"
        ready,ready_sha=ready_for[cid]
        raw=put(stem+".nc",b"synthetic-hdf")
        log_sha=put(stem+".log",b"log");put(sub+"run.log",b"log")
        bind=dict(case_id=cid,wall_id=w,seed=seed,source_sha256=ready["source_sha256"],mesh_sha256=ready["mesh_sha256"],
                  wall_sha256=ready["walls"][w]["wall_sha256"],pfc_sha256=ready["walls"][w]["pfc_sha256"],ready_sha256=ready_sha)
        put(sub+"request.json",bind)
        native=dict(bind,raw_sha256=raw,raw_file=f"{seed}.nc",histories=4096,ranks=4,attempt="001",runtime={},
                    runtime_sha256=put(sub+"native.runtime.json",{}),log_sha256=log_sha)
        receipt=put(stem+".json",native);put(sub+"publish.json",native)
        attempt=f"dispatch-attempts/{cid}/{w}/{seed}/001"
        launch=put(attempt+"/launch.json",dict(job=[cid,w,seed],execution_sha256=seal))
        end=put(attempt+"/completion.json",{"returncode":0});put(attempt+"/stdout.log")
        put(f"receipts/{cid}/{w}/{seed}.json",dict(job=[cid,w,seed],execution_sha256=seal,returncode=0,
            attempt=attempt,launch_sha256=launch,completion_sha256=end,raw_sha256=raw,raw_path=stem+".nc",
            native_receipt_path=stem+".json",native_receipt_sha256=receipt))
        outputs[cid,w]["seeds"].append(dict(seed=seed,raw_sha256=raw,receipt_sha256=receipt))
    smoke=put("native/historical_smoke/historical.pending.nc")
    put("native/historical_smoke/audit.json",dict(raw_sha256=smoke,historical_raw_sha256=smoke))
    put("metrics.json",dict(schema="nf-stage5-whole-wall-metrics.v1",eligible_files=7,role_count=51,distinct_case_walls=46,
                            native_seeds=list(pack.SEEDS),registration_sha256=reg_sha,panel_sha256=panel_hash,
                            analyzer_sha256=driver,cases=metrics_cases))
    monkeypatch.setattr(pack,"collect",lambda study,stage4,repo:catalog)
    out=tmp_path/"bundle"
    sha=pack.build(src,out,mode)
    shutil.rmtree(src)
    with pack.verified(out,sha) as (manifest,reader):
        assert len(reader.files)<=len(manifest["files"])-7
        assert all("/mesh/mmesh.nc" not in name for name in reader.files)
        assert len(reader.blocks)<len(reader.files)
        name="native/results/"+"/".join(map(str,tasks[0]))+".nc"
        assert manifest["files"][name]["sha256"]==pack.codec.hashed(b"synthetic-hdf")
        if mode=="full":
            assert reader.content(name)==b"synthetic-hdf"
        else:
            assert name not in reader.files
        missing=dict(manifest["files"])
        del missing["publication/firefly-forecast.csv"]
        with pytest.raises(ValueError,match="evidence"):
            pack.validate(pack.Evidence(missing,reader.content))
