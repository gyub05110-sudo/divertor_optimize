"""Package and verify the completed Stage 5 native study using Python 3.11+ standard library."""
import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
from tempfile import TemporaryDirectory


HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
NAME="nf-stage5-firefly-20261110"
MANIFEST="stage5-manifest.json"
SEEDS=(51031,51047,51059,51071,51089,51109)
REGISTRATION="d0627f34407ee7cb3fd054a6c8eecef17b1e20026c28ce068fdc3f0d224b7102"
CODEC_SHA="1cac0d886d7aa125253a085e96f6c8639bcefd3940e0e12a00b09267b7e1f742"


def load_codec():
    path=HERE/"payload/package_evidence.py"
    if not path.is_file():
        path=ROOT/"results/nf-stage3-20261108/package_evidence.py"
    if hashlib.sha256(path.read_bytes().replace(b"\r\n",b"\n")).hexdigest()!=CODEC_SHA:
        raise ValueError("changed split-ZIP codec")
    spec=importlib.util.spec_from_file_location("stage5_zip_codec",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


codec=load_codec()
require=codec.require


def job_inventory(rows):
    jobs=[tuple(j) for j in rows]
    require(len(jobs)==len(set(jobs))==276,"incomplete or duplicate native job inventory")
    pairs={(case,wall) for case,wall,seed in jobs}
    require(len(pairs)==46 and len({case for case,wall in pairs})==7 and
            all(re.fullmatch(r"[a-z0-9]{10}",case) and re.fullmatch(r"[0-9a-f]{64}",wall) for case,wall in pairs)
            and set(jobs)=={(case,wall,seed) for case,wall in pairs for seed in SEEDS},
            "native job case, wall or paired seed differs")
    return jobs


def check_job(e,job,execution_sha,bindings):
    case,wall,seed=job
    stem=f"native/results/{case}/{wall}/{seed}"
    raw,receipt=stem+".nc",stem+".json"
    launch_receipt=e.read(f"receipts/{case}/{wall}/{seed}.json")
    require(launch_receipt["job"]==list(job) and launch_receipt["execution_sha256"]==execution_sha and
            launch_receipt["returncode"]==0 and launch_receipt["raw_path"]==raw and
            launch_receipt["native_receipt_path"]==receipt,"dispatch identity or artifact path differs")
    raw_sha=e.check(raw,launch_receipt["raw_sha256"])
    receipt_sha=e.check(receipt,launch_receipt["native_receipt_sha256"])
    attempt=codec.safe_path(launch_receipt["attempt"])
    require(re.fullmatch(f"dispatch-attempts/{case}/{wall}/{seed}/[0-9]+",attempt),"dispatch attempt path differs")
    e.check(attempt+"/launch.json",launch_receipt["launch_sha256"])
    e.check(attempt+"/completion.json",launch_receipt["completion_sha256"])
    e.check(attempt+"/stdout.log")
    launch=e.read(attempt+"/launch.json")
    require(launch["job"]==list(job) and launch["execution_sha256"]==execution_sha and
            e.read(attempt+"/completion.json")["returncode"]==0,"launch or completion provenance differs")
    native=e.read(receipt)
    expected=dict(bindings,case_id=case,wall_id=wall,seed=seed,histories=4096,ranks=4,raw_sha256=raw_sha)
    require(all(native.get(k)==v for k,v in expected.items()),"native receipt or raw SHA differs")
    require(native["raw_file"]==f"{seed}.nc" and re.fullmatch(r"[0-9]+",native["attempt"]),
            "native attempt or raw filename differs")
    native_attempt=f"native/results/{case}/{wall}/attempts/{seed}/"+native["attempt"]+"/"
    e.check(stem+".log",native["log_sha256"])
    e.check(native_attempt+"run.log",native["log_sha256"])
    e.check(native_attempt+"native.runtime.json",native["runtime_sha256"])
    request=e.read(native_attempt+"request.json")
    require(e.read(native_attempt+"native.runtime.json")==native["runtime"] and
            e.read(native_attempt+"publish.json")==native and all(native.get(k)==v for k,v in request.items()),
            "original native job request, runtime or final audit differs")
    return dict(seed=seed,raw_sha256=raw_sha,receipt_sha256=receipt_sha)


class Evidence:
    def __init__(self,files,content):
        self.files,self.content=files,content

    def check(self,name,sha=None):
        key=codec.safe_path(name)
        require(key in self.files and (sha is None or self.files[key]["sha256"]==sha),"missing/changed evidence: "+key)
        return self.files[key]["sha256"]

    def read(self,name):
        self.check(name)
        return codec.as_json(self.content(name))


def validate(e):
    reg_sha=e.check("wall-registration.json",REGISTRATION)
    reg=e.read("wall-registration.json")
    jobs=job_inventory(reg["jobs"])
    require(reg["seeds"]==list(SEEDS) and reg["histories_per_seed"]==4096 and reg["mpi_ranks"]==4,
            "registered heat settings differ")
    plan=e.read("execution.json")
    seal=e.check("execution.json")
    require(seal==e.content("execution.sha256").decode().strip() and plan["jobs"]==reg["jobs"] and
            plan["registration_sha256"]==reg_sha and plan["protocol_sha256"]==reg["protocol_sha256"],
            "native execution seal differs")
    e.check("publication/protocol.ko.md",reg["protocol_sha256"])
    e.check("publication/wall-registration.json",reg_sha)
    for name in ("firefly-forecast.csv","firefly-products.csv","firefly-comparisons.pdf","firefly-comparisons.png",
                 "REPORT.ko.md","summarize.py","test_summarize.py"):
        e.check("publication/"+name)
    for key,sha in plan["code"].items():
        e.check("source-code/"+codec.safe_path(key),sha)
    for key,sha in plan["files"].items():
        e.check(key,sha)
    panel=e.read("inputs/panel.json")
    for name,sha in panel["stage4_source_sha256"].items():
        e.check("stage4/"+codec.safe_path(name),sha)
    e.check("stage4/stage4_report.json",reg["stage4_report_sha256"])
    metrics=e.read("metrics.json")
    require(metrics["schema"]=="nf-stage5-whole-wall-metrics.v1" and metrics["eligible_files"]==7 and
            metrics["role_count"]==51 and metrics["distinct_case_walls"]==46 and metrics["native_seeds"]==list(SEEDS)
            and metrics["registration_sha256"]==reg_sha and metrics["panel_sha256"]==e.check("inputs/panel.json"),
            "complete native metrics and current panel required")
    e.check("analysis-code/verification/stage5_metrics.py",metrics["analyzer_sha256"])
    derived={v["case_id"]:v for v in metrics["cases"]}
    posted={(v["case_id"],v["role"],v["optimizer_seed"]):v for v in reg["roles"]}
    seen_roles=set()
    checked_jobs=set()
    for c in panel["cases"]:
        case=c["case_id"]
        inp=f"inputs/{case}/"
        base=f"native/cases/{case}/"
        case_sha=e.check(inp+"case.json",panel["case_record_sha256"][case])
        require(e.read(inp+"case.json")==c,"case record differs from panel")
        ready_sha=e.check(base+"ready.json")
        ready=e.read(base+"ready.json")
        require(ready["case_id"]==case and ready["case_sha256"]==case_sha and
                ready["source_sha256"]==e.check(inp+"OMP.source") and ready["source_csv_sha256"]==e.check(inp+"OMP.csv")
                and ready["material_inputs"]==c["files_sha256"],"native readiness bindings differ")
        for name,sha in c["files_sha256"].items():
            e.check(inp+codec.safe_path(name),sha)
        for name,sha in ready["native_libraries"].items():
            e.check("historical-runtime/"+codec.safe_path(name),sha)
        e.check(base+"mesh.json",ready["mesh_receipt_sha256"])
        mesh=e.read(base+"mesh.json")
        require(mesh["mesh_sha256"]==ready["mesh_sha256"]==e.check(base+"mesh/mmesh.nc")
                and mesh["mesh_bytes"]==e.files[base+"mesh/mmesh.nc"]["bytes"],"omitted mesh reference differs")
        e.check(base+"mesh/mmesh.ctr",mesh["settings"]["mesh_control_sha256"])
        e.check(base+"model/wall.txt",mesh["settings"]["casing_sha256"])
        e.check(base+"field_native.txt")
        require(ready["field_audit"]==e.read(base+"field_audit.json"),"field comparison differs")
        e.check("stage4/inputs/raw/d3d_shot_"+case+".parquet",c["raw_sha256"])
        for role in c["roles"]:
            rk=(case,role["role"],None if role["role"]=="baseline" else role["seed"])
            fixed=posted[rk]
            require(rk not in seen_roles and fixed["wall_id"]==role["wall_id"] and
                    fixed["measurement_sha256"]==role["files_sha256"],"registered role or wall differs")
            seen_roles.add(rk)
            directory="stage4/"+codec.safe_path(fixed["stage4_directory"])
            for name,sha in role["files_sha256"].items():
                e.check(directory+"/"+codec.safe_path(name),sha)
            run=directory.rsplit("/",1)[0]+"/"
            for name,key in (("result.json","result_sha256"),("events.jsonl","events_sha256"),("receipt.json","receipt_sha256")):
                e.check(run+name,fixed[key])
        walls={r["wall_id"] for r in c["roles"]}
        require(set(ready["walls"])==set(derived[case]["geometries"])==walls,"wrong native or metric wall panel")
        for wall in sorted(walls):
            asset=ready["walls"][wall]
            e.check(base+wall+".pfc.nc",asset["pfc_sha256"])
            e.check(base+wall+".open.nc",asset["open_pfc_sha256"])
            e.check(inp+"walls/"+wall+".wall",asset["wall_sha256"])
            strikes=asset["strikes"]
            require(e.read(base+wall+".strikes.json")==strikes and strikes["pass"],"strike audit differs")
            e.check(base+wall+".strikes.dat",strikes["strikes_sha256"])
            e.check(inp+"walls/"+wall+".hits",strikes["reference_sha256"])
            geometry=derived[case]["geometries"][wall]
            profile=f"profiles/{case}/{wall}.csv.gz"
            require(geometry["profile"]==profile,"portable ring profile path differs")
            e.check(profile,geometry["profile_sha256"])
            bindings=dict(source_sha256=ready["source_sha256"],wall_sha256=asset["wall_sha256"],
                          mesh_sha256=ready["mesh_sha256"],pfc_sha256=asset["pfc_sha256"],ready_sha256=ready_sha)
            raw_ledgers=[]
            for seed in SEEDS:
                job=(case,wall,seed)
                checked_jobs.add(job)
                raw_ledgers.append(check_job(e,job,seal,bindings))
            require(raw_ledgers==[{k:s[k] for k in ("seed","raw_sha256","receipt_sha256")} for s in geometry["seeds"]],
                    "pooled metrics do not bind every raw native seed")
    require(checked_jobs==set(jobs) and len(seen_roles)==len(posted)==51 and len(derived)==7,
            "registered case, geometry or seed missing")
    raw_names={p for p in e.files if p.startswith("native/results/") and p.endswith(".nc") and "/attempts/" not in p}
    require(raw_names=={f"native/results/{c}/{w}/{s}.nc" for c,w,s in jobs},"unexpected extra or missing raw native heat file")
    smoke=e.read("native/historical_smoke/audit.json")
    require(smoke["raw_sha256"]==smoke["historical_raw_sha256"]==e.check("native/historical_smoke/historical.pending.nc"),
            "historical replay digest differs")


def collect(study,stage4,repo):
    require((study/"metrics.json").is_file(),"complete metrics required before packaging")
    files={}
    permitted={"inputs","source-code","dispatch-attempts","receipts","profiles","native","input-cdf-correction"}
    for p in sorted(study.rglob("*")):
        if not p.is_file():
            continue
        rel=p.relative_to(study)
        if len(rel.parts)>1 and rel.parts[0] not in permitted:
            continue
        if any(v in ("__pycache__","runtime") for v in rel.parts) or p.suffix in (".pyc",".lock",".partial"):
            continue
        if ".pending" in p.suffixes and rel.parent.as_posix()!="native/historical_smoke":
            continue
        if "/attempts/" in rel.as_posix() and p.suffix==".nc":
            continue
        files[rel.as_posix()]=p
    reg=codec.as_json(files["wall-registration.json"].read_bytes())
    panel=codec.as_json(files["inputs/panel.json"].read_bytes())
    paper=repo/"results"/NAME
    for p in paper.iterdir():
        if p.is_file() and not p.name.startswith("evidence") and p.suffix.lower() not in (".zip",".pending"):
            files["publication/"+p.name]=p
    for p in (repo/"verification").glob("*stage5*"):
        if p.is_file() and p.suffix in (".py",".sh",".ps1"):
            files["analysis-code/verification/"+p.name]=p
    for name in panel["stage4_source_sha256"]:
        files["stage4/"+codec.safe_path(name)]=stage4/name
    for r in reg["roles"]:
        rel=codec.safe_path(r["stage4_directory"])
        for name in r["measurement_sha256"]:
            files["stage4/"+rel+"/"+codec.safe_path(name)]=stage4/rel/name
        folder=rel.rsplit("/",1)[0]
        for name in ("events.jsonl","result.json","receipt.json"):
            files["stage4/"+folder+"/"+name]=stage4/folder/name
    historical=study.parent/"firefly-baseline-20261105"
    for c in panel["cases"]:
        case=c["case_id"]
        name="inputs/raw/d3d_shot_"+case+".parquet"
        files["stage4/"+name]=stage4/name
        ready=codec.as_json(files[f"native/cases/{case}/ready.json"].read_bytes())
        for name in ready["native_libraries"]:
            files["historical-runtime/"+codec.safe_path(name)]=historical/name
    for name in ("firefly-midplane.patch","enter.sh","build_compat.sh","REPRODUCE.md"):
        files["historical-runtime/"+name]=historical/name
    for name in ("d3d_optimization/data/DIIID_geom.json","baseline_review/baseline_review.json",
                 "baseline_review/optimization_region_policy.json"):
        files["fixed-reference/"+name]=repo/name
    codec.check_paths(files)
    return files


def stored(name,mode):
    if re.fullmatch(r"native/cases/[0-9a-z]{10}/mesh/mmesh\.nc",name):
        return False
    if mode=="full":
        return True
    return (name.startswith(("publication/","profiles/","source-code/","analysis-code/")) or
            name.endswith((".json",".sha256",".log",".md",".py",".sh",".ctr")) or
            ("/" not in name and name.endswith(".csv")))


def build(study,bundle,mode,stage4=None,repo=ROOT):
    study,bundle,repo=map(lambda p:Path(p).resolve(),(study,bundle,repo))
    require(mode in ("compact","full"),"invalid package mode")
    require(not bundle.exists(),"bundle destination exists")
    catalog=collect(study,stage4 or study.parent/"nf-stage4-20261109",repo)
    files={k:dict(bytes=p.stat().st_size,sha256=codec.digest(p)) for k,p in sorted(catalog.items())}
    validate(Evidence(files,lambda k:catalog[k].read_bytes()))
    subset={k:p for k,p in catalog.items() if stored(k,mode)}
    bundle.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=bundle.parent,prefix=".stage5-evidence-") as temp:
        stage=Path(temp)/"bundle"
        payload=stage/"payload"
        payload.mkdir(parents=True)
        codec.write_payload(payload,subset,{})
        shutil.copyfile(__file__,stage/"stage5_evidence.py")
        manifest=dict(schema="nf-stage5-portable-evidence.v1",mode=mode,registration_sha256=REGISTRATION,
                      inventories=dict(cases=7,roles=51,case_walls=46,paired_seeds=6,native_runs=276),
                      inventory_sha256=codec.canonical(files),files=files,payload_sha256=codec.digest(payload/codec.MANIFEST),
                      verifier_sha256=codec.digest(stage/"stage5_evidence.py"),
                      omissions=dict(mesh="All seven native mesh/mmesh.nc files remain in the original D-drive study. "
                                     "Rebuild with the included case model, mmesh.ctr and stage5_native_worker.py mesh action.",
                                     runtime="The historical OCI/chroot image, compiler, MPI installation and transient caches are external. "
                                     "Executable hashes and FIREFLY adapter patch accompany the archive.",
                                     compact="Compact mode additionally omits bulk native NC, Stage 4 arrays, raw parquet and mesh diagnostics. "
                                     "Match inventory_sha256 to the full compressed bundle for every input, log, profile and heat file."))
        (stage/MANIFEST).write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf8",newline="\n")
        sha=codec.digest(stage/MANIFEST)
        (stage/(MANIFEST+".sha256")).write_text(sha+"\n",encoding="ascii",newline="\n")
        with verified(stage,sha):
            pass
        require(not bundle.exists(),"bundle destination exists")
        stage.rename(bundle)
    return sha


@contextmanager
def verified(bundle,expected=None):
    bundle=Path(bundle)
    sha=codec.digest(bundle/MANIFEST)
    require((bundle/(MANIFEST+".sha256")).read_text().strip()==sha and (expected is None or expected==sha),
            "portable manifest SHA-256 differs")
    m=codec.as_json((bundle/MANIFEST).read_bytes())
    require(m["schema"]=="nf-stage5-portable-evidence.v1" and m["mode"] in ("full","compact")
            and m["registration_sha256"]==REGISTRATION and m["verifier_sha256"]==codec.digest(bundle/"stage5_evidence.py"),
            "portable verifier or registration differs")
    files=m["files"]
    codec.check_paths(files)
    require(codec.canonical(files)==m["inventory_sha256"] and all(v["bytes"]>=0 and re.fullmatch(r"[0-9a-f]{64}",v["sha256"])
            for v in files.values()),"invalid portable file inventory")
    with codec.open_payload(bundle/"payload",m["payload_sha256"]) as reader:
        require(set(reader.files)=={k for k in files if stored(k,m["mode"])} and not reader.external,
                "payload inventory differs from selected package mode")
        for k,v in reader.files.items():
            require(v["bytes"]==files[k]["bytes"] and v["sha256"]==files[k]["sha256"],"payload SHA-256 differs: "+k)
        validate(Evidence(files,reader.content))
        yield m,reader


def restore(bundle,out,expected=None):
    bundle,out=Path(bundle).resolve(),Path(out).resolve()
    codec.disjoint(bundle,out)
    require(not out.exists(),"restore destination exists")
    out.parent.mkdir(parents=True,exist_ok=True)
    with verified(bundle,expected) as (_,reader), TemporaryDirectory(dir=out.parent,prefix=".stage5-restore-") as temp:
        stage=Path(temp)/"study"
        stage.mkdir()
        for name,record in reader.files.items():
            dest=stage/name
            dest.parent.mkdir(parents=True,exist_ok=True)
            with dest.open("xb") as h:
                for block in reader.chunks(name):
                    h.write(block)
            require(codec.digest(dest)==record["sha256"],"restored file differs: "+name)
        stage.rename(out)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest="action",required=True)
    build_cli=commands.add_parser("build")
    build_cli.add_argument("--study",type=Path,required=True)
    build_cli.add_argument("--bundle",type=Path,required=True)
    build_cli.add_argument("--mode",choices=("full","compact"),required=True)
    build_cli.add_argument("--stage4",type=Path)
    for action in ("verify","restore"):
        cli=commands.add_parser(action)
        cli.add_argument("--bundle",required=True,type=Path)
        cli.add_argument("--manifest-sha256")
        if action=="restore":
            cli.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    if args.action=="build":
        print(build(args.study,args.bundle,args.mode,args.stage4))
    elif args.action=="restore":
        restore(args.bundle,args.out,args.manifest_sha256)
    else:
        with verified(args.bundle,args.manifest_sha256) as (manifest,raw):
            print("Verified",manifest["mode"],"package and 276 native receipts.",len(raw.files),"stored files.",
                  "Inventory SHA-256",manifest["inventory_sha256"])
