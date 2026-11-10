"""Replay the execution-sealed Stage 5 ring/seed analyzer on restored evidence, at any path."""
import argparse
from contextlib import contextmanager
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path, PurePosixPath
import sys
from tempfile import TemporaryDirectory


REGISTRATION="d0627f34407ee7cb3fd054a6c8eecef17b1e20026c28ce068fdc3f0d224b7102"
EXECUTION="3151ee0a78584acd70714847b3957a16730c13347053d3ef2d2add070770e023"
ENGINE="verification/stage5_metrics.py"


def require(condition,message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def relative_folder(parent,location):
    key=PurePosixPath(location)
    require(not key.is_absolute() and key.parts and ".." not in key.parts and ":" not in location
            and "\\" not in location,"invalid registered relative Stage 4 path")
    path=(parent/key).resolve()
    require(path.is_relative_to(parent.resolve()) and path.is_dir(),"restored Stage 4 measurement missing: "+str(path))
    return path


@contextmanager
def analyzer(restored,stage,plan,panel,registration):
    """Load a private copy of the sealed module; adapt file locations only."""
    source=restored/"source-code"/ENGINE
    require(digest(source)==plan["code"][ENGINE],"sealed analyzer code differs")
    spec=importlib.util.spec_from_file_location("_stage5_restored_analyzer",source)
    module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module
    try:
        # Keep the immutable source snapshot read-only, including on filesystems
        # where the interpreter could otherwise create a __pycache__ directory.
        exec(compile(source.read_bytes(),str(source),"exec"),module.__dict__)
        by_role={(r["case_id"],r["role"],r["optimizer_seed"]):r for r in registration["roles"]}
        require(len(by_role)==len(registration["roles"]),"duplicate registered role")
        locations={}
        meshes={}
        for case in panel["cases"]:
            cid=case["case_id"]
            folder=restored/"native/cases"/cid
            ready_path=folder/"ready.json"
            require(digest(ready_path)==plan["files"][f"native/cases/{cid}/ready.json"],"changed sealed ready dossier")
            ready=read(ready_path)
            require(digest(folder/"mesh.json")==ready["mesh_receipt_sha256"] and
                    ready["case_sha256"]==panel["case_record_sha256"][cid],"mesh receipt/case binding differs")
            mesh_sha=read(folder/"mesh.json")["mesh_sha256"]
            require(mesh_sha==ready["mesh_sha256"],"mesh hash differs across ready and mesh dossiers")
            mesh=folder/"mesh/mmesh.nc"
            if mesh.is_file():
                require(digest(mesh)==mesh_sha,"restored mesh bytes changed")
            else:
                meshes[mesh]=mesh_sha
            for role in case["roles"]:
                name=role["role"]
                fixed=by_role.pop((cid,name,None if name=="baseline" else role["seed"]))
                require(role["wall_id"]==fixed["wall_id"] and role["files_sha256"]==fixed["measurement_sha256"],
                        "registered Stage 4 wall or measurement digest differs")
                # Use the role in the hash-checked panel solely as a lookup key.
                # Its original host directory is never used as a filesystem path.
                key=id(role)
                locations[key]=relative_folder(restored/"stage4",fixed["stage4_directory"])
        require(not by_role,"registered role absent from input panel")
        original_read,original_write,original_sha=module.read_model,module.write_profile,module.sha
        profiles=restored/"profiles"

        def read_model(role,*args):
            return original_read(dict(role,stage4_directory=str(locations[id(role)])),*args)

        def write_profile(path,*args):
            original_write(stage/"profiles"/path.relative_to(profiles),*args)

        def sha(path):
            path=Path(path)
            if path in meshes:
                return meshes[path]
            if path.is_relative_to(profiles):
                return original_sha(stage/"profiles"/path.relative_to(profiles))
            return original_sha(path)

        module.read_model,module.write_profile,module.sha=read_model,write_profile,sha
        yield module,sorted(path.parent.parent.name for path in meshes)
    finally:
        sys.modules.pop(spec.name,None)


def compare_values(expected,actual,location="case"):
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and set(actual)==set(expected),"different fields at "+location)
        for k in expected:
            compare_values(expected[k],actual[k],location+"/"+k)
    elif isinstance(expected,list):
        require(isinstance(actual,list) and len(actual)==len(expected),"different row count at "+location)
        for i,(left,right) in enumerate(zip(expected,actual,strict=True)):
            compare_values(left,right,f"{location}/{i}")
    elif type(expected) is float:
        require(type(actual) in (int,float) and math.isclose(expected,actual,rel_tol=5e-12,abs_tol=5e-12),
                "different floating-point result at "+location)
    else:
        require(expected==actual,"different categorical or integral result at "+location)


def compare_profile(original,replayed):
    """When compressor or arithmetic differs, compare every numeric column and label."""
    with gzip.open(original,"rt",newline="",encoding="utf8") as a,gzip.open(replayed,"rt",newline="",encoding="utf8") as b:
        left,right=csv.reader(a),csv.reader(b)
        fields=next(left)
        require(next(right)==fields,"different ring CSV header")
        for i,(u,v) in enumerate(zip(left,right,strict=True)):
            for field,x,y in zip(fields,u,v,strict=True):
                if x==y:
                    continue
                require(field!="region" and math.isclose(float(x),float(y),rel_tol=5e-12,abs_tol=5e-12),
                        f"different ring cell at row {i}, {field}")


def replay(restored,out,*,registration_sha=REGISTRATION,seal_sha=EXECUTION):
    restored,out=Path(restored).resolve(),Path(out).resolve()
    require(not out.exists(),"output already exists")
    require(not out.is_relative_to(restored) and not restored.is_relative_to(out),"output and restored evidence overlap")
    require(digest(restored/"execution.json")==seal_sha==(restored/"execution.sha256").read_text().strip()
            and digest(restored/"wall-registration.json")==registration_sha,"execution or wall-registration seal differs")
    plan=read(restored/"execution.json")
    registration=read(restored/"wall-registration.json")
    panel=read(restored/"inputs/panel.json")
    reference=read(restored/"metrics.json")
    require(reference["analyzer_sha256"]==plan["code"][ENGINE] and
            reference["registration_sha256"]==plan["registration_sha256"]==registration_sha and
            reference["panel_sha256"]==plan["files"]["inputs/panel.json"]==digest(restored/"inputs/panel.json"),
            "metrics, execution or panel input digest differs")
    require(len(reference["cases"])==len(panel["cases"]),"incomplete archived metrics")
    out.parent.mkdir(exist_ok=True,parents=True)
    with TemporaryDirectory(dir=out.parent,prefix=".portable-metrics-") as temp:
        stage=Path(temp)/"replay"
        stage.mkdir()
        actual=[]
        byte_identical=0
        jobs=0
        with analyzer(restored,stage,plan,panel,registration) as (engine,omitted):
            weights,frequency=engine.bootstrap_counts(len(engine.SEEDS))
            for c,expected in zip(panel["cases"],reference["cases"],strict=True):
                cid=c["case_id"]
                case=engine.analyze_case(restored,c,panel["case_record_sha256"][cid],weights,frequency)
                actual.append(case)
                # Reconstruct the result before comparing it to the archived answer.
                checked=json.loads(json.dumps(case))
                for wall,original in expected["geometries"].items():
                    name=f"profiles/{cid}/{wall}.csv.gz"
                    profile=case["geometries"][wall]
                    require(original["profile"]==profile["profile"]==name and digest(restored/name)==original["profile_sha256"]
                            and digest(stage/name)==profile["profile_sha256"],"original or replay ring file differs")
                    jobs+=len(profile["seeds"])
                    if original["profile_sha256"]==profile["profile_sha256"]:
                        byte_identical+=1
                    else:
                        compare_profile(restored/name,stage/name)
                        checked["geometries"][wall]["profile_sha256"]=original["profile_sha256"]
                compare_values(expected,checked)
                print("Replayed",cid,len(case["geometries"]),"material walls",flush=True)
        wall_count=sum(len(c["geometries"]) for c in actual)
        audit=dict(execution_sha256=seal_sha,registration_sha256=registration_sha,
                   source_metrics_sha256=digest(restored/"metrics.json"),analyzer_sha256=plan["code"][ENGINE],
                   case_count=len(actual),case_wall_count=wall_count,seed_jobs=jobs,
                   exact_cases=(reference["cases"]==actual),exact_profile_bytes=byte_identical,
                   mesh_digests_from_receipts=omitted,floating_relative_tolerance=5e-12,floating_absolute_tolerance=5e-12)
        # All numerical and categorical case rows were generated by the sealed
        # analyzer. Preserve the archived fixed protocol metadata verbatim.
        (stage/"metrics.json").write_text(json.dumps(dict(reference,cases=actual),indent=2,allow_nan=False)+"\n",encoding="utf8")
        (stage/"replay.json").write_text(json.dumps(audit,indent=2,allow_nan=False)+"\n",encoding="utf8")
        stage.rename(out)
    return audit


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--restored",type=Path,required=True)
    cli.add_argument("--out",type=Path,required=True)
    args=cli.parse_args()
    print(json.dumps(replay(args.restored,args.out),indent=2))
