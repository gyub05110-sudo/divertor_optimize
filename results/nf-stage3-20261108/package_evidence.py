"""Build, verify or restore the completed Stage 3 evidence using only the Python standard library."""
import argparse
from contextlib import contextmanager, ExitStack
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import re
import shutil
import stat
import tempfile
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


HERE = Path(__file__).resolve().parent
LIMIT = 45 * 1024**2
PIECE = 4 * 1024**2
MANIFEST = "evidence-manifest.json"
SOURCE_REV = "1e280905b85f2a6fdde7e06fca8cf3a1edf447cb"
SOURCE_BASE = f"https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/resolve/{SOURCE_REV}/"
ARMS = ("floor-screened", "full-screened", "full-direct")
SEEDS = (20711, 20729, 20747)
HISTORY = {
    "nf-prospective-inputs-20261107": {"hub-inventory.json","prior-sources.json","prior-frame-index.json",
                                      "downloads.json","dataset-README.md"},
    "unused-equilibrium-20261106": {"source-audit.json"},
}
EXPORTS = {"geometry.csv", "donors.csv.gz", "profile.csv.gz", "boundary.npz"}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def hashed(data):
    return hashlib.sha256(data).hexdigest()


def canonical(obj):
    return hashed(json.dumps(obj,sort_keys=True,separators=(",",":"),allow_nan=False).encode())


def as_json(data):
    def pairs(rows):
        obj = dict(rows)
        require(len(obj) == len(rows), "duplicate JSON field")
        return obj
    return json.loads(data, object_pairs_hook=pairs)


def safe_path(name):
    require(isinstance(name,str) and name and len(name)<4096, "invalid relative path")
    for part in name.split("/"):
        require(part and part not in (".","..") and not part.endswith((".", " ")) and
                not re.search(r'[<>:"\\|?*\x00-\x1f]',part) and
                not re.fullmatch(r"(CON|PRN|AUX|NUL|CONIN\$|CONOUT\$|COM[1-9¹²³]|LPT[1-9¹²³])(\..*)?",
                                 part, flags=re.I), f"unsafe path: {name}")
    return name


def check_paths(names):
    files, directories = set(), {}
    for name in names:
        safe_path(name)
        key = name.casefold()
        require(key not in files and key not in directories, f"path collision: {name}")
        files.add(key)
        parts = name.split("/")
        for i in range(1,len(parts)):
            parent = "/".join(parts[:i])
            folded = parent.casefold()
            require(folded not in files and directories.setdefault(folded, parent) == parent,
                    f"path collision: {name}")


def check_link(path):
    info = path.lstat()
    require(not stat.S_ISLNK(info.st_mode) and not (getattr(info,"st_file_attributes",0) & 0x400),
            f"linked file or directory: {path}")
    return info


def walk(root):
    files = {}
    def visit(folder):
        for path in sorted(folder.iterdir()):
            rel = safe_path(path.relative_to(root).as_posix())
            info = check_link(path)
            if rel == "runtime":
                require(stat.S_ISDIR(info.st_mode), "runtime must be a cache directory")
            elif stat.S_ISDIR(info.st_mode):
                visit(path)
            else:
                require(stat.S_ISREG(info.st_mode), f"non-regular evidence: {rel}")
                files[rel] = path
    visit(root)
    return files


def write_payload(out, catalog, external, *, piece=PIECE, limit=LIMIT):
    """One compressed representative per byte block. Small files have human-readable ZIP member paths."""
    require(0 < piece <= PIECE and 0 < limit <= LIMIT, "invalid archive size settings")
    check_paths(catalog)
    known = {v["sha256"]: k for k,v in external.items()}
    blocks, files, volumes, duplicates = {}, {}, {}, {}
    archive = None
    volume = None
    central = 0

    def finish():
        archive.close()
        path=out/volume
        require(path.stat().st_size < limit, f"volume exceeds configured limit: {volume}")
        volumes[volume] = dict(bytes=path.stat().st_size, sha256=digest(path))
    try:
        for name,path in sorted(catalog.items()):
            size, sha = path.stat().st_size, digest(path)
            row = dict(bytes=size, sha256=sha)
            files[name] = row
            if name in external:
                require(external[name]["sha256"] == sha and external[name]["bytes"] == size,
                        f"external source differs from registration: {name}")
            if sha in known:
                row["external"] = known[sha]
                continue
            if sha in duplicates:
                row["blocks"] = duplicates[sha]
                continue
            parts = []
            checked = hashlib.sha256()
            with path.open("rb") as stream:
                while data := stream.read(piece):
                    checked.update(data)
                    key = hashed(data)
                    parts.append(key)
                    if key in blocks:
                        continue
                    member = "files/"+name if len(data)==size else "chunks/"+key
                    name_bytes = len(member.encode())
                    bound = len(data) + len(data)//500 + 280 + 2*name_bytes
                    require(bound+22 < limit, "piece exceeds ZIP volume limit")
                    if archive is not None and archive.fp.tell()+central+bound+22 >= limit:
                        finish()
                        archive = None
                    if archive is None:
                        volume=f"evidence-{len(volumes)+1:04d}.zip"
                        archive=ZipFile(out/volume, "x")
                        central = 0
                    info=ZipInfo(member, (1980,1,1,0,0,0))
                    archive.writestr(info, data, compress_type=ZIP_DEFLATED, compresslevel=3)
                    central += 46+name_bytes
                    blocks[key] = dict(volume=volume, member=member, bytes=len(data))
            require(checked.hexdigest() == sha, f"source changed while reading: {name}")
            row["blocks"] = duplicates[sha] = parts
        if archive is not None:
            finish()
            archive = None
    finally:
        if archive is not None:
            archive.close()
    shutil.copyfile(__file__,out/"package_evidence.py")
    manifest = dict(schema="nf-stage3-portable-evidence.v1", omitted_directories=["runtime"],
                    verifier_sha256=digest(out/"package_evidence.py"),
                    external=external, volumes=volumes, blocks=blocks, files=files)
    (out/MANIFEST).write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    (out/(MANIFEST+".sha256")).write_text(digest(out/MANIFEST)+"\n", encoding="ascii")


class Reader:
    def __init__(self, manifest, archives):
        self.files = manifest["files"]
        self.blocks = manifest["blocks"]
        self.external = manifest["external"]
        self.archives = archives

    def chunks(self,name):
        for key in self.files[name]["blocks"]:
            item = self.blocks[key]
            yield self.archives[item["volume"]].read(item["member"])

    def content(self,name):
        return b"".join(self.chunks(name))

    def read(self,name):
        return as_json(self.content(name))

    def check(self,name,sha=None):
        safe_path(name)
        require(name in self.files and (sha is None or self.files[name]["sha256"]==sha), f"missing/changed evidence: {name}")
        return self.files[name]["sha256"]


@contextmanager
def open_payload(out, expected=None):
    for name in (MANIFEST, MANIFEST+".sha256","package_evidence.py"):
        check_link(out/name)
    sha=digest(out/MANIFEST)
    require((out/(MANIFEST+".sha256")).read_text().strip()==sha and (expected is None or sha==expected),
            "manifest digest differs")
    m=as_json((out/MANIFEST).read_bytes())
    require(m["schema"]=="nf-stage3-portable-evidence.v1" and m["omitted_directories"]==["runtime"], "unexpected archive schema")
    require(digest(out/"package_evidence.py")==m["verifier_sha256"],"portable verifier hash differs")
    check_paths(m["files"])
    require(set(m["volumes"])=={p.name for p in out.glob("*.zip")}, "missing or additional ZIP volume")
    by_volume={n:{} for n in m["volumes"]}
    for key,record in m["blocks"].items():
        require(re.fullmatch("[0-9a-f]{64}",key) and 0 < record["bytes"] <= PIECE, "invalid block")
        name=safe_path(record["member"])
        require(record["volume"] in by_volume and name not in by_volume[record["volume"]], "duplicate archive member")
        by_volume[record["volume"]][name]=(key,record["bytes"])
    with ExitStack() as stack:
        archives={}
        for name,row in m["volumes"].items():
            require(re.fullmatch(r"evidence-\d{4,}\.zip",name), "invalid volume path")
            p=out/name
            check_link(p)
            require(0 < p.stat().st_size == row["bytes"] < LIMIT and digest(p)==row["sha256"], f"volume digest differs: {name}")
            archive=archives[name]=stack.enter_context(ZipFile(p))
            members=by_volume[name]
            require(len(archive.infolist())==len(members) and set(archive.namelist())==set(members), f"volume members differ: {name}")
            for info in archive.infolist():
                key, size = members[info.filename]
                require(info.file_size==size and info.compress_type==ZIP_DEFLATED and not info.is_dir() and
                        not stat.S_ISLNK(info.external_attr>>16), "invalid ZIP member metadata")
                require(hashed(archive.read(info))==key, f"block digest differs: {key}")
        reader=Reader(m,archives)
        checked, used=set(),set()
        for name,row in reader.files.items():
            sha=row["sha256"]
            require(re.fullmatch("[0-9a-f]{64}",sha) and type(row["bytes"]) is int and row["bytes"]>=0,
                    f"invalid file metadata: {name}")
            if "external" in row:
                ext=reader.external[row["external"]]
                require("blocks" not in row and (ext["sha256"],ext["bytes"])==(sha,row["bytes"]), "invalid external reference")
                continue
            parts=tuple(row["blocks"])
            require(all(k in reader.blocks for k in parts) and sum(reader.blocks[k]["bytes"] for k in parts)==row["bytes"],
                    "missing block or inconsistent file size")
            used.update(parts)
            if (sha,parts) not in checked:
                content=hashlib.sha256()
                for data in reader.chunks(name):
                    content.update(data)
                require(content.hexdigest()==sha,f"file digest differs: {name}")
                checked.add((sha,parts))
        require(used==set(reader.blocks), "unreferenced archive block")
        yield reader


def external_sources(registration):
    require(registration["revision"] == SOURCE_REV and len(registration["cases"])==12, "expected twelve pinned sources")
    records={}
    for row in registration["cases"]:
        ident=row["id"]
        require(re.fullmatch("[a-z0-9]{10}",ident), "invalid registered source ID")
        base="d3d_shot_"+ident+".parquet"
        path="data/diii_d_train/"+base
        require(row["dataset_path"]==path and row["url"]==SOURCE_BASE+path and
                re.fullmatch("[0-9a-f]{64}",row["sha256"]) and row["bytes"]>0,"invalid pinned parquet provenance")
        records["inputs/raw/"+base]={k:row[k] for k in ("id","dataset_path","url","bytes","sha256")}
    require(len(records)==len({v["sha256"] for v in records.values()})==12,"duplicate registered source")
    return records


def history_record(name):
    parts=PureWindowsPath(name).parts
    require(len(parts)>=2 and parts[-1] in HISTORY.get(parts[-2],set()), "unexpected historical metadata reference")
    return parts[-2]+"/"+parts[-1]


def validate_study(r):
    """Follow the study's published digest chains entirely through archive-relative names."""
    expected=r.read("expected.json")
    plan_sha=r.check("expected.json")
    require(r.content("expected.sha256").decode().strip()==plan_sha, "plan seal differs")
    reg_sha=r.check("inputs/registration.json",expected["registration_sha256"])
    reg=r.read("inputs/registration.json")
    require(r.external==external_sources(reg),"external source manifest differs from registration")
    for name,row in r.external.items():
        r.check(name,row["sha256"])
        require(r.files[name].get("external")==name, "registered parquet must remain external")
    r.check("inputs/selection.json",reg["selection_sha256"])
    r.check("inputs/selection-seal.json",reg["selection_seal_sha256"])
    r.check("inputs/protocol.ko.md",reg["protocol_sha256"])
    selection=r.read("inputs/selection.json")
    require(len(selection["selected"])==12 and all(all(x[k]==c[k] for k in x)
            for x,c in zip(selection["selected"],reg["cases"])), "registered population differs from selection")
    for p,v in selection["historical_records"].items():
        r.check("provenance/"+history_record(p),v)
    for p,v in reg["files_sha256"].items():
        r.check("inputs/"+safe_path(p),v)
    for p,v in expected["files"].items():
        r.check("source-code/"+safe_path(p),v)
    cases={c["id"]:c for c in reg["cases"]}
    for ident,c in cases.items():
        name=f"inputs/receipts/{ident}.json"
        r.check(name,c["receipt_sha256"])
        rec=r.read(name)
        require(all(rec[k]==c[k] for k in ("id","rank","url","bytes","sha256")),"different download receipt")
        for ending in ("started","finished"):
            r.check(f"inputs/requests/{ident}/{rec['success_attempt']}.{ending}.json")

    report=r.read("stage3_report.json")
    report_sha=r.check("stage3_report.json")
    require(report["complete"] is True and report["expected_sha256"]==plan_sha and report["registration_sha256"]==reg_sha,
            "complete report with matching plan and registration required")
    r.check("support/verification/stage3_report.py",report["report_driver_sha256"])
    for p,v in report["companion_sha256"].items():
        r.check(safe_path(p),v)
    require(set(report["companion_sha256"])=={"stage3_runs.csv","stage3_pairs.csv"}, "incomplete tabular report")
    def task_key(row):
        return row["case_id"],row["seed"],row["arm"]
    tasks={task_key(t):t for t in expected["tasks"]}
    runs={task_key(t):t for t in report["runs"]}
    keys={(ident,seed,arm) for ident in cases for seed in SEEDS for arm in ARMS}
    require(len(expected["tasks"])==len(report["runs"])==108 and set(tasks)==set(runs)==keys,"incomplete scheduled runs")
    canonical_plan_sha=canonical(expected)
    artifact_names={"result_sha256":"result.json","events_sha256":"events.jsonl","stdout_sha256":"stdout.log"}
    for (ident,seed,arm),task in tasks.items():
        folder=f"runs/{ident}/{seed}/{arm}"
        row=runs[ident,seed,arm]
        require(task["relative_output"]==row["relative_output"]==folder and task["sha256"]==cases[ident]["sha256"],
                "incorrect run-to-case binding")
        key=folder+"/receipt.json"
        r.check(key,row["receipt_sha256"])
        receipt=r.read(key)
        require(task_key(receipt)==task_key(task) and receipt["registration_sha256"]==reg_sha and
                receipt["job_sha256"]==canonical(task) and receipt["plan_sha256"]==canonical_plan_sha and
                all(row[k]==receipt[k] for k in ("events_sha256","result_sha256")) and
                receipt["status"]==row["receipt_status"] and receipt["status"] in
                ("registration_excluded","completed","execution_failure"),"unexpected run receipt")
        for field,name in artifact_names.items():
            p=folder+"/"+name
            if receipt[field] is not None:
                r.check(p,receipt[field])
            else:
                require(p not in r.files,"unreceipted endpoint file: "+p)
        if receipt["status"]=="registration_excluded":
            require(cases[ident]["status"]!="ready","ready case missing run")
        else:
            attempt=safe_path(receipt["source_attempt"])
            require(re.fullmatch(r"attempts/"+ident+"/"+str(seed)+"/"+arm+r"/\d+",attempt), "wrong attempt binding")
            r.check(attempt+"/launch.json",receipt["launch_sha256"])
            r.check(attempt+"/command.json")
            r.check(attempt+"/stdout.log",receipt["stdout_sha256"])
            if receipt["completion_sha256"] is not None:
                r.check(attempt+"/completion.json",receipt["completion_sha256"])
        if receipt["status"]=="completed":
            require(all(receipt[k] for k in artifact_names),"completed run lacks its raw evidence")
            check_exports(r,folder,r.read(folder+"/result.json"))
        observation=f"observations/{ident}/{seed}/{arm}/"
        require(any(p.startswith(observation) and p.endswith(".started.json") for p in r.files),
                "missing parent observation: "+folder)

    summary=r.read("replays/summary.json")
    require(summary["complete"] is True and summary["search_report_sha256"]==report_sha and
            len(summary["cases"])==12 and {x["case_id"] for x in summary["cases"]}==set(cases),"incomplete replay summary")
    for row in summary["cases"]:
        ident=row["case_id"]
        folder="replays/"+ident
        r.check(folder+"/report.json",row["report_sha256"])
        replay=r.read(folder+"/report.json")
        require(replay["complete"] is True and replay["case_id"]==ident,"incomplete replay: "+ident)
        r.check(folder+"/protocol.json",replay["protocol_sha256"])
        proto=r.read(folder+"/protocol.json")
        require(proto["search_report_sha256"]==report_sha and proto["case_id"]==ident,"replay input mismatch")
        r.check("support/verification/stage3_replay.py",proto["replay_driver_sha256"])
        declared=replay["measurement_receipts_sha256"]
        needed=set()
        if proto["common"]:
            require(set(replay["results"])==set(proto["levels"]), "missing replay numerical level")
            for level,grid in replay["results"].items():
                r.check(folder+"/"+safe_path(level)+"/setup.json")
                if grid["status"]=="measured":
                    needed.update(level+"/"+role+"/measurement.json" for role in ("baseline",*proto["candidates"]))
                else:
                    require(grid["status"]=="calculation_error","unresolved replay setup")
        require(needed<=set(declared), "replay measurement receipts incomplete")
        for p,v in declared.items():
            relative=folder+"/"+safe_path(p)
            r.check(relative,v)
            measurement=r.read(relative)
            if "sha256" in measurement:
                check_export(r,relative.rsplit("/",1)[0],measurement)
    commission="commissioning/direct_commissioning_summary.json"
    r.check(commission,expected["files"]["results/nf-stage3-20261108/commissioning.json"])
    for p,v in r.read(commission)["artifacts_sha256"].items():
        r.check("commissioning/direct/"+safe_path(p),v)
    check_exports(r,"commissioning/direct",r.read("commissioning/direct/result.json"))
    for p in ("support/results/nf-stage3-20261108/run-study.ps1","workflow.out.log","workflow.err.log",
              "preflight-tests.log","test-suite.log","progress.json"):
        r.check(p)
    require(r.read("progress.json")["completed"]==108,"incomplete campaign progress")


def check_export(r,folder,row):
    require(row["directory"]==folder.split("/")[-1] and set(row["sha256"])==EXPORTS, "incomplete native export")
    for p,v in row["sha256"].items():
        r.check(folder+"/"+p,v)


def check_exports(r,folder,result):
    for name in ("baseline","selected"):
        row=result.get(name) or {}
        if "directory" in row or "sha256" in row:
            check_export(r,folder+"/"+name,row)


def disjoint(a,b):
    a,b=a.resolve(),b.resolve()
    require(not a.is_relative_to(b) and not b.is_relative_to(a), "source and destination paths overlap")


def collect(study,repo):
    catalog=walk(study)
    selection=as_json((study/"inputs/selection.json").read_bytes())
    def include(name,source,checksum=None):
        require(name not in catalog and source.is_file(),"missing or conflicting supplementary file: "+name)
        check_link(source)
        if checksum is not None:
            require(digest(source)==checksum,"historical checksum differs: "+name)
        catalog[name]=source
    for src,sha in selection["historical_records"].items():
        rel=history_record(src)
        include("provenance/"+rel,study.parent/rel,sha)
    for name in ("verification/stage3_report.py","verification/stage3_replay.py",
                 "results/nf-stage3-20261108/run-study.ps1"):
        include("support/"+name,repo/name)
    return catalog


def build(study,out,repo=None):
    study=study.resolve()
    repo=repo.resolve() if repo else HERE.parents[1]
    disjoint(study,out)
    require(not out.exists(),"bundle destination exists")
    for name in ("stage3_report.json","replays/summary.json"):
        require((study/name).is_file() and as_json((study/name).read_bytes()).get("complete") is True,
                "complete timed report and replays required: "+name)
    catalog=collect(study,repo)
    stamps=lambda items: {n:(p.stat().st_size,p.stat().st_mtime_ns) for n,p in items.items()}
    before=stamps(catalog)
    ext=external_sources(as_json((study/"inputs/registration.json").read_bytes()))
    require(set(ext)<=set(catalog),"registered raw files missing during build")
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent,prefix=".stage3-build-") as tmp:
        stage=Path(tmp)/"bundle"
        stage.mkdir()
        write_payload(stage,catalog,ext)
        with open_payload(stage) as r:
            validate_study(r)
        require(stamps(collect(study,repo))==before,"study or workflow traces changed during packaging")
        publish_dir(stage,out)
    return digest(out/MANIFEST)


def publish_dir(stage,out):
    # Reserving the destination is atomic even if another process is publishing.
    out.mkdir(exist_ok=False)
    for item in sorted(stage.iterdir(),key=lambda p:p.name==MANIFEST):
        item.rename(out/item.name)


def restore(bundle,out,raw_dir=None,expected=None):
    disjoint(bundle,out)
    if raw_dir:
        disjoint(raw_dir,out)
    require(not out.exists(),"restore destination exists")
    out.parent.mkdir(parents=True,exist_ok=True)
    with open_payload(bundle,expected) as r:
        validate_study(r)
        with tempfile.TemporaryDirectory(dir=out.parent,prefix=".stage3-restore-") as tmp:
            stage=Path(tmp)/"tree"
            stage.mkdir()
            for name,row in r.files.items():
                dest=stage/name
                if "external" in row and raw_dir is None:
                    continue
                dest.parent.mkdir(parents=True,exist_ok=True)
                with dest.open("xb") as stream:
                    if "external" in row:
                        src=raw_dir/Path(row["external"]).name
                        check_link(src)
                        with src.open("rb") as source:
                            shutil.copyfileobj(source,stream)
                    else:
                        for data in r.chunks(name):
                            stream.write(data)
                require(dest.stat().st_size==row["bytes"] and digest(dest)==row["sha256"], f"restored file differs: {name}")
            publish_dir(stage,out)
    return len(r.external)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    commands=p.add_subparsers(dest="action",required=True)
    b=commands.add_parser("build",help="verify and pack a finished study into a new directory")
    b.add_argument("--study",type=Path,required=True)
    b.add_argument("--bundle",type=Path,required=True)
    b.add_argument("--repo",type=Path,help="repository containing the matching analysis and workflow drivers")
    for action in ("verify","restore"):
        v=commands.add_parser(action)
        v.add_argument("--bundle",type=Path,required=True)
        v.add_argument("--manifest-sha256",help="independently published SHA-256 of the manifest")
        if action=="restore":
            v.add_argument("--out",type=Path,required=True)
            v.add_argument("--raw-dir",type=Path,help="optional directory of the 12 already-downloaded pinned parquet files")
    args=p.parse_args()
    if args.action=="build":
        print("Manifest SHA-256",build(args.study,args.bundle,args.repo))
    elif args.action=="verify":
        with open_payload(args.bundle,args.manifest_sha256) as r:
            validate_study(r)
            print(f"Verified {len(r.files)} logical files, {len(r.external)} external inputs and {len(r.archives)} volumes.")
    else:
        n=restore(args.bundle,args.out,args.raw_dir,args.manifest_sha256)
        print("Restored",args.out,"with" if args.raw_dir else "with external references for",n,"registered parquet inputs.")


if __name__=="__main__":
    main()
