"""Build, verify and restore compact or full Stage 4 evidence with Python 3.11+."""
import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
from tempfile import TemporaryDirectory


HERE = Path(__file__).resolve().parent
MANIFEST = "stage4-manifest.json"
TABLES = {"registered-runs.csv", "policy-predictions.csv", "product-pairs.csv"}
SEEDS = (30713, 30731, 30757)
PRODUCT = ("hotspot", "full")
ARMS = (*PRODUCT, "floor-only", "floor-entry", "floor-shelf", "floor-connector")
CODEC_SHA = "1cac0d886d7aa125253a085e96f6c8639bcefd3940e0e12a00b09267b7e1f742"


def _codec():
    """Reuse the verified, deterministic split-ZIP format from the Stage 3 archive."""
    path = HERE / "payload/package_evidence.py"
    if not path.is_file():
        path = HERE.parent / "nf-stage3-20261108/package_evidence.py"
    if hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest() != CODEC_SHA:
        raise ValueError("split-ZIP codec differs from the published version")
    spec = importlib.util.spec_from_file_location("stage4_zip_codec", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


codec = _codec()
require = codec.require


def stored(name, mode):
    # The compact archive retains complete candidate journals, receipts and input
    # fingerprints. Large physical measurements and raw equilibria live on D:.
    return mode == "full" or (not name.lower().endswith((".parquet", ".npz", ".gz")) and
                              Path(name).name != "geometry.csv")


class Evidence:
    def __init__(self, files, content):
        self.files = files
        self.content = content

    def check(self, name, sha=None):
        codec.safe_path(name)
        require(name in self.files and (sha is None or self.files[name]["sha256"] == sha),
                "missing or changed evidence: "+name)
        return self.files[name]["sha256"]

    def read(self, name):
        self.check(name)
        return codec.as_json(self.content(name))


def validate(e):
    plan = e.read("expected.json")
    plan_sha = e.check("expected.json")
    require(plan["schema"] == "nf-stage4-hotspot.v1" and
            e.content("expected.sha256").decode().strip() == plan_sha, "plan seal differs")
    reg = e.read("inputs/registration.json")
    reg_sha = e.check("inputs/registration.json", plan["registration_sha256"])
    for prefix, names in (("source-code/", plan["files"]), ("inputs/", reg["files_sha256"])):
        for name, sha in names.items():
            e.check(prefix+codec.safe_path(name), sha)
    for name, field in (("protocol.ko.md","protocol_sha256"), ("selection.json","selection_sha256"),
                        ("selection-seal.json","selection_seal_sha256")):
        e.check("inputs/"+name, reg[field])
    selection = e.read("inputs/selection.json")
    seal = e.read("inputs/selection-seal.json")
    require(seal["selection_sha256"] == reg["selection_sha256"] and
            seal["protocol_sha256"] == reg["protocol_sha256"] == selection["protocol_sha256"],
            "input selection seal differs")
    cases = {r["id"]:r for r in reg["cases"]}
    require(len(cases) == len(reg["cases"]) == len(selection["selected"]) == 24 and
            all(all(c[k] == v for k,v in s.items()) for c,s in zip(reg["cases"],selection["selected"])),
            "expected all 24 registered cases in frozen selection order")
    for cid, case in cases.items():
        require(re.fullmatch("[a-z0-9]{10}", cid), "invalid registered case ID")
        if case.get("receipt_sha256"):
            e.check(f"inputs/receipts/{cid}.json", case["receipt_sha256"])
            raw = "inputs/raw/"+Path(case["dataset_path"]).name
            e.check(raw, case["sha256"])
            require(e.files[raw]["bytes"] == case["bytes"], "registered source size differs")
        else:
            require(case["status"] == "acquisition_failed", "missing registered source receipt")
        for i, attempt in enumerate(case.get("attempts", []), 1):
            for phase in ("started", "finished"):
                if phase+"_sha256" in attempt:
                    e.check(f"inputs/requests/{cid}/{i}.{phase}.json", attempt[phase+"_sha256"])
    report = e.read("stage4_report.json")
    report_sha = e.check("stage4_report.json")
    require(report["complete"] is True and report["plan_sha256"] == plan_sha and
            report["registration_sha256"] == reg_sha, "complete timed report with matching seal required")
    e.check("source-code/verification/stage4_report.py", report["report_driver_sha256"])
    require(set(report["companion_sha256"]) == TABLES, "missing published tables")
    for name, sha in report["companion_sha256"].items():
        e.check(name, sha)
    validate_jobs(e, plan, cases, report, reg_sha)
    validate_replays(e, cases, report_sha)


def validate_jobs(e, plan, cases, report, reg_sha):
    key = lambda r: (r["case_id"], r["seed"], r["arm"])
    jobs = {key(t):t for t in plan["tasks"]}
    rows = {key(r):r for r in report["runs"]}
    schedule = {(cid, seed, arm) for cid in cases for seed in SEEDS for arm in (ARMS if seed == SEEDS[0] else PRODUCT)}
    require(len(plan["tasks"]) == len(report["runs"]) == 240 and set(jobs) == set(rows) == schedule,
            "incomplete Stage 4 schedule or result inventory")
    plan_id = codec.canonical(plan)
    for k, job in jobs.items():
        cid, seed, arm = k
        folder = f"runs/{cid}/{seed}/{arm}"
        r = rows[k]
        require(job["relative_output"] == r["relative_output"] == folder and job["sha256"] == cases[cid]["sha256"],
                "job binding differs")
        e.check(folder+"/receipt.json", r["receipt_sha256"])
        receipt = e.read(folder+"/receipt.json")
        require(key(receipt) == k and receipt["plan_sha256"] == plan_id and receipt["job_sha256"] == codec.canonical(job)
                and receipt["registration_sha256"] == reg_sha and receipt["status"] == r["receipt_status"]
                and receipt["status"] in ("completed","registration_excluded","execution_failure")
                and (receipt["status"] == "registration_excluded") == (cases[cid]["status"] != "ready"),
                "run receipt binding differs: "+folder)
        for label, name in (("result","result.json"), ("events","events.jsonl"), ("stdout","stdout.log")):
            field = label+"_sha256"
            if label != "stdout":
                require(r[field] == receipt[field], "report disagrees with run receipt")
            if receipt[field] is None:
                require(folder+"/"+name not in e.files and receipt["status"] != "completed", "incomplete run evidence")
            else:
                e.check(folder+"/"+name, receipt[field])
        if receipt["status"] != "registration_excluded":
            attempt = codec.safe_path(receipt["source_attempt"])
            require(re.fullmatch(f"attempts/{cid}/{seed}/{arm}/[0-9]+", attempt), "attempt binding differs")
            e.check(attempt+"/launch.json", receipt["launch_sha256"])
            e.check(attempt+"/command.json")
            if receipt.get("completion_sha256"):
                e.check(attempt+"/completion.json", receipt["completion_sha256"])
        if receipt["status"] == "completed":
            codec.check_exports(e, folder, e.read(folder+"/result.json"))


def validate_replays(e, cases, report_sha):
    summary = e.read("replays/summary.json")
    require(summary["complete"] is True and summary["search_report_sha256"] == report_sha and
            len(summary["cases"]) == 24 and {r["case_id"] for r in summary["cases"]} == set(cases),
            "complete replay decision for every case required")
    e.check("stage4_replay.json", e.check("replays/summary.json"))
    for row in summary["cases"]:
        cid = row["case_id"]
        folder = f"replays/{cid}/"
        e.check(folder+"report.json", row["report_sha256"])
        replay = e.read(folder+"report.json")
        e.check(folder+"protocol.json", replay["protocol_sha256"])
        proto = e.read(folder+"protocol.json")
        require(replay["complete"] is True and replay["case_id"] == proto["case_id"] == cid and
                replay["search_report_sha256"] == proto["search_report_sha256"] == report_sha, "replay binding differs")
        e.check("source-code/verification/stage4_replay.py", proto["replay_driver_sha256"])
        e.check("source-code/verification/stage3_replay.py", proto["fixed_wall_driver_sha256"])
        expected_grids = {"cubic1280","cubic640","quintic1280"} if proto["common"] else set()
        require(set(replay["grid_checkpoints_sha256"]) == set(replay["results"]) == expected_grids, "replay grids differ")
        for level, sha in replay["grid_checkpoints_sha256"].items():
            e.check(folder+level+".json", sha)
            grid = e.read(folder+level+".json")
            require(grid["protocol_sha256"] == replay["protocol_sha256"] and grid["grid"] == replay["results"][level],
                    "replay checkpoint differs")
            records = grid["files_sha256"]
            require("setup.json" in records, "missing replay setup")
            if grid["grid"]["status"] == "measured":
                require(all(role+"/measurement.json" in records for role in ("baseline", *proto["candidates"])),
                        "missing replay measurement")
            for relative, checksum in records.items():
                name = folder+level+"/"+codec.safe_path(relative)
                e.check(name, checksum)
                if name.endswith("/measurement.json"):
                    measurement = e.read(name)
                    if measurement.get("sha256"):
                        codec.check_export(e, name.rsplit("/",1)[0], measurement)


def build(study, out, mode):
    study, out = Path(study).resolve(), Path(out).resolve()
    codec.disjoint(study, out)
    require(mode in ("compact","full") and not out.exists(), "unknown mode or bundle destination exists")
    for name in ("stage4_report.json","replays/summary.json"):
        require((study/name).is_file() and codec.as_json((study/name).read_bytes()).get("complete") is True,
                "complete timed report and replay summary required")
    catalog = codec.walk(study)
    def stamps(paths):
        return {n:(p.stat().st_size,p.stat().st_mtime_ns) for n,p in paths.items()}
    before = stamps(catalog)
    files = {n:dict(bytes=p.stat().st_size,sha256=codec.digest(p)) for n,p in catalog.items()}
    codec.check_paths(files)
    validate(Evidence(files, lambda name: catalog[name].read_bytes()))
    out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=out.parent, prefix=".stage4-pack-") as scratch:
        stage = Path(scratch)/"bundle"
        payload = stage/"payload"
        payload.mkdir(parents=True)
        codec.write_payload(payload, {n:p for n,p in catalog.items() if stored(n,mode)}, {})
        shutil.copyfile(__file__, stage/"package_evidence.py")
        manifest = dict(schema="nf-stage4-evidence.v1", mode=mode, source_study=str(study),
                        omitted_directories=["runtime"], files=files, inventory_sha256=codec.canonical(files),
                        payload_sha256=codec.digest(payload/codec.MANIFEST),
                        packager_sha256=codec.digest(stage/"package_evidence.py"))
        (stage/MANIFEST).write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
        sha = codec.digest(stage/MANIFEST)
        (stage/(MANIFEST+".sha256")).write_text(sha+"\n",encoding="ascii")
        with verified(stage, sha):
            pass
        require(stamps(codec.walk(study)) == before, "study changed during packaging")
        require(not out.exists(), "bundle destination exists")
        stage.rename(out)
    return sha


@contextmanager
def verified(bundle, expected=None):
    bundle = Path(bundle).resolve()
    for name in (MANIFEST, MANIFEST+".sha256", "package_evidence.py", "payload"):
        codec.check_link(bundle/name)
    sha = codec.digest(bundle/MANIFEST)
    require((bundle/(MANIFEST+".sha256")).read_text().strip() == sha and (expected is None or expected == sha),
            "Stage 4 manifest digest differs")
    manifest = codec.as_json((bundle/MANIFEST).read_bytes())
    mode = manifest["mode"]
    files = manifest["files"]
    require(manifest["schema"] == "nf-stage4-evidence.v1" and mode in ("compact","full") and
            manifest["omitted_directories"] == ["runtime"] and manifest["inventory_sha256"] == codec.canonical(files)
            and manifest["packager_sha256"] == codec.digest(bundle/"package_evidence.py"), "Stage 4 manifest differs")
    codec.check_paths(files)
    for row in files.values():
        require(type(row["bytes"]) is int and row["bytes"] >= 0 and re.fullmatch("[0-9a-f]{64}",row["sha256"]),
                "invalid file inventory")
    with codec.open_payload(bundle/"payload",manifest["payload_sha256"]) as reader:
        require(set(reader.files) == {n for n in files if stored(n,mode)} and not reader.external,
                "stored inventory differs from declared mode")
        for name, row in reader.files.items():
            require(all(row[k] == files[name][k] for k in ("bytes","sha256")), "stored evidence differs: "+name)
        validate(Evidence(files,reader.content))
        yield manifest, reader


def restore(bundle, out, expected=None, study=None):
    bundle, out = Path(bundle).resolve(), Path(out).resolve()
    codec.disjoint(bundle,out)
    if study is not None:
        study = Path(study).resolve()
        codec.disjoint(study,out)
    require(not out.exists(), "restore destination exists")
    out.parent.mkdir(parents=True, exist_ok=True)
    with verified(bundle,expected) as (manifest,reader), TemporaryDirectory(dir=out.parent,prefix=".stage4-restore-") as tmp:
        stage = Path(tmp)/"study"
        stage.mkdir()
        for name,row in manifest["files"].items():
            if name not in reader.files and study is None:
                continue
            dest = stage/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            with dest.open("xb") as output:
                if name in reader.files:
                    for block in reader.chunks(name):
                        output.write(block)
                else:
                    source = study/name
                    require(source.resolve().is_relative_to(study), "external file outside original study")
                    codec.check_link(source)
                    with source.open("rb") as stream:
                        shutil.copyfileobj(stream,output)
            require(dest.stat().st_size == row["bytes"] and codec.digest(dest) == row["sha256"],
                    "restored file differs: "+name)
        require(not out.exists(), "restore destination exists")
        stage.rename(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    build_cli = commands.add_parser("build",help="pack a finished study")
    build_cli.add_argument("--study",type=Path,required=True)
    build_cli.add_argument("--bundle",type=Path,required=True)
    build_cli.add_argument("--mode",choices=("compact","full"),required=True)
    for command in ("verify","restore"):
        p = commands.add_parser(command)
        p.add_argument("--bundle",type=Path,required=True)
        p.add_argument("--manifest-sha256",help="independently retained SHA-256 of the Stage 4 manifest")
        if command=="restore":
            p.add_argument("--out",type=Path,required=True)
            p.add_argument("--study",type=Path,help="supply omitted compact-mode files from the original study")
    a = parser.parse_args()
    if a.command=="build":
        print("Stage 4 manifest SHA-256",build(a.study,a.bundle,a.mode))
    elif a.command=="restore":
        restore(a.bundle,a.out,a.manifest_sha256,a.study)
        print("Restored",a.out)
    else:
        with verified(a.bundle,a.manifest_sha256) as (manifest,reader):
            print(f"Verified {len(reader.files)} stored files of {len(manifest['files'])} inventoried "
                  f"in {manifest['mode']} mode. Inventory SHA-256 {manifest['inventory_sha256']}")


if __name__ == "__main__":
    main()
