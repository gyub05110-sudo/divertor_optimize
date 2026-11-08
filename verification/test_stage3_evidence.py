import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import pytest


ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT/"results/nf-stage3-20261108/package_evidence.py"


@pytest.fixture
def pack():
    assert PACK.is_file(), "portable Stage 3 evidence packager required"
    spec = importlib.util.spec_from_file_location("stage3_evidence_pack", PACK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_chunked_dedup_and_roundtrip_without_absolute_source(pack, tmp_path):
    src = tmp_path/"archive-inputs"
    src.mkdir()
    source = src/"journal.jsonl"
    source.write_bytes(bytes(range(251))*35)
    empty = src/"stdout.log"
    empty.write_bytes(b"")
    catalog = {"runs/job/events.jsonl": source, "copies/journal.jsonl": source, "runs/job/stdout.log": empty}
    bundle = tmp_path/"bundle"
    bundle.mkdir()
    pack.write_payload(bundle, catalog, {}, piece=1024, limit=2500)
    source.unlink()
    with pack.open_payload(bundle) as reader:
        assert reader.content("runs/job/events.jsonl") == bytes(range(251))*35
        assert reader.files["copies/journal.jsonl"] == reader.files["runs/job/events.jsonl"]
        assert reader.content("runs/job/stdout.log") == b""
        assert len(list(bundle.glob("*.zip"))) > 1
        assert all(p.stat().st_size < 2500 for p in bundle.glob("*.zip"))
        assert len(reader.blocks) == 9


def test_rejects_traversal_windows_aliases_and_corrupt_archive(pack, tmp_path):
    for path in ("../x", "a/../x", "C:/x", "//server/share", "a\\b", "file:stream", "NUL", "nice."):
        with pytest.raises(ValueError, match="path"):
            pack.safe_path(path)
    f=tmp_path/"e"; f.write_bytes(b"evidence")
    out=tmp_path/"bundle"; out.mkdir()
    with pytest.raises(ValueError, match="collision"):
        pack.write_payload(out, {"A/a":f,"a/A":f}, {})
    pack.write_payload(out, {"valid":f}, {})
    volume=next(out.glob("*.zip"))
    volume.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="volume"):
        with pack.open_payload(out):
            pass


def test_refuses_incomplete_study_before_archiving(pack, tmp_path):
    study=tmp_path/"study"; study.mkdir()
    (study/"stage3_report.json").write_text('{"complete":false}')
    bundle=tmp_path/"evidence"
    with pytest.raises(ValueError, match="complete"):
        pack.build(study, bundle)
    assert not bundle.exists()


def synthetic_study(pack,study,repo):
    def save(root,name,value):
        path=root/name
        path.parent.mkdir(parents=True,exist_ok=True)
        data=value if isinstance(value,bytes) else json.dumps(value).encode()
        path.write_bytes(data)
        return hashlib.sha256(data).hexdigest()
    def put(name,value):
        return save(study,name,value)
    def text(name):
        return put(name,b"tiny synthetic log\n")
    support="results/nf-stage3-20261108/"
    save(repo,support+"run-study.ps1",b"# synthetic workflow")
    report_driver=save(repo,"verification/stage3_report.py",b"# synthetic report")
    replay_driver=save(repo,"verification/stage3_replay.py",b"# synthetic replay")
    commission=dict(artifacts_sha256={p:put("commissioning/direct/"+p,v) for p,v in
                                      (("events.jsonl",b"{}\n"),("result.json",{"status":"fixture"}))})
    cs=put("commissioning/direct_commissioning_summary.json",commission)
    put("source-code/"+support+"commissioning.json",commission)
    for p in ("workflow.out.log","workflow.err.log","test-suite.log","preflight-tests.log"):
        text(p)
    raw_directory=study/"inputs"/"raw"
    records=[]
    for i in range(12):
        ident=f"{i:010x}"
        filename=f"d3d_shot_{ident}.parquet"
        rel="data/diii_d_train/"+filename
        raw=("fake-"+ident).encode()
        sha=put("inputs/raw/"+filename,raw)
        rec=dict(id=ident,rank=i+1,sha256=sha,bytes=len(raw),dataset_path=rel,url=pack.SOURCE_BASE+rel,
                 status="ready" if not i else "missing_time",path="Q:/old/host/"+filename,time_ms=1400.)
        receipt={**rec,"success_attempt":1}
        rec["receipt_sha256"]=put(f"inputs/receipts/{ident}.json",receipt)
        for ending in ("started","finished"):
            put(f"inputs/requests/{ident}/1.{ending}.json",{})
        records.append(rec)
    historical=save(study.parent/"nf-prospective-inputs-20261107","prior-sources.json",{"samples":[]})
    sel=put("inputs/selection.json",dict(selected=[{"id":v["id"]} for v in records],
            historical_records={r"Q:\irrelevant\nf-prospective-inputs-20261107\prior-sources.json":historical}))
    raw_support=text("inputs/prior-frame-index.json")
    reg=put("inputs/registration.json",dict(revision=pack.SOURCE_REV,cases=records,selection_sha256=sel,
            protocol_sha256=text("inputs/protocol.ko.md"),selection_seal_sha256=put("inputs/selection-seal.json",{}),
            files_sha256={"prior-frame-index.json":raw_support}))
    tasks=[dict(case_id=c["id"],seed=s,arm=a,sha256=c["sha256"],relative_output=f"runs/{c['id']}/{s}/{a}")
           for c in records for s in pack.SEEDS for a in pack.ARMS]
    plan=dict(registration_sha256=reg,tasks=tasks,files={support+"commissioning.json":cs})
    plan_sha=put("expected.json",plan)
    put("expected.sha256",(plan_sha+"\n").encode())
    report_rows=[]
    for t in tasks:
        job=t["relative_output"]
        status="completed" if t["case_id"]==records[0]["id"] else "registration_excluded"
        receipt=dict(**t,registration_sha256=reg,plan_sha256=pack.canonical(plan),job_sha256=pack.canonical(t),
                     status=status,stdout_sha256=None,events_sha256=None,result_sha256=None)
        if status=="completed":
            receipt["result_sha256"]=put(job+"/result.json",dict(status="input_error"))
            receipt["events_sha256"]=put(job+"/events.jsonl",b"{}\n")
            receipt["stdout_sha256"]=put(job+"/stdout.log",b"invalid input\n")
            attempt="attempts/"+job.removeprefix("runs/")+"/001"
            receipt.update(source_attempt=attempt,launch_sha256=put(attempt+"/launch.json",{}),
                           completion_sha256=put(attempt+"/completion.json",{}))
            put(attempt+"/stdout.log",b"invalid input\n")
            put(attempt+"/command.json",[])
        put("observations/"+job.removeprefix("runs/")+"/001.started.json",{})
        report_rows.append(dict(**t,receipt_sha256=put(job+"/receipt.json",receipt),receipt_status=status,
                                result_sha256=receipt["result_sha256"],events_sha256=receipt["events_sha256"]))
    rs=put("stage3_report.json",dict(complete=True,expected_sha256=plan_sha,registration_sha256=reg,
           runs=report_rows,report_driver_sha256=report_driver,companion_sha256={k:text(k) for k in
           ("stage3_runs.csv","stage3_pairs.csv")}))
    replays=[]
    for c in records:
        ident=c["id"]
        f="replays/"+ident
        proto=put(f+"/protocol.json",dict(case_id=ident,search_report_sha256=rs,common=None,replay_driver_sha256=replay_driver))
        checksum=put(f+"/report.json",dict(case_id=ident,complete=True,protocol_sha256=proto,measurement_receipts_sha256={}))
        replays.append(dict(case_id=ident,report_sha256=checksum))
    put("replays/summary.json",dict(complete=True,search_report_sha256=rs,cases=replays))
    put("progress.json",dict(completed=108))
    text("runtime/temp/scratch")
    put("unknown-trace.txt",b"captured unknown trace\n")
    return raw_directory


def test_full_build_portable_verify_restore_external_and_reject_missing_export(pack,tmp_path):
    study=tmp_path/"studies"/"fixture"
    repo=tmp_path/"repo"
    raw=synthetic_study(pack,study,repo)
    bundle=tmp_path/"bundle"
    sha=pack.build(study,bundle,repo)
    other=tmp_path/"foreign"
    bundle.rename(other)
    detached=tmp_path/"public-inputs"
    raw.rename(detached)
    shutil.rmtree(study)
    shutil.rmtree(repo)
    with pack.open_payload(other,sha) as r:
        pack.validate_study(r)
        assert "unknown-trace.txt" in r.files
        assert "runtime/temp/scratch" not in r.files
        assert r.files["inputs/raw/d3d_shot_0000000000.parquet"]["external"].endswith("0000000000.parquet")
    out=tmp_path/"portable-copy"
    pack.restore(other,out,detached,sha)
    assert (out/"inputs/raw/d3d_shot_0000000000.parquet").read_bytes()==b"fake-0000000000"
    with pytest.raises(ValueError, match="exists"):
        pack.restore(other,out,detached,sha)
    next(detached.iterdir()).write_bytes(b"bad")
    with pytest.raises(ValueError, match="restored file differs"):
        pack.restore(other,tmp_path/"bad-inputs",detached,sha)
    assert not (tmp_path/"bad-inputs").exists()
    damaged=out/"runs/0000000000/20711/full-screened/events.jsonl"
    damaged.unlink()
    # The portable restore contains support and prior metadata; provide their published source locations for rebuilding.
    with pytest.raises(ValueError, match="missing/changed evidence"):
        with pack.open_payload(other) as r:
            del r.files[damaged.relative_to(out).as_posix()]
            pack.validate_study(r)
