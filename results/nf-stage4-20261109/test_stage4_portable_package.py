"""Exercise portable Stage 4 archives exclusively with synthetic source bytes."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.fixture
def pack():
    path = Path(__file__).with_name("package_evidence.py")
    assert path.is_file(), "Stage 4 packager is missing"
    spec = importlib.util.spec_from_file_location("stage4_archiver", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def study_fixture(root):
    def put(name, obj):
        raw = obj if isinstance(obj, bytes) else json.dumps(obj, indent=2).encode()
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        return hashlib.sha256(raw).hexdigest()
    def canonical(obj):
        return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    cases = []
    protocol = put("inputs/protocol.ko.md", b"synthetic protocol")
    report_driver = put("source-code/verification/stage4_report.py", b"# synthetic report")
    replay_driver = put("source-code/verification/stage4_replay.py", b"# synthetic replay")
    fixed_driver = put("source-code/verification/stage3_replay.py", b"# synthetic fixed evaluator")
    for rank in range(1, 25):
        cid = f"{rank:010}"
        filename = f"d3d_shot_{cid}.parquet"
        raw = f"synthetic source {cid}".encode()
        case = dict(id=cid, rank=rank, dataset_path="data/diii_d_train/"+filename,
                    sha256=put("inputs/raw/"+filename, raw), bytes=len(raw), time_ms=1400,
                    status="ready" if rank == 1 else "missing_time", path=str(root/"inputs/raw"/filename))
        case["receipt_sha256"] = put(f"inputs/receipts/{cid}.json", dict(case))
        cases.append(case)
    selection = put("inputs/selection.json", dict(selected=[{k:r[k] for k in ("id","rank","bytes","sha256")}
                                                            for r in cases], protocol_sha256=protocol))
    sel_seal = put("inputs/selection-seal.json", dict(selection_sha256=selection, protocol_sha256=protocol))
    reg = put("inputs/registration.json", dict(cases=cases, selection_sha256=selection, protocol_sha256=protocol,
                                              files_sha256={"protocol.ko.md":protocol},selection_seal_sha256=sel_seal))
    seeds = (30713, 30731, 30757)
    arms = ("hotspot","full","floor-only","floor-connector","floor-shelf","floor-entry")
    jobs = [dict(case_id=c["id"], seed=seed, arm=arm, registration_status=c["status"], sha256=c["sha256"],
                 relative_output=f"runs/{c['id']}/{seed}/{arm}") for c in cases for seed in seeds
            for arm in (arms if seed == seeds[0] else arms[:2])]
    files = {"verification/stage4_report.py":report_driver, "verification/stage4_replay.py":replay_driver,
             "verification/stage3_replay.py":fixed_driver}
    plan = dict(schema="nf-stage4-hotspot.v1", registration_sha256=reg, tasks=jobs, files=files)
    plan_sha = put("expected.json", plan)
    put("expected.sha256", (plan_sha+"\n").encode())
    rows = []
    for job in jobs:
        out = job["relative_output"]
        ident = {k:job[k] for k in ("case_id","seed","arm")}
        receipt = dict(ident, status="registration_excluded", job_sha256=canonical(job), plan_sha256=canonical(plan),
                       registration_sha256=reg, events_sha256=None, stdout_sha256=None, result_sha256=None)
        if job["registration_status"] == "ready":
            attempt = out.replace("runs/", "attempts/")+"/001"
            exports = {name:put(out+"/baseline/"+name, (name+out).encode())
                       for name in ("geometry.csv", "donors.csv.gz", "profile.csv.gz", "boundary.npz")}
            receipt.update(status="completed", events_sha256=put(out+"/events.jsonl", b'{"test":"journal"}\n'),
                           result_sha256=put(out+"/result.json", dict(status="no_design",
                                             baseline=dict(directory="baseline",sha256=exports))),
                           stdout_sha256=put(out+"/stdout.log",b"ok"), source_attempt=attempt,
                           completion_sha256=put(attempt+"/completion.json",{}),
                           launch_sha256=put(attempt+"/launch.json",{}))
            put(attempt+"/command.json", [])
        rec = put(out+"/receipt.json",receipt)
        rows.append(dict(ident, relative_output=out, receipt_status=receipt["status"], receipt_sha256=rec,
                         events_sha256=receipt["events_sha256"],result_sha256=receipt["result_sha256"]))
    companion = {name:put(name, b"synthetic,column\n") for name in
                 ("policy-predictions.csv","product-pairs.csv","registered-runs.csv")}
    rep = put("stage4_report.json", dict(complete=True, plan_sha256=plan_sha, registration_sha256=reg, runs=rows,
                                        report_driver_sha256=report_driver,companion_sha256=companion))
    summaries = []
    for c in cases:
        cid = c["id"]
        folder = "replays/"+cid+"/"
        active = c["rank"] == 1
        psha = put(folder+"protocol.json",dict(case_id=cid, search_report_sha256=rep, common={"test":True} if active else None,
                                              candidates={"testcandidate":{}} if active else {},
                                              replay_driver_sha256=replay_driver,fixed_wall_driver_sha256=fixed_driver))
        grids, checkpoints = {}, {}
        for level in ("cubic1280","cubic640","quintic1280") if active else ():
            files = {"setup.json":put(folder+level+"/setup.json",{"status":"measured"})}
            for role in ("baseline","testcandidate"):
                where = folder+level+"/"+role+"/"
                exports = {n:put(where+n, ("synthetic "+role+level+n).encode())
                           for n in ("geometry.csv","boundary.npz","donors.csv.gz","profile.csv.gz")}
                files[role+"/measurement.json"] = put(where+"measurement.json",dict(directory=role,sha256=exports))
            grids[level] = dict(status="measured")
            checkpoints[level] = put(folder+level+".json",dict(protocol_sha256=psha,files_sha256=files,grid=grids[level]))
        rsha = put(folder+"report.json",dict(case_id=cid,complete=True,protocol_sha256=psha,search_report_sha256=rep,
                                             grid_checkpoints_sha256=checkpoints,results=grids))
        summaries.append(dict(case_id=cid,report_sha256=rsha))
    summary = dict(complete=True,search_report_sha256=rep,cases=summaries)
    put("replays/summary.json",summary)
    put("stage4_replay.json",summary)
    put("runtime/cache.tmp",b"machine cache")
    return root


def test_portable_compact_and_full_roundtrip(pack, tmp_path):
    study = study_fixture(tmp_path/"study")
    digests = []
    for mode in ("compact", "full"):
        out = tmp_path/mode
        sha = pack.build(study,out,mode)
        portable = tmp_path/("detached-"+mode)
        out.rename(portable)
        result = subprocess.run([sys.executable,"-B",str(portable/"package_evidence.py"),"verify","--bundle",
                                 str(portable),"--manifest-sha256",sha],capture_output=True,text=True)
        assert result.returncode==0, result.stderr
        manifest = json.loads((portable/"stage4-manifest.json").read_text())
        digests.append(manifest["inventory_sha256"])
        restore = tmp_path/("restored-"+mode)
        pack.restore(portable,restore,sha)
        event = Path("runs/0000000001/30713/hotspot/events.jsonl")
        assert (restore/event).read_bytes() == (study/event).read_bytes()
        export = event.parent/"baseline/boundary.npz"
        assert (restore/export).is_file() == (mode == "full")
        raw = Path("inputs/raw/d3d_shot_0000000001.parquet")
        assert (restore/raw).is_file() == (mode == "full")
        replay = Path("replays/0000000001/quintic1280/testcandidate/boundary.npz")
        assert (restore/replay).is_file() == (mode == "full")
        assert not (restore/"runtime").exists()
        with pytest.raises(ValueError, match="destination exists"):
            pack.restore(portable,restore,sha)
        if mode == "compact":
            full = tmp_path/"restored-with-study"
            pack.restore(portable,full,sha,study)
            assert (full/export).read_bytes() == (study/export).read_bytes()
            (study/raw).write_bytes(b"altered")
            with pytest.raises(ValueError, match="differ"):
                pack.restore(portable,tmp_path/"corrupt",sha,study)
            assert not (tmp_path/"corrupt").exists()
            (study/raw).write_bytes((full/raw).read_bytes())
    assert digests[0] == digests[1]


def test_incomplete_campaign_and_missing_declared_export_do_not_publish(pack, tmp_path):
    study = study_fixture(tmp_path/"study")
    replay = study/"replays/summary.json"
    saved = replay.read_bytes()
    replay.write_bytes(b'{"complete":false}')
    with pytest.raises(ValueError, match="complete"):
        pack.build(study,tmp_path/"incomplete","compact")
    assert not (tmp_path/"incomplete").exists()
    replay.write_bytes(saved)
    next((study/"runs").rglob("boundary.npz")).unlink()
    with pytest.raises(ValueError, match="missing|changed"):
        pack.build(study,tmp_path/"missing","compact")
    assert not (tmp_path/"missing").exists()


def test_replay_export_corruption_is_detected_before_build(pack, tmp_path):
    study = study_fixture(tmp_path/"study")
    (study/"replays/0000000001/cubic1280/testcandidate/boundary.npz").write_bytes(b"changed")
    with pytest.raises(ValueError, match="changed"):
        pack.build(study,tmp_path/"bad-replay","full")
    assert not (tmp_path/"bad-replay").exists()
