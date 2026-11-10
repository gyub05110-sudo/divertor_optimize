import hashlib
from pathlib import Path
import pytest

from verification import stage5_dispatch as dispatch


def test_native_result_record_requires_both_original_files(tmp_path):
    job=("abc","d"*64,51031)
    raw, ledger=dispatch.result_paths(tmp_path,job)
    raw.parent.mkdir(parents=True)
    raw.write_bytes(b"sample native bytes")
    with pytest.raises((ValueError,FileNotFoundError)):
        dispatch.fingerprints(tmp_path,job)
    ledger.write_bytes(b'{"wall_W":1000000}')
    record=dispatch.fingerprints(tmp_path,job)
    assert record["raw_sha256"]==hashlib.sha256(raw.read_bytes()).hexdigest()
    assert record["native_receipt_sha256"]==hashlib.sha256(ledger.read_bytes()).hexdigest()


def test_existing_slot_receipt_must_match_seal_job_and_live_files(tmp_path):
    job=("case","e"*64,51031)
    raw, ledger=dispatch.result_paths(tmp_path,job)
    raw.parent.mkdir(parents=True)
    raw.write_bytes(b"raw"); ledger.write_bytes(b"{}")
    r={"job":list(job),"execution_sha256":"f"*64,**dispatch.fingerprints(tmp_path,job)}
    dispatch.validate_receipt(tmp_path,job,"f"*64,r)
    with pytest.raises(ValueError): dispatch.validate_receipt(tmp_path,job,"a"*64,r)
    raw.write_bytes(b"modified")
    with pytest.raises(ValueError): dispatch.validate_receipt(tmp_path,job,"f"*64,r)


def test_run_uses_preexisting_receipt_and_journals_last_job(tmp_path,monkeypatch):
    jobs=[("x","a"*64,51031),("x","b"*64,51031)]
    seal="f"*64
    monkeypatch.setattr(dispatch,"execution",lambda study:({"jobs":jobs},seal))
    def publish(study,job,identity):
        raw,ledger=dispatch.result_paths(study,job)
        raw.parent.mkdir(parents=True)
        raw.write_bytes(b"native")
        ledger.write_text("{}")
        dispatch.public_save(dispatch.receipt_path(study,job),dict(job=list(job),execution_sha256=identity,
                            **dispatch.fingerprints(study,job)))
    publish(tmp_path,jobs[0],seal)
    monkeypatch.setattr(dispatch,"execute",publish)
    dispatch.run(tmp_path,workers=2)
    progress=dispatch.read(tmp_path/"progress.json")
    assert progress["last"]==list(jobs[1]) and progress["completed"]==2
    dispatch.run(tmp_path,workers=2)
