"""Stage 4 cohort selection, historical source and EFIT-only registration tests."""
import hashlib
from importlib import import_module, util
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest


@pytest.fixture
def driver():
    name = "verification.stage4_inputs"
    assert util.find_spec(name), "Stage 4 cohort registration is required"
    return import_module(name)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def test_salted_draw_uses_24_unique_objects_after_all_path_checksum_and_id_aliases(driver):
    prefix = "data/diii_d_train/d3d_shot_"
    paths = [prefix + f"{i:010x}.parquet" for i in range(40)]
    hashes = [sha(path.encode()) for path in paths]
    hashes[2] = hashes[1]
    hashes[4] = hashes[3]
    hashes[6] = hashes[5]
    hashes[9] = hashes[8]
    entries = [dict(rfilename=p, size=9, lfs=dict(size=9, sha256=h)) for p, h in zip(paths, hashes)]
    sources = [dict(dataset_path=paths[1], sha256="a" * 64), dict(source_id="d3d_shot_0000000003", sha256="b" * 64),
               dict(published_shot_number="0000000005", sha256="c" * 64), dict(sha256=hashes[7])]
    selected, excluded = driver.draw(entries[::-1], sources)
    order = sorted(paths, key=lambda p: sha(("nf-stage4-hotspot-v1|" + p).encode()))
    seen = set(hashes[1:8])
    golden = []
    for p in order:
        h = hashes[paths.index(p)]
        if h not in seen:
            golden.append(p)
        seen.add(h)
    assert [row["dataset_path"] for row in selected] == golden[:24]
    assert set(paths[1:8]) <= set(excluded)
    assert len({r["sha256"] for r in selected}) == 24
    assert [r["rank"] for r in selected] == list(range(1, 25))
    assert all(r["bytes"] == 9 and r["sha256"] == hashes[paths.index(r["dataset_path"])] for r in selected)


def test_inventory_audit_hashes_renamed_inputs_and_fails_on_unregistered_real_archives(driver, tmp_path):
    source = tmp_path / "d3d_shot_known.parquet"
    source.write_bytes(b"known file")
    renamed = tmp_path / "renamed.parquet"
    renamed.write_bytes(source.read_bytes())
    fixture = tmp_path / "test_fixture0" / "case.parquet"
    fixture.parent.mkdir()
    fixture.write_bytes(b"dummy")
    audit = driver.audit_paths([source, renamed, fixture], {sha(b"known file")}, [tmp_path])
    assert audit == {sha(b"known file"): sorted([str(source), str(renamed)])}
    novel = tmp_path / "unknown.parquet"
    novel.write_bytes(b"unregistered file")
    with pytest.raises(ValueError, match="unregistered"):
        driver.audit_paths([novel], {sha(b"known file")}, [tmp_path])


def frame_file(path, times, frames):
    pq.write_table(pa.table({"source": ["DIII-D"], "efit_times": [times], "efit_psirz": [frames]}), path)


def test_non_target_frames_of_previous_missing_time_file_are_excluded_and_case_counts_preserved(driver, tmp_path):
    previous = tmp_path / "previous.parquet"
    frame = np.arange(6.).reshape(2, 3)
    frame_file(previous, [1450., 1500.], [frame.tolist(), (frame**2).tolist()])
    recs, index = driver.stage3.historical_frames([dict(local=previous, source_id="old", sha256=sha(previous.read_bytes()))])
    assert len(recs[0]["frames"]) == 2
    case = tmp_path / "new.parquet"
    frame_file(case, [1400.], [(frame**2 * 2. + 8.).tolist()])
    result = driver.assess(case, index)
    assert result["status"] == "duplicate_frame"
    assert result["exact_frame_count"] == 1
    assert result["duplicate_matches"]["normalized_sha256"]["time_ms"] == 1500.
    frame_file(case, [1349., 1450.], [frame.tolist()] * 2)
    assert driver.assess(case, index)["exact_frame_count"] == 0
    frame_file(case, [1400., 1400.], [frame.tolist()] * 2)
    assert driver.assess(case, index)["exact_frame_count"] == 2
    frame_file(case, [1400.], [(frame**3).tolist()])
    assert driver.assess(case, index)["status"] == "ready"
def test_new_files_also_index_off_target_frames_including_missing_time(driver, tmp_path):
    index = driver.shared.FrameIndex()
    frame = np.arange(6.).reshape(2, 3)
    old = tmp_path / "other.parquet"
    frame_file(old, [1300., 1500.], [frame.tolist(), (frame**2).tolist()])
    outcome = driver.register_frames(old, index, dict(id="new-old", sha256=sha(old.read_bytes())))
    assert outcome["status"] == "missing_time" and outcome["file_frame_count"] == 2
    new = tmp_path / "last.parquet"
    frame_file(new, [1400.], [(frame**2).tolist()])
    outcome = driver.register_frames(new, index, dict(id="new-last", sha256=sha(new.read_bytes())))
    assert outcome["status"] == "duplicate_frame"
    assert outcome["duplicate_matches"]["dimensional_sha256"] == dict(
        id="new-old", sha256=sha(old.read_bytes()), frame_index=1, time_ms=1500.)


def test_missing_download_retained_without_replacement_and_anchored_to_seal(driver, tmp_path):
    row = dict(rank=1, id="x", dataset_path="data/diii_d_train/d3d_shot_x.parquet", sha256="d"*64,
               bytes=100, url="https://example.test/pinned")
    driver.shared.publish(tmp_path / "selection-seal.json", dict(frozen_utc="2026-11-09T00:00:00+00:00",
                                                                selection_sha256="a"*64))
    log = tmp_path / "requests/x"
    start = dict(number=1, started_utc="2026-11-09T00:01:00+00:00", selection_sha256="a"*64,
                 requested_url=row["url"])
    driver.shared.publish(log / "1.started.json", start)
    driver.shared.publish(log / "1.finished.json", dict(**start, finished_utc="2026-11-09T00:02:00+00:00",
                                                        error="HTTP 404", http_status=404))
    record = driver.case_record(tmp_path, row, "a"*64, driver.shared.FrameIndex())
    assert (record["rank"],record["id"],record["status"]) == (1,"x","acquisition_failed")
    with pytest.raises(ValueError, match="selection"):
        driver.case_record(tmp_path, row, "b"*64, driver.shared.FrameIndex())


def test_local_train_filename_aliases_are_excluded_by_identifier_and_public_lfs_alias(driver):
    paths = [f"data/diii_d_train/d3d_shot_{i:010x}.parquet" for i in range(35)]
    hashes = [sha(s.encode()) for s in paths]
    hashes[3] = hashes[2]
    entries = [dict(rfilename=p, size=9, lfs=dict(size=9, sha256=h)) for p,h in zip(paths, hashes)]
    local = {hashes[0]: [str(Path("D:/archive") / Path(paths[2]).name)]}
    aliases = driver.local_aliases(local)
    _, excluded = driver.draw(entries, [dict(source_id="0000000000", sha256=hashes[0]), *aliases])
    assert paths[2] in excluded and paths[3] in excluded


def verify_published_study(inputs):
    """Independent byte-level audit of the frozen cohort, with no production driver imports."""
    import json
    import re

    def read(p):
        return json.loads(p.read_text(encoding="utf-8"))

    def rows(raw):
        table = pq.read_table(raw, columns=["source", "efit_times", "efit_psirz"])
        assert table["source"].to_pylist() == ["DIII-D"]
        times = table["efit_times"][0].as_py()
        grids = table["efit_psirz"][0].values
        assert len(times) == len(grids)
        for t, grid in zip(times, grids):
            lines = grid.values
            widths = np.diff(lines.offsets.to_numpy())
            assert np.all(widths == widths[0])
            psi = lines.flatten().to_numpy().astype("<f8").reshape(len(lines), int(widths[0]))
            assert np.isfinite(psi).all() and psi.max() > psi.min()
            psi[psi == 0] = 0
            scaled = (psi - psi.min()) / np.ptp(psi)
            shape = np.asarray(psi.shape, dtype="<u8").tobytes()
            yield dict(time_ms=t, shape=list(psi.shape), dimensional_sha256=sha(shape + psi.tobytes()),
                       normalized_sha256=sha(shape + scaled.astype("<f8").tobytes()))

    selection = read(inputs / "selection.json")
    seal = read(inputs / "selection-seal.json")
    reg = read(inputs / "registration.json")
    assert sha((inputs / "selection.json").read_bytes()) == seal["selection_sha256"] == reg["selection_sha256"]
    assert sha((inputs / "selection-seal.json").read_bytes()) == reg["selection_seal_sha256"]
    for name, expected in reg["files_sha256"].items():
        assert sha((inputs / name).read_bytes()) == expected, name
    previous = inputs.parent.parent / "nf-stage3-20261108" / "inputs"
    stage3 = read(previous / "registration.json")
    old_selection = read(previous / "selection.json")
    sources = selection["historical_sources"]
    assert sources[:-12] == old_selection["historical_sources"]
    assert [s["sha256"] for s in sources[-12:]] == [c["sha256"] for c in stage3["cases"]]
    digests = {s["sha256"] for s in sources}
    ids = {s["source_id"].removeprefix("d3d_shot_") for s in sources}
    ids.update(str(s["published_shot_number"]) for s in sources if "published_shot_number" in s)
    audit = read(inputs / "local-source-audit.json")["verified_copies"]
    assert set(audit) == digests
    for h, paths in audit.items():
        for path in map(Path, paths):
            assert sha(path.read_bytes()) == h, path
            m = re.fullmatch(r"d3d_shot_(.+)\.parquet", path.name)
            if m:
                ids.add(m.group(1))
    inventory = read(inputs / "hub-inventory.json")
    training = {e["rfilename"]: e["lfs"] for e in inventory["siblings"]
                if re.fullmatch(r"data/diii_d_train/d3d_shot_[A-Za-z0-9_-]+\.parquet", e["rfilename"])}
    banned = {s["dataset_path"] for s in sources if "dataset_path" in s}
    banned.update(p for p in training if Path(p).stem[9:] in ids)
    digests.update(training[p]["sha256"] for p in banned)
    ordinal = lambda p: sha(("nf-stage4-hotspot-v1|" + p).encode("utf-8"))
    order = []
    for path in sorted(training, key=lambda p: (ordinal(p), p)):
        if path not in banned and training[path]["sha256"] not in digests:
            order.append(path)
        digests.add(training[path]["sha256"])
    assert [c["dataset_path"] for c in reg["cases"]] == order[:24] == [s["dataset_path"] for s in selection["selected"]]

    index = read(inputs / "prior-frame-index.json")
    dimensional, normalized = set(), set()
    n_hist = 0
    for s, expected in zip(sources, index, strict=True):
        raw = Path(s["local"])
        assert sha(raw.read_bytes()) == s["sha256"] == expected["sha256"] and raw.stat().st_size == s["bytes"]
        for i, (row, actual) in enumerate(zip(expected["frames"], rows(raw), strict=True)):
            assert row == dict(frame_index=i, **actual)
            dimensional.add(actual["dimensional_sha256"])
            normalized.add(actual["normalized_sha256"])
            n_hist += 1
    assert n_hist == reg["historical_frames"] == reg["valid_historical_frames"]
    counts = []
    for rank, case in enumerate(reg["cases"], 1):
        raw = Path(case["path"])
        lfs = training[case["dataset_path"]]
        assert case["rank"] == rank and case["time_ms"] == 1400.
        assert raw.stat().st_size == case["bytes"] == lfs["size"]
        assert sha(raw.read_bytes()) == case["sha256"] == lfs["sha256"]
        assert case["order_sha256"] == ordinal(case["dataset_path"])
        receipt = inputs / "receipts" / (case["id"] + ".json")
        assert sha(receipt.read_bytes()) == case["receipt_sha256"]
        assert read(receipt)["selection_sha256"] == reg["selection_sha256"]
        log = inputs / "requests" / case["id"]
        for n, attempt in enumerate(case["attempts"], 1):
            start, end = log / f"{n}.started.json", log / f"{n}.finished.json"
            assert sha(start.read_bytes()) == attempt["started_sha256"]
            assert sha(end.read_bytes()) == attempt["finished_sha256"]
            assert read(start)["selection_sha256"] == reg["selection_sha256"]
            assert seal["frozen_utc"] <= read(start)["started_utc"] <= read(end)["finished_utc"]
        file_rows = list(rows(raw))
        target = [i for i, r in enumerate(file_rows) if r["time_ms"] == 1400.]
        assert len(target) == case["exact_frame_count"] and len(file_rows) == case["file_frame_count"]
        assert len(file_rows) == case["valid_file_frames"]
        if len(target) == 1:
            i = target[0]
            fp = file_rows[i]
            assert i == case["frame_index"]
            assert all(fp[k] == v for k,v in case["fingerprint"].items())
            duplicate = fp["dimensional_sha256"] in dimensional or fp["normalized_sha256"] in normalized
            assert case["status"] == ("duplicate_frame" if duplicate else "ready")
        else:
            assert case["status"] == ("missing_time" if not target else "ambiguous_time")
        dimensional.update(r["dimensional_sha256"] for r in file_rows)
        normalized.update(r["normalized_sha256"] for r in file_rows)
        counts.append({k: case[k] for k in ("rank", "id", "file_frame_count", "exact_frame_count", "status")})
    return dict(registration_sha256=sha((inputs / "registration.json").read_bytes()),
                selection_sha256=reg["selection_sha256"], seal_sha256=reg["selection_seal_sha256"],
                historical_sources=len(sources), historical_copies=sum(map(len, audit.values())),
                historical_frames=n_hist, per_file_counts=counts)


if __name__ == "__main__":
    import json
    import sys

    print(json.dumps(verify_published_study(Path(sys.argv[1])), indent=2))
