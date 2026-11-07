import hashlib
from importlib import import_module, util

import numpy as np
import pytest


@pytest.fixture
def driver():
    name = "verification.nf_prospective_inputs"
    assert util.find_spec(name), "prospective acquisition driver is required"
    return import_module(name)


def test_order_uses_full_posix_path_and_excludes_before_ranking(driver):
    folder = "data/diii_d_train/"
    paths = [folder + f"d3d_shot_{i:010x}.parquet" for i in range(90)]
    banned = {paths[3], paths[47]}
    golden = sorted(set(paths) - banned, key=lambda p: (
        hashlib.sha256(("nf-prospective-20261107-v1/" + p).encode()).hexdigest(), p))
    assert driver.ordered_paths(paths[::-1] + ["data/other/d3d_shot_000.parquet"], banned) == golden
    with pytest.raises(ValueError, match="duplicat"):
        driver.ordered_paths(paths + [paths[3]], banned)


def test_shape_signed_zero_and_affine_duplicate_policy(driver):
    frame = np.array([[-0., 1.], [2., 3.]], dtype=">f4")
    first = driver.frame_fingerprint(frame)
    frame[0, 0] = +0.
    assert first == driver.frame_fingerprint(frame)
    changed_units = driver.frame_fingerprint(frame * 8. + 13.)
    assert first["dimensional_sha256"] != changed_units["dimensional_sha256"]
    assert first["normalized_sha256"] == changed_units["normalized_sha256"]
    assert first["dimensional_sha256"] != driver.frame_fingerprint(frame.reshape(1, 4))["dimensional_sha256"]
    index = driver.FrameIndex()
    index.add(first, {"source_id": "old", "frame_index": 17})
    assert index.matches(changed_units) == {"normalized_sha256": {"source_id": "old", "frame_index": 17}}
    with pytest.raises(ValueError, match="finite|constant"):
        driver.frame_fingerprint([[2., 2.], [2., 2.]])


def test_staged_metadata_has_only_geometry_and_no_historical_warm_shapes(driver):
    source = {"wall": [[1, 2]], "optimization_policy": {"a": True},
              "optimization_footprints": [{"source": "IMP", "segment": 69, "points": [[1, 2]],
                                            "sampled_hit_count": 22}], "sources": {"OMP": {"heat": 13}},
              "candidate": [5]}
    minimal = driver.wall_definition(source)
    assert set(minimal) == {"wall", "optimization_policy", "optimization_footprints"}
    assert minimal["optimization_footprints"] == [{"source": "IMP", "segment": 69, "points": [[1, 2]]}]
    assert driver.input_manifest(1400., "new-shot", "old-geometry") == {
        "time_ms": 1400., "shot_sha256": "new-shot", "geometry": {"sha256": "old-geometry"}}


def test_verified_receipt_rejects_replacement_bytes_and_path_escape(driver, tmp_path):
    path = "data/diii_d_train/d3d_shot_0123456789.parquet"
    raw = tmp_path / "raw" / "d3d_shot_0123456789.parquet"
    raw.parent.mkdir()
    raw.write_bytes(b"source")
    checksum = hashlib.sha256(raw.read_bytes()).hexdigest()
    entry = {"rfilename": path, "size": 6, "lfs": {"size": 6, "sha256": checksum}}
    receipt = {"path": path, "url": driver.raw_url(path), "local": str(raw),
               "bytes": 6, "sha256": checksum, "selection_sha256": "registered"}
    driver.verify_receipt(tmp_path, entry, "registered", receipt)
    outside = tmp_path / "also-valid.parquet"
    outside.write_bytes(raw.read_bytes())
    with pytest.raises(ValueError, match="path"):
        driver.verify_receipt(tmp_path, entry, "registered", {**receipt, "local": str(outside)})
    raw.write_bytes(b"tamper")
    with pytest.raises(ValueError, match="hash|checksum"):
        driver.verify_receipt(tmp_path, entry, "registered", receipt)
    with pytest.raises(ValueError, match="selection"):
        driver.verify_receipt(tmp_path, entry, "different", receipt)


def test_failed_receipt_requires_three_finished_unsuccessful_requests(driver, tmp_path):
    path = "data/diii_d_train/d3d_shot_0123456789.parquet"
    entry = {"rfilename": path, "size": 6, "lfs": {"size": 6, "sha256": "a" * 64}}
    receipt = {"path": path, "url": driver.raw_url(path), "error": "all attempts failed",
               "attempts": [], "selection_sha256": "registered"}
    with pytest.raises(ValueError, match="failed attempts"):
        driver.verify_receipt(tmp_path, entry, "registered", receipt)


def test_unrelated_pytest_fixture_cannot_become_an_excluded_dataset_input(driver, tmp_path):
    fixture = tmp_path / "d3d_shot_0123456789.parquet"
    fixture.write_bytes(b"fake test data")
    extra = driver.additional_archives([fixture], set(), {})
    assert extra == []


def test_mid_request_resume_counts_interruption_and_preserves_pinned_bytes(driver, tmp_path, monkeypatch):
    source = b"new input"
    sha = hashlib.sha256(source).hexdigest()
    path = "data/diii_d_train/d3d_shot_456789abcd.parquet"
    entry = {"rfilename": path, "size": len(source), "lfs": {"size": len(source), "sha256": sha}}
    log = tmp_path / "requests" / "d3d_shot_456789abcd"
    driver.publish(log / "1.started.json", {"number": 1, "started_utc": "before interruption",
                                          "url": driver.raw_url(path), "selection_sha256": "registration"})

    class Response:
        status_code = 200
        content = source
        def __enter__(self):
            return self
        def __exit__(self, *exception):
            return False

    monkeypatch.setattr(driver.requests, "get", lambda *args, **kwargs: Response())
    record = driver.download(tmp_path, entry, "registration")
    assert len(record["attempts"]) == 2 and "interrupted" in record["attempts"][0]["error"]
    assert "error" not in record["attempts"][1]
    driver.audit_requests(tmp_path, record)
    (log / "1.finished.json").write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="request checkpoint"):
        driver.audit_requests(tmp_path, record)
