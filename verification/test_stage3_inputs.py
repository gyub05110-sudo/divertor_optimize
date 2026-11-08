import hashlib
from importlib import import_module, util
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest


@pytest.fixture
def driver():
    name = "verification.stage3_inputs"
    assert util.find_spec(name), "stage-three input registration is required"
    return import_module(name)


def test_deterministic_draw_excludes_historical_public_aliases_and_duplicates(driver):
    prefix = "data/diii_d_train/d3d_shot_"
    paths = [prefix + f"{i:010x}.parquet" for i in range(24)]
    sums = [hashlib.sha256(p.encode()).hexdigest() for p in paths]
    sums[2] = sums[1]  # public alias of a historical local file with different bytes
    sums[4] = sums[3]  # two eligible paths point to the same LFS object
    entries = [dict(rfilename=p, size=8, lfs=dict(size=8, sha256=s)) for p,s in zip(paths,sums)]
    historical = [dict(dataset_path=paths[1], sha256="f" * 64)]
    rows, excluded = driver.draw(entries[::-1], historical, 12)
    order = sorted(paths, key=lambda p: hashlib.sha256(("nf-stage3-v1|" + p).encode()).digest())
    seen = {sums[1]}
    golden = []
    for path in order:
        checksum = sums[paths.index(path)]
        if checksum not in seen:
            golden.append(path)
        seen.add(checksum)
    assert [row["dataset_path"] for row in rows] == golden[:12]
    assert len({row["sha256"] for row in rows}) == 12
    assert paths[1] in excluded and paths[2] in excluded


def write_frame_file(raw, times, frames):
    pq.write_table(pa.table({"source": ["DIII-D"], "efit_times": [times], "efit_psirz": [frames]}), raw)


def test_only_exact_unique_time_and_all_prior_frames_can_decide_readiness(driver, tmp_path):
    index = driver.shared.FrameIndex()
    frame = np.arange(6.).reshape(2, 3)
    index.add(driver.shared.frame_fingerprint(frame), {"id": "earlier", "time_ms": 1575.0})
    raw = tmp_path / "case.parquet"
    write_frame_file(raw, [1399., 1401.], [frame.tolist()] * 2)
    assert driver.assess(raw, index)["status"] == "missing_time"
    write_frame_file(raw, [1400., 1400.], [frame.tolist()] * 2)
    assert driver.assess(raw, index)["status"] == "ambiguous_time"
    write_frame_file(raw, [1400.], [(frame * 8. + 16.).tolist()])
    duplicate = driver.assess(raw, index)
    assert duplicate["status"] == "duplicate_frame"
    assert duplicate["duplicate_matches"]["normalized_sha256"]["id"] == "earlier"
    write_frame_file(raw, [1400.], [(frame ** 2).tolist()])
    assert driver.assess(raw, index)["status"] == "ready"
    write_frame_file(raw, [1400.], [np.zeros_like(frame).tolist()])
    assert driver.assess(raw, index)["status"] == "invalid_input"


def test_no_substitution_for_unavailable_or_tampered_pinned_bytes(driver, tmp_path):
    raw = tmp_path / "case.parquet"
    raw.write_bytes(b"same size")
    expected = dict(sha256=hashlib.sha256(b"diff size").hexdigest(), bytes=9)
    with pytest.raises(ValueError, match="checksum"):
        driver.check_raw(raw, expected)
    raw.write_bytes(b"diff size")
    driver.check_raw(raw, expected)

