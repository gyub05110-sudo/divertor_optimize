"""Preregister, acquire and screen fresh EFIT input frames using geometry alone."""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import struct
import subprocess
import tempfile

import numpy as np
import pyarrow.parquet as pq
import requests

from d3d_optimization.broadening_study import digest, load_json as read
from d3d_optimization.case import load_case
from verification import nf_unused_inputs as prior
from verification import nf_unused_live as live


ROOT = live.CODE_ROOT
STUDIES = Path("D:/DivertorOptimization/studies")
STUDY = STUDIES / "nf-prospective-inputs-20261107"
PREVIOUS = STUDIES / "unused-equilibrium-20261106"
EXAMPLES = STUDIES / "firefly-baseline-20261105/common_rank"
DATASET = "Sophelio/fusion-equilibrium-challenge"
REVISION = "1e280905b85f2a6fdde7e06fca8cf3a1edf447cb"
PREFIX = "nf-prospective-20261107-v1/"
TRAIN = "data/diii_d_train/"
SHOT_NAME = re.compile(r"d3d_shot_([a-zA-Z0-9_-]+)\.parquet")
HUB = f"https://huggingface.co/datasets/{DATASET}"
API = f"https://huggingface.co/api/datasets/{DATASET}/revision/{REVISION}?blobs=true"
SOURCE_RECORDS = ("protocol.ko.md", "selection.json", "downloads.json", "preflight.json")
SCAN_ROOTS = (Path("D:/DivertorOptimization"), Path("C:/Users/FPL-user/Documents/Divertor Optimization"))
REFERENCE_RECORDS = {"selection.json": "9ceffedce5633cbf7d39b2ab928a6f8bb2f53d4cf9b18295e5e4f318c19df3ee",
                     "downloads.json": "d2db95748ae57a2e49953d0c77a4f03d5435da6291d8541ce3f64273bad4d562"}
DRIVERS = ("verification/nf_prospective_inputs.py", "verification/nf_unused_inputs.py",
           "verification/nf_equilibrium_times.py", "verification/nf_unused_live.py",
           "d3d_optimization/case.py", "d3d_optimization/frozen_field.py",
           "d3d_optimization/frozen_mapping.py", "d3d_optimization/continuous_boundary.py",
           "d3d_optimization/independent_sources.py", "d3d_optimization/broadening_study.py")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def raw_url(path):
    return f"{HUB}/resolve/{REVISION}/{path}"


def publish(path, data):
    """Publish complete immutable bytes. Identical content supports idempotent repair."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not isinstance(data, bytes):
        data = (json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode()
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".writing-", delete=False) as temp:
        tmp = Path(temp.name)
        try:
            temp.write(data)
            temp.flush()
            os.fsync(temp.fileno())
        except BaseException:
            temp.close()
            tmp.unlink()
            raise
    try:
        try:
            os.link(tmp, path)
        except FileExistsError:
            require(path.read_bytes() == data, f"immutable record differs: {path}")
    finally:
        tmp.unlink()


@contextmanager
def study_lock(study):
    """OS-released process lock also permits recovery after an interrupted worker."""
    import msvcrt
    with (study / "acquisition.lock").open("a+b") as stream:
        if stream.tell() == 0:
            stream.write(b"0")
            stream.flush()
        stream.seek(0)
        msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        try:
            yield
        finally:
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)


def ordered_paths(paths, excluded):
    pool = [p for p in paths if p.startswith(TRAIN) and SHOT_NAME.fullmatch(p[len(TRAIN):])]
    require(len(set(pool)) == len(pool), "duplicate training path")
    return sorted(set(pool) - set(excluded),
                  key=lambda p: (hashlib.sha256((PREFIX + p).encode()).hexdigest(), p))


def driver_hashes():
    return {p: digest(ROOT / p) for p in DRIVERS}


def source_file(path, expected=None):
    actual = digest(path)
    if expected:
        require(actual == expected, f"historical file hash mismatch: {path}")
    return dict(local=str(path.resolve()), sha256=actual, bytes=path.stat().st_size)


def additional_archives(paths, known_names, entries):
    grouped = {}
    for p in paths:
        dataset_path = TRAIN + p.name
        if p.name in known_names or dataset_path not in entries or not SHOT_NAME.fullmatch(p.name):
            continue
        expected = entries[dataset_path]["lfs"]["sha256"]
        if p.stat().st_size == entries[dataset_path]["size"] and digest(p) == expected:
            grouped.setdefault(p.name, []).append(p)
    return [dict(source_id=Path(name).stem, dataset_path=TRAIN + name, **source_file(min(copies)))
            for name, copies in sorted(grouped.items())]


def previous_sources(study, entries):
    for name, value in REFERENCE_RECORDS.items():
        source_file(PREVIOUS / name, value)
    old = read(PREVIOUS / "selection.json")
    downloads = read(PREVIOUS / "downloads.json")
    require(old["revision"] == downloads["revision"] == REVISION and len(old["selected"]) == 12
            and old["selected"] == [r["path"] for r in downloads["files"]], "historical selection changed")
    original_shot = ROOT / live.SHOT
    original_hash = read(PREVIOUS / "preflight.json")["original_sha256"][str(Path(live.SHOT))]
    rows = [dict(source_id=original_shot.stem, dataset_path=TRAIN + original_shot.name,
                 **source_file(original_shot, original_hash))]
    for r in downloads["files"]:
        rows.append(dict(source_id=Path(r["path"]).stem, dataset_path=r["path"],
                         **source_file(Path(r["local"]), r["sha256"])))
    for r in read(EXAMPLES / "public_shot_comparison.json")["examples"]:
        p = EXAMPLES / "public_examples" / f"{r['shot']}.parquet"
        rows.append(dict(source_id=f"starter_{r['shot']}", published_shot_number=r["shot"],
                         source_url=r["url"], **source_file(p, r["lfs_sha256"])))

    # Aliases of the original filename also occur in old staged and temporary roots.
    # Identify other raw input filenames while the new input study is still empty.
    paths = subprocess.run(["rg", "--files", "--hidden", "--no-ignore", "-g", "*.parquet",
                            "-g", "!**/.git/**", "-g", "!**/.local_dependencies/**",
                            "-g", "!**/node_modules/**", *map(str, SCAN_ROOTS)],
                           check=True, capture_output=True, text=True).stdout.splitlines()
    paths = [Path(p).resolve() for p in paths if not Path(p).resolve().is_relative_to(study)]
    names = {Path(r["local"]).name for r in rows}
    extra = additional_archives(paths, names, entries)
    rows.extend(extra)
    audit = dict(created_utc=utc(), roots=list(map(str, SCAN_ROOTS)), matching_filenames=sorted({p.name for p in paths}),
                 excluded_directories=[".git", ".local_dependencies", "node_modules", str(study)],
                 additional_anonymous_paths=[r["dataset_path"] for r in extra], files=rows,
                 referenced_records={str(p): digest(p) for p in
                     [*(PREVIOUS / n for n in (*REFERENCE_RECORDS, "preflight.json")),
                      EXAMPLES / "public_shot_comparison.json"]})
    publish(study / "prior-sources.json", audit)
    return audit


def inventory(study):
    path = study / "hub-inventory.json"
    if not path.exists():
        response = requests.get(API, timeout=60)
        response.raise_for_status()
        publish(path, response.content)
    data = read(path)
    require(data["id"] == DATASET and data["sha"] == REVISION and data["cardData"]["license"] == "cc-by-4.0",
            "pinned dataset inventory or license differs")
    siblings = data["siblings"]
    require(len(siblings) == len({r["rfilename"] for r in siblings}), "duplicate Hub path")
    return {r["rfilename"]: r for r in siblings}


def register(study):
    if (study / "registration.json").exists():
        return check_registration(study)
    require((study / "protocol.ko.md").is_file(), "copy the exclusive input protocol before registration")
    require(not any((study / p).exists() for p in ("cases", "raw", "requests", "receipts")),
            "candidate content exists before registration")
    entries = inventory(study)
    historical = (read(study / "prior-sources.json") if (study / "prior-sources.json").exists()
                  else previous_sources(study, entries))
    banned = {row["dataset_path"] for row in historical["files"] if "dataset_path" in row}
    checksums = {row["sha256"] for row in historical["files"]}
    banned.update(p for p, e in entries.items() if p.startswith(TRAIN)
                  and e.get("lfs", {}).get("sha256") in checksums)
    require(banned <= set(entries), "unknown historical training paths")
    order = ordered_paths(entries, banned)
    require(len(order) >= 72, "need 72 preregistered eligible filenames")
    publish(study / "ranked-paths.json", order)
    card = study / "dataset-README.md"
    if not card.exists():
        response = requests.get(raw_url("README.md"), timeout=60)
        response.raise_for_status()
        data = response.content
        blob = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
        require(blob == entries["README.md"]["blobId"], "dataset-card git blob differs")
        publish(card, data)
    if not (study / "selection.json").exists():
        selection = dict(schema="nf-prospective-input-selection.v1", created_utc=utc(),
                         dataset=DATASET, revision=REVISION, license="cc-by-4.0", api_url=API,
                         dataset_card_url=f"{HUB}/blob/{REVISION}/README.md",
                         hub_inventory_saved_utc=datetime.fromtimestamp((study / "hub-inventory.json").stat().st_mtime, timezone.utc).isoformat(),
                         protocol_saved_utc=datetime.fromtimestamp((study / "protocol.ko.md").stat().st_mtime, timezone.utc).isoformat(),
                         repository_file_count=len(entries), training_file_count=len(order) + len(banned),
                         population=len(order), salt=PREFIX, exclude=sorted(banned),
                         block_size=24, target_passes=3, selected=order[:24],
                         reserve_blocks=[order[24:48], order[48:72]],
                         protocol_sha256=digest(study / "protocol.ko.md"), drivers=driver_hashes(),
                         original_sha256={p: digest(ROOT / p) for p in live.INPUTS},
                         files_sha256={n: digest(study / n) for n in
                             ("hub-inventory.json", "prior-sources.json", "ranked-paths.json", "dataset-README.md")})
        publish(study / "selection.json", selection)
    publish(study / "registration.json", dict(selection_sha256=digest(study / "selection.json"),
            registered_utc=utc(), staged_or_downloaded_candidates=0))
    return check_registration(study)


def check_registration(study):
    selection = read(study / "selection.json")
    require(read(study / "registration.json")["selection_sha256"] == digest(study / "selection.json"),
            "selection checkpoint differs")
    require(selection["protocol_sha256"] == digest(study / "protocol.ko.md"), "protocol checksum differs")
    require(selection["drivers"] == driver_hashes(), "registered acquisition or geometry driver changed")
    require(selection["original_sha256"] == {p: digest(ROOT / p) for p in live.INPUTS}, "original inputs changed")
    for p, value in selection["files_sha256"].items():
        require(digest(study / p) == value, f"source record checksum differs: {p}")
    order = ordered_paths(inventory(study), selection["exclude"])
    require(order == read(study / "ranked-paths.json") and selection["selected"] == order[:24]
            and selection["reserve_blocks"] == [order[24:48], order[48:72]], "preregistered order differs")
    return selection


def frame_fingerprint(frame):
    arr = np.array(frame, dtype="<f8", order="C", copy=True)
    require(arr.ndim == 2 and arr.size and np.isfinite(arr).all(), "need finite 2D EFIT grid")
    low, high = float(arr.min()), float(arr.max())
    require(high > low, "constant EFIT grid")
    arr[arr == 0.] = 0.
    normalized = (arr - low) / (high - low)
    shape = struct.pack("<QQ", *arr.shape)
    return dict(shape=list(arr.shape),
                dimensional_sha256=hashlib.sha256(shape + arr.tobytes()).hexdigest(),
                normalized_sha256=hashlib.sha256(shape + normalized.astype("<f8", copy=False).tobytes()).hexdigest())


class FrameIndex:
    keys = ("dimensional_sha256", "normalized_sha256")

    def __init__(self):
        self.by_digest = {k: {} for k in self.keys}

    def add(self, fingerprint, location):
        for k in self.keys:
            self.by_digest[k].setdefault(fingerprint[k], location)

    def matches(self, fingerprint):
        return {k: self.by_digest[k][fingerprint[k]] for k in self.keys if fingerprint[k] in self.by_digest[k]}


def read_frames(path):
    t = pq.read_table(path, columns=["source", "efit_times", "efit_psirz"], use_threads=False)
    require(t.num_rows == 1 and t["source"][0].as_py() == "DIII-D", "expected one DIII-D record")
    times, frames = t["efit_times"][0].values.to_numpy(zero_copy_only=False), t["efit_psirz"][0]
    require(len(times) == len(frames), "EFIT time/frame lengths differ")
    return times, frames


def array_from_frame(frame):
    rows = frame.values
    sizes = np.diff(rows.offsets.to_numpy(zero_copy_only=False))
    require(len(sizes) and sizes[0] > 0 and np.all(sizes == sizes[0]), "irregular EFIT frame dimensions")
    return rows.flatten().to_numpy(zero_copy_only=False).reshape(len(sizes), sizes[0])


def archive_index(study):
    path = study / "prior-frame-index.json"
    old = read(study / "prior-sources.json")
    if not path.exists():
        records = []
        for source in old["files"]:
            p = Path(source["local"])
            source_file(p, source["sha256"])
            times, frames = read_frames(p)
            rows = []
            for i, time in enumerate(times):
                r = dict(frame_index=i, time_ms=float(time) if np.isfinite(time) else None)
                try:
                    r.update(frame_fingerprint(array_from_frame(frames[i])))
                except ValueError as exc:
                    r["invalid"] = str(exc)
                rows.append(r)
            records.append(dict(source_id=source["source_id"], file_sha256=source["sha256"], frames=rows))
        publish(path, dict(created_utc=utc(), prior_sources_sha256=digest(study / "prior-sources.json"),
                           selection_sha256=digest(study / "selection.json"), files=records))
    data = read(path)
    require(data["prior_sources_sha256"] == digest(study / "prior-sources.json")
            and data["selection_sha256"] == digest(study / "selection.json"), "archive index provenance differs")
    require([(r["source_id"], r["file_sha256"]) for r in data["files"]]
            == [(s["source_id"], s["sha256"]) for s in old["files"]], "prior frame source order differs")
    index = FrameIndex()
    for source in data["files"]:
        for row in source["frames"]:
            if "invalid" not in row:
                index.add(row, dict(source_id=source["source_id"], frame_index=row["frame_index"]))
    return index


def verify_receipt(study, entry, selection_hash, row):
    path = entry["rfilename"]
    require(row["selection_sha256"] == selection_hash, "receipt selection differs")
    require(row["path"] == path and row["url"] == raw_url(path), "receipt source path/URL differs")
    if "error" in row:
        require(len(row["attempts"]) == 3 and all(a.get("error") for a in row["attempts"]), "missing failed attempts")
        return
    raw = study / "raw" / PurePosixPath(path).name
    require(Path(row["local"]).resolve() == raw.resolve(), "unexpected raw file path")
    require(raw.stat().st_size == row["bytes"] == entry["size"] == entry["lfs"]["size"], "source size mismatch")
    require(digest(raw) == row["sha256"] == entry["lfs"]["sha256"], "source hash mismatch")


def download(study, entry, selection_hash):
    path = entry["rfilename"]
    name = PurePosixPath(path).name
    stem = Path(name).stem
    raw = study / "raw" / name
    receipt = study / "receipts" / (stem + ".json")
    if receipt.exists():
        row = read(receipt)
        verify_receipt(study, entry, selection_hash, row)
        return row
    row = dict(path=path, url=raw_url(path), selection_sha256=selection_hash,
               expected_bytes=entry["size"], expected_sha256=entry["lfs"]["sha256"], attempts=[])
    log = study / "requests" / stem
    for number in range(1, 4):
        start, finish = log / f"{number}.started.json", log / f"{number}.finished.json"
        if not start.exists():
            publish(start, dict(number=number, started_utc=utc(), url=row["url"], selection_sha256=selection_hash))
            fresh = True
        else:
            fresh = False
        if not finish.exists():
            attempt = dict(started_sha256=digest(start), **read(start))
            try:
                if raw.exists():
                    require(digest(raw) == row["expected_sha256"] and raw.stat().st_size == row["expected_bytes"],
                            "published source differs from pinned bytes")
                    attempt["recovered_published_bytes"] = True
                elif not fresh:
                    raise ValueError("interrupted request without verified published bytes")
                else:
                    with requests.get(row["url"], timeout=(20, 90)) as response:
                        attempt["status_code"] = response.status_code
                        require(response.status_code == 200, f"HTTP {response.status_code} at pinned source")
                        data = response.content
                        require(len(data) == row["expected_bytes"], "download size differs from Hub")
                        require(hashlib.sha256(data).hexdigest() == row["expected_sha256"], "download hash differs from Hub")
                        publish(raw, data)
            except (ValueError, requests.RequestException) as exc:
                attempt["error"] = str(exc) if isinstance(exc, ValueError) else type(exc).__name__
            attempt["finished_utc"] = utc()
            publish(finish, attempt)
        attempt = read(finish)
        require(attempt["started_sha256"] == digest(start) and attempt["selection_sha256"] == selection_hash,
                "request checkpoint changed")
        row["attempts"].append(dict(**attempt, finished_sha256=digest(finish)))
        if "error" not in attempt:
            row.update(local=str(raw.resolve()), bytes=raw.stat().st_size, sha256=digest(raw))
            break
    else:
        row["error"] = "all three pinned source requests failed"
    row["recorded_utc"] = utc()
    verify_receipt(study, entry, selection_hash, row)
    publish(receipt, row)
    return row


def audit_requests(study, row):
    log = study / "requests" / Path(row["path"]).stem
    attempts = row["attempts"]
    require(1 <= len(attempts) <= 3, "request checkpoint count differs")
    for n, attempt in enumerate(attempts, 1):
        start, finish = log / f"{n}.started.json", log / f"{n}.finished.json"
        require(attempt["number"] == n and attempt["started_sha256"] == digest(start)
                and attempt["finished_sha256"] == digest(finish)
                and {k: v for k, v in attempt.items() if k != "finished_sha256"} == read(finish)
                and all(attempt.get(k) == v for k, v in read(start).items())
                and attempt["selection_sha256"] == row["selection_sha256"] and attempt["url"] == row["url"],
                "request checkpoint differs")
    require(all("error" in a for a in attempts[:-1]) and ("error" in attempts[-1]) == ("error" in row),
            "request disposition differs")


def wall_definition(original):
    return dict(wall=original["wall"], optimization_policy=original["optimization_policy"],
                optimization_footprints=[{k: row[k] for k in ("source", "segment", "points")}
                                         for row in original["optimization_footprints"]])


def input_manifest(time, shot_hash, geometry_hash):
    return dict(time_ms=time, shot_sha256=shot_hash, geometry={"sha256": geometry_hash})


def check_staged(root):
    names = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    require(names == set(live.INPUTS), "unexpected files under prospective input root")
    review = read(root / live.REVIEW)
    require(set(review) == {"wall","optimization_policy","optimization_footprints","time_ms","legs","sources"}
            and review["sources"] == {"OMP": {}, "IMP": {}}, "unexpected staged review fields")
    require(all(set(row) == {"source", "segment", "points"} for row in review["optimization_footprints"]),
            "unexpected footprint metadata")
    record = live.input_record(root)
    require(read(root / live.MANIFEST) == input_manifest(record["time_ms"], record["hashes"][live.SHOT],
                                                         record["hashes"][live.GEOMETRY]), "unexpected manifest fields")
    return record


def stage_geometry(root, source, record):
    for name, raw in ((live.SHOT, Path(source["local"])), (live.GEOMETRY, ROOT / live.GEOMETRY)):
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(raw, root / name)
    original = wall_definition(read(ROOT / live.REVIEW))
    publish(root / live.REVIEW, original)
    time = record["time_ms"]
    publish(root / live.MANIFEST, input_manifest(time, source["sha256"], digest(ROOT / live.GEOMETRY)))
    old_root = prior.geometry.ROOT
    try:
        prior.geometry.ROOT = root
        details, review, reason = prior.geometry.geometry_at(time)
    finally:
        prior.geometry.ROOT = old_root
    record.update(details)
    if reason:
        raise ValueError(reason)
    # Only this private temporary tree is updated. Published roots remain immutable.
    (root / live.REVIEW).unlink()
    publish(root / live.REVIEW, review)
    verified = check_staged(root)
    field, lcfs, axis = load_case(root / live.SHOT, time)
    axis_flux = float(field.values([axis])[0])
    xp_flux = float(field.values([details["X_RZ_m"]])[0])
    with np.errstate(divide="raise", invalid="raise"):
        lcfs_flux = (field.values(lcfs) - axis_flux) / (xp_flux - axis_flux)
    require(np.isfinite(lcfs_flux).all(), "nonfinite normalized LCFS flux")
    record.update(inputs=verified["hashes"], tf_current_kA=field.tf_current_ka,
                  grid_shape=list(field.psi.shape), lower_X_flux=xp_flux,
                  lcfs_normalized_flux_quantiles=np.quantile(lcfs_flux, [0, .05, .5, .95, 1]).tolist(),
                  supplied_wall_equal=True, supplied_geometry_equal=True, omp_endpoints_equal=True,
                  engineering_policy_equal=True)


def screen_one(study, rank, source, index, parent_hash):
    name = Path(source["path"]).stem.removeprefix("d3d_shot_")
    folder = study / "cases" / name
    root = folder / "input-root"
    checkpoint = folder / "preflight.json"
    if not checkpoint.exists():
        row = dict(case_id=name, rank=rank, source=source, status="excluded", started_utc=utc(),
                   previous_checkpoint_sha256=parent_hash, root=str(root))
        staging = study / "staging"
        staging.mkdir(exist_ok=True)
        folder.parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=staging) as temporary:
            bundle = Path(temporary) / "case"
            staged = bundle / "input-root"
            try:
                require("error" not in source, source.get("error", ""))
                times, frames = read_frames(source["local"])
                i, time = prior.nearest_time(times)
                row.update(frame_index=i, time_ms=time)
                signature = frame_fingerprint(array_from_frame(frames[i]))
                matches = index.matches(signature)
                row.update(fingerprint=signature, duplicate_matches=matches)
                require(not matches, "identical dimensional or normalized EFIT frame")
                stage_geometry(staged, source, row)
                row.update(status="passed", geometry_passed=True)
            except (ValueError, RuntimeError, FloatingPointError, AssertionError) as exc:
                row["reason"] = f"{type(exc).__name__}: {exc}"
                if staged.exists():
                    shutil.rmtree(staged)
            row["completed_utc"] = utc()
            publish(bundle / "preflight.json", row)
            bundle.rename(folder)
    row = read(checkpoint)
    require(row["case_id"] == name and row["rank"] == rank and row["source"] == source
            and row["previous_checkpoint_sha256"] == parent_hash, "case checkpoint provenance differs")
    if "fingerprint" in row:
        require(row["duplicate_matches"] == index.matches(row["fingerprint"]), "duplicate history differs")
        index.add(row["fingerprint"], dict(source_id=name, frame_index=row["frame_index"]))
    if row["status"] == "passed":
        require(check_staged(root)["hashes"] == row["inputs"], "staged hashes differ")
    else:
        require(row["status"] == "excluded" and not root.exists(), "excluded case unexpectedly staged")
    return row, digest(checkpoint)


def screen(study):
    selection = check_registration(study)
    entries = inventory(study)
    index = archive_index(study)
    parent_hash = digest(study / "prior-frame-index.json")
    selection_hash = digest(study / "selection.json")
    register_hash = digest(study / "registration.json")
    selected = selection["selected"] + sum(selection["reserve_blocks"], [])
    rows, checks = [], []
    for rank, path in enumerate(selected, 1):
        receipt = download(study, entries[path], selection_hash)
        audit_requests(study, receipt)
        row, parent_hash = screen_one(study, rank, receipt, index, parent_hash)
        rows.append(row)
        checks.append(parent_hash)
        print(rank, row["case_id"], row["status"], row.get("time_ms"), row.get("reason", ""), flush=True)
        if sum(r["status"] == "passed" for r in rows) == selection["target_passes"]:
            break
    files = [r["source"] for r in rows]
    downloads_path = study / "downloads.json"
    downloads = dict(revision=REVISION, protocol_sha256=selection["protocol_sha256"],
                     selection_sha256=selection_hash, registration_sha256=register_hash, files=files)
    publish(downloads_path, downloads)
    preflight_path = study / "preflight.json"
    created = read(preflight_path)["created_utc"] if preflight_path.exists() else utc()
    preflight = dict(created_utc=created, schema="nf-prospective-input-preflight.v1", cases=rows,
                     protocol_sha256=selection["protocol_sha256"], selection_sha256=selection_hash,
                     downloads_sha256=digest(downloads_path), registration_sha256=register_hash,
                     prior_sources_sha256=digest(study / "prior-sources.json"),
                     prior_frame_index_sha256=digest(study / "prior-frame-index.json"),
                     case_checkpoint_sha256=checks, driver_sha256=selection["drivers"][DRIVERS[0]],
                     geometry_driver_sha256=selection["drivers"]["verification/nf_equilibrium_times.py"],
                     original_sha256=selection["original_sha256"],
                     unexamined_preregistered=selected[len(rows):])
    publish(preflight_path, preflight)
    cases = live.load_cases(study)
    require([r["case_id"] for r in cases] == [r["case_id"] for r in rows if r["status"] == "passed"],
            "live loader disagrees with preregistered cases")
    return cases


def verify(study):
    """Independently reread all archived/new EFIT bytes and published input roots."""
    selection = check_registration(study)
    index = FrameIndex()
    metadata = read(study / "prior-sources.json")
    archived = read(study / "prior-frame-index.json")
    require(archived["prior_sources_sha256"] == digest(study / "prior-sources.json")
            and archived["selection_sha256"] == digest(study / "selection.json"), "archive linkage differs")
    for path, value in metadata["referenced_records"].items():
        require(digest(path) == value, f"historical metadata changed: {path}")
    count, valid = 0, 0
    for original, stored in zip(metadata["files"], archived["files"], strict=True):
        raw = Path(original["local"])
        require(raw.stat().st_size == original["bytes"] and digest(raw) == original["sha256"] == stored["file_sha256"]
                and original["source_id"] == stored["source_id"], "historical parquet provenance differs")
        times, frames = read_frames(raw)
        require(len(times) == len(stored["frames"]), "historical frame count differs")
        for i, row in enumerate(stored["frames"]):
            require(row["frame_index"] == i and row["time_ms"] == (float(times[i]) if np.isfinite(times[i]) else None),
                    "historical frame/time differs")
            count += 1
            try:
                signature = frame_fingerprint(array_from_frame(frames[i]))
            except ValueError as exc:
                require(row.get("invalid") == str(exc), "historical invalid frame differs")
                continue
            require(all(row[k] == v for k, v in signature.items()), "historical EFIT fingerprint differs")
            index.add(signature, dict(source_id=original["source_id"], frame_index=i))
            valid += 1

    entries = inventory(study)
    report, downloads = read(study / "preflight.json"), read(study / "downloads.json")
    for n, field in (("protocol.ko.md", "protocol_sha256"), ("selection.json", "selection_sha256"),
                     ("downloads.json", "downloads_sha256"), ("registration.json", "registration_sha256"),
                     ("prior-sources.json", "prior_sources_sha256"), ("prior-frame-index.json", "prior_frame_index_sha256")):
        require(report[field] == digest(study / n), f"preflight linkage differs: {n}")
    for n in ("protocol_sha256", "selection_sha256", "registration_sha256"):
        require(downloads[n] == report[n], "download record linkage differs")
    require(downloads["revision"] == REVISION, "download revision differs")
    order = read(study / "ranked-paths.json")
    rows = report["cases"]
    require([r["source"]["path"] for r in rows] == order[:len(rows)] and len(rows) <= 72
            and report["unexamined_preregistered"] == order[len(rows):72], "evaluated prefix differs")
    registration = read(study / "registration.json")
    prior_hash = report["prior_frame_index_sha256"]
    for rank, (r, source, checksum) in enumerate(zip(rows, downloads["files"], report["case_checkpoint_sha256"], strict=True), 1):
        folder = study / "cases" / r["case_id"]
        require(r["source"] == source and r["rank"] == rank
                and read(study / "receipts" / (Path(source["path"]).stem + ".json")) == source
                and r == read(folder / "preflight.json") and digest(folder / "preflight.json") == checksum
                and r["previous_checkpoint_sha256"] == prior_hash, "case checkpoint chain differs")
        prior_hash = checksum
        verify_receipt(study, entries[source["path"]], report["selection_sha256"], source)
        audit_requests(study, source)
        require(selection["protocol_saved_utc"] <= selection["created_utc"] <= registration["registered_utc"]
                <= source["attempts"][0]["started_utc"] <= source["recorded_utc"] <= r["started_utc"]
                <= r["completed_utc"] <= report["created_utc"], "prospective chronology differs")
        if "error" not in source:
            times, frames = read_frames(source["local"])
            try:
                i, time = prior.nearest_time(times)
                signature = frame_fingerprint(array_from_frame(frames[i]))
            except ValueError as exc:
                require(r["status"] == "excluded" and r["reason"] == f"ValueError: {exc}", "invalid frame disposition differs")
            else:
                duplicates = index.matches(signature)
                require(r["frame_index"] == i and r["time_ms"] == time and signature == r["fingerprint"]
                        and duplicates == r["duplicate_matches"], "new EFIT frame or duplicate status differs")
                if duplicates:
                    require(r["status"] == "excluded", "duplicate input accepted")
                index.add(signature, dict(source_id=r["case_id"], frame_index=i))
        if r["status"] == "passed":
            require(check_staged(folder / "input-root")["hashes"] == r["inputs"], "accepted input changed")
        else:
            require(r["status"] == "excluded" and "reason" in r and not (folder / "input-root").exists(),
                    "failed input disposition differs")
    cases = live.load_cases(study)
    require([r["case_id"] for r in cases] == [r["case_id"] for r in rows if r["status"] == "passed"]
            and ((len(cases) == 3 and rows[-1]["status"] == "passed") or (len(cases) < 3 and len(rows) == 72)),
            "first-three stopping rule differs")
    return dict(archived_sources=len(metadata["files"]), archived_frames=count, valid_archived_frames=valid,
                examined=len(rows), cases=[r["case_id"] for r in cases], source_records={n: digest(study / n) for n in SOURCE_RECORDS})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=STUDY)
    parser.add_argument("--phase", required=True, choices=("register", "screen", "verify"))
    args = parser.parse_args()
    study = args.study.resolve()
    require(study.is_relative_to(STUDIES.resolve()) and study != PREVIOUS.resolve(), "use a new D: study root")
    with study_lock(study):
        if args.phase == "register":
            selection = register(study)
            print("registered", len(selection["selected"]), "plus", [len(b) for b in selection["reserve_blocks"]],
                  "population", selection["population"], "sha256", digest(study / "selection.json"))
        elif args.phase == "screen":
            cases = screen(study)
            print("load_cases", [(c["case_id"], c["time_ms"]) for c in cases])
        else:
            print(json.dumps(verify(study), indent=2))


if __name__ == "__main__":
    main()
