"""Freeze and verify the Stage 3 public DIII-D input cohort without evaluating heat."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path, PurePosixPath
import time

import numpy as np
import pyarrow as pa
import requests

from verification import nf_prospective_inputs as shared


STUDY = shared.STUDIES / "nf-stage3-20261108" / "inputs"
PREVIOUS = shared.STUDY
PROTOCOL = shared.ROOT / "results/nf-stage3-20261108/protocol.ko.md"
SALT = "nf-stage3-v1|"
TIME_MS = 1400.0
SOURCE_HASHES = {
    "hub-inventory.json": "d868f9da1c59d0ca54050c8ba180c10cf8351422dc9f7e82c39395de51757eca",
    "prior-sources.json": "0c68ae9069ea5d2204bd6c825d00c39b6d3264201655aaccc35e7676a4ebc997",
    "prior-frame-index.json": "4fd99847de03645bac4894bda5f149fda16bb774e3b2344f4c15ecdd897f8474",
    "downloads.json": "e0d4bb8682be42d07f7f7280d40e68b70aa25fb1f28bf8c7223a749910c57eec",
    "dataset-README.md": "a54d197b1ed20f77843a42222ae1ce2ba17ab34cc22a12f5e1f65d4c11b639d7",
}
AUDIT = shared.PREVIOUS / "source-audit.json"
AUDIT_SHA = "5f08ad5012fdf4004ea2c391799b84e865cfa4b3de4c1e0f853672d7564863aa"
PROTOCOL_SHA = "9dc07e514aa99c744b72fd26bacfb5d081c3a13302e28c52893fe9216e6263b5"
require, read, digest, publish = shared.require, shared.read, shared.digest, shared.publish


def draw(entries, prior_sources, count):
    train = {e["rfilename"]: e for e in entries if e["rfilename"].startswith(shared.TRAIN)
             and shared.SHOT_NAME.fullmatch(e["rfilename"][len(shared.TRAIN):])}
    require(len({e["rfilename"] for e in entries}) == len(entries), "duplicate inventory path")
    paths = {s["dataset_path"] for s in prior_sources if s.get("dataset_path")}
    require(paths <= train.keys(), "unknown historical training path")
    seen = {s["sha256"] for s in prior_sources}
    seen.update(train[p]["lfs"]["sha256"] for p in paths)
    seen_ids = {s["source_id"].removeprefix("d3d_shot_") for s in prior_sources if "source_id" in s}
    ranked = sorted((hashlib.sha256((SALT + p).encode("utf-8")).hexdigest(), p) for p in train)
    selected, excluded = [], {}
    for key, path in ranked:
        entry = train[path]
        sha = entry["lfs"]["sha256"]
        shot = shared.SHOT_NAME.fullmatch(PurePosixPath(path).name).group(1)
        require(len(sha) == 64 and entry["size"] == entry["lfs"]["size"] > 0, "invalid LFS metadata")
        if path in paths or shot in seen_ids or sha in seen:
            excluded[path] = "known_path_or_id" if path in paths or shot in seen_ids else "same_file_sha256"
            continue
        seen.add(sha)
        seen_ids.add(shot)
        if len(selected) < count:
            selected.append(dict(rank=len(selected) + 1, id=shot, dataset_path=path, order_sha256=key,
                                 url=shared.raw_url(path), sha256=sha, bytes=entry["size"]))
    require(len(selected) == count, "insufficient unique public objects")
    return selected, excluded


def historical_sources():
    for name, sha in SOURCE_HASHES.items():
        require(digest(PREVIOUS / name) == sha, f"historical record differs: {name}")
    require(digest(AUDIT) == AUDIT_SHA, "source audit changed")
    sources = list(read(PREVIOUS / "prior-sources.json")["files"])
    for r in read(PREVIOUS / "downloads.json")["files"]:
        sources.append(dict(source_id=Path(r["path"]).stem, dataset_path=r["path"], local=r["local"],
                            sha256=r["sha256"], bytes=r["bytes"]))
    # The public, diverted-only edition has different bytes from the bundled baseline.
    public = read(AUDIT)["original_hf_revision_identity"]
    path = sources[0]["dataset_path"]
    require(public["url"] == shared.raw_url(path), "public baseline alias path differs")
    raw = Path(public["downloaded_local"])
    sources.append(dict(source_id=Path(path).stem, dataset_path=path, local=str(raw),
                        sha256=public["downloaded_sha256"], bytes=raw.stat().st_size))
    for source in sources:
        check_raw(Path(source["local"]), source)
    return sources


def check_raw(raw, expected):
    require(raw.is_file() and raw.stat().st_size == expected["bytes"], f"file size mismatch: {raw}")
    require(digest(raw) == expected["sha256"], f"file checksum mismatch: {raw}")


def freeze(out):
    sources = historical_sources()
    require(digest(PROTOCOL) == PROTOCOL_SHA, "frozen protocol changed")
    inventory = read(PREVIOUS / "hub-inventory.json")
    require(inventory["id"] == shared.DATASET and inventory["sha"] == shared.REVISION and
            inventory["cardData"]["license"] == "cc-by-4.0", "dataset provenance differs")
    selected, excluded = draw(inventory["siblings"], sources, 12)
    checks = {str(PREVIOUS / n): v for n,v in SOURCE_HASHES.items()}
    checks[str(AUDIT)] = AUDIT_SHA
    selection = dict(schema="nf-stage3-input-selection.v1", protocol_sha256=PROTOCOL_SHA,
                     protocol_commit="b75c2a9", dataset=shared.DATASET, revision=shared.REVISION,
                     license=inventory["cardData"]["license"], inventory_url=shared.API,
                     dataset_card_url=shared.raw_url("README.md"), salt=SALT, time_ms=TIME_MS,
                     historical_records=checks, historical_sources=sources, excluded=excluded, selected=selected)
    if not (out / "selection.json").exists():
        require(not any((out / name).exists() for name in ("raw", "requests", "registration.json")),
                "candidate records exist before selection freeze")
    publish(out / "selection.json", selection)
    seal = out / "selection-seal.json"
    sha = digest(out / "selection.json")
    if not seal.exists():
        publish(seal, dict(frozen_utc=shared.utc(), selection_sha256=sha, protocol_sha256=PROTOCOL_SHA))
    require(read(seal)["selection_sha256"] == sha, "selection seal differs")
    return selection


def download(out, selected, selection_sha):
    raw = out / "raw" / PurePosixPath(selected["dataset_path"]).name
    receipt = out / "receipts" / (selected["id"] + ".json")
    log = out / "requests" / selected["id"]
    if receipt.is_file():
        item = read(receipt)
        check_raw(raw, selected)
        require(item["selection_sha256"] == selection_sha and item["sha256"] == selected["sha256"],
                "download receipt provenance differs")
        return
    attempts = len(list(log.glob("*.started.json")))
    for i in range(1, 5):
        number = attempts + i
        url = selected["url"]
        if number > 1:
            # A fresh resolve prevents a web cache from returning an expired CDN redirect.
            url += f"?download=true&stage3_retry={number}"
        start = dict(started_utc=shared.utc(), selection_sha256=selection_sha, requested_url=url, number=number)
        publish(log / f"{number}.started.json", start)
        result = dict(start, finished_utc=None)
        try:
            if not raw.exists():
                with requests.get(url, timeout=(20, 120)) as response:
                    result["http_status"] = response.status_code
                    require(response.status_code == 200, f"HTTP {response.status_code}")
                    data = response.content
                    require(len(data) == selected["bytes"] and hashlib.sha256(data).hexdigest() == selected["sha256"],
                            "response differs from pinned file size or checksum")
                    publish(raw, data)
            check_raw(raw, selected)
        except (requests.RequestException, ValueError) as exc:
            result["error"] = type(exc).__name__ if isinstance(exc, requests.RequestException) else str(exc)
        result["finished_utc"] = shared.utc()
        publish(log / f"{number}.finished.json", result)
        if "error" not in result:
            publish(receipt, dict(rank=selected["rank"], id=selected["id"], url=selected["url"],
                                  selection_sha256=selection_sha, sha256=selected["sha256"], bytes=selected["bytes"],
                                  success_attempt=number, finished_utc=result["finished_utc"]))
            print("verified raw", selected["rank"], selected["id"], flush=True)
            return
        time.sleep(min(8, i * 2))
    raise RuntimeError(f"source pending for {selected['id']}; retry the same frozen source")


def historical_frames(sources):
    records, index = [], shared.FrameIndex()
    for source in sources:
        times, frames = shared.read_frames(source["local"])
        rows = []
        for i, timestamp in enumerate(times):
            row = dict(frame_index=i, time_ms=float(timestamp) if np.isfinite(timestamp) else None)
            try:
                row.update(shared.frame_fingerprint(shared.array_from_frame(frames[i])))
            except ValueError as exc:
                row["invalid"] = str(exc)
            else:
                index.add(row, dict(id=source["source_id"], sha256=source["sha256"],
                                    frame_index=i, time_ms=row["time_ms"]))
            rows.append(row)
        records.append(dict(id=source["source_id"], sha256=source["sha256"], frames=rows))
    return records, index


def assess(raw, index):
    try:
        times, frames = shared.read_frames(raw)
        indices = np.flatnonzero(times == TIME_MS)
        if len(indices) != 1:
            return dict(status="ambiguous_time" if len(indices) else "missing_time",
                        reason=f"found {len(indices)} exact frames at {TIME_MS:g} ms", exact_frame_count=len(indices))
        i = int(indices[0])
        fingerprint = shared.frame_fingerprint(shared.array_from_frame(frames[i]))
    except (ValueError, TypeError, AttributeError, OverflowError, pa.ArrowException) as exc:
        return dict(status="invalid_input", reason=f"{type(exc).__name__}: {exc}")
    matches = index.matches(fingerprint)
    outcome = dict(frame_index=i, fingerprint=fingerprint, duplicate_matches=matches,
                   status="duplicate_frame" if matches else "ready")
    if matches:
        outcome["reason"] = "identical dimensional or normalized EFIT frame"
    return outcome


def register(out, selection, verify=False):
    selection_sha = digest(out / "selection.json")
    if not verify:
        with ThreadPoolExecutor(max_workers=3) as pool:
            futures = [pool.submit(download, out, r, selection_sha) for r in selection["selected"]]
            for future in futures:
                future.result()
    records, index = historical_frames(selection["historical_sources"])
    files = {}
    for name, content in (("prior-frame-index.json", records),
                          ("protocol.ko.md", PROTOCOL.read_bytes()),
                          ("dataset-README.md", (PREVIOUS / "dataset-README.md").read_bytes()),
                          ("hub-inventory.json", (PREVIOUS / "hub-inventory.json").read_bytes())):
        publish(out / name, content)
        files[name] = digest(out / name)
    drivers = ("verification/stage3_inputs.py", "verification/nf_prospective_inputs.py",
               "d3d_optimization/broadening_study.py")
    for name in drivers:
        publish(out / "source-code" / name, (shared.ROOT / name).read_bytes())
        files["source-code/" + name] = digest(out / "source-code" / name)
    cases = []
    for row in selection["selected"]:
        raw = out / "raw" / PurePosixPath(row["dataset_path"]).name
        check_raw(raw, row)
        receipt_path = out / "receipts" / (row["id"] + ".json")
        receipt = read(receipt_path)
        require(all(receipt[k] == row[k] for k in ("id", "rank", "sha256", "url", "bytes")) and
                receipt["selection_sha256"] == selection_sha, "download receipt differs")
        attempt = read(out / "requests" / row["id"] / f"{receipt['success_attempt']}.finished.json")
        require("error" not in attempt and attempt["selection_sha256"] == selection_sha and
                attempt["finished_utc"] == receipt["finished_utc"], "successful request provenance differs")
        case = dict(**row, path=str(raw.resolve()), time_ms=TIME_MS,
                    receipt_sha256=digest(receipt_path), **assess(raw, index))
        if "fingerprint" in case:
            index.add(case["fingerprint"], dict(id=row["id"], sha256=row["sha256"],
                                                frame_index=case["frame_index"], time_ms=TIME_MS))
        cases.append(case)
    result = dict(schema="nf-stage3-input-registration.v1", protocol_sha256=PROTOCOL_SHA,
                  selection_sha256=selection_sha, selection_seal_sha256=digest(out / "selection-seal.json"),
                  dataset=shared.DATASET, revision=shared.REVISION, license=selection["license"],
                  historical_source_objects=len(records), historical_frames=sum(len(s["frames"]) for s in records),
                  valid_historical_frames=sum("invalid" not in f for s in records for f in s["frames"]),
                  files_sha256=files, cases=cases)
    if verify:
        require(read(out / "registration.json") == result, "archived registration differs from recomputation")
    else:
        publish(out / "registration.json", result)
    print("registration", digest(out / "registration.json"), [(r["id"], r["status"]) for r in cases], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("freeze", "acquire", "verify"))
    parser.add_argument("--output", type=Path, default=STUDY)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with shared.study_lock(args.output):
        selection = freeze(args.output)
        print("selection", digest(args.output / "selection.json"), flush=True)
        if args.phase != "freeze":
            register(args.output, selection, verify=args.phase == "verify")


if __name__ == "__main__":
    main()
