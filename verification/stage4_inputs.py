"""Freeze, acquire and verify the Stage 4 cohort using only public input metadata and EFIT."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import re
import subprocess

import numpy as np
import pyarrow as pa

from verification import nf_prospective_inputs as shared
from verification import stage3_inputs as stage3


STUDY = shared.STUDIES / "nf-stage4-20261109" / "inputs"
PROTOCOL = shared.ROOT / "results/nf-stage4-20261109/protocol.ko.md"
PROTOCOL_SHA = "506285e6e9d7c36aead7f5a38a096bea23aeda836d43444de76f0e975fdcf47d"
PROTOCOL_COMMIT = "cb9358fb52e38b83e3915708f766a32472aa26b1"
STAGE3_SHA = "35f77695555964d5bf6e6a8da71c94d6442624d25477802bc70f8fcb48709a83"
SALT = "nf-stage4-hotspot-v1|"
COUNT = 24
TIME_MS = stage3.TIME_MS
DRIVERS = ("verification/stage4_inputs.py", "verification/stage3_inputs.py",
           "verification/nf_prospective_inputs.py", "d3d_optimization/broadening_study.py")
require, read, digest, publish = shared.require, shared.read, shared.digest, shared.publish


def draw(entries, sources):
    train = {e["rfilename"]: e for e in entries if e["rfilename"].startswith(shared.TRAIN)
             and shared.SHOT_NAME.fullmatch(e["rfilename"][len(shared.TRAIN):])}
    require(len({e["rfilename"] for e in entries}) == len(entries), "duplicate inventory path")
    paths = {s["dataset_path"] for s in sources if s.get("dataset_path")}
    require(paths <= train.keys(), "unknown historical training path")
    ids = {s["source_id"].removeprefix("d3d_shot_") for s in sources if "source_id" in s}
    ids.update(str(s["published_shot_number"]) for s in sources if "published_shot_number" in s)
    seen = {s["sha256"] for s in sources}
    for path, e in train.items():
        file_id = shared.SHOT_NAME.fullmatch(path[len(shared.TRAIN):]).group(1)
        if path in paths or file_id in ids:
            paths.add(path)
            seen.add(e["lfs"]["sha256"])
    ranked = sorted((hashlib.sha256((SALT + path).encode("utf-8")).hexdigest(), path) for path in train)
    selected, excluded = [], {}
    for key, path in ranked:
        e = train[path]
        sha = e["lfs"]["sha256"]
        require(re.fullmatch("[0-9a-f]{64}", sha) and e["size"] == e["lfs"]["size"] > 0, "invalid LFS metadata")
        if path in paths or sha in seen:
            excluded[path] = "known_path_or_id" if path in paths else "same_file_sha256"
            continue
        seen.add(sha)
        if len(selected) < COUNT:
            file_id = shared.SHOT_NAME.fullmatch(path[len(shared.TRAIN):]).group(1)
            selected.append(dict(rank=len(selected)+1, id=file_id, dataset_path=path, order_sha256=key,
                                 sha256=sha, bytes=e["size"], url=shared.raw_url(path)))
    require(len(selected) == COUNT, "insufficient distinct public files")
    return selected, excluded


def audit_paths(paths, known, roots=shared.SCAN_ROOTS):
    """Confirm all local study inputs, including renamed staged files, match historical bytes."""
    aliases = {}
    for path in sorted(set(paths)):
        sha = digest(path)
        if sha not in known:
            # Repository unit tests create tiny synthetic parquet fixtures. All public
            # train objects at this revision exceed 2 MB. A study file always fails closed.
            relative = next(path.relative_to(root) for root in roots if path.is_relative_to(root))
            test_fixture = any(part.startswith(("test_", "pytest-")) for part in relative.parts[:-1])
            require(test_fixture and path.stat().st_size < 10_000, f"unregistered historical input: {path}")
            continue
        aliases.setdefault(sha, []).append(str(path))
    return {sha: sorted(names) for sha, names in sorted(aliases.items())}


def local_aliases(copies):
    pairs = {(m.group(1), sha) for sha, paths in copies.items() for path in paths
             if (m := shared.SHOT_NAME.fullmatch(Path(path).name))}
    return [dict(source_id=file_id, sha256=sha) for file_id, sha in sorted(pairs)]


def historical_sources():
    reg_path = stage3.STUDY / "registration.json"
    require(digest(reg_path) == STAGE3_SHA, "stage 3 registration differs")
    reg = read(reg_path)
    select_path = stage3.STUDY / "selection.json"
    require(digest(select_path) == reg["selection_sha256"], "stage 3 selection differs")
    select = read(select_path)
    sources = stage3.historical_sources()
    require(sources == select["historical_sources"], "stage 3 historical population differs")
    require(len(reg["cases"]) == len(select["selected"]) == 12, "incomplete stage 3 history")
    for case, row in zip(reg["cases"], select["selected"], strict=True):
        require(all(case[k] == v for k,v in row.items()), "stage 3 source metadata differs")
        stage3.check_raw(Path(case["path"]), row)
        sources.append(dict(source_id=row["id"], dataset_path=row["dataset_path"], local=case["path"],
                            sha256=row["sha256"], bytes=row["bytes"]))
    return sources


def freeze(out):
    require(digest(PROTOCOL) == PROTOCOL_SHA, "stage 4 protocol differs")
    sources = historical_sources()
    audit = out / "local-source-audit.json"
    if not (out / "selection.json").is_file():
        require(not any((out / p).exists() for p in ("raw", "requests", "receipts", "registration.json")),
                "prospective candidate data precedes selection")
    if not audit.is_file():
        scan = subprocess.run(["rg", "--files", "--hidden", "--no-ignore", "-g", "*.parquet",
                               "-g", "!**/.git/**", "-g", "!**/.local_dependencies/**", "-g", "!**/node_modules/**",
                               *map(str, shared.SCAN_ROOTS)], check=True, capture_output=True, text=True)
        aliases = audit_paths((Path(s).resolve() for s in scan.stdout.splitlines()),
                              {s["sha256"] for s in sources})
        require(set(aliases) == {s["sha256"] for s in sources}, "historical raw omitted from disk audit")
        publish(audit, dict(roots=list(map(str, shared.SCAN_ROOTS)), verified_copies=aliases))
    inventory_path = stage3.PREVIOUS / "hub-inventory.json"
    inv = read(inventory_path)
    require(inv["id"] == shared.DATASET and inv["sha"] == shared.REVISION and
            inv["cardData"]["license"] == "cc-by-4.0", "public data provenance differs")
    selected, excluded = draw(inv["siblings"], sources + local_aliases(read(audit)["verified_copies"]))
    files = {}
    for name, src in [("protocol.ko.md", PROTOCOL), ("hub-inventory.json", inventory_path),
                      ("dataset-README.md", stage3.PREVIOUS / "dataset-README.md"),
                      *((f"source-code/{d}", shared.ROOT / d) for d in DRIVERS)]:
        publish(out / name, src.read_bytes())
        files[name] = digest(out / name)
    files[audit.name] = digest(audit)
    selection = dict(schema="nf-stage4-input-selection.v1", protocol_sha256=PROTOCOL_SHA,
                     protocol_commit=PROTOCOL_COMMIT, salt=SALT, time_ms=TIME_MS, dataset=shared.DATASET,
                     revision=shared.REVISION, license=inv["cardData"]["license"], inventory_url=shared.API,
                     dataset_card_url=shared.raw_url("README.md"), files_sha256=files,
                     stage3_registration_sha256=STAGE3_SHA, historical_sources=sources, excluded=excluded,
                     selected=selected)
    selection_path = out / "selection.json"
    seal_path = out / "selection-seal.json"
    publish(selection_path, selection)
    sha = digest(selection_path)
    if not seal_path.is_file():
        require(not any((out / p).exists() for p in ("raw", "requests", "receipts", "registration.json")),
                "prospective candidate data precedes the selection seal")
        publish(seal_path, dict(frozen_utc=shared.utc(), selection_sha256=sha, protocol_sha256=PROTOCOL_SHA))
    seal = read(seal_path)
    require(seal["selection_sha256"] == sha and seal["protocol_sha256"] == PROTOCOL_SHA, "selection seal differs")
    return selection


def assess(path, index):
    outcome = stage3.assess(path, index)
    outcome.setdefault("exact_frame_count", 1 if "frame_index" in outcome else None)
    return outcome


def register_frames(path, index, row):
    result = assess(path, index)
    try:
        times, frames = shared.read_frames(path)
    except (ValueError, TypeError, AttributeError, OverflowError, pa.ArrowException):
        require(result["status"] == "invalid_input", "inconsistent frame read")
        return dict(result, file_frame_count=None, valid_file_frames=None)
    result["exact_frame_count"] = int(np.count_nonzero(times == TIME_MS))
    valid = 0
    for i, t in enumerate(times):
        try:
            fp = shared.frame_fingerprint(shared.array_from_frame(frames[i]))
        except (ValueError, TypeError, AttributeError, OverflowError, pa.ArrowException):
            continue
        index.add(fp, dict(id=row["id"], sha256=row["sha256"], frame_index=i,
                           time_ms=float(t) if np.isfinite(t) else None))
        valid += 1
    return dict(result, file_frame_count=len(times), valid_file_frames=valid)


def case_record(out, row, selection_sha, index):
    seal = read(out / "selection-seal.json")
    require(seal["selection_sha256"] == selection_sha, "selection linkage differs")
    log = out / "requests" / row["id"]
    attempts = []
    for path in sorted(log.glob("*.started.json"), key=lambda p: int(p.name.split(".")[0])):
        start = read(path)
        n = start["number"]
        require(n == len(attempts)+1 and path.name == f"{n}.started.json" and
                start["selection_sha256"] == selection_sha and start["started_utc"] >= seal["frozen_utc"] and
                start["requested_url"].split("?")[0] == row["url"], "request selection or chronology differs")
        finish = log / f"{n}.finished.json"
        attempt = dict(started_sha256=digest(path))
        if finish.is_file():
            done = read(finish)
            require(all(done.get(k) == v for k,v in start.items()) and
                    done["finished_utc"] >= start["started_utc"], "request completion differs")
            attempt.update(finished_sha256=digest(finish), error=done.get("error"))
        else:
            attempt["error"] = "interrupted"
        attempts.append(attempt)
    require(attempts, "no pinned object request")
    receipt_path = out / "receipts" / (row["id"] + ".json")
    result = dict(**row, time_ms=TIME_MS, attempts=attempts)
    raw = out / "raw" / Path(row["dataset_path"]).name
    if not receipt_path.exists():
        require(all(a["error"] for a in attempts) and not raw.exists(), "unresolved pinned download")
        return dict(result, status="acquisition_failed", reason=attempts[-1]["error"], exact_frame_count=None)
    receipt = read(receipt_path)
    n = receipt["success_attempt"]
    require(n == len(attempts) and attempts[-1].get("error") is None and receipt["selection_sha256"] == selection_sha
            and receipt["finished_utc"] == read(log / f"{n}.finished.json")["finished_utc"] and
            all(receipt[k] == row[k] for k in ("id", "rank", "url", "bytes", "sha256")), "download receipt differs")
    stage3.check_raw(raw, row)
    return dict(result, path=str(raw.resolve()), receipt_sha256=digest(receipt_path), **register_frames(raw, index, row))


def register(out, selection, verify):
    sha = digest(out / "selection.json")
    registration_path = out / "registration.json"
    if not verify and not registration_path.is_file():
        with ThreadPoolExecutor(max_workers=3) as pool:
            jobs = [pool.submit(stage3.download, out, r, sha) for r in selection["selected"]]
            for job in jobs:
                try:
                    job.result()
                except RuntimeError as exc:
                    # Retry budget exhausted for this member. Preserve its evidence and its rank.
                    print(str(exc), flush=True)
    records, index = stage3.historical_frames(selection["historical_sources"])
    previous_index = stage3.STUDY / "prior-frame-index.json"
    require(digest(previous_index) == read(stage3.STUDY / "registration.json")["files_sha256"][previous_index.name]
            and records[:-12] == read(previous_index), "recomputed pre-stage-3 frame history differs")
    index_path = out / "prior-frame-index.json"
    publish(index_path, records)
    cases = [case_record(out, r, sha, index) for r in selection["selected"]]
    result = dict(schema="nf-stage4-input-registration.v1", protocol_sha256=PROTOCOL_SHA,
                  selection_sha256=sha, selection_seal_sha256=digest(out / "selection-seal.json"),
                  dataset=shared.DATASET, revision=shared.REVISION, license=selection["license"],
                  files_sha256={**selection["files_sha256"], index_path.name: digest(index_path)},
                  historical_source_objects=len(records), historical_frames=sum(len(s["frames"]) for s in records),
                  valid_historical_frames=sum("invalid" not in f for s in records for f in s["frames"]),
                  cases=cases)
    if verify:
        require(read(registration_path) == result, "registration differs from input-byte recomputation")
    else:
        publish(registration_path, result)
    print("registration", digest(registration_path), [(r["id"], r["status"]) for r in cases], flush=True)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("phase", choices=("freeze", "acquire", "verify"))
    cli.add_argument("--output", type=Path, default=STUDY)
    args = cli.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with shared.study_lock(args.output):
        selection = freeze(args.output)
        print("selection", digest(args.output / "selection.json"), flush=True)
        if args.phase != "freeze":
            register(args.output, selection, verify=args.phase == "verify")


if __name__ == "__main__":
    main()
