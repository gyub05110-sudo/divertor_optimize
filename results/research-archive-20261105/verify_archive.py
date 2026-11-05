"""Verify or restore the November 5 research ZIPs using Python's standard library."""
import argparse
import hashlib
import json
import os
import shutil
import tempfile
import zipfile
from pathlib import Path, PurePosixPath


def check_stream(stream, expected, output=None):
    digest = hashlib.sha256()
    total = 0
    while chunk := stream.read(1024 * 1024):
        digest.update(chunk)
        total += len(chunk)
        if output is not None:
            output.write(chunk)
    if total != expected["bytes"] or digest.hexdigest() != expected["sha256"]:
        raise ValueError(f"Size or SHA-256 mismatch for {expected.get('name', expected.get('path'))}")


def check_path(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or ".." in path.parts or ":" in name or "\\" in name:
        raise ValueError(f"Unsafe archive path: {name}")
    return path


def verify(assets_dir, inventory, stage):
    seen = set()
    total = 0
    for asset in inventory["assets"]:
        name = asset["name"]
        asset_path = check_path(name)
        if len(asset_path.parts) != 1:
            raise ValueError(f"Invalid asset name: {name}")
        records = inventory["members"][name]
        rows = {row["path"]: row for row in records}
        if len(rows) != len(records) or len(rows) != asset["files"] or not seen.isdisjoint(rows):
            raise ValueError(f"Duplicate or missing inventory member in {name}")

        with (assets_dir / name).open("rb") as file:
            check_stream(file, asset)
            file.seek(0)
            with zipfile.ZipFile(file) as archive:
                if len(archive.infolist()) != len(rows) or set(archive.namelist()) != rows.keys():
                    raise ValueError(f"ZIP member list differs from inventory for {name}")
                for member in archive.infolist():
                    relative = check_path(member.filename)
                    record = rows[member.filename]
                    if member.file_size != record["bytes"]:
                        raise ValueError(f"Member size mismatch: {member.filename}")
                    with archive.open(member) as source:
                        if stage is None:
                            check_stream(source, record)
                        else:
                            dest = stage.joinpath(*relative.parts)
                            if not dest.resolve().is_relative_to(stage):
                                raise ValueError(f"Member outside restoration directory: {relative}")
                            dest.parent.mkdir(parents=True, exist_ok=True)
                            with dest.open("xb") as output:
                                check_stream(source, record, output)
                            if os.name == "posix":
                                dest.chmod(0o755 if member.filename.endswith(".sh") else 0o644)
            seen.update(rows)
            total += sum(row["bytes"] for row in rows.values())
        print(f'{name}: {len(rows)} files verified', flush=True)
    return len(seen), total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assets", type=Path, help="Directory containing all 13 ZIPs")
    parser.add_argument("--inventory", type=Path, default=Path(__file__).with_name("inventory.json"))
    parser.add_argument("--restore", type=Path, help="New output directory, which must not exist yet")
    args = parser.parse_args()
    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    if inventory["schema"] != 1 or inventory["tag"] != "research-archive-20261105":
        parser.error("unexpected inventory schema or tag")
    stage = None
    try:
        if args.restore is not None:
            dest = args.restore.resolve()
            if dest.exists():
                parser.error(f"restore path already exists: {dest}")
            dest.parent.mkdir(parents=True, exist_ok=True)
            stage = Path(tempfile.mkdtemp(prefix=".restore-research-", dir=dest.parent)).resolve()
        count, total = verify(args.assets, inventory, stage)
        if stage is not None:
            stage.rename(dest)
            stage = None
            print(f"Restored to {dest}")
        print(f"Verified {count:,} files, {total:,} uncompressed bytes.")
    finally:
        if stage is not None:
            shutil.rmtree(stage)


if __name__ == "__main__":
    main()
