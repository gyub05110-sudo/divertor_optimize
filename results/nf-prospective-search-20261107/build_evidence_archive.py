"""Archive the closed experiment, selected input trees and acquisition receipts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


def sha(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def walk_files(root):
    for folder, directories, names in os.walk(root):
        directories[:] = sorted(d for d in directories if d not in {"cache", "posttrial-runtime", "__pycache__"})
        for name in sorted(names):
            if not name.endswith(".lock"):
                yield Path(folder)/name


def add(bundle, name, data):
    entry = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
    entry.external_attr = 0o100644 << 16
    bundle.writestr(entry, data, compress_type=ZIP_DEFLATED, compresslevel=6)


def archive(study, output):
    from d3d_optimization.campaign_state import process_lock
    with process_lock(study/"campaign.lock"), process_lock(study/"posttrial.lock"):
        registration = json.loads((study/"frozen.json").read_text(encoding="utf-8"))
        inputs = Path(registration["source_study"]).resolve()
        if not (study/"numerics/summary.json").is_file() or not (output/"summary.json").is_file():
            raise ValueError("Finish the scores and independent report before packaging")
        if any(output.is_relative_to(root) for root in (study, inputs)):
            raise ValueError("Store the archive outside both source directories")
        out_zip = output/"evidence.zip"
        manifest_file = output/"evidence-manifest.json"
        sum_file = output/"evidence.zip.sha256"
        if any(p.exists() for p in (out_zip, manifest_file, sum_file)):
            raise FileExistsError("Choose a new destination to preserve the previous evidence")
        public_sources = json.loads((inputs/"downloads.json").read_text(encoding="utf-8"))["files"]
        published = {Path(row["local"]).resolve():row for row in public_sources}
        included, external = {}, {}
        with ZipFile(out_zip, "x", allowZip64=True) as bundle:
            for root in (study, inputs):
                for path in walk_files(root):
                    name = (Path("studies")/root.name/path.relative_to(root)).as_posix()
                    data = path.read_bytes()
                    row = {"bytes":len(data), "sha256":sha(data)}
                    if root == inputs and path.resolve() in published:
                        source = published[path.resolve()]
                        if (source["sha256"], source["bytes"]) != (row["sha256"], row["bytes"]):
                            raise ValueError("Public source changed: " + name)
                        external[name] = dict(row, url=source["url"])
                    else:
                        included[name] = row
                        add(bundle, name, data)
            manifest = json_bytes(dict(schema="nf-prospective-archive.v1", registration_sha256=sha((study/"frozen.json").read_bytes()),
                                       members=included, public_raw_files=external,
                                       excluded_directories=["cache", "posttrial-runtime", "__pycache__"],
                                       excluded_file_suffixes=[".lock"]))
            add(bundle, "evidence-manifest.json", manifest)
        with ZipFile(out_zip) as bundle:
            if set(bundle.namelist()) != set(included) | {"evidence-manifest.json"}:
                raise ValueError("Archive member set differs")
            for name, expected in included.items():
                data = bundle.read(name)
                if (len(data), sha(data)) != (expected["bytes"], expected["sha256"]):
                    raise ValueError("Archive member round trip differs: " + name)
            if bundle.read("evidence-manifest.json") != manifest:
                raise ValueError("Archive manifest differs")
        manifest_file.write_bytes(manifest)
        with out_zip.open("rb") as source:
            out_hash = hashlib.file_digest(source, "sha256").hexdigest()
        sum_file.write_text(out_hash+"  evidence.zip\n", encoding="utf-8")
        return dict(members=len(included), externally_downloadable=len(external), bytes=out_zip.stat().st_size, sha256=out_hash)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    study, output = args.study.resolve(), args.output.resolve()
    if study.drive.lower() != "d:" or output.drive.lower() != "d:":
        raise ValueError("Keep the study and the archive on D")
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    print(json.dumps(archive(study, output)))
