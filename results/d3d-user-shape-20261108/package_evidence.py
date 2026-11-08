"""Pack complete search artifacts and replay summaries; verify checksums without extraction."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

HERE = Path(__file__).resolve().parent
CASES = ("e806f47778", "1a3155370c")
ENTRY = ("entry-demo", "entry-332a6d97fc", "entry-51d83022ca", "entry-release")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def pack(study, name, files):
    with ZipFile(HERE / (name + ".zip"), "x", compression=ZIP_DEFLATED, compresslevel=6) as archive:
        members = {}
        for path in sorted(files):
            key = path.relative_to(study).as_posix()
            data = path.read_bytes()
            archive.writestr(key, data)
            members[key] = dict(bytes=len(data), sha256=digest(data))
    return members


def build(study):
    groups = {"entry-checks": [p for run in ENTRY for p in (study/run).rglob("*") if p.is_file()]}
    for run in ENTRY:
        if not (study/run/"result.json").is_file():
            raise ValueError(f"missing entry check: {run}")
    for case in CASES:
        source = study / case
        paths = []
        for family in ("floor-only", "full"):
            path = source/family/"result.json"
            if not path.is_file():
                raise ValueError(f"missing search: {path}")
            paths.extend(p for p in path.parent.rglob("*") if p.is_file())
        replay = json.loads((source/"replay"/"report.json").read_text(encoding="utf-8"))
        if (not replay["complete"] or len(replay["results"]) != 3 or
                not replay["results"]["cubic1280"]["reproduction"]["passed"]):
            raise ValueError(f"incomplete fixed-shape verification: {case}")
        for filename in ("protocol.json", "report.json"):
            paths.append(source/"replay"/filename)
        groups[case] = paths
    manifest = {name+".zip": dict(members=pack(study, name, paths)) for name, paths in groups.items()}
    for name, group in manifest.items():
        group["sha256"] = digest((HERE/name).read_bytes())
    (HERE/"evidence-manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def verify():
    manifest = json.loads((HERE/"evidence-manifest.json").read_text(encoding="utf-8"))
    checked = 0
    for name, group in manifest.items():
        path = HERE/name
        if digest(path.read_bytes()) != group["sha256"]:
            raise ValueError(f"archive hash mismatch: {name}")
        with ZipFile(path) as archive:
            if sorted(archive.namelist()) != sorted(group["members"]):
                raise ValueError(f"archive member mismatch: {name}")
            for member, expected in group["members"].items():
                data = archive.read(member)
                if len(data) != expected["bytes"] or digest(data) != expected["sha256"]:
                    raise ValueError(f"member mismatch: {member}")
                checked += 1
    print(f"Verified {checked} members in {len(manifest)} archives.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-from", type=Path, help="create archives from completed local studies")
    args = parser.parse_args()
    if args.build_from:
        build(args.build_from)
    verify()
