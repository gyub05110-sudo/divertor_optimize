"""Publication fingerprints resolve only explicitly pinned historical copies."""
import hashlib
import importlib
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
OLD_AGENTS = "7d76769a11ebe52acc2da743dd5290d2139634640f0c1f22cc3f0cbab2cb2067"
OLD_BOUNDARY = "1207eb10590fe3a38942bfc90e0864c3f05a4153fd0186d0eeb6aee685d4144d"
OLD_README = "f0a57aba864a3156de9fa5d18ea43a53776a88bfe60bb380bd28caae6aa8b40b"


def api():
    return importlib.import_module("tools.verify_publication")


def write(root, name, payload):
    path = root/name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def historical_file(root, name):
    archived = "docs/publication/legacy-20261029/"+name
    write(root, name, (ROOT/name).read_bytes())
    write(root, archived, (ROOT/archived).read_bytes())
    return root/archived


@pytest.mark.parametrize("name,expected", [("AGENTS.md", OLD_AGENTS),
                                           ("d3d_optimization/boundary_evaluator.py", OLD_BOUNDARY)])
def test_known_historical_copy_matches_original_expected_hash_without_replacing_current(name, expected):
    module = api()
    actual = module.resolve_recorded_file(ROOT, name, expected)
    assert actual == ROOT/"docs/publication/legacy-20261029"/name
    assert module.digest(actual) == expected
    assert module.digest(ROOT/name) != expected


def test_existing_readme_legacy_mapping_remains_hash_bound(tmp_path):
    module = api()
    write(tmp_path, "README.md", b"current entry point")
    write(tmp_path, "docs/README_LEGACY.md", (ROOT/"docs/README_LEGACY.md").read_bytes())
    assert module.resolve_recorded_file(tmp_path, "README.md", OLD_README) == tmp_path/"docs/README_LEGACY.md"
    with pytest.raises(ValueError, match="fingerprint"):
        module.resolve_recorded_file(tmp_path, "README.md", "0"*64)


@pytest.mark.parametrize("name,expected", [("AGENTS.md", OLD_AGENTS),
                                           ("d3d_optimization/boundary_evaluator.py", OLD_BOUNDARY)])
def test_historical_mapping_does_not_hide_an_unknown_current_change(tmp_path, name, expected):
    module = api()
    historical_file(tmp_path, name)
    write(tmp_path, name, b"an additional unregistered current change")
    with pytest.raises(ValueError, match="unregistered current"):
        module.resolve_recorded_file(tmp_path, name, expected)


def test_historical_copy_tampering_or_changed_expected_hash_is_not_accepted(tmp_path):
    module = api()
    archived = historical_file(tmp_path, "d3d_optimization/boundary_evaluator.py")
    archived.write_bytes(b"changed historical bytes")
    with pytest.raises(ValueError, match="fingerprint"):
        module.resolve_recorded_file(tmp_path, "d3d_optimization/boundary_evaluator.py", OLD_BOUNDARY)
    with pytest.raises(ValueError, match="fingerprint"):
        module.resolve_recorded_file(tmp_path, "d3d_optimization/boundary_evaluator.py", module.digest(archived))


def test_unlisted_archive_never_masks_changed_current_production(tmp_path):
    module = api()
    expected = write(tmp_path, "docs/publication/legacy-20261029/d3d_optimization/other.py", b"original")
    write(tmp_path, "d3d_optimization/other.py", b"changed")
    with pytest.raises(ValueError, match="fingerprint"):
        module.resolve_recorded_file(tmp_path, "d3d_optimization/other.py", expected)


def historical_study(root):
    module = api()
    preserved = {"dependency.py": write(root, "dependency.py", b"fixed dependency")}
    write(root, module.PRESERVATION, json.dumps(preserved).encode())
    files = dict(preserved)
    for name in set(module.REQUIRED_INPUTS) | set(module.MODEL_FILES):
        if name == "d3d_optimization/boundary_evaluator.py":
            historical_file(root, name)
            files[name] = OLD_BOUNDARY
        elif name == module.PRESERVATION:
            files[name] = module.digest(root/name)
        else:
            files[name] = write(root, name, ("recorded "+name).encode())
    files["extra/seed.txt"] = write(root, "extra/seed.txt", b"additional recorded input")
    folder = root/"saved-study"
    record = dict(schema="boundary-study-provenance.v1", files=files)
    write(root, "saved-study/provenance.json", json.dumps(record).encode())
    return folder, record


def test_historical_study_checks_schema_required_inputs_and_all_recorded_hashes(tmp_path):
    module = api()
    folder, record = historical_study(tmp_path)
    assert module.verify_historical_study(tmp_path, folder)["files_verified"] == len(record["files"])
    write(tmp_path, "extra/seed.txt", b"extra input changed")
    with pytest.raises(ValueError, match="fingerprint"):
        module.verify_historical_study(tmp_path, folder)


@pytest.mark.parametrize("change", ["schema", "required_input", "preservation_input"])
def test_historical_study_rejects_incomplete_provenance(tmp_path, change):
    module = api()
    folder, record = historical_study(tmp_path)
    if change == "schema":
        record["schema"] = "other.v1"
    else:
        record["files"].pop("dependency.py" if change == "preservation_input" else module.REQUIRED_INPUTS[0])
    (folder/"provenance.json").write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="incomplete boundary provenance"):
        module.verify_historical_study(tmp_path, folder)


def test_historical_study_rejects_replacing_original_hash_with_known_current_hash(tmp_path):
    module = api()
    folder, record = historical_study(tmp_path)
    name = "d3d_optimization/boundary_evaluator.py"
    record["files"][name] = module.HISTORICAL_COPIES[name]["known_current"]
    (folder/"provenance.json").write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="historical expected fingerprint"):
        module.verify_historical_study(tmp_path, folder)


def test_preserved_dependency_change_still_fails(tmp_path):
    module = api()
    folder, _ = historical_study(tmp_path)
    write(tmp_path, "dependency.py", b"changed dependency")
    with pytest.raises(ValueError, match="fingerprint"):
        module.verify_historical_study(tmp_path, folder)
