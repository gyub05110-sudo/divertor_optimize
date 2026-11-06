"""Verify archived publication fingerprints without changing current production code."""
import hashlib
import json
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]
PRESERVATION = "docs/superpowers/plans/2026-10-28-boundary-preservation.json"
# Requirements of the archived schema; no imports or patches to live evaluators.
REQUIRED_INPUTS = (
    "d3d_optimization/__init__.py", "d3d_optimization/case.py",
    "d3d_optimization/frozen_field.py", "d3d_optimization/frozen_mapping.py",
    "d3d_optimization/multisurface.py", "d3d_optimization/shape.py",
    "d3d_optimization/independent_sources.py", "d3d_optimization/continuous_boundary.py",
    "optimize_independent_sources.py", "baseline_review/baseline_review.json",
    "input/D3D/d3d/d3d_shot_1a3155370c.parquet",
    "input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json",
    "input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json",
)
MODEL_FILES = (
    "d3d_optimization/boundary_transport.py", "d3d_optimization/boundary_evaluator.py",
    "optimize_boundary_transport.py", "verification/validate_boundary_transport.py", PRESERVATION,
)
# Each relocation is bound to the original recorded hash, never just a filename.
# Current hashes identify intentional changes in 7a889e8; original bytes are from
# df239e3 and match the unchanged publication/preservation/provenance manifests.
HISTORICAL_COPIES = {
    "README.md": {
        "expected": "f0a57aba864a3156de9fa5d18ea43a53776a88bfe60bb380bd28caae6aa8b40b",
        "path": "docs/README_LEGACY.md",
    },
    "AGENTS.md": {
        "expected": "7d76769a11ebe52acc2da743dd5290d2139634640f0c1f22cc3f0cbab2cb2067",
        "path": "docs/publication/legacy-20261029/AGENTS.md",
        "known_current": "d17a8bd0a47db4b96af2bbd81af59aeaf20624bfd797b8340be93e4b80343571",
    },
    "d3d_optimization/boundary_evaluator.py": {
        "expected": "1207eb10590fe3a38942bfc90e0864c3f05a4153fd0186d0eeb6aee685d4144d",
        "path": "docs/publication/legacy-20261029/d3d_optimization/boundary_evaluator.py",
        "known_current": "33ae04f9a95245cd74525c719aa35b2ee6f96f45cd6811a63b62973c8dc6ced0",
    },
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def recorded_path(root, name):
    require(isinstance(name, str) and name and "\\" not in name and ":" not in name,
            "invalid recorded relative path")
    relative = PurePosixPath(name)
    require(not relative.is_absolute() and ".." not in relative.parts, "invalid recorded relative path")
    path = (Path(root)/name).resolve()
    require(path.is_relative_to(Path(root).resolve()), "recorded path leaves publication root")
    return path


def resolve_recorded_file(root, name, expected):
    """Return only the current file or one explicit, hash-bound historical copy."""
    require(isinstance(expected, str) and len(expected) == 64 and
            all(c in "0123456789abcdef" for c in expected), "invalid recorded SHA-256 fingerprint")
    path = recorded_path(root, name)
    historical = HISTORICAL_COPIES.get(name)
    if historical is not None:
        require(expected == historical["expected"], "historical expected fingerprint differs: "+name)
        if "known_current" in historical:
            require(path.is_file() and digest(path) in (expected, historical["known_current"]),
                    "unregistered current change behind historical copy: "+name)
        path = recorded_path(root, historical["path"])
    require(path.is_file() and digest(path) == expected, "recorded fingerprint differs: "+name)
    return path


def verify_files(root, files):
    require(isinstance(files, dict) and files, "missing recorded file fingerprints")
    return {name: resolve_recorded_file(root, name, expected) for name, expected in files.items()}


def verify_historical_study(root, folder):
    """Validate the historical schema, required inputs and every recorded hash."""
    preserved = load(recorded_path(root, PRESERVATION))
    verify_files(root, preserved)
    record = load(Path(folder)/"provenance.json")
    files = record.get("files")
    required = set(REQUIRED_INPUTS) | set(MODEL_FILES) | set(preserved)
    require(record.get("schema") == "boundary-study-provenance.v1" and
            isinstance(files, dict) and required.issubset(files), "incomplete boundary provenance")
    resolved = verify_files(root, files)
    return dict(files_verified=len(resolved), required_files=len(required))


def verify_publication(root=ROOT):
    root = Path(root).resolve()
    copied = load(root/"docs/publication/source_copy_hashes.json")
    resolved = verify_files(root, copied["files"])
    preserved = load(root/"results/regional-20261028/preservation.json")
    verify_files(root, preserved)
    for source in ("IMP", "OMP"):
        verify_historical_study(root, root/"results/regional-20261028"/source)
        folder = root/"results/regional-20261028/qualified"/source
        verify_historical_study(root, folder)
        require(load(folder/"validation.json")["verdict"]["improvement_confirmed"] is True,
                "saved qualified verdict did not confirm improvement: "+source)
        require(load(folder/"regional_validation.json")["passed"] is True,
                "saved regional verdict did not pass: "+source)
    relocated = {name: dict(path=path.relative_to(root).as_posix(), sha256=copied["files"][name])
                 for name, path in resolved.items() if path != root/name}
    return dict(copied_project_files_verified=len(copied["files"]),
                preserved_python_files_verified=len(preserved), study_provenance_verified=4,
                saved_qualified_verdicts_passed=2, historical_resolutions=relocated,
                scope="historical fingerprints and saved verdicts; not current-model identity or a fresh numerical refinement")


if __name__ == "__main__":
    print(json.dumps(verify_publication()))
