"""Replay matched OMP search winners and fixed geometric controls on common physical inputs."""
import argparse
import json
from math import isclose, isfinite
from pathlib import Path
from time import perf_counter

import numpy as np

from optimize_divertor import FINE, serial, software
from d3d_optimization.campaign_geometry import Candidate, evaluator
from d3d_optimization.campaign_state import public_save
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.user_input import GEOMETRY, POLICY, load_user_input
from d3d_optimization.user_output import export_measurement, file_hash, segment_regions
from d3d_optimization.user_physics import evaluate_user_candidate, prepare_user, shelf_limits


LEVELS = dict(cubic1280=FINE, quintic1280=("quintic", *FINE[1:]),
              cubic640=("cubic", 640, .00025, 321, .0005))
DEFAULT_LEVELS = ("cubic1280", "quintic1280")
INPUT_KEYS = ("equilibrium_sha256", "time_ms", "geometry_sha256", "design_policy_sha256",
              "design_reference_sha256", "wall_sha256", "omp_endpoints_sha256")
PHYSICS_MODULES = {
    "user_input.py", "user_physics.py", "campaign_geometry.py", "tangent_geometry.py", "tangent_direct.py",
    "tangent_evaluator.py", "case.py", "frozen_field.py", "frozen_mapping.py", "shape.py", "multisurface.py",
    "continuous_boundary.py", "independent_sources.py", "reference_transport.py", "broadening.py",
    "boundary_evaluator.py", "boundary_transport.py", "accelerated_transport.py", "gpu_footprint.py",
}


def _check(ok, reason):
    if not ok:
        raise ValueError(reason)


def _manifest(path, family):
    path = Path(path).resolve()
    if path.is_dir():
        path = path/"result.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    args = data["arguments"]
    _check(data["status"] == "qualified" and args["family"] == family,
           f"expected qualified {family} search at {path}")
    _check(data["search"].get("termination", "deadline") == "deadline", "search ended on an error")
    _check(data["baseline"]["reference_only"] and data["baseline"]["physics_measurable"]
           and data["selected"]["candidate_eligible"], "missing reference or eligible replacement")
    _check(tuple(data["levels"]["fine"]) == FINE, "input search used a different canonical fine level")
    winner = Candidate.from_record(data["candidate"])
    _check(winner.identity == data["candidate"]["id"], "selected candidate identity mismatch")
    stream = path.parent/"events.jsonl"
    events = [json.loads(line) for line in stream.read_text(encoding="utf-8").splitlines() if line.strip()]
    admitted = []
    for row in events:
        if row["stage"] != "fine" or not row["within_budget"] or not row["eligible"]:
            continue
        _check(row["status"] == "completed" and 0 <= row["start_s"] <= row["end_s"] <= args["budget_seconds"],
               "fine event admission or deadline mismatch")
        _check(isfinite(row["peak_q"]) and row["peak_q"] > 0., "invalid admitted fine peak")
        candidate = Candidate.from_record(row["candidate"])
        _check(row["id"] == candidate.identity and candidate.source == "OMP"
               and candidate.radius_mm == args["radius_mm"], "fine evidence geometry mismatch")
        _check(family != "floor-only" or not any(candidate.vector[:12]), "nonfloor displacement in floor arm")
        admitted.append((row, candidate))
    _check(bool(admitted), "no eligible on-time fine evidence")
    best, geometry = min(admitted, key=lambda pair: pair[0]["peak_q"])
    _check(geometry == winner, "returned wall differs from the strict eligible fine minimum")
    expected_peak = best["peak_q"]*data["reference_SI"]["power_W"]/1e6
    _check(isclose(expected_peak, data["selected"]["peak_MW_m2"], rel_tol=2e-14),
           "selected export and fine event peaks differ")
    for key in ("baseline", "selected"):
        record = data[key]
        for name, checksum in record.get("sha256", {}).items():
            _check(Path(name).name == name and record["directory"] == key, "invalid artifact location")
            _check(file_hash(path.parent/key/name) == checksum, f"changed {key} artifact {name}")
    source = dict(result=str(path), result_sha256=file_hash(path), events_sha256=file_hash(stream),
                  baseline={k: data["baseline"][k] for k in ("peak_MW_m2", "A2_m2")},
                  selected={k: data["selected"][k] for k in ("peak_MW_m2", "A2_m2")})
    return data, source, events, admitted, winner


def prepare_protocol(full, floor, levels=DEFAULT_LEVELS, backend="cpu"):
    """Seal geometric attribution and control selection before replaying either winner."""
    _check(backend in ("cpu", "gpu"), "unsupported replay backend")
    _check(levels and len(set(levels)) == len(levels) and all(k in LEVELS for k in levels),
           "select distinct supported replay levels")
    a, source_a, events_a, fine_a, winner_a = _manifest(full, "full")
    b, source_b, events_b, fine_b, winner_b = _manifest(floor, "floor-only")
    for key in INPUT_KEYS:
        _check(a["inputs"][key] == b["inputs"][key], f"unmatched input {key}")
    for key in ("reference_SI", "levels"):
        _check(a[key] == b[key], f"unmatched {key}")
    args = {k: a["arguments"][k] for k in ("seed", "budget_seconds", "radius_mm", "backend")}
    for key, value in args.items():
        _check(b["arguments"][key] == value, f"unmatched search setting {key}")
    _check(a["software"]["code_sha256"] == b["software"]["code_sha256"], "unmatched search code hashes")
    _check(source_a["baseline"] == source_b["baseline"], "unmatched fine reference metrics")

    roles, candidates = {}, {}
    def add(role, candidate, reason):
        roles[role] = dict(id=candidate.identity if candidate else None, selection=reason)
        if candidate is not None:
            candidates[candidate.identity] = candidate.record()

    def ablate(candidate):
        return Candidate("OMP", candidate.radius_mm, (0.,)*12 + candidate.vector[12:])

    add("full_winner", winner_a, "Minimum completed eligible full-arm fine peak")
    add("floor_winner", winner_b, "Minimum completed eligible floor-only fine peak")
    add("full_winner_floor_ablation", ablate(winner_a), "Full winner with the first twelve coefficients zeroed")
    if not any(winner_a.vector[:12]):
        nonfloor = [(row, c) for row, c in fine_a if any(c.vector[:12])]
        alternate = min(nonfloor, key=lambda pair: pair[0]["peak_q"])[1] if nonfloor else None
        add("full_nonfloor_finalist", alternate, "Lowest admitted full-arm fine peak with nonzero nonfloor coefficients")
        add("full_nonfloor_finalist_floor_ablation", ablate(alternate) if alternate else None,
            "Nonfloor finalist with the first twelve coefficients zeroed")
    origins = []
    for events in (events_a, events_b):
        origin = next((e for e in events if e["stage"] == "proposal" and e.get("proposal_kind") == "projected_origin"
                       and e["within_budget"] and e.get("candidate")), None)
        c = Candidate.from_record(origin["candidate"]) if origin else None
        _check(c is None or c.identity == origin["id"], "projected origin identity mismatch")
        origins.append(c)
    _check(origins[0] == origins[1], "unmatched projected cold origins")
    add("projected_origin", origins[0], "Common projected cold origin from both proposal journals")
    restart = next(((c, arm) for arm, rows in (("full", fine_a), ("floor-only", fine_b))
                    for event, c in rows if event["proposal_kind"] == "projected_restart"), None)
    add("restart_control", restart[0] if restart else None,
        ("First admitted eligible independent projected restart in the "+restart[1]+" arm")
        if restart else "No independent restart received an eligible on-time fine score")
    zero = Candidate("OMP", args["radius_mm"], (0.,)*16)
    add("g1_zero", zero, "All free coefficients zero with mandatory G1 floor endpoint tangent")
    for mm in (-40, -20, 0, 20, 40):
        c = Candidate("OMP", args["radius_mm"], (0.,)*12 + (mm/1000.,)*4)
        add(f"one_amplitude_{mm:+d}mm", c, f"Four equal floor coefficients of {mm} mm with all constraints active")
    return dict(schema="user-omp-shape-replay.v1", sources=dict(full=source_a, floor=source_b),
                inputs=a["inputs"], reference_SI=a["reference_SI"], matched_search=args,
                search_levels=a["levels"], search_code_sha256=a["software"]["code_sha256"],
                backend=backend, levels={name: LEVELS[name] for name in levels},
                roles=roles, candidates=candidates)


def _measurement(folder, e, raw):
    row = export_measurement(folder, e, raw)
    region = segment_regions(raw)
    hits = raw["hits"].segment
    regions = np.full(len(hits), -1, int)
    regions[hits >= 0] = region[hits[hits >= 0]]
    weights = e.source.weights
    for k, name in enumerate(raw["region_names"]):
        keep = (regions == k) & (weights > 0)
        for key, array in (("hit_weighted_physical_width_m", raw["deposition"]["sigma_physical_m"]),
                           ("hit_weighted_parallel_length_m", raw["thermal"]["length_m"])):
            row["surfaces"][name][key] = float(np.average(array[keep], weights=weights[keep])) if keep.any() else None
    row["separatrix_witness"] = raw.get("boundary_witness")
    return row


def _reproduction(protocol, grid):
    comparisons = {}
    for arm, role in (("full", "full_winner"), ("floor", "floor_winner")):
        rows = dict(baseline=grid["baseline"], selected=grid["candidates"][protocol["roles"][role]["id"]])
        for label, row in rows.items():
            expected = protocol["sources"][arm][label]
            comparison = {}
            for key in ("peak_MW_m2", "A2_m2"):
                value = row.get(key)
                ratio = value / expected[key] - 1. if value is not None else None
                comparison[key+"_relative_difference"] = ratio
            valid = row.get("candidate_eligible") if label == "selected" else row["physics_measurable"]
            comparison["passed"] = bool(valid and all(v is not None and abs(v) <= 1e-9
                                                      for v in comparison.values()))
            comparisons[arm+"_"+label] = comparison
    return dict(passed=all(v["passed"] for v in comparisons.values()), tolerance_relative=1e-9,
                comparisons=comparisons)


def _replay_input_paths(source, *, equilibrium=None, geometry=None, design_policy=None):
    paths = {}
    for key, override, default in (("equilibrium", equilibrium, None), ("geometry", geometry, GEOMETRY),
                                   ("design_policy", design_policy, POLICY)):
        path = Path(override) if override is not None else Path(source[key+"_path"])
        if override is None and not path.is_file() and default is not None:
            path = default
        if not path.is_file():
            raise FileNotFoundError(f"missing replay input {path}; supply --{key.replace('_', '-')} with a matching copy")
        _check(file_hash(path) == source[key+"_sha256"], f"{key} checksum differs from the archived input at {path}")
        paths[key] = str(path.resolve())
    return paths


def run(full, floor, *, output, levels=DEFAULT_LEVELS, backend="cpu", equilibrium=None,
        geometry=None, design_policy=None):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"replay output already exists: {output}")
    protocol = prepare_protocol(full, floor, levels, backend)
    source = protocol["inputs"]
    paths = _replay_input_paths(source, equilibrium=equilibrium, geometry=geometry, design_policy=design_policy)
    protocol["replay_input_paths"] = paths
    user = load_user_input(paths["equilibrium"], source["time_ms"], paths["geometry"], paths["design_policy"])
    for key in INPUT_KEYS:
        _check(source[key] == user.provenance[key], f"replay source changed for {key}")
    reference = ReferencePoint(**protocol["reference_SI"])
    protocol["replay_software"] = software()
    for module, checksum in protocol["search_code_sha256"].items():
        if module in PHYSICS_MODULES:
            _check(protocol["replay_software"]["code_sha256"].get(module) == checksum,
                   f"replay solver source changed for {module}")
    protocol["driver_sha256"] = file_hash(__file__)
    output.mkdir(parents=True, exist_ok=False)
    public_save(output/"protocol.json", serial(protocol))
    report = dict(complete=False, protocol_sha256=file_hash(output/"protocol.json"), results={})
    public_save(output/"report.json", report)
    for name, resolution in protocol["levels"].items():
        folder = output/name
        folder.mkdir()
        start = perf_counter()
        original = prepare_user(user, reference, resolution)
        solver = evaluator(original, protocol["matched_search"]["radius_mm"], backend)
        grid = dict(resolution=resolution, prepare_s=perf_counter()-start, fixed_shelf=shelf_limits(user, resolution[0]),
                    candidates={})
        base = original.baseline
        grid["baseline"] = (_measurement(folder/"baseline", original, base) if "deposition" in base
                            else dict(physics_measurable=False, errors=base.get("errors", [])))
        report["results"][name] = grid
        public_save(output/"report.json", serial(report))
        for identity, record in protocol["candidates"].items():
            candidate = Candidate.from_record(record)
            t = perf_counter()
            try:
                raw = evaluate_user_candidate(solver, candidate)
            except (ValueError, ArithmeticError) as exc:
                raw = dict(feasible=False, errors=[f"{type(exc).__name__}: {exc}"])
            compute = perf_counter()-t
            if "deposition" in raw:
                row = _measurement(folder/identity, solver, raw)
            else:
                row = dict(physics_measurable=False, candidate_eligible=False, errors=raw.get("errors", []),
                           geometry_screen=raw.get("geometry_screen"), shape_constraints=raw.get("shape_constraints"),
                           precise_length_ratios=raw.get("precise_length_ratios"))
            grid["candidates"][identity] = dict(row, evaluation_s=compute, id=identity)
            public_save(output/"report.json", serial(report))
            del raw
        if tuple(resolution) == FINE:
            grid["reproduction"] = _reproduction(protocol, grid)
        grid["elapsed_s"] = perf_counter()-start
        public_save(output/"report.json", serial(report))
        del original, solver, base
    report["complete"] = True
    public_save(output/"report.json", serial(report))
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", type=Path, required=True, help="full search directory or result.json")
    parser.add_argument("--floor", type=Path, required=True, help="floor-only search directory or result.json")
    parser.add_argument("--output", type=Path, required=True, help="new exclusive replay directory")
    parser.add_argument("--levels", choices=LEVELS, nargs="+", default=DEFAULT_LEVELS)
    parser.add_argument("--backend", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--equilibrium", type=Path, help="relocated DIII-D parquet matching the archived SHA-256")
    parser.add_argument("--geometry", type=Path, help="relocated material geometry matching the archived SHA-256")
    parser.add_argument("--design-policy", type=Path, help="relocated design policy matching the archived SHA-256")
    args = parser.parse_args(argv)
    report = run(**vars(args))
    check = report["results"].get("cubic1280", {}).get("reproduction")
    return 1 if check and not check["passed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
