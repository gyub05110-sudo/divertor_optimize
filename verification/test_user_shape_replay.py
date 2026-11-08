import importlib
import importlib.util
import json
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

import pytest

from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.reference_transport import ReferencePoint


def replay():
    assert importlib.util.find_spec("verification.user_shape_replay"), "matched-run replay is missing"
    return importlib.import_module("verification.user_shape_replay")


def wall(first=0., floor=0.):
    return Candidate("OMP", 20, (first,)*12 + (floor,)*4)


def write(path, data):
    path.write_text(json.dumps(data), encoding="utf-8")


@pytest.fixture
def pair(tmp_path):
    inputs = dict(equilibrium_path="shot", geometry_path="wall", design_policy_path="policy",
                  time_ms=1400., **{key+"_sha256": key for key in (
                      "equilibrium", "geometry", "design_policy", "wall", "omp_endpoints", "design_reference")})
    for key in ("equilibrium", "geometry", "design_policy"):
        path = tmp_path/inputs[key+"_path"]
        path.write_bytes(key.encode())
        inputs[key+"_path"] = str(path)
        inputs[key+"_sha256"] = sha256(path.read_bytes()).hexdigest()
    for arm in ("full", "floor-only"):
        folder = tmp_path / arm
        folder.mkdir()
        origin, winner, other = wall(), wall(0., .01), wall(.001 if arm == "full" else 0., .02)
        journal = []
        for i, (candidate, peak, kind) in enumerate(zip(
                [origin, winner, other], [8., 4., 5.], ["projected_origin", "local_projected", "projected_restart"])):
            row = dict(id=candidate.identity, index=i, candidate=candidate.record(), proposal_kind=kind,
                       within_budget=True, status="completed", start_s=i*3., end_s=i*3.+1)
            journal.extend([dict(row, stage="proposal"), dict(row, stage="coarse", eligible=True, peak_q=peak),
                            dict(row, stage="fine", end_s=i*3.+2, eligible=True, peak_q=peak)])
        journal.append(dict(journal[-1], stage="fine", within_budget=False, end_s=61., status="late", peak_q=.1))
        (folder/"events.jsonl").write_text("".join(json.dumps(e)+"\n" for e in journal), encoding="utf-8")
        report = dict(schema_version=1, status="qualified", arguments=dict(family=arm, budget_seconds=60.,
                    seed=3, radius_mm=20, backend="cpu"), inputs=inputs, reference_SI=asdict(ReferencePoint()),
                    levels=dict(coarse=["cubic",320,.0005,161,.001], fine=["cubic",1280,.000125,641,.00025]),
                    search=dict(termination="deadline"), software=dict(code_sha256={"user_physics.py": "123"}),
                    candidate=dict(id=winner.identity, **winner.record()),
                    selected=dict(peak_MW_m2=4., A2_m2=.8, candidate_eligible=True),
                    baseline=dict(peak_MW_m2=8., A2_m2=.4, physics_measurable=True, reference_only=True))
        write(folder/"result.json", report)
    return tmp_path/"full", tmp_path/"floor-only"


def test_protocol_freezes_matched_arms_and_controls_deduplicated_by_identity(pair):
    api = replay()
    protocol = api.prepare_protocol(pair[0]/"result.json", pair[1], ["cubic1280","quintic1280"], "cpu")
    roles, geometries = protocol["roles"], protocol["candidates"]
    assert roles["full_winner"]["id"] == roles["floor_winner"]["id"] == roles["full_winner_floor_ablation"]["id"]
    assert roles["full_nonfloor_finalist"]["id"] == roles["restart_control"]["id"] == wall(.001,.02).identity
    assert roles["full_nonfloor_finalist_floor_ablation"]["id"] == wall(0.,.02).identity
    assert roles["g1_zero"]["id"] == roles["projected_origin"]["id"] == roles["one_amplitude_+0mm"]["id"]
    assert len([k for k in roles if k.startswith("one_amplitude_")]) == 5
    assert set(geometries) == {role["id"] for role in roles.values() if role["id"] is not None}
    assert all(Candidate.from_record(row).identity == key for key,row in geometries.items())


@pytest.mark.parametrize("tamper", ["seed", "geometry_sha256", "source_code", "selected_score", "winner", "budget"])
def test_mismatched_or_impossible_study_evidence_is_rejected(pair, tamper):
    api = replay()
    p = pair[0]/"result.json"
    r = json.loads(p.read_text())
    if tamper == "seed":
        r["arguments"]["seed"] += 1
    elif tamper == "geometry_sha256":
        r["inputs"]["geometry_sha256"] += "other"
    elif tamper == "source_code":
        r["software"]["code_sha256"]["user_physics.py"] += "different"
    elif tamper == "selected_score":
        r["selected"]["peak_MW_m2"] = 1.
    elif tamper == "winner":
        r["candidate"] = dict(id=wall().identity, **wall().record())
    else:
        r["arguments"]["budget_seconds"] = 1.
    write(p, r)
    with pytest.raises(ValueError):
        api.prepare_protocol(pair[0], pair[1], ["cubic1280"], "cpu")


def test_replay_refuses_to_write_an_existing_directory_before_seeking_heat(pair, tmp_path):
    api = replay()
    out = tmp_path/"full"
    before = (out/"result.json").read_bytes()
    with pytest.raises(FileExistsError):
        api.run(*pair, output=out)
    assert before == (out/"result.json").read_bytes()


def test_one_calculation_per_geometry_per_level_and_preservation_of_rejections(pair, tmp_path, monkeypatch):
    api = replay()
    protocol = api.prepare_protocol(*pair, ["cubic1280"], "cpu")
    canonical = protocol["inputs"]
    monkeypatch.setattr(api, "software", lambda: {"code_sha256":protocol["search_code_sha256"]})
    monkeypatch.setattr(api, "load_user_input", lambda *args: SimpleNamespace(provenance=canonical))
    monkeypatch.setattr(api, "shelf_limits", lambda *a: {"admissible": True})
    original = SimpleNamespace(baseline=dict(physics_measurable=True, deposition={}, peak_q=8.))
    monkeypatch.setattr(api, "prepare_user", lambda *args: original)
    monkeypatch.setattr(api, "evaluator", lambda *args: object())
    calls = []
    def evaluate(e, c):
        calls.append(c.identity)
        if c.vector[-1] < 0:
            return dict(feasible=False, errors=["continuous_shape_constraint_failed"],
                        shape_constraints={"passed": False})
        return dict(physics_measurable=True, feasible=True, deposition={},
                    peak_q=4. if c==wall(0,.01) else 5.)
    monkeypatch.setattr(api, "evaluate_user_candidate", evaluate)
    def measure(out, e, r):
        return dict(candidate_eligible=bool(r.get("feasible")), physics_measurable=True,
                    peak_MW_m2=r["peak_q"], A2_m2=.4 if e is original else .8)
    monkeypatch.setattr(api, "_measurement", measure)
    out = tmp_path/"replayed"
    relocated = tmp_path/"relocated-shot.parquet"
    relocated.write_bytes(Path(canonical["equilibrium_path"]).read_bytes())
    report = api.run(*pair, output=out, backend="cpu", levels=["cubic1280"], equilibrium=relocated)
    assert len(set(calls)) == len(calls) == len(protocol["candidates"])
    grid = report["results"]["cubic1280"]
    assert grid["reproduction"]["passed"]
    assert grid["candidates"][wall(0.,-.04).identity]["shape_constraints"] == {"passed":False}
    assert (out/"protocol.json").is_file() and (out/"report.json").is_file()
    saved = json.loads((out/"protocol.json").read_text())
    assert saved["inputs"]["equilibrium_path"] == canonical["equilibrium_path"]
    assert saved["replay_input_paths"] == dict(equilibrium=str(relocated),
        geometry=canonical["geometry_path"], design_policy=canonical["design_policy_path"])


def test_replay_detects_changes_to_the_underlying_solver_before_new_heat_is_evaluated(pair, tmp_path, monkeypatch):
    api = replay()
    protocol = api.prepare_protocol(*pair)
    monkeypatch.setattr(api, "load_user_input", lambda *args: SimpleNamespace(provenance=protocol["inputs"]))
    monkeypatch.setattr(api, "software", lambda: {"code_sha256":{"user_physics.py":"changed"}})
    monkeypatch.setattr(api, "prepare_user", lambda *a: pytest.fail("heat setup started with changed physics code"))
    with pytest.raises(ValueError, match="replay solver source changed"):
        api.run(*pair, output=tmp_path/"changed")


def test_missing_archived_material_and_policy_use_checksum_verified_checkout_files(pair, tmp_path, monkeypatch):
    api = replay()
    source = api.prepare_protocol(*pair)["inputs"]
    overrides = {}
    for key, variable in (("geometry", "GEOMETRY"), ("design_policy", "POLICY")):
        moved = tmp_path/("repo-"+key)
        Path(source[key+"_path"]).rename(moved)
        monkeypatch.setattr(api, variable, moved)
        overrides[key] = str(moved)
    paths = api._replay_input_paths(source)
    assert paths == dict(equilibrium=source["equilibrium_path"], **overrides)
    missing_shot = Path(source["equilibrium_path"])
    missing_shot.unlink()
    with pytest.raises(FileNotFoundError, match="--equilibrium"):
        api._replay_input_paths(source)


@pytest.mark.parametrize("key", ["equilibrium", "geometry", "design_policy"])
def test_explicit_overrides_take_precedence_and_reject_a_mismatched_copy(pair, tmp_path, key):
    api = replay()
    source = api.prepare_protocol(*pair)["inputs"]
    copy = tmp_path/"override"
    copy.write_bytes(b"other")
    with pytest.raises(ValueError, match=key):
        api._replay_input_paths(source, **{key:copy})
    copy.write_bytes(Path(source[key+"_path"]).read_bytes())
    assert api._replay_input_paths(source, **{key:copy})[key] == str(copy)
