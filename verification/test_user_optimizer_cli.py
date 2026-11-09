import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "optimize_divertor.py"


def invoke(out, source, *args):
    return subprocess.run([sys.executable, str(CLI), "--output", str(out),
                           "--equilibrium", str(source), "--time-ms", "1400", *args],
                          capture_output=True, text=True, timeout=45)


def test_missing_file_has_machine_readable_error_and_preserves_it(tmp_path):
    out = tmp_path / "result"
    run = invoke(out, tmp_path / "missing.parquet")
    assert run.returncode == 2
    result = out / "result.json"
    assert json.loads(result.read_text())["status"] == "input_error"
    original = result.read_bytes()
    assert invoke(out, tmp_path / "different.parquet").returncode != 0
    assert result.read_bytes() == original


def test_expired_budget_cannot_return_a_shape(tmp_path):
    out = tmp_path / "budget"
    source = ROOT / "input/D3D/d3d/d3d_shot_1a3155370c.parquet"
    result = invoke(out, source, "--budget-seconds", "0.000001")
    assert result.returncode == 0, result.stderr
    row = json.loads((out / "result.json").read_text())
    assert row["status"] == "budget_exhausted"
    assert row["baseline"] is None and row["selected"] is None


def test_invalid_operating_point_is_rejected_before_run(tmp_path):
    out = tmp_path / "invalid"
    assert invoke(out, "example.parquet", "--lambda-q-mm", "nan").returncode == 2
    assert not out.exists()


@pytest.fixture
def synthetic_preparation(tmp_path, monkeypatch):
    import optimize_divertor as cli
    from d3d_optimization import campaign_geometry, user_input, user_output, user_physics, user_precheck, user_search

    path = tmp_path / "named.parquet"
    path.touch()
    user = SimpleNamespace(provenance={}, review={k: {} for k in ("axis_RZ_m", "xpoint_RZ_m", "legs")})
    fixed = dict(excluded=False, classification="no_fixed_wall_lcfs_intersection")
    calls = []
    def precheck(item):
        assert item is user
        calls.append("fixed_wall")
        return fixed
    monkeypatch.setattr(user_precheck, "fixed_wall_precheck", precheck)
    monkeypatch.setattr(user_input, "load_user_input", lambda *a: user)
    def limits(*args):
        calls.append("shelf_limits")
        return {"admissible": True}
    monkeypatch.setattr(user_physics, "shelf_limits", limits)
    names = ["fixed", "omp_entry", "shelf", "omp_connector", "floor", "imp_sloped", "imp_vertical"]
    segments = np.arange(len(names))
    native = dict(physics_measurable=True, reference_only=True, wall=np.zeros((8, 2)),
                  region_names=names, surface_segments={name: [i] for i, name in enumerate(names[1:], 1)},
                  hits=SimpleNamespace(segment=np.array([4])),
                  deposition=dict(cell_parent=segments, sample_parent=segments,
                                  cell_power=np.array([.1, .1, .1, .3, .2, .1, .1]),
                                  q=np.array([1., 2., 1., 4., 3., 0., 1.])))
    original = SimpleNamespace(baseline=native, source=SimpleNamespace(weights=np.ones(1)),
                               reference=SimpleNamespace(power_W=1e6))
    def prepare(_user, _reference, level):
        calls.append("fine" if level == cli.FINE else "coarse")
        if level == cli.COARSE:
            return SimpleNamespace(baseline={"physics_measurable": False})
        return original
    monkeypatch.setattr(user_physics, "prepare_user", prepare)
    monkeypatch.setattr(campaign_geometry, "evaluator", lambda owner, *args: owner)
    monkeypatch.setattr(user_output, "export_measurement", lambda *a: {"physics_measurable": True})
    monkeypatch.setattr(cli, "software", lambda: {})
    out = tmp_path / "result"
    def search(coarse, fine, *, family, **kwargs):
        assert fine is original and coarse is not original
        calls.append(family)
        events = [json.loads(line) for line in (out / "events.jsonl").read_text().splitlines()]
        decision, = [e for e in events if e["stage"] == "surface_policy_decision" and e["status"] == "recorded"]
        assert decision["policy"]["selected_family"] == family
        return dict(status="no_design", selected=None, winner=None)
    monkeypatch.setattr(user_search, "search_peak", search)
    def run(*args):
        code = cli.main(["--equilibrium", str(path), "--output", str(out), "--time-ms", "1400", *args])
        row = json.loads((out / "result.json").read_text())
        return code, row
    return SimpleNamespace(run=run, fixed=fixed, calls=calls, original=original, prepare=prepare)


def test_completed_reference_survives_coarse_preparation_failure(synthetic_preparation, monkeypatch):
    import optimize_divertor as cli
    from d3d_optimization import user_physics
    setup = synthetic_preparation
    def prepare(user, reference, level):
        if level == cli.COARSE:
            raise ValueError("coarse trace unresolved")
        return setup.prepare(user, reference, level)
    monkeypatch.setattr(user_physics, "prepare_user", prepare)
    code, row = setup.run()
    assert code == 1
    assert row["baseline"]["physics_measurable"] and row["selected"] is None
    assert row["status"] == "calculation_error" and "coarse trace" in row["reason"]


@pytest.mark.parametrize("family", ["hotspot", "full", "floor-only", "floor-entry", "floor-shelf", "floor-connector"])
def test_cli_freezes_original_fine_prediction_and_executes_requested_arm(synthetic_preparation, family):
    setup = synthetic_preparation
    code, row = setup.run("--family", family)
    selected = "floor-connector" if family == "hotspot" else family
    assert code == 0 and row["status"] == "no_design"
    assert setup.calls == ["fixed_wall", "shelf_limits", "fine", "coarse", selected]
    policy = row["surface_policy"]
    assert policy["requested_family"] == family and policy["selected_family"] == selected
    assert policy["predicted_family"] == "floor-connector" and policy["predicted_floor_sufficient"] is False
    assert policy["regional_metrics"]["omp_connector"] == {"peak_MW_m2": 4., "deposited_W": 3e5}
    stages = [t["stage"] for t in row["timings"]["preparation"]]
    assert stages.index("surface_policy") < stages.index("prepare_fine_design")


def test_fixed_wall_exclusion_precedes_any_shelf_or_thermal_evaluation(synthetic_preparation):
    setup = synthetic_preparation
    setup.fixed.update(excluded=True, classification="fixed_wall_lcfs_intersection",
                       fixed_segment_intersections=[{"segment": 11}])
    code, row = setup.run()
    assert code == 0 and row["status"] == "geometry_excluded"
    assert "fixed" in row["reason"] and row["fixed_wall"] == setup.fixed
    assert row["baseline"] is row["selected"] is row["surface_policy"] is None
    assert setup.calls == ["fixed_wall"]


def test_late_presearch_prediction_cannot_start_design_or_enter_record(synthetic_preparation, monkeypatch):
    import optimize_divertor as cli
    from d3d_optimization import user_output
    now = [100.]
    monkeypatch.setattr(cli, "perf_counter", lambda: now[0])
    original = user_output.surface_metrics
    def delayed(*args):
        now[0] += 11.
        return original(*args)
    monkeypatch.setattr(user_output, "surface_metrics", delayed)
    code, row = synthetic_preparation.run("--budget-seconds", "10")
    assert code == 0 and row["status"] == "budget_exhausted"
    assert row["surface_policy"] is row["candidate"] is None
    assert synthetic_preparation.calls == ["fixed_wall", "shelf_limits", "fine"]
    assert row["timings"]["preparation"][-1]["within_budget"] is False


def test_late_coarse_reference_cannot_start_design_setup(synthetic_preparation, monkeypatch):
    import optimize_divertor as cli
    from d3d_optimization import campaign_geometry, user_physics
    setup = synthetic_preparation
    now = [100.]
    monkeypatch.setattr(cli, "perf_counter", lambda: now[0])
    def prepare(user, reference, level):
        if level == cli.COARSE:
            now[0] += 11.
        return setup.prepare(user, reference, level)
    design_setups = []
    def evaluator(owner, *args):
        design_setups.append(owner)
        return owner
    monkeypatch.setattr(user_physics, "prepare_user", prepare)
    monkeypatch.setattr(campaign_geometry, "evaluator", evaluator)
    code, row = setup.run("--budget-seconds", "10")
    assert code == 0 and row["status"] == "budget_exhausted" and row["candidate"] is None
    assert design_setups == [setup.original]


def test_late_fixed_wall_exclusion_is_never_admitted(synthetic_preparation, monkeypatch):
    import optimize_divertor as cli
    from d3d_optimization import user_precheck
    now = [100.]
    monkeypatch.setattr(cli, "perf_counter", lambda: now[0])
    def delayed(user):
        now[0] += 11.
        return {"excluded": True, "classification": "fixed_wall_lcfs_intersection"}
    monkeypatch.setattr(user_precheck, "fixed_wall_precheck", delayed)
    code, row = synthetic_preparation.run("--budget-seconds", "10")
    assert code == 0 and row["status"] == "budget_exhausted"
    assert row["fixed_wall"] is row["surface_policy"] is row["candidate"] is None
    assert not synthetic_preparation.calls


def test_screened_measurement_without_regional_evidence_fails_before_search(synthetic_preparation):
    setup = synthetic_preparation
    del setup.original.baseline["hits"]
    code, row = setup.run()
    assert code == 1 and row["status"] == "calculation_error"
    assert "hits" in row["reason"] and row["surface_policy"] is None
    assert setup.calls == ["fixed_wall", "shelf_limits", "fine"]


@pytest.mark.parametrize("family", ["hotspot", "floor-only", "floor-entry", "floor-shelf", "floor-connector"])
def test_cli_direct_still_requires_full(tmp_path, family):
    proc = invoke(tmp_path / "out", tmp_path / "unused", "--method", "direct", "--family", family)
    assert proc.returncode == 2 and "requires --family full" in proc.stderr


def test_reference_corner_is_assessed_on_original_straight_floor():
    from d3d_optimization.user_output import reference_junction
    wall = np.array([[-1., 0.], [0., 0.], [1., 0.], [1., 1.]])
    assert reference_junction(wall, 1)["tangent_error_deg"] == pytest.approx(90.)
    assert reference_junction(wall, 1)["passed"] is False
    wall[-1] = [2., 0.]
    assert reference_junction(wall, 1)["passed"] is True


def test_deposited_and_first_hit_power_distinguish_unhit_and_transfer():
    from types import SimpleNamespace
    from d3d_optimization.user_output import surface_metrics
    native = dict(wall=np.zeros((4, 2)), region_names=["fixed", "floor"],
                  surface_segments={"floor": [1, 2]},
                  hits=SimpleNamespace(segment=np.array([1, 0, -1])),
                  deposition=dict(cell_parent=np.array([0, 1, 2]),
                                  cell_power=np.array([.25, .30, .35]),
                                  sample_parent=np.array([0, 1, 2]),
                                  q=np.array([2., 4., 3.])))
    rows = surface_metrics(native, np.array([.75, .15, .10]), 2e6)
    assert rows["floor"]["first_hit_W"] == 1.5e6
    assert rows["floor"]["deposited_W"] == pytest.approx(1.3e6)
    assert rows["fixed"]["first_hit_W"] == .3e6
    assert rows["floor"]["peak_MW_m2"] == 8
