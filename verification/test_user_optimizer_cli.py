import json
from pathlib import Path
import subprocess
import sys

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
