"""Standalone equilibrium inputs feed the same fixed-wall physics as archived runs."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from d3d_optimization.case import load_case
from d3d_optimization.independent_sources import split_design_wall


ROOT = Path(__file__).resolve().parents[1]
SHOT = ROOT / "input/D3D/d3d/d3d_shot_1a3155370c.parquet"
GEOMETRY = ROOT / "d3d_optimization/data/DIIID_geom.json"
POLICY = ROOT / "baseline_review/optimization_region_policy.json"
REVIEW = ROOT / "baseline_review/baseline_review.json"


def loader():
    assert importlib.util.find_spec("d3d_optimization.user_input"), "user-equilibrium loader is missing"
    from d3d_optimization.user_input import load_user_input
    return load_user_input


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def omp_endpoints(review):
    return {s["name"]: s["endpoints_RZ_m"] for s in split_design_wall(review)[1]
            if s["name"] in {"shelf", "floor", "omp_entry", "omp_connector"}}


def test_renamed_equilibrium_uses_requested_field_and_fixed_design(tmp_path):
    path = tmp_path / "my-discharge.data"
    shutil.copyfile(SHOT, path)
    # Equivalent, user-named geometry/policy with different byte encoding.
    geom, policy = tmp_path / "material.json", tmp_path / "policy.json"
    for dest, src in ((geom, GEOMETRY), (policy, POLICY)):
        dest.write_text(json.dumps(json.loads(src.read_bytes())), encoding="utf8")
    model = loader()(path, 1200., geom, policy)
    archived = json.loads(REVIEW.read_bytes())
    assert model.time_ms == 1200 and model.equilibrium_sha256 == file_hash(path)
    assert model.geometry == json.loads(GEOMETRY.read_bytes())
    assert model.review["wall"] == archived["wall"]
    assert omp_endpoints(model.review) == omp_endpoints(archived)
    assert model.review["legs"]["outer"]["strike_RZ_m"] != archived["legs"]["outer"]["strike_RZ_m"]
    for leg in model.review["legs"].values():
        assert 0 <= leg["refinement_m"] <= 1e-5
        wall = np.array(model.review["wall"])
        a, b = wall[leg["segment"]:leg["segment"] + 2]
        np.testing.assert_allclose(np.linalg.norm(np.subtract(leg["strike_RZ_m"], a)) +
                                   np.linalg.norm(np.subtract(leg["strike_RZ_m"], b)), np.linalg.norm(b-a))
    for method in ("cubic", "quintic"):
        field, lcfs, axis = model.field(method)
        ref, ref_lcfs, ref_axis = load_case(path, 1200, method)
        np.testing.assert_array_equal(lcfs, ref_lcfs)
        np.testing.assert_array_equal(axis, ref_axis)
        np.testing.assert_array_equal(field.magnetic([[1.44, -1.24], [1.83, -0.49]]),
                                      ref.magnetic([[1.44, -1.24], [1.83, -0.49]]))
    assert model.provenance["geometry_sha256"] == file_hash(geom)
    assert model.provenance["design_policy_sha256"] == file_hash(policy)
    assert model.provenance["wall_sha256"] == sha256(np.array(archived["wall"], dtype="<f8").tobytes()).hexdigest()
    assert model.provenance["frame_index"] == 50
    assert model.provenance["baseline_outer_strike_design"] == {"location":"adjustable_interior", "surfaces":["floor"]}
    assert model.provenance["baseline_lower_strike_legs"] == model.review["legs"]


def test_supported_input_requires_requested_time(tmp_path):
    with pytest.raises(ValueError, match="exactly once"):
        loader()(SHOT, 1400.0000001)
    with pytest.raises(ValueError, match="finite"):
        loader()(SHOT, float("nan"))


def test_wrong_material_and_unsupported_policy_rejected(tmp_path):
    geom = json.loads(GEOMETRY.read_bytes())
    geom["limiter"][0][0] += .001
    output = tmp_path / "custom.json"
    output.write_text(json.dumps(geom))
    with pytest.raises(ValueError, match="material geometry"):
        loader()(SHOT, 1400, geometry=output)
    geom = json.loads(GEOMETRY.read_bytes())
    geom["vv"][0]["pts"][0][0] += .001
    output.write_text(json.dumps(geom))
    with pytest.raises(ValueError, match="material geometry"):
        loader()(SHOT, 1400, geometry=output)
    policy = json.loads(POLICY.read_bytes())
    policy["excluded_omp_segments"].append(59)
    output.write_text(json.dumps(policy))
    with pytest.raises(ValueError, match="design policy"):
        loader()(SHOT, 1400, design_policy=output)


def sample_record():
    grid = np.linspace(1, 3, 6).tolist()
    return dict(source="DIII-D", efit_times=[1000.], efit_grid_R=grid, efit_grid_Z=grid,
                efit_psirz=[np.ones((6, 6)).tolist()], efit_lcfs_n=[4.],
                efit_lcfs_r=[[1.3, 1.7, 1.6, 1.3]], efit_lcfs_z=[[1.3, 1.3, 1.7, 1.3]],
                efit_r_axis=[1.5], efit_z_axis=[1.4], magnetics_bcoil=[20., 21.],
                magnetics_time=[998., 1002.])


@pytest.mark.parametrize("update,error", [
    ({"source":"other device"}, "DIII-D"),
    ({"magnetics_time":[999., 999.9]}, "cover"),
    ({"magnetics_time":[1001., 999.]}, "increasing"),
    ({"magnetics_bcoil":[20., float("nan")]}, "finite"),
    ({"efit_times":[1000., 1000.]}, "exactly once"),
    ({"efit_lcfs_n":[20.]}, "LCFS"),
    ({"efit_r_axis":[float("nan")]}, "axis"),
])
def test_invalid_parquet_reports_physical_input_error(tmp_path, update, error):
    path = tmp_path / "invalid.parquet"
    pq.write_table(pa.Table.from_pylist([dict(sample_record(), **update)]), path)
    with pytest.raises(ValueError, match=error):
        loader()(path, 1000)


def test_multishot_file_and_missing_schema_are_rejected(tmp_path):
    path = tmp_path / "not-a-discharge.parquet"
    pq.write_table(pa.Table.from_pylist([sample_record(), sample_record()]), path)
    with pytest.raises(ValueError, match="one DIII-D"):
        loader()(path, 1000)
    pq.write_table(pa.Table.from_pylist([{"source":"DIII-D"}]), path)
    with pytest.raises(ValueError, match="missing.*efit_times"):
        loader()(path, 1000)

@pytest.mark.parametrize("update,error", [
    ({"magnetics_time":[[998.], [1002.]]}, "magnetics_time"),
    ({"efit_psirz":[[None]]}, "efit_psirz"),
    ({"efit_grid_Z":[0., 1., 2.]}, "grid"),
    ({"efit_grid_R":[1., 1., 1., 2., 3., 4.]}, "grid|increasing"),
    ({"efit_lcfs_n":["4"]}, "efit_lcfs_n"),
    ({"efit_lcfs_r":[["x", "y", "a", "x"]]}, "efit_lcfs_r"),
])
def test_bad_column_types_are_actionable_value_errors(tmp_path, update, error):
    path = tmp_path / "wrong-shape.parquet"
    pq.write_table(pa.Table.from_pylist([dict(sample_record(), **update)]), path)
    with pytest.raises(ValueError, match=error):
        loader()(path, 1000)


def test_modified_packaged_reference_cannot_silently_move_engineering_endpoints(tmp_path, monkeypatch):
    loader()
    from d3d_optimization import user_input
    path = tmp_path / "design.json"
    ref = json.loads(REVIEW.read_bytes())
    next(f for f in ref["optimization_footprints"] if f["segment"] == 58)["points"][0][0] += .001
    path.write_text(json.dumps(ref))
    monkeypatch.setattr(user_input, "REFERENCE", path)
    with pytest.raises(ValueError, match="engineering.*reference"):
        user_input.load_user_input(SHOT, 1400)

def test_outer_strike_on_fixed_wall_is_supported_without_retargeting_design():
    from d3d_optimization import user_input
    classify = getattr(user_input, "_outer_strike_location", None)
    assert classify is not None, "report the strike separately from the engineering policy"
    endpoints = omp_endpoints(json.loads(REVIEW.read_bytes()))
    shelf = np.asarray(endpoints["shelf"])
    fixed = classify([1.90, -1.05], endpoints)
    assert fixed == {"location":"fixed_wall", "surfaces":[]}
    assert classify(shelf[1], endpoints) == {"location":"adjustable_endpoint", "surfaces":["shelf"]}
    assert classify(shelf.mean(axis=0), endpoints) == {"location":"adjustable_interior", "surfaces":["shelf"]}

def test_numeric_null_lcfs_is_an_input_error(tmp_path):
    record = dict(sample_record(), efit_lcfs_r=[[1.3, None, 1.6, 1.3]])
    path = tmp_path / "missing-coordinate.parquet"
    pq.write_table(pa.Table.from_pylist([record]), path)
    with pytest.raises(ValueError, match="LCFS"):
        loader()(path, 1000)


def test_zero_tf_current_is_preserved_when_input_supplier_sets_it():
    from d3d_optimization.user_input import _parquet
    record = dict(sample_record(), magnetics_bcoil=[0., 0.])
    stream = pa.BufferOutputStream()
    pq.write_table(pa.Table.from_pylist([record]), stream)
    field, _, _, _ = _parquet(stream.getvalue(), 1000)
    assert field.tf_current_ka == 0 and field.magnetic([[1.5,1.5]])[0,2] == 0

def test_unrelated_missing_coil_sample_does_not_invalidate_named_frame():
    from d3d_optimization.user_input import _parquet
    # The finite samples at 998 and 1002 ms bracket the requested 1000 ms.
    record = dict(sample_record(), magnetics_time=[300., 998., 1002., 1700.],
                  magnetics_bcoil=[float("nan"), 20., 26., float("nan")])
    stream = pa.BufferOutputStream()
    pq.write_table(pa.Table.from_pylist([record]), stream)
    field, _, _, meta = _parquet(stream.getvalue(), 1000)
    assert field.tf_current_ka == 23.
    assert meta["tf_interpolation_times_ms"] == [998., 1002.]
    assert meta["tf_interpolation_currents_kA"] == [20., 26.]
