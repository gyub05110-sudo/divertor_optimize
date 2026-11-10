import h5py
import numpy as np
import pytest

from verification import stage5_native as native


@pytest.mark.skipif(not (native.HISTORICAL/"blind_screen/results/ref_OMP_32027.nc").is_file(),
                    reason="requires the retained native FIREFLY reference NetCDF")
def test_ledger_checks_real_sector_tally_and_rejects_censor_corruption(tmp_path):
    case = native.HISTORICAL / "blind_screen"
    raw = case / "results/ref_OMP_32027.nc"
    wall = np.loadtxt(case / "walls/ref.wall")
    tally = native.audit_tally(raw, wall, 1024)
    assert tally["wall_W"] > 900000
    assert abs(tally["closure_residual_W"]) < 0.3
    assert tally["max_wall_distance_m"] < 1e-10
    with pytest.raises(ValueError, match="histories"):
        native.audit_tally(raw, wall, 4096)
    import shutil
    bad = tmp_path / "bad.nc"
    shutil.copyfile(raw, bad)
    with h5py.File(bad, "r+") as f:
        f["parameters/censored_fraction"][()] += 0.01
    with pytest.raises(ValueError, match="closure"):
        native.audit_tally(bad, wall, 1024)


def test_casing_encloses_two_unchanged_walls_and_retains_source_midplane():
    wall = np.array([[1., -1.], [2., -1.], [2., 1.], [1., 1.], [1., -1.]])
    altered = wall.copy()
    altered[1, 0] += .020
    casing, clearance = native.enclose_walls([wall, altered])
    assert clearance > .010
    assert casing[:, 0].max() >= 2.030
    assert np.array_equal(altered[1], [2.02, -1.])


def test_source_rejects_permutation_of_heat_weights(tmp_path):
    power = np.array([200000., 800000.-1e-8, 1e-8])
    radii, z = np.array([2.25, 2.255, 2.275]), np.zeros(3)
    np.savetxt(tmp_path/"OMP.csv", np.c_[radii, z, -np.ones(3), power/1e6, power],
               header="R_m,Z_m,sign_along_B,power_fraction,power_W", delimiter=",", comments="")
    cdf = np.c_[radii, z, np.cumsum(power)/1e6]
    cdf[-1,2] = 1.
    np.savetxt(tmp_path / "OMP.source", cdf, header=str(len(cdf)), comments="")
    source = native.audit_source(tmp_path)
    assert source["source_power_W"] == pytest.approx(1e6)
    cdf[-2,2] = 1.0000000000000002
    np.savetxt(tmp_path / "OMP.source", cdf, header=str(len(cdf)), comments="")
    with pytest.raises(ValueError, match="CDF"):
        native.audit_source(tmp_path)
    cdf[-2,2] = 1.
    cdf[1, 2] = cdf[0, 2]
    np.savetxt(tmp_path / "OMP.source", cdf, header=str(len(cdf)), comments="")
    with pytest.raises(ValueError, match="CDF"):
        native.audit_source(tmp_path)


@pytest.mark.parametrize("toroidal_sign,direction", [(-1, 1), (1, -1)])
def test_lower_leg_uses_signed_toroidal_field(tmp_path, toroidal_sign, direction):
    import json
    inp = tmp_path/"inputs"/"eq"
    base = tmp_path/"native/cases"/"eq"
    (inp/"walls").mkdir(parents=True)
    base.mkdir(parents=True)
    (inp/"case.json").write_text(json.dumps({"topology":{"F_T_m":toroidal_sign}}))
    (inp/"OMP.csv").write_text("Z_m,sign_along_B,power_W\n0,-1,250000\n0,-1,750000\n")
    np.savetxt(inp/"walls/fixed.hits", [[1.5,-1.,25.],[1.6,-1.,25.]])
    native_trace = np.array([[25.,25.,0,1.,1.,1.,1.],[25.,25.,0,1.,1.,1.,1.]])
    col = 3 if direction<0 else 5
    native_trace[:,col:col+2] = [[1.5002,-1.],[1.6002,-1.]]
    np.savetxt(base/"fixed.strikes.dat", native_trace)
    result = native.strike_audit(tmp_path,"eq","fixed")
    assert result["pass"] and result["toroidal_direction"] == direction
    assert result["weighted_p95_hit_error_m"] == pytest.approx(.0002)
    native_trace[1,col+1] = 0.2
    np.savetxt(base/"fixed.strikes.dat", native_trace)
    result = native.strike_audit(tmp_path,"eq","fixed")
    assert not result["pass"] and result["above_midplane_W"] == 750000
