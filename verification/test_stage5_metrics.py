from itertools import product

import h5py
import numpy as np
import pytest


WALL = np.array([[1., 0.], [1.006, 0.], [1.006, .004], [1., .004], [1., 0.]])
REGIONS = ["floor", "shelf", "other", "other"]


def test_conservative_area_ring_projection_and_surface_identity():
    from verification.stage5_metrics import common_rings, project, summary
    grid = common_rings(WALL, REGIONS)
    areas = np.pi * (WALL[1:, 0] + WALL[:-1, 0]) * np.linalg.norm(np.diff(WALL, axis=0), axis=1)
    watts = areas * [3e6, 2e6, 1e6, 1e6]
    model = project(WALL, watts, grid)
    assert len(model) == 10
    assert model.sum() == pytest.approx(watts.sum())
    assert model[0] == pytest.approx(np.pi*(1. + 1.002)*.002*3e6)
    assert summary(model, grid)["surface_W"]["shelf"] == pytest.approx(watts[1])
    assert summary(model, grid)["peak_surface"] == "floor"
    assert summary(model, grid)["A2_m2"] == pytest.approx(watts.sum()**2 / np.sum(watts**2/areas))
    # Same perimeter, area, and seam, different poloidal orientation.
    with pytest.raises(ValueError, match="geometry"):
        project(WALL[::-1], watts[::-1], grid)


def test_exact_paired_bootstrap_pools_bins_then_reselects_minimum_on_every_draw():
    from verification.stage5_metrics import bootstrap_counts, distribution, contrast
    weights, freq = bootstrap_counts(3)
    assert sum(freq) == 27
    x = np.array([[3., 1.], [1., 3.], [3., 1.]])*1e6
    y = np.array([[2., 2.], [6., 0.], [2., 2.]])*1e6
    a = distribution(x, np.ones(2), weights)
    b = distribution(y, np.ones(2), weights)
    assert a["peak"] == pytest.approx(7/3)
    candidate = {"a": a, "b": b, "c": b, "d": b}
    regret = contrast(a, candidate, freq)
    brute = []
    for ids in product(range(3), repeat=3):
        q1 = x[list(ids)].mean(axis=0).max()
        q2 = y[list(ids)].mean(axis=0).max()
        brute.append(100*(q1/min(q1,q2)-1))
    assert regret["estimate_pct"] == 0
    assert regret["ci95_pct"] == pytest.approx(np.quantile(brute,[.025,.975]))
    assert regret["ci95_pct"][1] > 0
    assert regret["verdict_1pct"] == "unresolved"
    assert contrast(a, {"itself":a},freq)["ci95_pct"] == [0.,0.]


@pytest.mark.parametrize("pct,ci,cls,verdict",[
    (.5, [0.,1.], "supported","supported"),
    (3., [1.,5.], "unresolved","unresolved"),
    (5., [2.,8.], "contradicted","unresolved"),
    (3., [2.75,3.25], "contradicted","contradicted")])
def test_one_pp_precision_gate_is_separate_from_ci_only_decision(pct,ci,cls,verdict):
    from verification.stage5_metrics import inference
    result = inference(pct,ci)
    assert result["ci_only_1pct"] == cls
    assert result["verdict_1pct"] == verdict


def native_fixture(path):
    with h5py.File(path,"w") as f:
        f.attrs["kind"] = np.bytes_("heat_load_proxy")
        for part, p in (("wall", WALL), ("reflect_midplane",WALL[:2])):
            g = f.create_group("grid/domain/" + part)
            g["x"] = np.repeat(p[:,None,:],3,axis=1)
            u = np.array([0.,1.,360/64])
            g["domain/u"] = u
        parts = np.vstack((np.diff(WALL,axis=0),np.diff(WALL[:2],axis=0)))
        radii = np.r_[WALL[:-1,0]+WALL[1:,0],WALL[0,0]+WALL[1,0]]
        area = (np.pi*radii*np.linalg.norm(parts,axis=1))[:,None] * np.diff(u)[None,:]/360
        f["area"] = area.flatten()
        f["hload"] = np.repeat([.2/4,.2/4],5).reshape(-1)/area.flatten()/64
        f["hload"].attrs["units"] = np.bytes_("m**(-2)")
        f["area"].attrs["units"] = np.bytes_("m**2")
        for k,v in {"nparticles":4096,"n0":1e19,"T0":100.,"chi":.45,"tau":1e-7,
                    "max_history_seconds":.2,"domain_loss_fraction":.4,"censored_fraction":.1}.items():
            f["parameters/"+k] = v


def test_native_power_is_full_torus_and_includes_cap_and_losses(tmp_path):
    from verification.stage5_metrics import read_native
    file = tmp_path/"power.nc"
    native_fixture(file)
    nodes, watts, ledger = read_native(file, 1e6)
    np.testing.assert_array_equal(nodes,WALL)
    np.testing.assert_allclose(watts,1e5)
    assert ledger["wall_W"] == pytest.approx(4e5)
    assert ledger["cap_W"] == pytest.approx(1e5)
    assert ledger["domain_loss_W"] == pytest.approx(4e5)
    assert abs(ledger["closure_residual_W"]) < 1e-9
    with h5py.File(file,"r+") as f:
        f["area"][0] *= 1.02
    with pytest.raises(ValueError,match="cell area"):
        read_native(file,1e6)


def test_prediction_uses_forced_winner_and_preserves_extra_product_seed():
    from verification.stage5_metrics import comparisons, bootstrap_counts
    _,freq = bootstrap_counts(3)
    def result(role,seed,peak):
        return dict(role=role,seed=seed,wall_id=str(peak),stage4={"peak_MW_m2":peak},
                    model={"peak_MW_m2":peak},firefly={"peak_MW_m2":peak})
    peaks={"baseline":10.,"hotspot":2.,"full":2.5,"floor-only":3.,"floor-entry":6.,"floor-shelf":4.,
           "floor-connector":5.}
    roles=[result(role,30713,peak) for role,peak in peaks.items()]
    roles.extend([result("hotspot",30757,3.),result("full",30757,2.)])
    dist={r["wall_id"]:{"peak":r["firefly"]["peak_MW_m2"],"boot":np.full(len(freq),r["firefly"]["peak_MW_m2"])}
          for r in roles}
    primary,extra=comparisons(roles,"floor-shelf",dist,freq)
    assert primary["predicted_forced_vs_best_forced"]["firefly"]["estimate_pct"] == pytest.approx(100/3)
    assert primary["floor_only_vs_best_forced"]["firefly"]["estimate_pct"] == 0
    assert primary["hotspot_vs_full"]["firefly"]["estimate_pct"] == pytest.approx(-20.)
    assert extra["hotspot_vs_full"]["firefly"]["estimate_pct"] == pytest.approx(50.)
    assert "predicted_forced_vs_best_forced" not in extra


def test_stage4_checks_native_sample_winner_and_area_units(tmp_path):
    from verification.stage5_metrics import read_model,sha
    area = np.pi * (WALL[1:,0]+WALL[:-1,0]) * np.linalg.norm(np.diff(WALL,axis=0),axis=1)
    fraction=np.array([.25,.25,.25,.25])
    rms=(fraction.sum())**2 / np.sum(fraction*fraction/area)
    def save(edges):
        np.savez(tmp_path/"boundary.npz",wall=WALL,assigned_power_W=1e6,receiving_wall=WALL,cell_parent=np.arange(4),
                 cell_power=fraction,area_edges=edges,q_MW_m2=[2.,7.,3.,1.],sample_parent=np.arange(4),
                 points=.5*(WALL[1:]+WALL[:-1]))
    save(np.r_[0.,np.cumsum(area)])
    with (tmp_path/"geometry.csv").open("w") as f:
        f.write("vertex,R_m,Z_m,outgoing_region\n")
        for i,(point,region) in enumerate(zip(WALL,[*REGIONS,""])):
            f.write(f"{i},{point[0]},{point[1]},{region}\n")
    role=dict(stage4_directory=str(tmp_path), files_sha256={k:sha(tmp_path/k) for k in ("geometry.csv","boundary.npz")},
              stage4_peak_MW_m2=7.,stage4_peak_region="shelf",stage4_A2_m2=rms)
    grid, watts, original = read_model(role,WALL,1e6)
    assert original["peak_surface"] == "shelf"
    assert watts.sum() == pytest.approx(1e6)
    with pytest.raises(ValueError,match="original peak"):
        read_model(dict(role,stage4_peak_region="floor"),WALL,1e6)
    save(np.r_[0.,np.cumsum(area)]*1.01)
    role["files_sha256"]["boundary.npz"]=sha(tmp_path/"boundary.npz")
    with pytest.raises(ValueError,match="model cell area"):
        read_model(role,WALL,1e6)


def test_native_receipt_binds_watts_and_exact_case_geometry_seed(tmp_path):
    from verification.stage5_metrics import read_seed,read_native,sha
    import json
    file=tmp_path/"51031.nc"
    native_fixture(file)
    _,_,ledger=read_native(file,1e6)
    expected=dict(case_id="case",wall_id="a",source_sha256="source",mesh_sha256="mesh",
                  wall_sha256="geometry",pfc_sha256="pfc",ready_sha256="ready")
    receipt=dict(expected,seed=51031,histories=4096,ranks=4,raw_sha256=sha(file),**ledger)
    def write(value):
        file.with_suffix(".json").write_text(json.dumps(value),encoding="utf-8")
    write(receipt)
    grid=__import__("verification.stage5_metrics",fromlist=["common_rings"]).common_rings(WALL,REGIONS)
    rows, checked=read_seed(file,51031,expected,1e6,grid)
    assert rows.sum() == pytest.approx(4e5)
    assert checked["raw_sha256"] == sha(file)
    write(dict(receipt,seed=51047))
    with pytest.raises(ValueError,match="receipt identity"):
        read_seed(file,51031,expected,1e6,grid)
    write(dict(receipt,wall_W=1e6))
    with pytest.raises(ValueError,match="receipt watts"):
        read_seed(file,51031,expected,1e6,grid)
