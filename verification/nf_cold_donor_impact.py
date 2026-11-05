"""Quantify the frozen two-point donor tail on the original OMP deposition model."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from d3d_optimization.boundary_transport import AreaBoundary, wrapped_density, wrapped_mass
from d3d_optimization.physical_profile import read_physical_profile


ROOT = Path(__file__).resolve().parents[1]
DATA = Path("D:/DivertorOptimization")
OUTPUT = ROOT / "results/nf-finalist-robustness-20261106"
LABELS = {"baseline": "baseline", "OMP20": "stage18", "OMP50": "stage34"}


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def scores(power, area):
    flux = power / area
    return dict(deposited_MW=power.sum(), peak_MW_m2=flux.max(),
                hottest_ring=int(np.argmax(flux)), A2_m2=power.sum()**2 / (power @ flux))


def evaluate(path, name, donor, study, prior):
    with np.load(path, allow_pickle=False) as archive:
        z = dict(archive)
    physical = read_physical_profile(path, operating_power_W=1.e6)
    original = z["wall"]
    step = 2 * np.linalg.norm(np.diff(original, axis=0), axis=1).max()
    mesh = AreaBoundary(original, step)
    period = z["area_edges"][-1]
    np.testing.assert_allclose(mesh.total, period, rtol=1e-12)
    weights = donor["power_W"] / 1.e6
    collapsed = z["donor_region_power"].sum(axis=1)
    np.testing.assert_allclose(weights, donor["power_fraction"], rtol=2e-15)
    np.testing.assert_allclose(weights.sum(), 1., rtol=0., atol=1e-12)
    np.testing.assert_allclose(z["source_x_m"], donor["x_mid_m"], rtol=0., atol=2e-15)
    np.testing.assert_allclose(weights, collapsed, rtol=1e-11, atol=1e-18)
    centers = z["donor_area_m2"]
    widths = z["sigma_area_m2"]
    np.testing.assert_allclose(mesh.hit_area(z["hit_points"], z["hit_segment"]), centers,
                               rtol=0., atol=1e-11)
    np.testing.assert_allclose(2*np.pi*z["hit_points"][:, 0]*z["sigma_physical_m"], widths)
    temperature = z["thermal_Tt_eV"]
    if not np.isfinite(temperature).all():
        raise ValueError("incomplete two-point output")
    tagged = temperature < 5.
    withheld = weights[tagged].sum()
    np.testing.assert_allclose(withheld, prior["thermal"]["cold_Tt_below5_power_fraction"], atol=2e-14)
    edges = {"native": z["area_edges"]}
    arclength, saved, hashes = {}, {}, {}
    for mm in (2, 1, 4):
        n = round(mesh.arc[-1] / (mm / 1000))
        cuts = np.linspace(0., mesh.arc[-1], n+1)
        area_edges = mesh.area_at_arc(cuts)
        area_edges[-1] = period
        edges[mm] = area_edges
        arclength[mm] = (cuts[:-1] + cuts[1:]) / 2
        ring_file = study / (name + "_fine") / f"omp_{mm}mm_0.npz"
        with np.load(ring_file, allow_pickle=False) as rings:
            np.testing.assert_allclose(np.diff(area_edges), rings["area"], atol=1e-11, rtol=1e-9)
            np.testing.assert_allclose(arclength[mm], rings["s"], atol=1e-12)
            saved[mm] = rings["original_MW"]
        hashes[mm] = digest(ring_file)

    # Evaluate production interval integrals on the native and all three ring grids.
    lo = np.concatenate([e[:-1] for e in edges.values()])
    hi = np.concatenate([e[1:] for e in edges.values()])
    split = np.cumsum([len(e)-1 for e in edges.values()])[:-1]
    warm, cold = np.zeros_like(lo), np.zeros_like(lo)
    for mu, sigma, power, is_cold in zip(centers, widths, weights, tagged, strict=True):
        out = cold if is_cold else warm
        out += power * wrapped_mass(lo, hi, mu, sigma, period)
    warm = dict(zip(edges, np.split(warm, split), strict=True))
    cold = dict(zip(edges, np.split(cold, split), strict=True))
    native_error = warm["native"] + cold["native"] - z["cell_power"]
    native_q_error = np.max(np.abs(native_error) / physical["area_m2"])
    if native_q_error > 1.e-9:
        raise ValueError(f"{name}: saved native-cell replay differs by {native_q_error}")
    for key in edges:
        np.testing.assert_allclose(cold[key].sum(), withheld, atol=2e-12, rtol=0.)
        np.testing.assert_allclose(warm[key].sum(), weights.sum()-withheld, atol=2e-12, rtol=0.)

    # Independently evaluate both components at the archived numerical continuous maximum.
    j = int(np.argmax(z["q_MW_m2"]))
    peak_point = z["points"][j]
    a, delta = original[:-1], np.diff(original, axis=0)
    frac = np.clip(np.sum((peak_point-a)*delta, axis=1) / np.sum(delta**2, axis=1), 0., 1.)
    segment = int(np.argmin(np.linalg.norm(a + frac[:,None]*delta - peak_point, axis=1)))
    peak_area, = mesh.hit_area(peak_point[None, :], [segment])
    point_terms = np.array([w*wrapped_density(peak_area, c, s, period)
                            for c, s, w in zip(centers, widths, weights, strict=True)])
    np.testing.assert_allclose(point_terms.sum(), z["q_MW_m2"][j], rtol=1e-10)

    cases = {}
    for mm in arclength:
        area = np.diff(edges[mm])
        full = warm[mm] + cold[mm]
        all_on, all_off = scores(full, area), scores(warm[mm], area)
        tagged_q = cold[mm]/area
        full_j = all_on["hottest_ring"]
        # Independent saved profiles linearly rebin cumulative native cell powers.
        # Replay that convention exactly, then measure its difference from the exact kernel integral.
        rebinned = np.diff(np.interp(edges[mm], edges["native"], np.r_[0., np.cumsum(z["cell_power"])]))
        np.testing.assert_allclose(rebinned, saved[mm], rtol=0., atol=1.e-10)
        d_full = full @ (full/area)
        d_warm = warm[mm] @ (warm[mm]/area)
        cases[mm] = dict(rings=len(area), all_retained=all_on, all_withheld=all_off,
            tagged_at_full_peak_MW_m2=tagged_q[full_j],
            tagged_at_full_peak_pct=100*cold[mm][full_j]/full[full_j],
            max_tagged_any_ring_MW_m2=tagged_q.max(),
            tagged_only_hottest_ring=int(np.argmax(tagged_q)),
            peak_change_pct=100*(all_off["peak_MW_m2"]/all_on["peak_MW_m2"]-1),
            A2_change_pct=100*(all_off["A2_m2"]/all_on["A2_m2"]-1),
            peak_bounds_MW_m2=[all_off["peak_MW_m2"], all_on["peak_MW_m2"]],
            conservative_A2_bounds_m2=[warm[mm].sum()**2/d_full, full.sum()**2/d_warm],
            native_rebin_vs_saved_max_abs_MW=np.max(abs(rebinned-saved[mm])),
            exact_vs_saved_rebin_max_abs_MW_m2=np.max(abs(full-saved[mm])/area),
            exact_vs_saved_rebin_peak_relative=all_on["peak_MW_m2"]/np.max(saved[mm]/area)-1.)
    profile = {"s_m": arclength[2], "area_m2": np.diff(edges[2]), "untagged_MW": warm[2], "tagged_MW": cold[2]}
    return dict(npz_sha256=digest(path), archived_original_ring_sha256=hashes,
        donors=len(weights), tagged_donors=int(tagged.sum()), min_Tt_eV=temperature.min(),
        max_nt_m3=z["thermal_nt_m3"].max(), max_source_weight_abs_error=np.max(abs(weights-collapsed)),
        original_prescribed_MW=weights.sum(), withheld_at_all_off_MW=withheld, withheld_at_all_on_MW=0.,
        native_max_cell_power_error_MW=np.max(abs(native_error)), native_max_q_error_MW_m2=native_q_error,
        saved_continuous_peak=dict(full_MW_m2=z["q_MW_m2"][j], tagged_at_peak_MW_m2=point_terms[tagged].sum(),
            tagged_at_peak_pct=100*point_terms[tagged].sum()/point_terms.sum(), RZ_m=peak_point.tolist()),
        grids=cases), profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA, help="restored archive root")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--write", action="store_true", help="create the result and 2 mm powers")
    args = parser.parse_args()
    evidence = args.data / "analyses/submission-adversarial-20261105/audit.json"
    table = args.data / "studies/firefly-baseline-20261105/OMP.csv"
    source = np.genfromtxt(table, delimiter=",", names=True)
    records = read(evidence)
    report = dict(evidence_sha256=digest(evidence), original_source_csv_sha256=digest(table), cases={})
    arrays = {}
    for name, label in LABELS.items():
        meta = records["inputs"][name]
        # Archived paths refer to the original Windows drive. Accept a relocated archive.
        relative = meta["path"].replace("\\", "/").split("/studies/", 1)[1]
        path = args.data / "studies" / relative
        if digest(path) != meta["sha256"]:
            raise ValueError(f"{label}: archived digest mismatch")
        case, profile = evaluate(path, name, source, args.data/"studies/preimpact-firstwall-20261105",
                                 records["diagnostics"][name])
        report["cases"][label] = case
        arrays.update({f"{label}_{k}": value for k,value in profile.items()})
        print(f"{label} checked", flush=True)
    out = args.output / "cold-impact.json"
    profiles = args.output / "cold-rings-2mm.npz"
    payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.write:
        if profiles.exists() or out.exists():
            raise FileExistsError("choose a new output directory")
        args.output.mkdir(exist_ok=True, parents=True)
        with profiles.open("xb") as stream:
            np.savez_compressed(stream, **arrays)
        out.write_text(payload, encoding="utf-8", newline="\n")
    else:
        if json.loads(payload) != read(out):
            raise ValueError("published summary differs from recomputation")
        with np.load(profiles, allow_pickle=False) as archive:
            if set(archive.files) != set(arrays):
                raise ValueError("2 mm profile set mismatch")
            for k,v in arrays.items():
                np.testing.assert_array_equal(v, archive[k])
    print("Verified the native and 1, 2, 4 mm ledgers for all three geometries.")


if __name__ == "__main__":
    main()
