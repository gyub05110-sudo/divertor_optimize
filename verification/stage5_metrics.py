"""Audit a fixed Stage 4 panel against independent FIREFLY whole-wall heat tallies."""
import argparse
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
from itertools import combinations_with_replacement
import json
from math import factorial, prod
from pathlib import Path

import h5py
import numpy as np


SEEDS = (51031, 51047, 51059, 51071, 51089, 51109)
FORCED = ("floor-only", "floor-entry", "floor-shelf", "floor-connector")
WIDTH_M = .002
GEOMETRY_TOL_M = 2e-6
AREA_TOL_M2 = 5e-6


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def axes(points, closed=True):
    points = np.asarray(points, dtype=float)
    require(points.ndim == 2 and points.shape[1] == 2 and len(points) >= (4 if closed else 2)
            and np.isfinite(points).all() and (points[:, 0] > 0).all(), "invalid surface coordinates")
    require(not closed or np.linalg.norm(points[-1] - points[0]) < 1e-10, "open wall geometry")
    ds = np.linalg.norm(np.diff(points, axis=0), axis=1)
    require((ds > 0).all(), "zero-length geometry segment")
    return np.r_[0., np.cumsum(ds)], np.r_[0., np.cumsum(np.pi*(points[1:,0]+points[:-1,0])*ds)]


def at_arc(points, arc, area, cuts):
    idx = np.clip(np.searchsorted(arc, cuts, side="right") - 1, 0, len(arc)-2)
    delta = np.diff(points, axis=0)[idx]
    ds = cuts-arc[idx]
    frac = ds/(arc[idx+1]-arc[idx])
    rz = points[idx] + frac[:,None]*delta
    ae = area[idx] + np.pi*ds*(points[idx,0] + rz[:,0])
    return rz, ae


@dataclass
class Rings:
    wall: np.ndarray
    wall_s: np.ndarray
    cuts: np.ndarray
    area: np.ndarray
    rz: np.ndarray
    region: np.ndarray


def common_rings(wall, labels):
    wall = np.asarray(wall, dtype=float)
    arc, ae = axes(wall)
    labels = np.asarray(labels, dtype=str)
    require(labels.shape == (len(wall)-1,) and (labels != "").all(), "invalid wall surface labels")
    joints = np.r_[0, np.flatnonzero(labels[:-1] != labels[1:])+1, len(labels)]
    cuts = [0.]
    names = []
    for lo, hi in zip(joints[:-1], joints[1:]):
        count = max(1, round((arc[hi]-arc[lo])/WIDTH_M))
        cuts.extend(np.linspace(arc[lo], arc[hi], count+1)[1:])
        names.extend([labels[lo]]*count)
    cuts = np.asarray(cuts)
    rz, _ = at_arc(wall, arc, ae, .5*(cuts[:-1]+cuts[1:]))
    _, edges = at_arc(wall, arc, ae, cuts)
    return Rings(wall, arc, cuts, np.diff(edges), rz, np.asarray(names))


def project(points, power_W, grid):
    """Integrate a cellwise uniform areal density onto region-aligned physical rings."""
    points = np.asarray(points, dtype=float)
    power_W = np.asarray(power_W, dtype=float)
    arc, ae = axes(points)
    require(power_W.shape == (len(points)-1,) and np.isfinite(power_W).all()
            and (power_W >= 0).all(), "invalid cell watts")
    require(abs(arc[-1]-grid.wall_s[-1]) <= GEOMETRY_TOL_M, "wall geometry perimeter mismatch")
    # Compare both sets of vertices. Between consecutive checked offsets both polylines are linear.
    offsets = np.unique(np.r_[arc/arc[-1], grid.wall_s/grid.wall_s[-1]])
    target_rz = np.column_stack([np.interp(offsets*grid.wall_s[-1], grid.wall_s, col) for col in grid.wall.T])
    native_rz = np.column_stack([np.interp(offsets*arc[-1], arc, col) for col in points.T])
    require(np.linalg.norm(target_rz-native_rz, axis=1).max() <= GEOMETRY_TOL_M, "wall geometry path mismatch")
    _, edges = at_arc(points, arc, ae, grid.cuts/grid.cuts[-1]*arc[-1])
    edges[0], edges[-1] = ae[0], ae[-1]
    require(abs(ae[-1]-grid.area.sum()) <= AREA_TOL_M2 and
            np.max(np.abs(np.diff(edges)-grid.area)) <= AREA_TOL_M2, "geometry ring area mismatch")
    bins = np.diff(np.interp(edges, ae, np.r_[0., np.cumsum(power_W)]))
    require((bins >= 0).all() and np.isclose(bins.sum(), power_W.sum(), atol=1e-7, rtol=2e-12),
            "ring integration failed power conservation")
    return bins


def summary(power_W, grid):
    q = power_W/grid.area/1e6
    i = int(np.argmax(q))
    require(q[i] > 0, "whole-wall peak is zero")
    surfaces = np.unique(grid.region).tolist()
    return dict(peak_MW_m2=float(q[i]), peak_surface=str(grid.region[i]), peak_ring=i,
                peak_RZ_m=grid.rz[i].tolist(), peak_s_m=float(.5*(grid.cuts[i]+grid.cuts[i+1])),
                A2_m2=float(power_W.sum()**2 / np.sum(power_W**2/grid.area)), wall_W=float(power_W.sum()),
                surface_W={k:float(power_W[grid.region==k].sum()) for k in surfaces},
                surface_peak_MW_m2={k:float(q[grid.region==k].max()) for k in surfaces})


def bootstrap_counts(n):
    counts = np.array([np.bincount(row,minlength=n) for row in combinations_with_replacement(range(n),n)])
    freq = np.array([factorial(n)//prod(factorial(int(k)) for k in row) for row in counts])
    require(freq.sum() == n**n, "incomplete paired bootstrap")
    return counts/n, freq


def distribution(rows, area, weights):
    """Every row is one seed at the fixed launched watts; use the same seed weights across walls."""
    def peak_a2(watts):
        q = watts/area
        peak = q.max(axis=-1)/1e6
        require(np.isfinite(peak).all() and (peak > 0).all(), "invalid heat peak")
        return peak, watts.sum(axis=-1)**2 / np.sum(watts*q, axis=-1)
    peak, a2 = peak_a2(rows.mean(axis=0))
    boot, boot_a2 = peak_a2(weights @ rows)
    return dict(peak=float(peak), A2=float(a2), boot=boot, boot_A2=boot_a2)


def ci(values, freq):
    return np.quantile(np.repeat(values, freq), [.025,.975]).tolist()


def inference(estimate, bounds):
    low, high = bounds
    uncertainty = max(estimate-low, high-estimate)
    # A decimal 1.00% comparison can accumulate a few ulps in division/subtraction.
    ci_only = "supported" if high <= 1+1e-12 else "contradicted" if low > 1+1e-12 else "unresolved"
    precise = bool(uncertainty <= 1+1e-12)
    return dict(estimate_pct=float(estimate), ci95_pct=bounds, uncertainty_pp=float(uncertainty),
                direction="lower" if high < 0 else "higher" if low > 0 else "unresolved",
                ci_only_1pct=ci_only, precision_1pct=precise, verdict_1pct=ci_only if precise else "unresolved")


def contrast(target, references, freq):
    ref_peak = min(r["peak"] for r in references.values())
    ref_boot = np.minimum.reduce([r["boot"] for r in references.values()])
    draws = 100*(target["boot"]/ref_boot - 1)
    estimate = 100*(target["peak"]/ref_peak - 1)
    result = inference(estimate, ci(draws, freq))
    result.update(delta_MW_m2=float(target["peak"]-ref_peak),
                  delta_ci95_MW_m2=ci(target["boot"]-ref_boot, freq),
                  tied_pct=float(100*freq[target["boot"]==ref_boot].sum()/freq.sum()),
                  best_references=[k for k,r in references.items() if r["peak"]==ref_peak],
                  best_or_tied_probability={k:float(freq[r["boot"]==ref_boot].sum()/freq.sum())
                                            for k,r in references.items()})
    return result


def read_native(path, input_W):
    require(np.isfinite(input_W) and input_W > 0, "invalid input watts")
    with h5py.File(path) as file:
        require(file.attrs["kind"] == b"heat_load_proxy" and file["area"].attrs["units"] == b"m**2"
                and file["hload"].attrs["units"] == b"m**(-2)", "native quantity/units mismatch")
        domain = file["grid/domain"]
        require(set(domain) == {"wall","reflect_midplane"}, "unexpected native surfaces")
        a = np.asarray(file["area"], dtype=float)
        h = np.asarray(file["hload"], dtype=float)
        require(a.ndim == 1 and a.shape == h.shape and np.isfinite(h).all() and (h >= 0).all()
                and np.isfinite(a).all() and (a > 0).all(), "invalid native heat or cell area")
        offset = 0
        ring_W = {}
        area_errors = {}
        for name in ("wall", "reflect_midplane"):
            geom = domain[name]
            rz = np.asarray(geom["x"])
            require(rz.ndim == 3 and rz.shape[2] == 2 and np.max(abs(rz-rz[:,:1,:])) < 1e-9,
                    "native geometry is not axisymmetric")
            nodes = rz[:,0,:]
            _, edges = axes(nodes, closed=(name=="wall"))
            dphi = np.diff(geom["domain/u"][:])
            require(rz.shape[1] == len(dphi)+1 and (dphi > 0).all() and abs(dphi.sum()-360/64) < 1e-8,
                    "unexpected toroidal sector")
            sector = 360/dphi.sum()
            exact = np.diff(edges)[:,None]*dphi[None,:]/360
            cells = slice(offset, offset+exact.size)
            require(a[cells].size == exact.size, "incomplete native cell array")
            area_errors[name] = float(np.max(np.abs(a[cells].reshape(exact.shape)/exact-1)))
            require(area_errors[name] <= .0002, "native cell area or flattened surface indexing mismatch: "+name)
            ring_W[name] = (a[cells]*h[cells]).reshape(exact.shape).sum(axis=1)*sector*input_W
            offset += exact.size
            if name == "wall":
                wall_nodes = nodes
        require(offset == len(a), "unmapped native power cells")
        params = {k:float(file["parameters/"+k][()]) for k in (
            "nparticles","n0","T0","chi","tau","max_history_seconds","domain_loss_fraction","censored_fraction")}
    require(params["nparticles"] == 4096 and all(np.isfinite(v) for v in params.values())
            and all(0 <= params[k] <= 1 for k in ("domain_loss_fraction","censored_fraction")),
            "unexpected particle count or invalid native parameters")
    wall_W = float(ring_W["wall"].sum())
    cap_W = float(ring_W["reflect_midplane"].sum())
    lost = params["domain_loss_fraction"]*input_W
    censored = params["censored_fraction"]*input_W
    residual = wall_W+cap_W+lost+censored-input_W
    require(abs(residual) <= .3*input_W/1e6, "FIREFLY power ledger does not close")
    return wall_nodes, ring_W["wall"], dict(source_power_W=input_W, wall_W=wall_W, cap_W=cap_W,
        domain_loss_W=lost, censored_W=censored, closure_residual_W=residual, sector_multiplier=sector,
        cell_area_max_relative_error=area_errors, parameters=params, wall_cells=len(ring_W["wall"]))


def read_model(role, pinned_wall, input_W):
    folder = Path(role["stage4_directory"])
    for name, expected in role["files_sha256"].items():
        require(sha(folder/name) == expected, f"Stage 4 source hash mismatch: {folder/name}")
    with (folder/"geometry.csv").open(encoding="utf-8-sig",newline="") as f:
        records = list(csv.DictReader(f))
    vertices = np.asarray([[float(r["R_m"]),float(r["Z_m"])] for r in records])
    names = np.asarray([r["outgoing_region"] for r in records[:-1]])
    require([int(r["vertex"]) for r in records] == list(range(len(records))) and
            np.array_equal(vertices,pinned_wall), "Stage 4 material vertex IDs or geometry changed")
    grid = common_rings(vertices,names)
    with np.load(folder/"boundary.npz", allow_pickle=False) as raw:
        require(float(raw["assigned_power_W"]) == input_W and np.array_equal(raw["wall"],vertices),
                "Stage 4 source watts or material wall changed")
        receiver, fractions, edges = raw["receiving_wall"], raw["cell_power"], raw["area_edges"]
        s, ae = axes(receiver)
        ds = np.linalg.norm(np.diff(receiver,axis=0),axis=1)
        area = np.pi*(receiver[:-1,0]+receiver[1:,0])*ds
        require(edges.shape == ae.shape and edges[0] == 0 and
                np.allclose(np.diff(edges),area,rtol=2e-8,atol=2e-13), "model cell area units or geometry mismatch")
        parents = raw["cell_parent"]
        parent_mid = np.searchsorted(grid.wall_s, .5*(s[1:]+s[:-1]),side="right")-1
        require(np.array_equal(parents,parent_mid), "incorrect model-cell material parent")
        watts = fractions*input_W
        require(abs(watts.sum()-input_W) <= 1e-5, "reduced model fixed input W closure")
        ring_W = project(receiver,watts,grid)
        q, sample_parent, point = raw["q_MW_m2"], raw["sample_parent"], raw["points"]
        require(q.ndim == 1 and np.isfinite(q).all() and (q>=0).all() and point.shape==(len(q),2)
                and sample_parent.shape == q.shape and (sample_parent>=0).all() and (sample_parent<len(names)).all(),
                "invalid original peak samples")
        i = int(q.argmax())
        original = dict(peak_MW_m2=float(q[i]),peak_surface=str(names[sample_parent[i]]),
                        peak_RZ_m=point[i].tolist(),A2_m2=float(watts.sum()**2/np.sum(watts**2/area)),
                        wall_W=float(watts.sum()), surface_W={str(k):float(watts[names[parents]==k].sum())
                                                             for k in np.unique(names)})
        require(original["peak_MW_m2"]==role["stage4_peak_MW_m2"] and
                original["peak_surface"]==role["stage4_peak_region"],"Stage 4 original peak or winning material differs")
        require(np.isclose(original["A2_m2"],role["stage4_A2_m2"],rtol=1e-9,atol=0),
                "Stage 4 A2 does not reproduce the registered winner")
    return grid, ring_W, original


def comparisons(roles, predicted, distributions, freq):
    """Keep the designated forced-family winner and the actual product winner distinct."""
    baseline, = [r for r in roles if r["role"]=="baseline"]
    reports = []
    for seed in sorted({r["seed"] for r in roles if r["role"]!="baseline"}):
        by_role = {r["role"]:r for r in roles if r["role"]!="baseline" and r["seed"]==seed}
        by_role["baseline"] = baseline
        definitions = {"hotspot_vs_baseline":("hotspot",["baseline"]),
                       "hotspot_vs_full":("hotspot",["full"]),
                       "full_vs_baseline":("full",["baseline"])}
        if seed == 30713:
            require(predicted in FORCED, "unregistered forced-family prediction")
            definitions.update(predicted_forced_vs_best_forced=(predicted,FORCED),
                               floor_only_vs_best_forced=("floor-only",FORCED),
                               predicted_forced_vs_floor_only=(predicted,["floor-only"]))
        row = {"optimizer_seed":seed}
        for key,(name,refs) in definitions.items():
            target = by_role[name]
            comparator = {k:distributions[by_role[k]["wall_id"]] for k in refs}
            ff = contrast(distributions[target["wall_id"]],comparator,freq)
            pair = dict(target=name, references=list(refs),firefly=ff,
                        shared_reference_walls=[k for k in refs if by_role[k]["wall_id"]==target["wall_id"]])
            for model in ("model","stage4"):
                peaks = {k:by_role[k][model]["peak_MW_m2"] for k in refs}
                best = min(peaks.values())
                pair[model+"_pct"] = 100*(target[model]["peak_MW_m2"]/best-1)
                pair[model+"_best_references"] = [k for k,p in peaks.items() if p==best]
            row[key] = pair
        reports.append(row)
    return reports


def read_seed(path, seed, expected, input_W, grid):
    receipt_path = path.with_suffix(".json")
    receipt = read_json(receipt_path)
    for key,value in dict(expected,seed=seed,histories=4096,ranks=4,raw_sha256=sha(path)).items():
        require(receipt.get(key)==value, f"FIREFLY receipt identity mismatch: {key} at {path}")
    points, bank_W, ledger = read_native(path,input_W)
    for key in ("source_power_W","wall_W","cap_W","domain_loss_W","censored_W","closure_residual_W"):
        require(np.isclose(ledger[key],receipt[key],rtol=2e-12,atol=1e-6),"FIREFLY receipt watts differ: "+key)
    bins = project(points,bank_W,grid)
    return bins,dict(ledger,seed=seed,status="audited",raw_sha256=receipt["raw_sha256"],
                     receipt_sha256=sha(receipt_path),wall_above_midplane_W=receipt.get("wall_above_midplane_W"))


def write_profile(path, grid, model_W, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".pending")
    mean = rows.mean(axis=0)
    with tmp.open("wb") as raw, gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as zipped:
        with io.TextIOWrapper(zipped,encoding="utf-8",newline="") as stream:
            wr = csv.writer(stream,lineterminator="\n")
            wr.writerow(["s_start_m","s_end_m","R_mid_m","Z_mid_m","region","area_m2","model_W",
                         "firefly_W","model_MW_m2","firefly_MW_m2",*[f"seed_{k}_W" for k in SEEDS]])
            wr.writerows(zip(grid.cuts[:-1],grid.cuts[1:],*grid.rz.T,grid.region,grid.area,model_W,
                             mean,model_W/grid.area/1e6,mean/grid.area/1e6,*rows))
    tmp.replace(path)


def pooled_case_wall(study, case, role, grid, model_W, expected, weights, freq):
    cid,wall = case["case_id"],role["wall_id"]
    source = case["reference_SI"]["power_W"]
    power_rows, ledgers = [], []
    for seed in SEEDS:
        path = study/"native/results"/cid/wall/(str(seed)+".nc")
        watts, ledger = read_seed(path,seed,expected,source,grid)
        power_rows.append(watts)
        ledgers.append(ledger)
    rows = np.asarray(power_rows)
    dist = distribution(rows,grid.area,weights)
    summary_ff = summary(rows.mean(axis=0),grid)
    require(summary_ff["peak_MW_m2"]==dist["peak"],"pooled peak mismatch")
    summary_ff.update(peak_ci95_MW_m2=ci(dist["boot"],freq),A2_ci95_m2=ci(dist["boot_A2"],freq))
    pooled_q = (weights @ rows)/grid.area
    best_q = np.max(pooled_q,axis=1)
    summary_ff["peak_surface_best_or_tied_probability"] = {
        str(k):float(freq[pooled_q[:,grid.region==k].max(axis=1)==best_q].sum()/freq.sum()) for k in np.unique(grid.region)}
    file = study/"profiles"/cid/(wall+".csv.gz")
    write_profile(file,grid,model_W,rows)
    audit = dict(wall_id=wall,source_power_W=source,seeds=ledgers,ring_count=len(grid.area),
                 min_ring_width_m=float(np.diff(grid.cuts).min()),max_ring_width_m=float(np.diff(grid.cuts).max()),
                 total_wall_area_m2=float(grid.area.sum()),profile=file.relative_to(study).as_posix(),
                 profile_sha256=sha(file),input_sha256=expected,
                 pooled_ledger_W={k:float(np.mean([r[k] for r in ledgers])) for k in
                                  ("wall_W","cap_W","domain_loss_W","censored_W","closure_residual_W")})
    return summary_ff,dist,audit


def analyze_case(study, case, panel_digest, weights, freq, model_only=False):
    cid = case["case_id"]
    inp = study/"inputs"/cid
    require(sha(inp/"case.json")==panel_digest and read_json(inp/"case.json")==case, "case digest/record differs")
    for file,digest in case["files_sha256"].items():
        require(sha(inp/file)==digest, "case input hash mismatch: "+str(inp/file))
    native = study/"native/cases"/cid
    locked = {}
    if not model_only:
        locked = dict(case_id=cid,source_sha256=case["files_sha256"]["OMP.source"],
                      mesh_sha256=sha(native/"mesh/mmesh.nc"), ready_sha256=sha(native/"ready.json"))
    roles, geometries, distributions, ff, checked = [], {}, {}, {}, {}
    for role in case["roles"]:
        w = role["wall_id"]
        pinned_wall = np.loadtxt(inp/role["wall_file"])
        require(hashlib.sha256(np.ascontiguousarray(pinned_wall,dtype="<f8").tobytes()).hexdigest()==w,
                "native material coordinates differ from registered wall")
        grid, model_W, original = read_model(role,pinned_wall,case["reference_SI"]["power_W"])
        reduced = summary(model_W,grid)
        entry = dict(role=role["role"],seed=role["seed"],wall_id=w,candidate_id=role["candidate_id"],
                     stage4=original,model=reduced)
        roles.append(entry)
        if w in checked:
            area, labels, previous_W = checked[w]
            require(np.array_equal(area,grid.area) and np.array_equal(labels,grid.region) and
                    np.allclose(previous_W,model_W,rtol=1e-11,atol=1e-6), "aliased wall differs in region or model profile")
        else:
            checked[w] = (grid.area,grid.region,model_W)
        if model_only:
            continue
        if w not in ff:
            expected = dict(locked,wall_id=w,wall_sha256=case["files_sha256"][role["wall_file"]],
                            pfc_sha256=sha(native/(w+".pfc.nc")))
            ff[w],distributions[w],geometries[w] = pooled_case_wall(study,case,role,grid,model_W,expected,weights,freq)
        entry["firefly"] = ff[w]
    if model_only:
        return dict(case_id=cid,role_count=len(roles),wall_count=len(checked))
    baseline, = [r for r in roles if r["role"]=="baseline"]
    require(case["baseline_peak_region"]==baseline["stage4"]["peak_surface"], "original baseline hotspot mismatch")
    hotspots = dict(predicted_region=case["predicted_region"],stage4=case["baseline_peak_region"],
                    model=baseline["model"]["peak_surface"],firefly=baseline["firefly"]["peak_surface"])
    hotspots["stage4_matches_firefly"] = hotspots["stage4"] == hotspots["firefly"]
    hotspots["model_matches_firefly"] = hotspots["model"] == hotspots["firefly"]
    return dict(case_id=cid,predicted_family=case["predicted_family"],baseline_hotspot=hotspots,
                roles=roles,geometries=geometries,comparisons=comparisons(roles,case["predicted_family"],distributions,freq))


def run(study, model_only=False):
    from verification.stage5_campaign import verify_panel
    registration_path = study/"wall-registration.json"
    require(sha(registration_path)=="d0627f34407ee7cb3fd054a6c8eecef17b1e20026c28ce068fdc3f0d224b7102",
            "fixed wall registration changed")
    plan = read_json(registration_path)
    panel_file = study/"inputs/panel.json"
    require(sha(panel_file)=="d78c153199e79ed7462bbd09c5f17cc0ea88d635233eaa240166542908c794a1",
            "frozen native input panel changed")
    panel = read_json(panel_file)
    verify_panel(plan["roles"],panel)
    require(list(SEEDS)==plan["seeds"] and plan["primary_ring_m"]==WIDTH_M and
            plan["histories_per_seed"]==4096 and plan["mpi_ranks"]==4,"metric protocol differs from registration")
    weights,freq = bootstrap_counts(len(SEEDS))
    cases = []
    for case in panel["cases"]:
        cid = case["case_id"]
        require(plan["predictions"][cid] == case["predicted_family"], "changed predicted family")
        cases.append(analyze_case(study,case,panel["case_record_sha256"][cid],weights,freq,model_only))
        print(cid,len(case["roles"]),"roles audited",flush=True)
    if model_only:
        return cases
    result = dict(schema="nf-stage5-whole-wall-metrics.v1",panel_sha256=sha(panel_file),
                  registration_sha256=sha(registration_path),analyzer_sha256=sha(__file__),native_seeds=SEEDS,
                  unique_bootstrap_counts=len(freq),ordered_bootstrap_draws=int(freq.sum()),ring_width_m=WIDTH_M,
                  eligible_files=len(cases),original_registration_files=24,role_count=sum(len(c["roles"]) for c in cases),
                  distinct_case_walls=sum(len(c["geometries"]) for c in cases),cases=cases)
    temp = study/"metrics.json.pending"
    temp.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n", encoding="utf-8")
    temp.replace(study/"metrics.json")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path)
    parser.add_argument("--model-check",action="store_true",help="audit fixed model profiles without opening FIREFLY output")
    args = parser.parse_args()
    run(args.study.resolve(),args.model_check)


if __name__ == "__main__":
    main()
