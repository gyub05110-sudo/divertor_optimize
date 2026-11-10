"""Material, launch and heat-ledger checks for the sealed Stage 5 native FIREFLY runs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import h5py
import numpy as np
from shapely.geometry import LineString, MultiPoint, Polygon

HISTORICAL = Path("D:/DivertorOptimization/studies/firefly-baseline-20261105")
STUDY = HISTORICAL.parent / "nf-stage5-firefly-20261110"
HERE = Path(__file__).resolve().parent
SEEDS = (51031, 51047, 51059, 51071, 51089, 51109)
HISTORIES = 4096
RANKS = 4
TAU = 4 / (2 * (2000 / 1.602176634e-19 * 100**2.5 / 1e19))
PINNED = {
    "install/lib/libFIREFLY.so": "d5ebfe0dc123f5bc44ad78a524f93ca71b7c74c48c9f3f14cb762a06bf5d6a68",
    "install/lib/python/firefly/f2py.cpython-313-x86_64-linux-gnu.so":
        "d30420eee99eb9d1442ab64d768c76f72e67a9629d07e6eab2f2b24bbeb27aea",
}


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".pending")
    temp.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temp, path)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def enclose_walls(walls):
    for wall in walls:
        require(np.isfinite(wall).all() and np.linalg.norm(wall[0]-wall[-1]) < 1e-10
                and Polygon(wall).is_valid, "invalid material wall")
    casing = MultiPoint(np.concatenate(walls)).convex_hull.buffer(.012, quad_segs=4)
    clearance = min(casing.boundary.distance(LineString(w)) for w in walls)
    require(clearance > .010 and all(casing.contains(Polygon(w)) for w in walls), "casing clearance")
    return np.asarray(casing.exterior.coords), clearance


def audit_source(inp):
    donors = np.atleast_1d(np.genfromtxt(inp / "OMP.csv", names=True, delimiter=","))
    path = inp / "OMP.source"
    rows = np.loadtxt(path, skiprows=1, ndmin=2)
    with path.open() as stream:
        require(int(stream.readline()) == len(donors) == len(rows), "source count")
    power = donors["power_W"]
    require(np.isfinite(rows).all() and np.all(power > 0) and
            abs(power.sum()-1e6) < 1e-5, "source power")
    require(rows[0,2] > 0 and rows[-1,2] == 1 and np.all(np.diff(rows[:,2]) >= 0) and
            np.max(np.abs(rows[:, 2]-np.cumsum(power)/1e6)) < 1e-12 and
            np.max(np.abs(donors["power_fraction"] - power/1e6)) < 1e-13, "source CDF")
    require(np.max(np.abs(rows[:, :2]-np.c_[donors["R_m"], donors["Z_m"]])) < 1e-12 and
            np.ptp(rows[:, 1]) < 1e-10 and np.all(np.diff(rows[:, 0]) > 0), "source locations")
    require(np.unique(donors["sign_along_B"]).size == 1 and
            abs(donors["sign_along_B"][0]) == 1, "source direction")
    return dict(donors=len(rows), source_power_W=float(power.sum()), z_m=float(rows[0,1]),
                first_R_m=float(rows[0,0]), direction=int(donors["sign_along_B"][0]))


def arclength(p):
    return np.r_[0., np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]


def geometry_error(reference, actual):
    a, b = arclength(reference), arclength(actual)
    s = np.union1d(a, b)
    return float(np.max(np.hypot(*(np.interp(s, a, reference[:, i]) -
                                   np.interp(s, b, actual[:, i]) for i in (0, 1)))))


def audit_tally(path, reference_wall, histories=HISTORIES):
    with h5py.File(path) as f:
        domains = f["grid/domain"]
        require(list(domains) == ["wall", "reflect_midplane"], "native tally domains")
        require(f.attrs["kind"] == b"heat_load_proxy" and f["hload"].attrs["units"] == b"m**(-2)",
                "native heat-load units")
        par = f["parameters"]
        require(par["nparticles"][()] == histories, "native histories")
        for key, val in dict(n0=1e19, T0=100., chi=.45, tau=TAU, max_history_seconds=.2).items():
            require(np.isclose(par[key][()], val, rtol=2e-7, atol=0), f"native parameter {key}")
        area, heat = f["area"][:], f["hload"][:]
        require(np.all(np.isfinite(heat)) and np.all(heat >= 0) and np.all(area > 0), "invalid heat/area")
        offset = 0
        total = {}
        max_area_error = 0.
        for name in ("wall", "reflect_midplane"):
            surface = domains[name]
            coords = surface["x"][:]
            p, dphi = coords[:, 0, :], np.radians(np.diff(surface["domain/u"][:]))
            require(np.max(np.abs(coords-p[:,None,:])) < 1e-9 and np.all(dphi > 0) and
                    abs(dphi.sum()*64-2*np.pi) < 1e-9, "axisymmetric sector")
            if name == "wall":
                wall_error = geometry_error(reference_wall, p)
                require(wall_error < 1e-8, f"material geometry differs by {wall_error:g} m")
            exact = np.linalg.norm(np.diff(p, axis=0), axis=1)[:,None] * (
                .5*(p[:-1,0]+p[1:,0]))[:,None] * dphi
            end = offset + exact.size
            an = area[offset:end].reshape(exact.shape)
            max_area_error = max(max_area_error, float(np.max(np.abs(an/exact-1))))
            total[name] = float((heat[offset:end] * area[offset:end]).sum()*64*1e6)
            if name == "wall":
                wall_cell_W = (heat[offset:end]*area[offset:end]).reshape(exact.shape).sum(axis=1)*64*1e6
            offset = end
        require(offset == len(area) == len(heat) and max_area_error < .0002, "native cell areas")
        loss = float(par["domain_loss_fraction"][()])*1e6
        censor = float(par["censored_fraction"][()])*1e6
        require(min(loss,censor) >= 0 and total["reflect_midplane"] == 0, "native loss or cap")
        residual = total["wall"]+total["reflect_midplane"]+loss+censor-1e6
        require(abs(residual) < .3, f"native closure {residual:g} W")
        zcut = float(domains["reflect_midplane/x"][0,0,1])
        wp = domains["wall/x"][:,0,1]
        above = float(wall_cell_W[(wp[:-1]+wp[1:])*.5 > zcut].sum())
    return dict(source_power_W=1e6, wall_W=total["wall"], cap_W=total["reflect_midplane"],
                domain_loss_W=loss, censored_W=censor, closure_residual_W=residual,
                wall_above_midplane_W=above, max_wall_distance_m=wall_error,
                max_cell_area_relative_error=max_area_error)


def wsl_path(path):
    p = Path(path).resolve().as_posix()
    return "/mnt/" + p[0].lower() + p[2:]


def run_native(study, phase, args, log, ranks=1):
    log.parent.mkdir(parents=True, exist_ok=True)
    command = ["wsl", "-d", "Ubuntu", "-u", "root", "--", "bash", wsl_path(HERE/"stage5_native_enter.sh"),
               wsl_path(study), wsl_path(HERE)]
    if phase == "heat":
        case, wall, seed, _ = args
        command += ["flock", "-n", f"/stage5/native/results/{case}/{wall}/{seed}.linux.lock"]
    if ranks > 1:
        command += ["mpirun", "--allow-run-as-root", "--oversubscribe", "--bind-to", "none", "-np", str(ranks)]
    command += ["python3", "/stage5code/stage5_native_worker.py", phase, *map(str,args)]
    with log.open("wb") as stream:
        rc = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT).returncode
    require(rc == 0, f"native {phase} exited {rc}, inspect {log}")


def inputs_for(study, case):
    inp = study / "inputs" / case
    record = read(inp / "case.json")
    require(record["case_id"] == case and record["time_ms"] == 1400, "case/time")
    for name, signature in record["files_sha256"].items():
        require(digest(inp/name) == signature, f"case input changed {name}")
    registration = study / "wall-registration.json"
    require(digest(registration) == "d0627f34407ee7cb3fd054a6c8eecef17b1e20026c28ce068fdc3f0d224b7102", "prospective thermal job registration")
    expected = read(registration)
    walls = list(dict.fromkeys(r["wall_id"] for r in record["roles"]))
    require(set(walls) == {r["wall_id"] for r in expected["roles"] if r["case_id"] == case},
            "wall set differs from registered cohort")
    require(all(digest(HISTORICAL/k)==v for k,v in PINNED.items()), "native executable pin")
    return inp, record, walls


def prepare_geometry(study, case):
    inp, record, walls = inputs_for(study, case)
    source = audit_source(inp)
    base = study/"native/cases"/case
    model, mesh = base/"model", base/"mesh"
    for d in (model, mesh):
        d.mkdir(parents=True, exist_ok=True)
    for name in (".bfield", ".boundary", ".equi2d", "fixed.geqdsk"):
        shutil.copyfile(inp/"model"/name, model/name)
    casing, clearance = enclose_walls([np.loadtxt(inp/"walls"/(wall+".wall")) for wall in walls])
    np.savetxt(model/"wall.txt", casing, fmt="%.16e")
    # A closed surface ~8 cm inboard of the OMP separatrix seeds the core exclusion,
    # matching the historical inner boundary's offset without borrowing its R or Z.
    r = record["topology"]["omp_separatrix_RZ_m"][0] - .080
    z = source["z_m"]
    ctr = (f"[FLARE]\nmodel: /stage5/native/cases/{case}/model\n[mmesh_parameters]\n"
           'layout: "unstructured"\nsymmetry: 64\nstellarator_symmetry: False\n'
           f'pcoordinates: "cylindrical"\np1 = {r:.16e}, {z:.16e}, 0.0\n'
           f'p2 = {r+.0015:.16e}, {z:.16e}, 0.0\nnt = 8\nnp = 1000\ndelta_r = 0.0015\n')
    (mesh/"mmesh.ctr").write_text(ctr, newline="\n")
    return dict(case_id=case, wall_ids=walls, source=source, minimum_casing_clearance_m=clearance,
                case_sha256=digest(inp/"case.json"), mesh_control_sha256=digest(mesh/"mmesh.ctr"),
                casing_sha256=digest(model/"wall.txt"), geqdsk_sha256=digest(model/"fixed.geqdsk"))


def prepare_field_mesh(study, case):
    settings = prepare_geometry(study, case)
    inp = study/"inputs"/case
    base = study/"native/cases"/case
    receipt = base/"mesh.json"
    if receipt.is_file():
        old = read(receipt)
        require(old["settings"] == settings and old["mesh_sha256"] == digest(base/"mesh/mmesh.nc"),
                "previous common mesh differs")
        return old
    run_native(study, "field", [case], base/"field.log")
    measured = np.loadtxt(base/"field_native.txt")
    reference = np.loadtxt(inp/"field_probes.txt")
    require(measured.shape == reference.shape and np.max(abs(measured[:,:2]-reference[:,:2])) < 1e-12,
            "FLARE field probe coordinates")
    delta = np.max(abs(measured[:,2:]-reference[:,2:]), axis=0)
    stats = dict(count=len(reference), max_abs_error_B_T=delta[:3].tolist(),
                 max_abs_error_psiN=float(delta[3]))
    save(base/"field_audit.json", stats)
    require(np.max(delta[:3]) < 3e-5 and delta[3] < 1e-4, f"signed magnetic field/psi mismatch {stats}")
    run_native(study, "mesh", [case], base/"mesh/build.log")
    mesh = base/"mesh/mmesh.nc"
    with h5py.File(mesh) as f:
        require(len(f.keys()) > 0, "empty native common mesh")
    result = dict(settings=settings, field_audit=stats, mesh_sha256=digest(mesh),
                  mesh_bytes=mesh.stat().st_size)
    save(receipt, result)
    return result


def strike_audit(study, case, wall):
    base, inp = study/"native/cases"/case, study/"inputs"/case
    reg = read(inp/"case.json")
    d = np.genfromtxt(inp/"OMP.csv", delimiter=",", names=True)
    direct = np.loadtxt(inp/"walls"/(wall+".hits"))
    raw = np.loadtxt(base/(wall+".strikes.dat"))
    require(raw.shape == (len(d), 7) and direct.shape == (len(d), 3), "strike table shape")
    direction = d["sign_along_B"] * np.sign(reg["topology"]["F_T_m"])
    require(np.all(direction == direction[0]) and abs(direction[0]) == 1, "native toroidal direction")
    col = int(direction[0] > 0)
    hit, length = raw[:,3+2*col:5+2*col], raw[:,col]
    captured = np.isfinite(hit).all(axis=1) & (hit[:,0] > 0)
    lower = captured & (hit[:,1] < d["Z_m"])
    require(np.all(direct[:,2]>0) and np.all(direct[:,1]<d["Z_m"]), "stage4 lower leg")
    w = d["power_W"]
    err = np.linalg.norm(hit[lower] - direct[lower,:2], axis=1)
    require(err.size, "no lower strikes")
    order = np.argsort(err)
    quantiles = np.interp([.50,.95], np.cumsum(w[lower][order])/w[lower].sum(), err[order])
    relative_L = float(np.average(np.abs(length[lower]/direct[lower,2]-1), weights=w[lower]))
    # Native ierr=1 means invalid source coordinates, ierr=2 a toroidal
    # surface-map error in either direction. Lower and opposite hits are separate.
    mapped = (raw[:,2].astype(int) & 1) == 0
    stats = dict(toroidal_direction=int(direction[0]), donors=len(d), unmapped=int((~mapped).sum()),
                 captured_lower_W=float(w[lower].sum()), uncaptured_W=float(w[~captured].sum()),
                 above_midplane_W=float(w[captured & ~lower].sum()), lower_count=int(lower.sum()),
                 weighted_median_hit_error_m=float(quantiles[0]), weighted_p95_hit_error_m=float(quantiles[1]),
                 max_hit_error_m=float(err.max()), weighted_mean_abs_relative_length_error=relative_L,
                 strikes_sha256=digest(base/(wall+".strikes.dat")), reference_sha256=digest(inp/"walls"/(wall+".hits")))
    stats["pass"] = bool(np.all(mapped) and w[lower].sum()>=.99e6 and quantiles[1]<=.005 and relative_L<=.05)
    return stats


def prepare_walls(study, case):
    inp, record, walls = inputs_for(study, case)
    base = study/"native/cases"/case
    common = read(base/"mesh.json")
    require(common["settings"]["case_sha256"] == digest(inp/"case.json")
            and digest(base/"mesh/mmesh.nc") == common["mesh_sha256"], "case mesh binding")
    assets = {}
    for wall in walls:
        pfc = base/(wall+".pfc.nc")
        open_wall = base/(wall+".open.nc")
        diag = base/(wall+".strikes.json")
        if not diag.is_file():
            run_native(study, "pfc", [case,wall], base/(wall+".pfc.log"))
            with h5py.File(pfc) as f, h5py.File(open_wall) as op:
                source_wall = np.loadtxt(inp/"walls"/(wall+".wall"))
                require(np.array_equal(f["wall/rz"][:],source_wall) and
                        np.array_equal(op["wall/rz"][:],source_wall), "changed native PFC vertices")
                cut = f["reflect_midplane/rz"][:]
                require(np.all(cut[:,1] == common["settings"]["source"]["z_m"])
                        and min(source_wall[:,0]) > cut[0,0] and max(source_wall[:,0]) < cut[-1,0],
                        "incomplete reflecting plane")
            run_native(study, "trace", [case,wall], base/(wall+".trace.log"))
            save(diag, strike_audit(study, case, wall))
        stats = read(diag)
        require(stats == strike_audit(study, case, wall), "strike receipt differs")
        assets[wall] = dict(pfc_sha256=digest(pfc), open_pfc_sha256=digest(open_wall),
                            wall_sha256=digest(inp/"walls"/(wall+".wall")), strikes=stats)
    require(all(s["strikes"]["pass"] for s in assets.values()), "lower-leg strike checks failed: " +
            ", ".join(k for k,v in assets.items() if not v["strikes"]["pass"]))
    result = dict(case_id=case, mesh_sha256=common["mesh_sha256"], source_sha256=digest(inp/"OMP.source"),
                  source_csv_sha256=digest(inp/"OMP.csv"), case_sha256=digest(inp/"case.json"),
                  mesh_receipt_sha256=digest(base/"mesh.json"), walls=assets, field_audit=common["field_audit"],
                  material_inputs=record["files_sha256"], native_libraries=PINNED)
    save(base/"ready.json", result)
    return result


def active_transport(study):
    marker = study/"native-physics/active-library.txt"
    if not marker.is_file():
        return dict(active_transport_library="/work/install/lib/libFIREFLY.so",
                    active_transport_sha256=PINNED["install/lib/libFIREFLY.so"])
    relative = marker.read_text(encoding="ascii").strip()
    parts = relative.split("/")
    require(len(parts)==3 and parts[0]=="native-physics" and parts[2]=="libFIREFLY.so"
            and parts[1] and all(c.isalnum() or c=="-" for c in parts[1]), "thermal library path")
    mapped = "/stage5/"+relative
    env = read(study/"execution-environment.json")
    signature = digest(study/relative)
    require(env["active_transport_library"]==mapped and env["native_files_sha256"].get(mapped)==signature,
            "changed active thermal library")
    return dict(active_transport_library=mapped, active_transport_sha256=signature)


def context(study, case, wall, seed):
    registration = study/"wall-registration.json"
    require(digest(registration) == "d0627f34407ee7cb3fd054a6c8eecef17b1e20026c28ce068fdc3f0d224b7102"
            and [case,wall,seed] in read(registration)["jobs"], "job absent from seal")
    inp, base = study/"inputs"/case, study/"native/cases"/case
    ready_path = base/"ready.json"
    ready = read(ready_path)
    plan = study/"execution.json"
    plan_sha = digest(plan)
    execution = read(plan)
    require(plan_sha == (study/"execution.sha256").read_text().strip() and
            execution["files"][ready_path.relative_to(study).as_posix()] == digest(ready_path) and
            execution["registration_sha256"] == digest(registration), "native execution seal")
    for worker in ("stage5_native.py","stage5_native_worker.py","stage5_native_enter.sh"):
        require(execution["code"]["verification/"+worker] == digest(HERE/worker), "native driver differs from seal")
    paths = dict(mesh_sha256=base/"mesh/mmesh.nc", pfc_sha256=base/(wall+".pfc.nc"),
                 source_sha256=inp/"OMP.source", wall_sha256=inp/"walls"/(wall+".wall"))
    hashes = {key:digest(path) for key,path in paths.items()}
    require(hashes["mesh_sha256"] == ready["mesh_sha256"] and hashes["pfc_sha256"] == ready["walls"][wall]["pfc_sha256"]
            and hashes["source_sha256"] == ready["source_sha256"] and hashes["wall_sha256"] == ready["walls"][wall]["wall_sha256"]
            and ready["walls"][wall]["strikes"]["pass"] and ready["case_sha256"] == digest(inp/"case.json"),
            "changed native inputs after preflight")
    require(all(digest(HISTORICAL/k)==v for k,v in PINNED.items()), "changed native runtime")
    return dict(case_id=case, wall_id=wall, seed=seed, histories=HISTORIES, ranks=RANKS,
                seed_scheme="FIREFLY init_workspace", source_power_W=1e6, n0=1e19, T0=100.,
                chi=.45, tau=TAU, max_history_seconds=.2, reflect_midplane=True, symmetry=64,
                ready_sha256=digest(ready_path), registration_sha256=digest(registration),
                preparation_libraries=PINNED, execution_sha256=plan_sha, **active_transport(study), **hashes)


def final_evidence(path, final, request, wall, raw):
    runtime = read(path/"native.runtime.json")
    require(runtime["ranks"] == RANKS and runtime["cpu_seconds"] > 0 and runtime["seconds"] > 0 and
            [p for p in runtime["library_mappings"] if p.endswith("/libFIREFLY.so")]==[request["active_transport_library"]]
            and runtime["active_transport_sha256"]==request["active_transport_sha256"], "compiled runtime receipt")
    ledger = audit_tally(raw, wall)
    require(ledger["wall_above_midplane_W"] < .3, "upper-leg contamination")
    signature = digest(raw)
    return {**request, **ledger, "raw_sha256": signature, "log_sha256": digest(path/"run.log"),
            "raw_file": final.name, "attempt": path.name, "runtime": runtime,
            "runtime_sha256": digest(path/"native.runtime.json")}


def publish_job(study, case, wall, seed, run):
    import msvcrt
    request = context(study, case, wall, seed)
    folder = study/"native/results"/case/wall
    require(run or folder.is_dir(), "native job incomplete")
    folder.mkdir(parents=True, exist_ok=True)
    final = folder/(str(seed)+".nc")
    receipt, log = final.with_suffix(".json"), final.with_suffix(".log")
    reference = np.loadtxt(study/"inputs"/case/"walls"/(wall+".wall"))
    with final.with_suffix(".lock").open("a+b") as handle:
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        attempts = folder/"attempts"/str(seed)
        attempts.mkdir(parents=True, exist_ok=True)
        if receipt.is_file():
            old = read(receipt)
            path = attempts/old["attempt"]
            evidence = final_evidence(path, final, read(path/"request.json"), reference, final)
            require(old == evidence and digest(log) == old["log_sha256"]
                    and all(old[k] == v for k,v in request.items()), "completed native job changed")
            return old
        # An audited attempt whose publication was interrupted can be finished exactly.
        # Incomplete attempts survive under their original names.
        for path in sorted(attempts.iterdir()):
            marker = path/"publish.json"
            if not marker.is_file():
                continue
            old = read(marker)
            if not all(old[k] == v for k,v in request.items()):
                continue
            raw = path/"native.nc" if (path/"native.nc").is_file() else final
            matched = final_evidence(path, final, read(path/"request.json"), reference, raw)
            require(matched == old, "attempt publication record differs")
            if raw != final:
                os.replace(raw, final)
            shutil.copyfile(path/"run.log", log)
            save(receipt, matched)
            return matched
        require(not final.exists(), "orphan native output without audit marker")
        require(run, "native job incomplete")
        path = attempts/f"{len(list(attempts.iterdir()))+1:03}"
        path.mkdir()
        request["worker_sha256"] = digest(HERE/"stage5_native_worker.py")
        request["runner_sha256"] = digest(HERE/"stage5_native.py")
        request["entry_sha256"] = digest(HERE/"stage5_native_enter.sh")
        save(path/"request.json", request)
        run_native(study,"heat",[case,wall,str(seed),path.name],path/"run.log",ranks=RANKS)
        raw = path/"native.nc"
        evidence = final_evidence(path, final, request, reference, raw)
        save(path/"publish.json", evidence)
        os.replace(raw, final)
        shutil.copyfile(path/"run.log", log)
        save(receipt, evidence)
        return evidence


def historical_smoke(study):
    folder = study/"native/historical_smoke"
    folder.mkdir(parents=True, exist_ok=True)
    run_native(study, "smoke", [], folder/"runtime.log", ranks=RANKS)
    raw = folder/"historical.pending.nc"
    historical = HISTORICAL/"blind_screen/results/ref_OMP_32027.nc"
    ledger = audit_tally(raw, np.loadtxt(HISTORICAL/"blind_screen/walls/ref.wall"), 1024)
    with h5py.File(raw) as a, h5py.File(historical) as b:
        maximum = float(np.max(abs(a["hload"][:]-b["hload"][:])))
    require(maximum == 0, "historical native snapshot mismatch")
    result = dict(**ledger, historical_raw_sha256=digest(historical), raw_sha256=digest(raw),
                  max_hload_discrepancy=maximum, runtime=read(folder/"historical.pending.runtime.json"))
    save(folder/"audit.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("smoke","mesh","walls","prepare","heat","verify"))
    parser.add_argument("--study", type=Path, default=STUDY)
    parser.add_argument("--case")
    parser.add_argument("--wall")
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()
    if args.phase == "smoke":
        result = historical_smoke(args.study)
    elif args.phase in ("heat", "verify"):
        result = publish_job(args.study,args.case,args.wall,args.seed,args.phase=="heat")
    else:
        if args.phase in ("mesh","prepare"):
            result = prepare_field_mesh(args.study,args.case)
        if args.phase in ("walls","prepare"):
            result = prepare_walls(args.study,args.case)
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
