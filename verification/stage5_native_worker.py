"""Inside the retained WSL image, call the pinned, compiled FIREFLY and FLARE APIs."""
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np
from mpi4py import MPI

ROOT = Path("/stage5")
NATIVE = ROOT / "native"
OLD = Path("/work")
COMM = MPI.COMM_WORLD


def field(case):
    from flare import model
    from flare.analysis import bfield, equi2d
    from moose.grids import R3grid
    inp = ROOT / "inputs" / case
    mesh = NATIVE / "cases" / case
    model.load(str(mesh / "model"))
    q = np.loadtxt(inp / "field_probes.txt")
    field = bfield.eval(np.c_[q[:, :2], np.zeros(len(q))])
    direct = equi2d.psiN(q[:,0], q[:,1])
    np.savetxt(mesh / "field_native.txt", np.c_[q[:,:2], field, direct], fmt="%.16e")
    donor = np.loadtxt(inp / "OMP.source", skiprows=1)
    R3grid.rzmesh(donor[:,0], [donor[0,1]], 0.).savetxt(str(mesh / "donors.grid"))


def pfc(case, wall):
    from moose.geometry import Axisurf, Targets3d
    inp = ROOT / "inputs" / case
    dest = NATIVE / "cases" / case
    geom = np.loadtxt(inp / "walls" / (wall + ".wall"))
    casing = np.loadtxt(dest / "model/wall.txt")
    z = np.loadtxt(inp / "OMP.source", skiprows=1, max_rows=1)[1]
    surface = Axisurf(geom[:,0], geom[:,1], units="m", label="physical wall")
    cut = Axisurf(np.array([casing[:,0].min(),casing[:,0].max()]), np.array([z,z]),
                  units="m", label="reflecting lower-leg midplane")
    Targets3d({"wall": surface, "reflect_midplane": cut}).savenc(str(dest / (wall + ".pfc.nc")))
    Targets3d({"wall": surface}).savenc(str(dest / (wall + ".open.nc")))


def trace(case, wall):
    from firefly import geometry, tasks
    base = NATIVE / "cases" / case
    geometry.init_workspace(str(base / "mesh/mmesh.nc"), seed=51031)
    geometry.set_pfc(str(base / (wall + ".open.nc")), dphi=2., dl=.002)
    tasks.connection_length(str(base / "donors.grid"), max_lc=700., output=str(base / (wall + ".strikes.dat")))


def solve(mesh, wall, source, dest, seed, histories):
    from firefly import geometry, tasks
    if COMM.size != 4:
        raise ValueError("native run requires four MPI ranks")
    os.chdir(dest.parent)
    start = time.monotonic()
    geometry.init_workspace(str(mesh), seed=int(seed))
    geometry.set_pfc(str(wall), dphi=2., dl=.002)
    os.environ["FIREFLY_MIDPLANE_SOURCE"] = str(source)
    os.environ["FIREFLY_HISTORY_SECONDS"] = "0.2"
    tau = 4/(2*(2000/1.602176634e-19*100**2.5/1e19))
    tasks.heat_load_proxy(1e19, 100., .45, histories, pload=-1., tau=tau, output=str(dest))
    COMM.Barrier()
    usage = resource.getrusage(resource.RUSAGE_SELF)
    cpu = COMM.gather(usage.ru_utime + usage.ru_stime, root=0)
    if COMM.rank == 0:
        mapped = sorted({line.split()[-1] for line in Path("/proc/self/maps").read_text().splitlines()
                         if any(lib in line for lib in ("libFIREFLY.so", "libFLARE.so", "firefly/f2py."))})
        library = os.environ.get("STAGE5_THERMAL_LIBRARY", "/work/install/lib/libFIREFLY.so")
        if [p for p in mapped if p.endswith("/libFIREFLY.so")] != [library]:
            raise ValueError("incorrect compiled FIREFLY runtime: " + str(mapped))
        with open(library, "rb") as shared:
            signature = hashlib.file_digest(shared, "sha256").hexdigest()
        dest.with_suffix(".runtime.json").write_text(json.dumps(dict(ranks=COMM.size,
            seconds=time.monotonic()-start, cpu_seconds=sum(cpu), library_mappings=mapped, active_transport_sha256=signature)))


def main():
    phase, *args = sys.argv[1:]
    if phase == "field":
        return field(*args)
    if phase == "pfc":
        return pfc(*args)
    if phase == "trace":
        return trace(*args)
    if phase == "mesh":
        import subprocess
        base = NATIVE / "cases" / args[0] / "mesh"
        for task in ("inner_boundary", "flux_tubes"):
            subprocess.run(["flare", "-n", "1", "mgen", "-t", task], cwd=base, check=True)
        return
    if phase == "smoke":
        return solve(OLD/"blind_screen/mesh/mmesh.nc", OLD/"blind_screen/pfc/ref.nc", OLD/"OMP.source",
                     NATIVE/"historical_smoke/historical.pending.nc", 32027, 1024)
    if phase == "heat":
        case, wall, seed, attempt = args
        base = NATIVE / "cases" / case
        result = NATIVE / "results" / case / wall / "attempts" / seed / attempt / "native.nc"
        return solve(base/"mesh/mmesh.nc", base/(wall+".pfc.nc"), ROOT/"inputs"/case/"OMP.source",
                     result, int(seed), 4096)
    raise ValueError("unknown action " + phase)


if __name__ == "__main__":
    main()
