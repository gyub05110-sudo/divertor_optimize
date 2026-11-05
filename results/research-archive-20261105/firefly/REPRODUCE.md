# Reproducing the DIII-D midplane comparison

The completed D drive study contains the executable environment, input files, and raw output. Run the Windows commands from `D:\DivertorOptimization\studies\firefly-baseline-20261105`. `enter.sh` mounts the retained D-drive ext4 file and enters the CI image using WSL Ubuntu as root. Numerical runs write into `/work`, which is bound to this D-drive directory. The original optimization worktree is a read-only input to these scripts. Run series serially when repeating them on a smaller machine.

```powershell
$enter = '/mnt/d/DivertorOptimization/studies/firefly-baseline-20261105/enter.sh'
python export.py
wsl -d Ubuntu -u root -- bash $enter trace.sh
wsl -d Ubuntu -u root -- bash $enter tight.sh
python analyze.py
python prepare_candidates.py
wsl -d Ubuntu -u root -- bash $enter candidates.sh
python geom_candidate_audit.py
python casing_check.py
python analytic.py build
wsl -d Ubuntu -u root -- bash $enter analytic.sh
python analytic.py audit
```

The stored field-line meshes are usable directly. To regenerate the three meshes, run `mesh.sh`, `refinemesh.sh mesh_fine 640 0.002 8`, and `refinemesh.sh mesh_finer 1000 0.0015 8` through the same WSL entry point. For each mesh, run `diag.sh mesh`, `diag.sh mesh_fine`, or `diag.sh mesh_finer`. These use FIREFLY's reconstructed magnetic field lines. The `*_mesh_lc.dat` outputs include the two independently traced terminal coordinates. `mesh_compare.py` reports missing hits separately from misclassified wall sections. Small direct-mapping differences at a vertex can change a subsegment id.

```powershell
foreach ($src in @('OMP','IMP')) {
    wsl -d Ubuntu -u root -- bash $enter bothlegs.sh $src
    wsl -d Ubuntu -u root -- bash $enter fullseries.sh $src
    wsl -d Ubuntu -u root -- bash $enter longseries.sh $src
    wsl -d Ubuntu -u root -- bash $enter replicates.sh $src
}
python mesh_compare.py
python thermal_compare.py
python summarize.py
python report.py
python verify_study.py
```

The two `bothlegs.sh` controls remove the internal reflecting midplane and retain the imported midplane heat source. The two physical endpoints of a field line need not receive equal power. `fullseries.sh` includes the 2 ms sensitivity and zero-transverse-diffusion control, `longseries.sh` uses 0.2 s per history for finite transport. `replicates.sh` adds four seeds at nominal temperature, spatial resolution, and time step. Numerical output for a seed also depends on the number of MPI ranks. The archived nominal runs used five ranks for the first two seeds and four ranks for the other four. The six-seed combined result is written to `{OMP,IMP}_pooled_profile.csv`. The independent ring calculation reads SI coordinates directly from NetCDF, checks the native cell areas against the analytic areas of revolution, integrates the archived reference cell powers into identical rings, and never normalizes by deposited power.

The unmodified FLARE application is from the official image `hfrerichs/flare:v1.4.7`, manifest `sha256:cc26d9250125d515eabc6d1a43eb5de98de11f1463770f36f66c3d2f5432b41c`. The OCI blobs are saved by digest under `oci`. The FLARE runtime independently prints `development-1.4.7-6-g3e06acfd-dirty`. FIREFLY source is pinned at `85a1e6c006e8af0c6a6cb6775ff78efee61af47e`. The retained installed shared library and Python wrapper in `install` implement the study. To rebuild separately from pristine FIREFLY, apply `firefly-midplane.patch` to this exact commit and build with `-DPARTICLE_EXHAUST=OFF` as in its CI. Install the Python wrappers and library together. The Linux build requires LF line endings for the f2py `autowrap` script. For NumPy f2py 2.2.4, `build_compat.sh` records the single generated wrapper import correction. Compare the resulting library and runs to the recorded outputs, allowing for compiler, rank and random-number differences.

The source table has full-torus watts, R, Z and the original sign of the downward leg. FIREFLY draws equal normalized sector-power histories by inverse CDF within deterministic cumulative-power strata and a uniform toroidal angle. The reflecting plane has its own named zero-load tally. Modified source launches are checked by a magnetic coordinate round trip. For DIII-D the adapter uses deuterium and `kappa0 = 2000`, and thermal coefficients are evaluated at the stated fixed density and temperature. The additional diagnostic columns do not change the shape geometry or the transport equation. `prelim_cap` holds superseded exploratory outputs, excluded from the report. When developing the cap, the boundary-event status on a rejected transversal step and the handling of a previously commented inner-boundary diagnostic were corrected before the reported cap runs. Exploratory build and diagnostic logs are kept separately from the results.

Environment requirements include Windows Python with NumPy, SciPy, Shapely, PyArrow, h5py and matplotlib, WSL Ubuntu, the saved image and checkout, the archived DIII-D source tables and the paths in `export.py`, `prepare_candidates.py` and `thermal_compare.py`. `checksums.json` records the selected executable, archived raw calculations, and final files. The reported bidirectional controls use the final binary, five ranks and the same 0.2 s cutoff as the corresponding reflecting-boundary cases. On completion, the chroot mounts can be released by calling `detach.sh` directly through WSL. It leaves the D-drive data intact.
