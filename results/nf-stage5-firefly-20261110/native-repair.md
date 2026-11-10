# Native reflecting-boundary correction

The retained FIREFLY source is the historical commit 85a1e6c006e8af0c6a6cb6775ff78efee61af47e with its archived midplane patch. The additional `cap-sync.patch` adds one executable statement in the heat transport task. At a transverse hit on the inner magnetic boundary FIREFLY clamps the magnetic coordinate. The correction synchronizes the cached physical position in a flagged tube when the reflecting cut is active. Diffusion coefficients, source, random-number generator, mesh, physical wall, trajectory cutoff, MPI layout and tally definitions are identical between the control and the patched executable.

The control source and the instrumented diagnostic reproduced the archived historical and the failed Stage 5 seed NetCDFs byte for byte, respectively. In MPI rank 2, history 763, iteration 376422, the unmodified routine recorded a cap hit whose magnetic position lay 0.390 nm above the actual plane. A second apparent hit reversed the remainder again and sent it upward. `detail.log` records both coordinate representations and the received signs. `replay.f90` and `replay.log` isolate this same parallel draw. Starting with the stale coordinates reproduces the original two contacts and final Z of +0.654560476 m. Synchronizing the coordinates gives one contact and final Z of −0.628222774 m.

From PowerShell in this directory, the following commands regenerate both shared libraries from the retained Linux build objects. The original installed library and the original study are inputs to these commands. The new library is loaded separately by process. The required native include directories, original link command and compiler flags are present in the archived chroot. `build.py` uses the original five object files for the other FIREFLY modules. Running the build again in a separate scratch directory gave the same corrected SHA-256.

```powershell
$entry = '/mnt/d/DivertorOptimization/studies/nf-stage5-firefly-cap-repair-20261110/enter.sh'
python make_sources.py control
python make_sources.py synchronized
wsl -d Ubuntu -u root -- bash $entry python3 /repair/build.py control
wsl -d Ubuntu -u root -- bash $entry python3 /repair/build.py synchronized
wsl -d Ubuntu -u root -- bash $entry bash /repair/run.sh control historical ref 32027
wsl -d Ubuntu -u root -- bash $entry bash /repair/run.sh synchronized historical ref 32027
wsl -d Ubuntu -u root -- bash $entry bash /repair/run.sh synchronized 3264fdd71f 021606f6597e71819ae87d1a57605978aae101e680c5a4c7e84de77b2a153318 51031
wsl -d Ubuntu -u root -- bash $entry bash /repair/run.sh synchronized 3264fdd71f 892b5e0939244002c986d16aa445f123f3f2ee3eeb024fbe63a125ceff5ef10c 51031
python audit_controls.py
```

To repeat the small deterministic native step test, compile the C-interoperable test procedure against the archived modules, then initialize the case geometry and call it twice. The original and synchronized paths run sequentially on rank zero with one shared read-only mesh. The test checks the original terminal Z and requires one cap contact and a lower terminal Z on the synchronized path.

```powershell
wsl -d Ubuntu -u root -- bash $entry gfortran -fPIC -shared -O3 -I/tmp/firefly-build/modules -I/usr/local/include/flare -I/usr/local/include/moose -I/usr/include -I/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/openmpi /repair/replay.f90 -L/work/install/lib -lFIREFLY -L/usr/local/lib -lFLARE -o /repair/replay.so
wsl -d Ubuntu -u root -- bash $entry mpirun --allow-run-as-root --oversubscribe --bind-to none -np 4 python3 /repair/replay.py
```

The control launcher writes each raw NetCDF and timing under its variant directory. Rerunning a command rewrites its scratch output. The corrected runtime is `synchronized/libFIREFLY.so` with SHA-256 `3f5da3b88f37ba6420297c1590fd4207f4b5f2c0635cf97a5fed53bd806b34c0`. The historical FIREFLY shared library retains `d5ebfe0dc123f5bc44ad78a524f93ca71b7c74c48c9f3f14cb762a06bf5d6a68`. The original FLARE and FIREFLY Python extension shared libraries retain their hashes in the Stage 5 execution environment. The corrected trajectories can consume different downstream random draws. Compare controls with identical MPI rank count, seed and native units. `control-results.json` reports cell-integrated full-torus watts, upper-wall watts, geometry, area, ledger closure and power-profile differences.
