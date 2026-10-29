# DIII-D / CPD-2 local 3D viewer

Open `http://127.0.0.1:8765/` while its local server is running. To restart, run `start_viewer.ps1` in PowerShell. The server binds only 127.0.0.1; it does not publish files online. Assets and data are local. Opening index.html directly is unsupported because ES modules and data fetch need HTTP.

Choose D3D or CPD-2 in **장치 선택**, then select an equilibrium/run and a recorded stage. Rotate with left drag, pan with right drag, zoom with wheel/pinch, or use camera presets and zoom buttons. Both devices use the common five lambda_q widths (1, 2, 3, 5, 10 mm). Four upstream equivalent diffusion assumptions (sigma/lambda_q=0, 0.25, 0.5, 1) are available through the run selector. Reduce rotation angle to inspect the interior. Hover over the heat surface to read its interpolated vertex heat, or use the visible summary and 2D comparison. PNG export saves the 3D canvas, not the surrounding labels. Device, case, stage, width and revolution angle are stored in the URL fragment.

This is an axisymmetric revolution of each recorded R–Z geometry, not a measured 3D CAD model. Heat is uniform in the toroidal direction. Grey is missing sampled support, not zero. The common color maximum is the greater of the baseline and selected-stage peaks at the selected width. Current common-model cases use q/P_in in m^-2: unit input power is assigned to the modeled lower-outboard branch. Source hashes are in each case bundle; pinned Three.js 0.180.0 provenance and license are in vendor/.

## Equilibrium and history selection

Choose a registered equilibrium/time/run from the list or select the original D3D `.parquet` or CPD `.geqdsk` file. Its SHA-256 selects only matching registered results. Different time slices remain separate entries. An unknown file clears the previous scene and history; no optimization is started by the file chooser.

Current common-model registry entries compare a numerically verified baseline and selected candidate. Original D3D history assets remain available separately. For the original D3D run, baseline + 27 improving coarse candidates + verified final = 29 snapshots. Coarse candidates are explicitly marked; case 6952 failed refined geometry checks and is retained as a historical coarse improvement. The final verified stage uses case 6758. A final verified score can be slightly worse than the raw coarse best.

New CPD-2 results are calculated from OMP radial exponential profiles and one lower-outboard direction, with modeled branch power 1 and a fixed 2 degree floor on positive-power first hits of the prescribed outer target. Other first-wall heat remains in the objective and power integral. `cpd2_lower_outboard_20261004` searches at least 1,500 fresh geometries, followed by independent refined peak, angle, receiver-power and coverage checks. Only verified comparisons are exported. The baseline is the previously optimized configured start. Historical native narrow-width and bidirectional 50:50 result files remain archived; they are not renamed as single-branch calculations.

LCFS is the stored EFIT boundary. Core flux surfaces are closed axis-enclosing cubic-EFIT contours at psi_N=.2/.4/.6/.8/.95, normalized using magnetic-axis psi and lower-X-point psi. SOL displays the lower outboard poloidal branch from upstream midplane to each stage's first material hit. It is not the entire SOL, an actual helical magnetic field line, or a plasma transport solution. Every stage regenerates heat and SOL hit endpoints on that stage's wall; equilibrium remains fixed.

For CPD-2, LCFS is the stored GEQDSK boundary and the five closed core contours use the saved equilibrium normalization. SOL displays the calculated lower-outboard paths, cut at each stage's first wall contact. Historical bidirectional bundles retain their own path and power labels. Camera positions and 2D bounds follow the selected device's dimensions. The 1.5 degree D3D label is not applied to CPD.

## Register another finished compatible run

From a scientific Python environment, run:

```powershell
$env:PYTHONPATH = "$env:TEMP\codex-d3d-inspect"
python "C:\Users\FPL-user\Desktop\LAB\NEW AutoFRD\freegs_divertor\input\D3D\optimized\d3d-divertor-floor-shelf\viewer3d\export_history.py" --case-dir "PATH_TO_FINISHED_RUN"
```

The run must have the same eight-variable floor/shelf format, manifest, evaluations.jsonl, vessel_geometry_audit.csv and candidate_verification.json, including a verified final. Existing archived math is loaded in an isolated namespace; `--code-dir` can select the archive explicitly. This writes only derived assets and upserts catalog.json in this viewer directory. Refresh Chrome to see new entries. Other formats require an exporter adapter; a raw equilibrium without registered results cannot provide an optimization history.

To regenerate CPD registration from the freegs_divertor root:

```powershell
python optimize_cpd_heat.py --verify-run input/Base_of_CPD2/optimized/cpd2_lower_outboard_20261004 --viewer input/D3D/optimized/d3d-divertor-floor-shelf/viewer3d
```

Tests in this directory: `node --test device.test.mjs geometry.test.mjs quality.test.mjs registry.test.mjs unified-data.test.mjs`; `python -B -m unittest test_history_data test_heat_quality test_quality_result test_source_domain`. From the root: `python -B -m pytest verification/test_cpd_viewer_export.py -q`. Browser QA evidence for the unified viewer is in `unified_browser_checks.json` and `unified-cpd2.png`. PNG download completion was not observable in the in-app browser used for QA; use the normal Chrome launcher for that control.
