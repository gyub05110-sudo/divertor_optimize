# Regional divertor optimization data

Repository: C:/Users/FPL-user/Desktop/LAB/NEW AutoFRD/freegs_divertor_no-transport

Two independent experiments at 1 MW per source; no simultaneous-source coupling assumed.
IMP varies imp_sloped and imp_vertical. OMP varies omp_entry, shelf, omp_connector and floor.
Other wall coordinates and the magnetic equilibrium remain fixed. The objective is global arriving peak reduction at conserved assigned power.

All large outputs, dependency files, snapshots and CUDA caches are on D. Original C project files are preserved. The optional CuPy runtime is at D:/DivertorOptimization/runtime/cupy, with environment override DIVERTOR_GPU_RUNTIME available.

## Commands used (from the repository directory)

Set process-local MKL_NUM_THREADS=1, OMP_NUM_THREADS=1, TEMP and TMP=D:/DivertorOptimization/tmp.

```powershell
python optimize_regional_transport.py --output D:/DivertorOptimization/studies/regional-20261028 --source IMP --backend gpu --generations 10
python optimize_regional_transport.py --output D:/DivertorOptimization/studies/regional-20261028 --source OMP --backend gpu --generations 10
python optimize_regional_transport.py --output D:/DivertorOptimization/studies/regional-20261028 --source both --mode verify
python qualify_regional_candidates.py --output D:/DivertorOptimization/studies/regional-20261028 --source IMP
python qualify_regional_candidates.py --output D:/DivertorOptimization/studies/regional-20261028 --source OMP
python optimize_regional_transport.py --output D:/DivertorOptimization/studies/regional-20261028/qualified --source both --mode verify
python visualize_regional_transport.py --output D:/DivertorOptimization/studies/regional-20261028
```

The two search commands ran in separate processes. Existing source-study directories cannot be overwritten: choose a new output root for another search. gpu_preflight.json contains the measured CPU/GPU timing and numerical parity of three actual OMP320 candidates; it is a prior calibration, not a whole-search speedup estimate.

## Main files

- REPORT.ko.md: final Korean results, assumptions, figures and links.
- heat_flux_comparison.png: common-scale baseline/selected IMP and OMP heat maps.
- summary.csv / summary.json: final comparison data.
- optimized_geometry.zip and final_geometry/: only the baseline/selected geometry CSVs and manifests for convenient reuse.
- completion_checks.json: independent delivery checks, conservation, exact CSV-array agreement and preserved source/snapshot hashes.
- verify_delivery.py: reproducible delivery checker; run from the project directory after rendering figures.
- IMP/ and OMP/: preserved original search results, including the lowest-peak candidates that failed the strict interpolation ratio threshold.
- qualified/IMP/ and qualified/OMP/: replacement candidates screened from the original history and subjected to full validation; consult validation.json for acceptance. Their copied search histories retain the ORIGINAL search minimum. manifest.candidates.selected identifies the replacement.
- qualification/: ranked-candidate interpolation screen; thresholds never relaxed.
- transport_comparison.png: parallel length including toroidal winding and footprint width versus source offset.
- */selected_<surface>_refined.csv: final named-surface R,Z geometry.
- */selected_wall_refined.csv: full wall with only that source region optimized.
- */selected_profile_refined.csv: arriving cell-average q and cell power/area.
- */selected_refined.npz: complete arrays, including sampled q and per-donor thermal state.
- */search/evaluations.jsonl: every unique proposal, including invalid proposals and constraint reasons.
- IMP/snapshot/ and OMP/snapshot/: exact production code and inputs captured before the original searches.
- qualified/*/qualification_evidence/: screening data and qualifier source, with hashes; see IMP provenance guard-fix record for the nonnumerical overwrite-protection correction.
- preservation.json: SHA256 of all 180 pre-existing Python files.

The optimization and receiving-boundary transport are reduced models. Numerical validation and fixed-region checks are required before interpreting improvements. Two independent optimized regions must be reevaluated together if later combined into one wall under simultaneous sources.
