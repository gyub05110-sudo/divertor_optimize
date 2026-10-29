# Runtime and archived evidence

The publication keeps the project source and input files byte-for-byte, except the root README, whose original is at `docs/README_LEGACY.md`. `.gitattributes` disables newline normalization because historical SHA256 checks include source bytes.

## CPU environment

The checked host uses Windows, Python 3.13.12, NumPy 2.4.3, SciPy 1.17.1, Shapely 2.1.2, Matplotlib 3.10.8 and pytest 9.1.1. See `requirements-repro.txt`.

Historical preservation inventories include 105 PyArrow Python files. Install **pyarrow==25.0.1** into `.local_dependencies/d3d-fixed-equilibrium` before running the preserved optimization/provenance paths. The package's native binaries and other installed files are intentionally not committed. The publication checks a fresh PyPI Windows CPython 3.13 wheel against the preservation hashes. Other operating systems/interpreter builds have not been validated for byte-identical restoration.

The old `verify_restoration.py` checks the entire historical copy, including third-party files and viewer data. It has a narrower historical purpose than the latest numerical validation. Use `tools/verify_publication.py` for the published source and study fingerprints, and the relevant numerical validator for a new calculation.

## Optional GPU environment

CPU operation does not require CUDA. The recorded GPU runtime is RTX 3090, CuPy 14.2.0, CUDA runtime 12.9, driver 610.88. Runtime packages were installed separately on D:

```powershell
python -m pip install --target D:/DivertorOptimization/runtime/cupy --no-deps cupy-cuda12x==14.2.0 fastrlock==0.8.4 nvidia-cuda-runtime-cu12==12.9.79 nvidia-cuda-nvrtc-cu12==12.9.86 cuda-pathfinder==1.8.4
$env:DIVERTOR_GPU_RUNTIME='D:/DivertorOptimization/runtime/cupy'
$env:CUPY_CACHE_DIR='D:/DivertorOptimization/runtime/cupy-cache'
```

Use `--backend gpu` only with a compatible GPU/driver and verified installation. The CPU mode is selected explicitly by `--backend cpu`; it is not an automatic fallback after GPU errors. CUDA packages/caches are excluded from Git. GPU tests require this runtime; CPU-only test instructions are in the README. [CuPy installation documentation](https://docs.cupy.dev/en/stable/install.html)

## Results and paths

`results/regional-20261028/` retains original search histories, both rejected lowest-peak candidates, qualified candidates, validation JSON, CSV/NPZ geometry and transport arrays, figures, source fingerprints and qualification evidence. The original `IMP/snapshot` and `OMP/snapshot` folders duplicate repository inputs/code and are excluded. The geometry ZIP is also omitted because its extracted `final_geometry` content is included. Full original study copies remain at D:/DivertorOptimization/studies/regional-20261028 on the authoring host.

Historical manifests, analysis scripts and publication inventories retain original absolute paths as evidence; they are not silently rewritten. Consequently an archived helper such as `verify_delivery.py`, `analyze.py`, or plotting a qualified historical manifest may refer to the original D/C paths. Static reports, CSVs, NPZs and images can be read directly anywhere. For a new computation, run the root CLIs from this clone with a new output directory. To re-run an archived helper elsewhere, make an explicitly separate path-adapted copy. The published provenance verifier resolves code/input fingerprints against this clone without editing original manifests.

Geometry and operating conditions are preserved, but platform/library changes can affect numerical results. Final candidate selection still requires the documented grid/interpolation checks; matching source hashes alone does not establish numerical convergence.
