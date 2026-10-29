# Independent regional optimization plan

Goal: independent IMP-source/IMP-region and OMP-source/OMP-region geometry optimization, fixed operating inputs, minimum arriving peak through spreading. All study outputs, snapshots, GPU dependencies and caches on D.

Reference: P=1 MW per independent source, n=1e19 m^-3, chi=0.3 m^2/s, lambda=3 mm, gamma=7, kappa0=2000. These are two independent experiments, not a 2 MW simultaneous discharge.

Regional variables: IMP imp_sloped + imp_vertical (8 coefficients); OMP omp_entry + shelf + omp_connector + floor (16 coefficients). All other coefficients exactly zero; all other wall coordinates independently checked unchanged. Global wall peak remains the objective to prevent dumping power onto an unoptimized surface. Preserve all old source files/hashes and all original constraints.

- [x] Add tests for coefficient expansion, non-region invariance, accelerated density parity and conserved deposition.
- [x] Implement additive regional adapter and optional GPU density backend. Keep field tracing, geometry intersections, CDF power and continuous peak refinement on CPU. Benchmark double-precision GPU vs CPU on the actual baseline; use GPU only after parity and timing checks.
- [x] Save CPU/GPU preflight and code/input provenance to D. Search both source/region pairs with 64 Sobol proposals, population 24, 10 generations and projected historical warm starts. Retain baseline; log every unique candidate.
- [x] Validate selected candidates with independent first hits/thermal quadrature, wall-only refinement and 320/640/1280 plus quintic640. Require fixed-region error zero; numerical acceptance inherited without weakening.
- [x] Export full-wall and named target geometry CSV/NPZ, heat-flux samples/cell power, donor thermal/path information, optimization history and final coefficient vectors.
- [x] Draw geometry overlays, heat maps, profiles, convergence history and width/connection-length comparisons. Inspect figures, report each source separately and limits.
- [x] Run regression tests, verify preserved hashes and D output inventory; record runtime, GPU use and verification status.

Review: regional restriction avoids the previous all-six-surface search. No simultaneous combined geometry optimum is claimed. GPU changes only evaluation arithmetic, not the transport model. No measurement-agreement gate. A bounded search cannot prove global optimality. If the top candidate fails refinement, evaluate ranked alternative candidates with the same fixed criteria, preserving every result.
