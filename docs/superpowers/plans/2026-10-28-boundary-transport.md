# Fast conservative boundary transport implementation plan

Execute inline with executing-plans. User explicitly authorized diagnosis, model repair, and optimization; no further permission gate is required.

**Goal:** Preserve the real baseline and hypothetical inputs while enabling fast OMP geometry optimization with a complete, explicit reduced deposition model. Experimental agreement is out of scope.

**Architecture:** Add diagnostics for the frozen legacy kernel, then a conservative boundary spreading operator with an explicit continuous receiving domain. Retain the existing thermal width calculation and fixed-field tracing where appropriate. Use independent kernel/conservation/geometry tests and a new CLI/output namespace. Review the receiving metric before implementation to prevent artificial shortcuts in optimization.

**Tech stack:** Python, NumPy/SciPy, Shapely, existing frozen-field evaluators, pytest, Matplotlib.

## Tasks

- [x] Inventory current Python hashes and reproduce legacy OMP missing power on the fixed baseline.
- [x] Add `verification/diagnose_boundary_transport.py`: per-donor branch-tail and blocked-power accounting, cold-tail attribution, destination/flux-range diagnostics, CSV/JSON/plot exports. Reconcile independently computed totals to the frozen ledger within 1e-10.
- [x] Complete `docs/superpowers/specs/2026-10-28-boundary-transport-design.md` after the diagnostic and adversarial review. Compare conservative receiving-domain approaches on locality, speed, broad/narrow limits and boundary topology. Do not silently fill unresolved probability into arbitrary walls.
- [x] Write failing tests in `verification/test_boundary_transport.py`: positive complete cell power, cyclic seam invariance, no shortcut across nonadjacent wall segments, near-isoflux wall handling, broad/narrow Gaussian limits, geometry-dependent physical width, invalid input errors.
- [x] Implement additive `d3d_optimization/boundary_transport.py` and evaluator; run the tests after observing RED. Preserve every old Python file.
- [x] Add a new optimization runner using fixed inputs, cached repeated proposals, absolute peak objective and independently checked plate-length constraints. Support baseline evaluation and bounded search with explicit budgets and exclusive output paths.
- [x] Benchmark against the old kernel on the same reference geometry, grid and input. Check OMP and IMP baseline conservation plus half/double input variants and wider numerical-width regimes.
- [x] Independently review implementation and repair supported findings with regressions before freezing the final run.
- [x] Run OMP optimization and fixed-candidate refinement with two trace/wall resolutions and an interpolation cross-check. Report qualified improvements or failed checks honestly, retaining a feasible baseline candidate.
- [x] Generate geometry, arriving-power and heat-load plots, concise Korean results and reproducible commands. Run full regression suite and verify original hashes.

Acceptance: fixed baseline/inputs, closed nonnegative power accounting without survivor renormalization, finite outputs for supported inputs, explicit receiver/locality assumptions, consistent arrival/heat distributions on refinement, existing geometric constraints, useful measured runtime. No measured-data agreement requirement.

Execution completed: see boundary-review.md, BOUNDARY_OPTIMIZATION.md, and input/D3D/optimized/boundary-optimization-20261028/OMP/REPORT.ko.md. Final OMP improvement confirmed by the preregistered numerical screens: 48.2519%; 234 search candidates in 280.6 s including setup. Original 169 Python files unchanged.
