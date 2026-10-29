# Connection broadening implementation plan

> Execute tasks in this session, with independent adversarial plan and code reviews. User has already authorized execution after review; no additional approval gate. Existing production files are preserved, so no Git/worktree operation applies.

**Goal:** Compare plate shapes under the existing direct model and an explicitly scoped L-dependent, power-accounted deposition surrogate.

**Architecture:** Add a numerical kernel, an evaluator wrapper, a dedicated comparison runner and verification/plotting modules. Keep all historical evaluator/optimizer files byte-identical. New schemas do not masquerade as historical verified direct runs.

**Tech stack:** Existing Python, NumPy, SciPy, Shapely, pytest and Matplotlib.

## Task 1 — plan and preservation
- [x] Inspect source, first-hit, continuous boundary guard, search and provenance code.
- [x] Establish baseline tests: 41 passed before implementation.
- [x] Save `2026-10-27-broadening-preservation.json` with all existing production hashes.
- [x] Write design and consult independent adversarial physics reviewer.
- [x] Resolve final plan findings and record them in `2026-10-27-broadening-reviews.md`.

## Task 2 — numerical kernel, tests first
Create `d3d_optimization/broadening.py` and `verification/test_broadening.py`.
- [x] RED: analytic constant helical field has L=5*s when Bp=1,Btor=sqrt(24), including fractional segment; invalid Bp/nonfinite lengths raise.
- [x] RED: standard Gaussian transfer equals CDF interval masses, nonnegative, no renormalization, lost tails explicit, weighted totals sum to one.
- [x] Implement `BroadeningConfig` with finite nonnegative sigma_ref_m, positive length_ref_m, exponent in {0.5,1}, fixed anchor and mesh controls.
- [x] Implement cached cumulative length with trapezoidal field metric, partial-hit interpolation and invalid/unhit handling. Prohibit using poloidal length directly as parallel length.
- [x] Implement receiving wall subdivision, checked monotone flux branches, donor branch association and vectorized Gaussian power/density kernels. Use stable Gaussian CDF differences in tails.
- [x] Add full-domain receiving-cell fan visibility without passing through LCFS/solids; blocked mass remains unresolved.
- [x] GREEN: analytic convolution/area, branch boundaries, blocking, invalid parameters and zero-weight rays pass.

## Task 3 — evaluator adapter, tests first
Create `d3d_optimization/broadened_evaluator.py` and `verification/test_broadened_evaluator.py`.
- [x] RED: sigma=0 yields exact old arrays/peak/feasibility; nonzero model preserves original source and cannot mutate direct baseline.
- [x] Wrap a prepared continuous evaluator. Use its geometry and first-hit evidence; compute broadened peaks separately from unsmeared boundary q. New receiving incidence checks cover powered shelf cells.
- [x] Return explicit deposition arrays, component ledger, length diagnostics, model config and validity flags. Reject nonfinite q/negative power/missing evidence and missing fraction >1e-3.
- [x] Confirm IMP and OMP unit powers separately and independent full-wall first-hit consistency.

## Task 4 — reproducible comparison
Create `compare_broadening.py`, `verification/validate_broadening.py`, `verification/plot_broadening.py`, tests for verdict and provenance.
- [x] Separate `preflight`, `optimize` and `verify` modes, finite argument validation, no overwrite. Capture all new/old model and actual input hashes BEFORE search; verification refuses changes.
- [x] Fixed seeds, matched Sobol/DE budget for each sigma, existing direct candidate as optional provenance-tracked warm start. Save all evaluations, coefficients and candidate checkpoints.
- [x] Build baseline/direct-best/broadening-best by model cross-evaluation, with per-width peaks and ratios, shape changes, power ledgers, length statistics and scoped validity status.
- [x] Independently recompute first hits and audit wall-area-integrated power, grid/interpolation sensitivity. No smoke search can claim a qualified optimum or confirmed physical improvement.
- [x] Render static geometry/profile comparison and save a readable report of model assumptions, findings, failed screens and reproducible commands.

## Task 5 — execution and independent review
- [x] Run targeted tests followed by all existing verification tests.
- [x] Run both-source preflight and bounded matched searches at declared sensitivity strengths. Increase numerical resolution only if evidence demands it.
- [x] Refine and cross-evaluate selected candidates; preserve failures rather than tuning acceptance to results.
- [x] Obtain independent spec/physics review, then code-quality review; repair supported findings and rerun affected checks.
- [x] Verify original production hashes and deliver code, reviewed plan, results and limitations.

## Commands and expected outcomes
`python -m pytest verification/test_broadening.py verification/test_broadened_evaluator.py -q`: initially import failure before implementation; afterwards all pass.

`python compare_broadening.py --mode preflight --source both --output input/D3D/optimized/broadening-preflight-20261027`: model/geometry/ledger diagnostics, never an optimization claim.

`python compare_broadening.py --mode optimize --source both --output input/D3D/optimized/broadening-comparison-20261027`: reproducible separate searches and raw cross-model table.

`python verification/validate_broadening.py --output input/D3D/optimized/broadening-comparison-20261027`: refined candidate comparisons with explicit pass/fail evidence.
