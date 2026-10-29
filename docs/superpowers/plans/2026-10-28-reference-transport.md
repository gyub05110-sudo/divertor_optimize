# Fixed reference transport implementation plan

User authorized execution of the preceding agreed method. No Git repository exists; additive files plus a SHA256 inventory preserve all existing Python files. Use test-first implementation and independent physics/code review. Scope: one-way geometry-dependent effective two-point closure, not coupled transverse transport.

## 1. Preserve and specify
- [x] Read existing source normalization, actual first-hit integration, conservative deposition and validation.
- [x] Record all starting Python file hashes in `2026-10-28-reference-preservation.json`.
- [x] Write `2026-10-28-reference-transport-design.md`, inspect primary sources, and incorporate physics reviewer corrections on heat capacity, n_t jump and width convention.

## 2. Thermal numerical module
- [x] Create `verification/test_reference_transport.py` importing `ReferencePoint`, `solve_two_point`, `ThermalPaths`; observe missing implementation failure.
- [x] Test analytic constant-B state and width with `S=a*sqrt(1.5*e*nu*Tu*chi/kappa)*2*L/(Tu**1.75+Tt**1.75)`; require `sigma=S/sqrt(2)` and independent sheath/conduction residuals <1e-10.
- [x] Test variable B (`q/B` constant), changed first hit, power-dependent sigma, chi sqrt scaling, invalid input rejection and actual partial-segment endpoints.
- [x] Implement immutable positive finite controls, log-bracket unique two-point root, cached field/path quantities and Gauss quadrature on actual cut poloidal paths. Use analytic cancellation of Bp in spreading metric; do not cap cold tails.
- [x] Export per-ray Tu,Tt,nt,sigma,L,q_u,q_t, residuals and weighted regime diagnostics.

## 3. Explicit-width deposition and evaluator
- [x] Test equality to old deposition for supplied sigma=old sigma(L), including cell power, q and all ledger terms.
- [x] Add `transport_deposition.py` with explicit-sigma adapter of the reviewed kernel; keep original module byte-identical.
- [x] Add `transport_evaluator.py` with `direct`, `frozen_transport`, `thermal` modes at nominal lambda=3mm. Reuse unchanged source/geometry/first-hit APIs; create new result records, never mutate the direct evaluator.
- [x] Test exact nominal-width direct result, fixed upstream input across geometry, thermal width response, fixed baseline-width control, independent integrated power and receiving-angle constraints.

## 4. Run/evidence interface
- [x] Add `compare_reference_transport.py` with preflight/optimize commands, exclusive outputs and before-search file/input provenance.
- [x] Use matched seeded bounded search with old direct candidates plus baseline as identical warm starts; optimize q_peak/P relative to the same-mode baseline. Single nominal source width replaces the old five-width minimax, explicitly in metadata.
- [x] Save all proposals, selected coefficients, per-ray thermal states/diagnostics, absolute q and q/P, cross-model comparisons and fixed-input metadata. Along-path temperature is used internally; the exports contain per-ray endpoint states, not full T(s) arrays.

## 5. Independent checks and study
- [x] Add `verification/validate_reference_transport.py`; independently recompute full first hits and thermal integrals with finer quadrature, power/source normalization and geometry constraints.
- [x] Run reference preflight before search; preserve any failed baseline rather than adjust reference to obtain a pass.
- [x] Search admissible source/models; cross-evaluate selected candidates at 320/640/1280 cubic and 640 quintic using the design's refined trace steps. Preserve original feasibility failures. IMP improvement confirmed; OMP thermal/frozen excluded by baseline and new direct candidate rejected by refined geometry.
- [x] Evaluate selected thermal/direct candidates at half/double power, density and chi (one at a time); retain unchanged shapes and label these sensitivity evaluations. Preserve failed parameter/anchor variants.
- [x] Independent code review; repair supported findings with regressions before freezing code for the final study.

## 6. Delivery
- [x] Generate static geometry/heat-load and state comparisons plus a Korean report, including invalid cases and one-way closure limitations.
- [x] Run `python -m pytest verification -q`, verify original hashes, and record exact numerical execution/remaining limitations. 132 passed, 1 skipped; 46 original Python files and both main provenance inventories match.
- [x] Deliver reference configuration, usage, implementation/review records and saved comparison results in `docs/REFERENCE_TRANSPORT.md` and `input/D3D/optimized/reference-comparison-20261028/{REPORT,INTERPRETATION}.md`.

Numerical acceptance: old first-hit/source/geometry audit plus independent thermal L/integral differences <=0.005; density-CDF power <=1e-4; absolute peak change <=0.005 and ratio change <=0.002; reduction >2*ratio_change+1e-4. All thresholds fixed before search. Cold/large free-streaming ratio and one-way assumptions are separate physical diagnostics, not experimental validation.
