# Continuous IMP Boundary Implementation Plan

> **For agentic workers:** Use executing-plans to execute in this already isolated restored directory. No Git repository is present. Follow the checklist in order.

**Goal:** Remove sampled IMP boundary gaps and false peak improvement, then rerun the authorized geometry with the preserved Sobol + DE method.

**Architecture:** Add a boundary-aware adapter around the unchanged archived heat/shape kernels. Version geometry policy and objective in each saved run. Add a zero-power separatrix-limit witness to the peak objective without changing power normalization. Keep all wall hits in the objective.

**Tech Stack:** Python, NumPy, SciPy, Shapely, pytest; existing local viewer.

## Ordered stages

- [x] Literature: `2026-10-24-boundary-literature.md`.
- [x] Plan self-review and adversarial review recorded before code edits.
- [x] Regression tests demonstrate the missing continuous boundary and endpoint peak guard.
- [x] Implement continuous IMP intervals and boundary witness; preserve legacy replay.
- [x] Integrate objective/export and multi-resolution acceptance.
- [x] Rerun searches, independently check first hits, and publish numerical findings with limitations.

## 1. Geometry and policy

Files: `d3d_optimization/independent_sources.py`, new `d3d_optimization/continuous_boundary.py`, new `verification/test_continuous_boundary.py`.

Compute ψX using the same interpolated equilibrium as the evaluator. Bracket all roots on original segment 69, require exactly one root consistent with the inner leg. Start the movable interval at this root; end at the actual 69/70 vertex. Start vertical 70 at that same vertex and end at Z=-1.15. Preserve all other selected and excluded regions. Keep `boundary_mode='legacy'` for replay. New runs explicitly select `continuous`.

Optional requested design extension: move the start along original 69 by exactly the authorized buffer distance toward PFR. No extension before the user's answer. The original geometry stays identical at zero coefficients. Record buffer in manifest, not in the historical baseline file.

Regression assertions:
```python
assert abs(field.values([strike])[0] - psi_sep) < 1e-12
assert np.linalg.norm(sloped_end - vertical_start) < 1e-12
assert vertical_end[1] == pytest.approx(-1.15)
assert {53, 54, 55}.isdisjoint(movable_original_segments)
```

Run `python -X utf8 -m pytest verification/test_continuous_boundary.py -q` before and after implementation. Retain the expected failing output in the execution record.

## 2. Peak at the continuous boundary

Trace each lower separatrix leg from the local X-point branch beyond the baseline wall, using DOP853 and ψ projection. Do not stop tracing at the old strike. Intersect against the complete candidate wall, then refine ψ=ψX on the hit surface. On a moving straight-reference B-spline use its analytic displacement and derivative; on fixed material use its exact straight tangent. Reject an unresolved boundary witness instead of treating it as zero.

Use `q_limit = abs(B_pol @ normal) / (abs(B_source) * source.normalization)` for each λ, and `peak_q = maximum(sampled_peak_q, q_limit)`. The witness carries no finite tube power. Check shelf angle at the SOL limit too. Preserve historical `maps` for power accounting and add explicit witness/peak fields to JSON and viewer data. A witness is a numerical guard, not a rigorous continuous maximum certificate.

Tests: baseline ratio exactly one; strict fixed strike has invariant q_limit for nonzero coefficients; cache agrees with full-wall intersection; source weights remain one; moving endpoints and witness normal agree with analytic spline derivatives. Test invalid/missing witness fails closed.

## 3. Search and validation

Files: `optimize_independent_sources.py`, new `verification/validate_continuous_run.py`.

Use the existing independent source normalization, coefficient bound, geometry constraints, Sobol seeds/counts and two DE seeds with 40 generations. Output to a fresh `input/D3D/optimized/continuous-boundary-20261024` tree. Old output and restoration archive must not change.

Validation levels: cubic (320,2 mm,161), (640,1 mm,321), (1280,0.5 mm,641), (2560,0.25 mm,1281), plus quintic at 1280. Full-wall first hits must match cached ones; fixed wall coordinates unchanged; no LCFS/vessel intersection; mapping/accounting constraints maintained. Refine the witness path with the validation level as well.

Acceptance defined BEFORE search: all checks feasible, worst ratio <=1+1e-10; fine-grid change in each of the five ratios <=0.002; cubic/quintic ratio difference <=0.002. Improvement must exceed twice the largest measured refinement/interpolation change plus 0.0001. This is a numerical screening tolerance, not a certified error bound. No verified-improvement flag merely because constraints pass.

Commands:
```powershell
python -X utf8 -m pytest verification -q
python -X utf8 verify_restoration.py
python -X utf8 verification/run_continuous_search.py --source both --boundary-mode continuous --mode optimize --output input/D3D/optimized/continuous-boundary-NEW-RUN
python -X utf8 verification/validate_continuous_run.py --output input/D3D/optimized/continuous-boundary-20261024
```
If an authorized buffer is supplied, pass it explicitly and use a separate result path. Compare strict and extended baseline geometry, not candidate coefficients on incompatible domains.

Post-implementation review amendment: new searches use the provenance launcher above. It captures all evaluator modules and the actual baseline-review geometry, not just the geometry reference file. Validation fails closed on absent/incomplete/mismatched fingerprints. The already-running 20261024 searches received sidecars labelled `during_existing_search`; they are not represented as pre-search snapshots. Original 898 archived files passed checks before the search. A historical OMP candidate is additionally replayed and independently validated because it scored better than this finite new stochastic search.

## 4. Reporting

Save manifest provenance, test/validation results, and a baseline/candidate boundary figure. Report physical assumptions, timings, cases, five-width ratios, minimum shelf angle, numerical uncertainty and acceptance. Experimental validation unavailable (user has no IR data); no global-optimality certificate. If strict IMP remains at ratio one, report that honestly rather than relaxing constraints or smoothing heat.
