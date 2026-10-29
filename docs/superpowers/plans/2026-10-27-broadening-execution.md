# Broadening execution evidence

## Preservation and tests

- No Git repository exists in the requested folder. Existing source was left in place; the new implementation is additive.
- Before implementation: 41 existing targeted tests passed.
- After implementation and review repairs: full `verification` suite 106 passed, 1 skipped (15.28 s). An additional resolved-nonzero-width limit test subsequently passed; a final full run follows the numerical study.
- `python verify_restoration.py`: 898 historical files verified, restored math verified. This legacy audit explicitly does not validate the new transport surrogate.
- The separate SHA256 preservation inventory confirms all 16 pre-existing production Python files are unchanged.

## Adversarial review

Plan/spec review selected a conservative flux-coordinate endpoint surrogate over wall-arclength blurring or an incomplete transport marcher. Implementation review reproduced and repaired thin-obstacle visibility leakage, zero-q grazing samples bypassing angle constraints, hidden paired flux extrema and exact-hit-vertex length rejection. Independent integration/code review additionally repaired erasure of original failed screens, validation overwrites and out-of-bound historical warm starts. Regression tests were run failing before each supported repair.

Physics/spec disposition: ready for exploratory execution, no remaining important blocker. Code-quality disposition: ready for bounded search and numerical validation. Neither disposition constitutes experimental or physical transport validation. Full review history is in `2026-10-27-broadening-reviews.md`.

## Runs

Smoke run: `input/D3D/optimized/broadening-smoke-20261027`, 48 rays and one DE generation; IMP 13.16 s, OMP 10.66 s. Its spiky profiles demonstrate why smoke output must not be treated as a physical or converged result.

Main run command:

```powershell
python compare_broadening.py --mode optimize --source both --tubes 320 --nodes 161 --step-mm 1 --wall-step-mm 2 --sobol-power 7 --population 24 --generations 8 --output input/D3D/optimized/broadening-comparison-20261027
```

Same search budget/seeds and historical warm-start eligibility are used for each model. sigma_ref=0,0.1,0.25,0.5 mm at L_ref=10m; exponent=1; anchor=5mm. Invalid baseline settings are explicitly skipped, not renormalized or rescued by a later grid.

## Evidence-driven verification refinement

Before final candidate validation, an independent baseline diagnostic at 320 paths / 2 mm tracing gave maximum cached-trapezoid versus midpoint-Simpson connection-length relative errors of 0.02182338 (IMP) and 0.01971968 (OMP). These FAIL the fixed 0.005 threshold. Density-vs-CDF errors were 6.0418e-5 and 1.2348e-6; maximum drift/sigma values were 0.0121872 and 0.00251652. The issue was trace/integration resolution close to the X-point, not power normalization.

At the same 320 paths with tracing refined to 0.5 mm, the connection-length discrepancies decreased to 0.00126325 and 0.00087175; both complete baseline audits passed. No acceptance threshold was relaxed, no physical parameter was changed, and no candidate was discarded.

After independent review, the entire verification matrix was uniformly refined for all saved candidates, models and sources:

| Field | Rays | Trace step (m) | Shape nodes | Receiving step (m) |
|---|---:|---:|---:|---:|
| cubic | 320 | 0.0005 | 161 | 0.002 |
| cubic | 640 | 0.00025 | 321 | 0.001 |
| cubic | 1280 | 0.000125 | 641 | 0.0005 |
| quintic | 640 | 0.00025 | 321 | 0.001 |

The search still used 320 rays / 1 mm tracing. It is not retroactively described as a fine-grid search. The exact verification matrix and initial failed diagnostic are saved in the run's `refinement_plan.json`. Existing hashed code and all thresholds remain unchanged; the verifier's existing `levels` parameter is used.

Exact per-source validation invocation (substitute OMP for IMP for the other source):

```powershell
@'
from pathlib import Path
from d3d_optimization.broadening_study import load_json
from verification.validate_broadening import validate_source
folder = Path('input/D3D/optimized/broadening-comparison-20261027')
validate_source(folder/'IMP', levels=load_json(folder/'refinement_plan.json')['levels'])
'@ | python -
```

Validation is deliberately non-overwriting. A repeated call on an existing complete/partial validation raises `FileExistsError`; retain the original evidence and use a new study destination for a new run.

## Result interpretation

The authoritative final per-candidate verdicts live in each source's `validation.json`. Both sources completed all four levels. No global optimum or experimentally validated transport prediction is claimed.

| Source/model | Finest worst peak ratio | Reduction | Maximum shape difference from direct selection | Interpretation |
|---|---:|---:|---:|---|
| IMP direct | 1.000000 | 0% | — | No improvement found within this search |
| IMP sigma_ref=0.1 mm | 0.941400 | 5.8600% | 6.1859 mm | Numerical checks passed; tested anchors both feasible |
| IMP sigma_ref=0.25 mm | 0.962512 | 3.7488% | 7.0562 mm | Numerical checks passed at default anchor; NOT robust to anchor |
| OMP direct | 0.642915 | 35.7085% | — | Numerical checks passed |
| OMP sigma_ref=0.1 mm | Not qualified | Not claimed | 42.0489 mm, exploratory only | Refined baseline unresolved power exceeds tolerance |

All reductions compare a candidate against the baseline geometry **within the same model**, taking the worst ratio over lambda_q=1,2,3,5,10 mm. Shape differences use the saved search geometry samples; these are not exact continuous Hausdorff distances.

For IMP sigma_ref=0.1 mm, the finest power-weighted mean L increases by 0.116–0.223 m across source widths. Re-evaluating this fixed candidate with exponent=0.5 at the 640-path level gives a worst ratio of 0.957898 (4.2102% reduction); this is a single-level sensitivity check, not another optimized/fully qualified run. For IMP sigma_ref=0.25 mm, the 2.5 mm anchor gives 4.6044% unresolved candidate power and fails the 0.1% gate. The default-anchor improvement must therefore remain conditional.

OMP sigma_ref=0.1 mm has baseline unresolved power 0.073871%, 0.260291%, and 0.163631% on the three cubic grids; the finer levels invalidate its apparent search improvement. No threshold was relaxed and no lost power was renormalized. IMP sigma_ref=0.5 mm and OMP sigma_ref=0.25/0.5 mm failed the original baseline screen and were excluded from optimization.

## Final reporting repair and provenance

The main search saved both numerical studies, then exited with a reporting-only KeyError because IMP validation was concurrently writing a partial checkpoint with an empty comparisons dictionary. Both numerical validations subsequently completed with exit code 0. Full original code/input provenance was checked after completion, and the unmodified report generator then ran successfully.

Only after numerical completion, the original plotting source was archived at `input/D3D/optimized/broadening-comparison-20261027/frozen_source/verification/plot_broadening.py`. The current plotting helper now tolerates partial validation, distinguishes search/refined power, and labels anchor-sensitive improvements. Two regressions failed before the repair and passed afterwards. A fresh independent review found no blockers. Old/new plotting hashes and scope are saved in `presentation_amendment.json`; historical provenance was NOT rewritten. All numerical code and inputs still match the original hashes. The strict historical verifier intentionally rejects the new plotting hash; exact historical reproduction requires the archived plotting source in a separate copy and a new output directory.

Final full suite: `python -m pytest verification -q` — **109 passed, 1 skipped in 14.50 s**. The skipped case is an OMP-specific shelf test on the IMP parametrization. Final preservation check: all **16 existing production files unchanged**. The earlier 898-file historical restoration audit also passed. Both final comparison figures were visually inspected. The completed report is `input/D3D/optimized/broadening-comparison-20261027/REPORT.md`.
