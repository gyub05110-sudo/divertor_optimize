# OMP qualification review and correction plan

User authorization: audit the prior work and correct or roll back unsupported changes. Implementation is authorized in this chat. Planning, debugging, test-driven development and completion verification skills apply proportionately; no new optimization or acceptance policy is requested.

## Fixed direction

Keep D3D geometry, operating inputs and assigned power fixed. The scientific objectives are lower peak and a broader heat footprint. Numerical checks measure confidence in those objectives and in local q; they are not interchangeable with the objectives. Argmax location is diagnostic for nearly equal peaks. Do not change the 1%/.5% numerical screen to get a desired candidate accepted. These are provisional engineering choices, not calibrated physical requirements or certified continuum error bounds.

## Decision

Full rollback would reinstate the inappropriate argmax-only rejection. Relaxing v3 to pass R20 would hide a reproduced local discrepancy. Preserve both historical policies and their outputs; audit the evidence and correct the presentation of what passing/failing proves. Do not promote a new final or change selection weights. Compare A2 with baseline as an explicit objective diagnostic, not merely A2 consistency between resolutions.

## Steps

- [x] Independently inspect policy, archived arrays and checksums; reproduce the largest discrepancy directly from the saved kernel, without interpolation or renormalization.
- [x] Add `audit_policy_review.py`, writing only to a new `D:/DivertorOptimization/studies/omp-policy-audit-20261104` directory. Verify sealed artifact manifests, protected code and original data, rerun the 20 archived v3 decisions, compute A2 and peak changes at each level, and record numerical variation without claiming it is a rigorous uncertainty bound.
- [x] Back up `dist/app.mjs`, `dist/index.html`, and `README.md`. Correct only wording: provisional screen, unresolved local numerical differences, and distinction from improvement. Add a linked static audit explanation with evidence; preserve existing historical results, profiles and pass flags.
- [x] Run relevant Python/JavaScript regressions, check actual browser OMP20/50/100 and IMP states, and verify old file hashes again. Record findings and limitations in a new audit report.
- [x] Obtain independent read-only code review and resolve any material findings before completion.

Scope exclusions: no new tolerance version, no new optimization, no changes to frozen numerical modules or saved candidate rankings, no claim of experimental validity or global optimality. Candidate ordering at one resolution is not proof that tiny differences are numerically resolved. Further refinement/calibration is a separate numerical investigation; retain current conservative full-distribution qualification until such evidence exists.
