# Archived OMP milestone qualification, 2026-11-04

User-approved scope: qualify the current geometry-valid OMP Rmin 20/50 mm
improvement milestones, retain the best passing candidate and expose it in the
local viewer. No new geometry search, no modification of the fixed D3D reference,
no change to the 1 MW operating reference, and no GitHub push.

## Selection and verification

- [x] Freeze the 19 + 34 candidate identities, coordinates and log provenance,
  reference, runtime, consumed source and protected input hashes on D:.
- [x] Evaluate every archived candidate at cubic 1280, with identical reference
  and geometry policy. Record failed geometry explicitly. Sort by fine peak.
- [x] In each radius case, fully qualify candidates in fine-peak order until
  the first passing improvement. Require cubic 320/640/1280, quintic 640,
  independent first-hit/thermal/geometry audits, wall-only mesh refinement,
  power accounting, and CPU/GPU q and cell-power agreement. A numerically valid
  candidate with no improvement must never be promoted.
- [x] Save immutable checksummed QualifiedStore bundles. Verify input/code
  preservation, and record all ranks and reasons for rejection or non-evaluation.
- [x] Export a hash-bound final overlay and finest-grid surface profiles to the
  viewer, preserving data.json, legs.json, all stages and earlier qualified files.
- [x] Verify selection adversarially and exercise both final cases, intermediate
  stages, heat colors, q/coordinate hover, legs and absent R100 final in-browser.
- [x] Write inspectable results with the explicit scope: best passing candidate
  among these archived milestones, not a proof of global optimality.

## Adversarial concerns and defenses

1. Coarse order can differ from fine order: score all 53 at 1280 before selection.
2. A low peak may fail physical geometry or convergence: never substitute coarse
   improvement or a numerical-only pass for complete qualification.
3. Fine ranking and validation could drift: bind candidate/context hashes and
   assert equal 1280 peaks. Verify running code and protected inputs unchanged.
4. Resume could mix settings or stale arrays: require exact context and artifact
   hashes before reuse. Refuse a second live writer via an OS file lock.
5. Best incumbent could be lost: publish append-only immutable bundles first;
   convenience indices are reconstructible. No writes into historical campaigns.
6. New viewer final could silently invalidate old profile or leg indices: use a
   separate overlay tied to original data hash; reuse the exact stage wall/ID.
7. A final could be misrepresented as an original 72h final: label the new
   qualification date and show its original archived stage separately.
8. Finding no passing candidate does not prove infeasibility. Stop only when a
   passing candidate is reached or the complete finite list is exhausted.

## Completion evidence

- Study: `D:/DivertorOptimization/studies/omp-stage-qualification-20261104`.
- All 53 scored at cubic 1280; 17 full qualifications attempted; 2 passed.
  OMP20 stage 4: 1.444160060 MW/m², 47.115011% reduction. OMP50 stage 34:
  1.786940489 MW/m², 34.562428% reduction. Remaining 36 candidates have no
  lower fine peak and were not fully qualified.
- The 15 lower-peak OMP20 candidates failed the unchanged 2 mm peak-location
  criterion; near-equal separated maxima exchange order across interpolation.
  This is not evidence that the entire footprint translated by that distance.
- `qualified-final-verification.json` verifies committed hashes, finite ranking,
  full qualification, identities, toroidal area, power and wall alignment.
- Python 15 tests and JavaScript 20 tests passed. Actual browser observations
  verified 4 q/coordinate values across R–Z, 3D and a coarse intermediate stage;
  3 regressions cover absent OMP100 final, preserved IMP20 final and baseline
  profile selection. Browser console contains no errors/warnings.
- Original study/files preserved; D3D worktree clean; no GitHub push.
- Detailed report and figures: study `report/REPORT.md`.
