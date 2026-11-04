# Leg overlay and qualification labels — 2026-11-04

Goal: show the inner/outer divertor separatrix legs from the fixed X-point to each displayed wall, while clearly separating current qualified results from archived refined results. This is a local viewer update, with no push, optimization or thermal solve.

Design: use the preserved D3D baseline review and its fixed psi grid. Trace the two poloidal separatrix branches once with the existing `extended_leg` geometry routine; clip each against each exported wall at its first intersection. Cache data in a separate `dist/legs.json`, keyed by exact candidate identity plus baseline. In 3D show R–Z leg curves on toroidal cuts, and optionally a transparent axisymmetric separatrix surface. This is a poloidal/separatrix representation, not a newly computed helical particle orbit or transport path.

Why this approach: drawing the old baseline legs unchanged would end above a lowered OMP plate. Straight extensions or interpolating strikes would not follow the fixed psi contour. New thermal or optimization calculations are unnecessary for this display.

- [x] Add first-hit clipping tests: two wall crossings chooses the first; a lower wall extends the branch; no hit returns an explicit unresolved result; baseline leg endpoints are reproduced.
- [x] Implement `export_legs.py`; validate fixed LCFS/baseline identity, finite coordinates, first-hit endpoints on the displayed wall, psi-contour drift, and unchanged original `dist/data.json` and stored optimization input hashes. Save per-identity clipped legs with source hashes and diagnostics for both traced branches.
- [x] Add a leg toggle, inner/outer color legend, X-point/strike markers to 2D and 3D; retain baseline versus selected distinction and all source/Rmin/run/step controls. Default legs on. Show a clear note if leg data fail to load without hiding otherwise valid shape data.
- [x] Rename final buttons/headings: current qualified results use final-validation language; older refined results explicitly say historical. State that the adopted numerical acceptance criteria and global optimality proof are distinct. Keep the original qualification files immutable.
- [x] Run unit/data checks and browser tests for lower OMP shape, IMP, baseline-only, empty R100 and leg toggle. Save an actual screenshot and update README/verification notes.

Adversarial checks: never extend a baseline polyline by a straight chord to a new strike; never select the farthest intersection; never reuse another candidate's strike; never treat clipping arc length as L_parallel; do not silently relabel the three historical refined candidates as latest qualified; preserve the run-wide code-scope audit exception in the written qualification explanation.

