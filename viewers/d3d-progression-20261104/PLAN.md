# D3D improvement progression viewer — 2026-11-04

Goal: a local 3D viewer of archived, baseline-improving IMP/OMP shapes across Rmin 20/50/100 mm. Preserve the original fixed D3D reference, independent 1 MW sources, all studies, and existing site. No GitHub push or hosting publication.

Reuse: the existing site's axisymmetric Three.js surface-of-revolution approach, cutaway, camera controls, 2D context and timeline. Copy only pinned vendor JavaScript and license. Do not read its old case/catalog/shape/heat JSON as research data.

Decision: a separate static D-drive viewer, with a new offline exporter. Adapting the old multi-device runtime would retain obsolete data assumptions; a hosted redesign adds unrequested deployment. Use the existing visual structure with a focused D3D interface instead.

Data contract:
- Separate experiment series by source, Rmin and run. Never concatenate runs as a fictitious optimization history.
- Baseline is step zero, even if it is outside the new candidate tangent family.
- Stream evaluations; select original feasible running-best records with objective/peak ratio < 1, then deduplicate and recheck exact current curvature/tangent/gap and full-wall collision policy. Exclude failed records without rewriting the original search history. This scope can omit later valid candidates that did not beat an invalid raw record; empty means no displayable recorded milestone, not no improving candidate exists. Cache screening by candidate identity.
- The campaign's compact feasible objective is peak/baseline_peak (source snapshot verified). Use that same coarse baseline to calculate recorded MW/m². Never compare this number to a refined baseline.
- Match final validation to immutable candidate identity. Final refined peak/heat is a separately labeled result at its own resolution; it is not a chronological search step.
- No invalid geometry, worse-than-baseline candidate, nonfinite value, or unproven historical promotion is presented as an improvement. An empty series is an explicit result, not a fabricated shape.
- Historical OMP20 and failed campaign OMP20 remain distinct. Latest requalification evidence is loaded for historical OMP50 and campaign IMP20/50.
- Reconstruct geometry from stored coefficients and the repository's shape classes / frozen fixture, without thermal solves. Save source paths and hashes, candidate identities, raw sequence and screening counts.
- Heat coloring only comes from matching saved NPZ arrays; intermediate shapes without profiles show solid material, with recorded peak but no fabricated heat. Shared source/baseline scale. Axisymmetric rendering, not a non-axisymmetric simulation.

UI: source tabs, Rmin buttons, run select, baseline/candidate/overlay, step selector and previous/play/next, peak reduction chart, 3D cutaway and focus presets, synchronized 2D R–Z overlay, final validated candidate shortcut, per-stage provenance, actual radius/gap and validation status, PNG export. Empty-state and WebGL fallback. No arbitrary Rmin slider interpolation between recorded studies.

Implementation sequence:
- [x] Write exporter tests for running best, rejected candidates, duplicate identities, finite ratios and geometry/heat separation; observe expected failure.
- [x] Implement export_data.py and stream the recorded studies into dist/data.json and optional profile files. Verify input hashes and preserved source files. Reconstruct curves with current geometry policy.
- [x] Copy pinned vendor assets only. Build dist/index.html, style.css, app.mjs, view-math.mjs; test axisymmetric coordinates and heat endpoint/missing data behavior.
- [x] Serve locally on 127.0.0.1. Verify source/Rmin/run navigation, timeline, final shortcut, empty cases, 2D/3D coherence, responsive layout and screenshots through the browser.
- [x] Write README.md, verification.json and launch script. Keep all new data on D. Open the working local viewer and provide its URL. No commit or push.

Acceptance: every displayed improvement has same-resolution baseline ratio < 1, recorded source provenance, current geometry pass, monotonic running-best sequence within its run; no old viewer research JSON in the data manifest. Final geometry can be compared with baseline with fixed camera and scale. Report unresolved thermal validation explicitly.
