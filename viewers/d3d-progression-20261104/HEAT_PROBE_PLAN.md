# Surface q probe and heat colors

Goal: point at either the 3D divertor or the R–Z wall to inspect the stored local heat flux, with one shared quantitative color scale. Local-only update; no push, optimization or heat solve.

Design: use the existing per-profile area-weighted display cells. Tooltip identifies q as a display-cell mean in MW/m² and reports R/Z in metres. No new interpolation of q across cells or unsaved candidates. Load the matching profile independently of the color checkbox so q inspection also works with colors hidden. At initial load choose a refined result when available; explicit timeline controls still select intermediate geometry. Baseline selections load the saved refined baseline and label its resolution.

3D picking: raycast the opaque displayed wall only, excluding legs, LCFS and ghost lines. Interpolated R/Z vertex attributes recover the poloidal coordinates without toroidal chord-radius error. Find the nearest profile cell only within the existing 20 µm shape-matching tolerance. 2D uses the exact displayed axis transform and closest wall point, with a 7 pixel hit threshold. A distinct baseline ghost line gets a baseline instruction rather than the selected candidate's q. Missing profile or uncovered wall is not q=0.

Colors: default on, applied to both views using the same fixed 0-to-baseline-refined-peak scale and the existing sequential heat palette. Add a larger numeric legend generated from that exact palette. Uncovered wall stays neutral. Keep the original peak metric separate from the display-cell average. Guard asynchronous loading on source/run/stage changes; clear stale tooltip on state/camera changes, pointer exit and drag.

- [x] Add failing tests for nearest projected segment/cell, zero vs missing data, distance tolerance, profile selection, and 3D R/Z interpolation.
- [x] Implement a small pure probe module and integrate shared profile loading, color legend and initial refined selection into app.mjs.
- [x] Add 2D pointer probing and 3D wall raycasting with a clamped tooltip; exclude annotation geometry and preserve orbit dragging.
- [x] Verify tests, source hashes and saved q values; browser-check nonzero/zero/missing readings, heat toggle, source changes and final heat colors. Save screenshots and concise documentation.

Adversarial review: a radial coordinate from hypot(x,z) on a faceted torus is biased, so interpolate supplied R/Z attributes; pointer distance measured in metres needs conversion from the same 2D pixel scale; do not pick a heat patch through an opaque wall; ignore leg/LCFS surfaces; distinguish unrecorded values from exactly zero; guard stale fetches and tooltips; no fake intermediate heat profiles.

