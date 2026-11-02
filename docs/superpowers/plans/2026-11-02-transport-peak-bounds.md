# Transport peak enclosure and coupled IMP pilot

User authorization: proceed with stronger transport-aware bounds and candidate error validation. Execute inline using test-driven-development; no additional approval gate.

## Design

Keep the original heat/geometry models unchanged. Introduce an independent verifier for the exact infinite periodic Gaussian mixture defined by frozen donor weights, arrival area fractions, diffusion widths and total wall area. Decimal directed arithmetic encloses arithmetic; correctly rounded exp/sqrt results are expanded by one adjacent representable number. Rational Machin pi bounds supply pi. An analytic positive image-tail bound covers omitted wraps. Branch the receiving coordinate interval, retaining all leaves, to bound its global maximum; use point lower enclosures, never numerical local optimizer outputs, as lower bounds.

Also allow explicit interval-valued kernel inputs. The pointwise lower enclosure then bounds every kernel's maximum, conditional on those input intervals. Such intervals are NOT automatically proven to contain the transport results of a geometry box. Do not infer a geometry-box enclosure from sampled extrema or claim full geometric optimality.

Run a finite 3-variable IMP20 pilot (offset indices 0,1,2, each +/-0.5 mm, 5 levels each). Preserve shared geometry and material constraints. Select numerically feasible candidates; re-evaluate the best and seed at 640/1280 tubes and report mesh sensitivity separately from rigorous frozen-kernel peak enclosures. Keep all output and source snapshots on D. Preserve old campaign and qualified incumbent files; do not promote a pilot candidate into the qualified store without all qualification gates.

## Tasks

- [x] Add failing tests for directed interval arithmetic, single Gaussian analytic peak, seam/broad mixtures, conservative truncation, partition coverage, tampered replay rejection, interval-input containment and invalid inputs.
- [x] Implement additive `d3d_optimization/transport_peak_bounds.py` with replayable certificates explicitly scoped to frozen/interval kernel inputs.
- [x] Add `verification/run_transport_peak_pilot.py` with immutable source/input provenance, 125-point coupled exploration, preserved best candidate, numerical resolution checks, fixed-kernel certification and explicit conditional robustness envelope.
- [x] Run bounded pilot, inspect outputs, and independently review the mathematical inequalities and scope labels. Repair issues and rerun affected verification.
- [x] Preserve raw D results, write compact repository report, test regression groups, verify old hashes, commit.

## Acceptance and limitations

The exact frozen kernel maximum must be enclosed to a stated absolute tolerance or return an explicit budget-exhausted bracket. Different resolutions are observations, not ODE/error certificates. Full transport parameter interval propagation, exact physical geometry feasibility and geometric global optimality remain unproven. Report these limits and actual improvements without promoting sampling to proof.
