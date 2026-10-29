# OMP tangent, curvature and opening constrained comparison

Accepted scope: preserve the original DIII-D baseline and fixed input power;
rerun OMP-only optimization with a tangent connection at the left floor end,
minimum radius 20, 50 and 100 mm, and the original lip-to-floor clearance.
The other region and all fixed wall segments remain unchanged. Results on D:.

## Design

Use an affine cubic normal-displacement spline: four free interior coefficients,
fixed endpoint positions, right endpoint slope zero, left endpoint slope derived
from the outgoing fixed wall tangent. This enforces oriented G1 continuity.
It does not enforce curvature continuity (G2). Keep old shape/evaluator files
unchanged; additive adapters must use this same spline for wall construction,
separatrix witness, length integration and continuous geometry constraints.

Apply absolute curvature bounds on all four OMP variable surfaces. Original
fixed engineering corners outside the requested junction are excluded. Calculate
curvature maxima and tip-to-floor minimum distance from polynomial stationary
points on every cubic span, including knots/endpoints. Clearance minimum is the
baseline point-to-floor distance (approximately 33.354336 mm), not the distance
between all wall segments, which vanishes at shared endpoints.

The comparison baseline is NOT a feasible candidate under the new tangent rule.
Keep its original evaluator and heat flux separate from the feasible candidate
pool. Search both displacement signs. Existing material-length, collision,
power, incidence and coefficient restrictions remain active. If a radius case
has no feasible/improving candidate, report that without weakening restrictions.

## Implementation / verification sequence

1. Tests first: exact oriented tangent, unchanged endpoints, old spline parity
   for zero imposed slope, both curvature signs, continuous extrema against
   independent dense samples, gap violation and invalid input rejection.
2. Add analytic shape/constraints and explicit evaluator/guard adapters. Preserve
   historical files byte-for-byte. Verify identical geometry in mapping/witness.
3. Find geometry-feasible seeds with constrained local solves from both signs;
   run bounded stochastic/local heat-flux search with independent starts for
   each radius. Log every heat evaluation and feasible pool. Share path caches.
4. Validate selected candidates at 320/640/1280 tubes and wall refinements, plus
   quintic field interpolation, independent thermal/first-hit/power checks and
   old convergence verdict. Refine/reselect if validation fails. Use GPU only
   after CPU parity; preserve CPU refined-peak and power integrals.
5. Save coefficients, geometry, heat flux, ledger, continuous constraint audits,
   comparison figures and Korean report on D:. Run regression/preservation
   checks and request independent code review before claiming completion.

No global optimality certificate or experimental calibration is claimed.

## Review-driven refinements

- Scalar maximum-curvature constraints made the first seed projection slow and
  unreliable. Restarted as a new v2 study with scaled variables and smooth
  sampled signed-curvature inequalities; continuous polynomial extrema still
  decide acceptance. The aborted study and its original snapshot are retained.
- Capture the historical warm-candidate manifest as an input; check CPU/GPU
  parity on the first feasible candidate as well as the selected candidate.
- A failed radius search reports no candidate and does not stop the other cases.
- IMP coefficients remain exactly zero. Coordinate invariance uses 1e-12 m
  tolerance because equivalent interpolation paths differ by about 2e-16 m.
- Keep numerical convergence distinct from improvement: a converged candidate
  with higher peak is a valid adverse result, not an improved geometry.
