# Restricted optimality bounds pilot

> **For agentic workers:** Execute inline with executing-plans and test-driven-development. The user has authorized implementation; no renewed approval gate.

**Goal:** Demonstrate a replayable lower-bound search on a one-coefficient D3D family, keeping numerical candidates separate from certified optimality claims.

**Architecture:** Freeze wall and spline constants as exact binary rational numbers. Bound toroidal receiving area by rational interval arithmetic and the peak by P/A. Subdivide the coefficient interval without discarding any region. Numerical transport probes use the existing shared geometry policy; their peaks are observations, not certified upper bounds.

**Tech stack:** Python, fractions, existing SciPy spline and transport evaluators, pytest. New output on D; expired campaign and incumbent bundles remain read-only.

## Scope and proof obligations

The lower bound applies to nonnegative, exactly power-conserving heat flux on the frozen continuous wall with positive R. It ignores shape constraints (a relaxation), so it remains valid on the feasible subset. The continuous model uses exact endpoints and frozen normal/spline constants; production floating-point construction and heat evaluation are not interval-enclosed. Therefore no certified gap or global optimum is asserted for the production transport solver.

## Tasks

- [x] Add failing tests in `verification/test_optimality_bounds.py`: interval spline containment, analytic cylinder area bound, exact partition replay, tamper rejection, refusal to turn a numerical peak into a certified gap.
- [x] Implement `d3d_optimization/optimality_bounds.py`: rational intervals, Machin pi bounds, checked square-root upper bounds, per-knot spline enclosures, whole-wall area upper bounds, deterministic partition history, replay.
- [x] Add `verification/run_optimality_bounds.py`: freeze an actual IMP20 candidate, vary coefficient 2 by +/-0.5 mm at fixed other coordinates, 16 bound subdivisions and 9 numerical probes at 320 tubes. Save inputs, bounds, candidate observations and limitations in a new D study directory; refuse existing output.
- [x] Execute the bounded pilot. Compare lower-bound history against numerical observations without treating the difference as certified. No original campaign restart.
- [x] Run targeted new tests and campaign regressions; review proof assumptions and report results. Preserve historical hashes and commit code plus compact summary.

## Acceptance

`python -m pytest verification/test_optimality_bounds.py -q` must pass. Replay must recompute bounds and reject altered histories. Every reported gap must explicitly distinguish a numerical diagnostic from a certified result. A loose lower bound is a valid result and must not be described as proof of near-optimality.
