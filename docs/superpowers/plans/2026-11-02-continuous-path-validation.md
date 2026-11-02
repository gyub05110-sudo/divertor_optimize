# Continuous paths and bounded geometry subdivision

## Research objective and approved scope

Fast divertor geometry optimization from the fixed D3D baseline and prescribed operating inputs, conserving assigned power while spreading/reducing peak heat flux. Experimental agreement is not a gate. Verification must help compare candidates, reuse fixed-field work, expose runtime, and preserve prior results. User approved continuing continuous-path error enclosure and unresolved geometry subdivision; no repeated design approval is required.

Keep the fast optimizer unchanged. Add a reusable continuous-path cache for the exact exported magnetic spline, then a verifier that consumes that cache for many geometry boxes. Never trade away power or geometry constraints to make a bound look successful. Retain unresolved results and separate proof scope from numerical evidence.

## Alternatives and decision

1. Interval time integration around RK4: general, but wrapping can grow across long tokamak legs and requires high-order/reconditioning work.
2. Validated graph charts plus flux invariant (chosen): axisymmetry gives an exact psi invariant; coordinate-plane steps with Picard existence and interval Newton contraction avoid transverse drift, preserve toroidal travel, and allow cache reuse.
3. Dense sampling: useful falsification only, insufficient for inclusion and excluded as a proof mechanism.

The nominal polyline guides coordinate-plane choices only. The new proof never assumes the ODE remains on that polyline. Exact exported spline/source constants remain a reduced model; field reconstruction error, continuous source integration, full feasibility and global optimality are not automatically certified.

## Mathematical checks before implementation

- Require C1 spline psi and local nonzero poloidal/chart gradient: the chart vector field is locally Lipschitz on each accepted rectangle.
- With axis A and other coordinate q, use A=A0+s*t, q'=s*b_other/b_axis, s=+/-1. A0 may be an interval at a chart switch; target coordinate is exact and terminal time H=s*(A_target-A0) is an interval with positive lower bound.
- Trial rectangle uses all t in [0,H.hi]. Strict Picard inclusion q0+[0,H.hi]*f(Y) inside the trial q range proves continuation/existence there. Verify direction*s*b_axis>0 throughout to prevent traversing the wrong branch.
- The endpoint q0+H*f(Y) includes all possible variable terminal times. Tighten using interval Newton of psi(q,A_target)=psi_source; only intersect proven enclosures. Never reset to a numerical guide point.
- Along the chart, L'=sqrt(g²+K²)/abs(grad_psi_other) and J'=(B/B0)*L'. Positive rectangle integrals enclose every continuous prefix. Mapping is g/grad_sep*B0/B.
- First-hit exclusion uses entire continuous tube rectangles. Uniform endpoint signs guarantee a crossing only along connected cells inside the same material surface parameter domain. Retain earlier possible crossings and tangencies. Partial-hit L/J uses the whole possible terminal cell, not linear interpolation of its endpoints.
- Replay continuous paths and geometry output; do not trust a mutable cache merely because a boolean says proven. Input/hash compatibility must include field, source point and verification controls.
- Geometry subdivision retains both children, exact common boundaries, unresolved leaves, bounded work. Normalized coordinate splitting avoids mixing metres with junction fractions. No objective pruning until a valid feasible incumbent and full lower-bound coverage are available.

## Tasks and acceptance

1. [x] Analytic tests first: constant/helical field, curved flux contour, chart switching, knot crossing, wrong direction, zero gradient, invalid C1 knots, deliberately insufficient Picard rectangle, budget exhaustion and tampered cache.
2. [x] Add continuous field/chart stepping and replayable path records. Positive continuous L/J/mapping/B enclosures. Record attempted/accepted cell budgets; never hide incomplete horizons.
3. [x] Add continuous tube first-hit and thermal/area propagation, preserving all donor power. Tests: earlier fixed wall, surface switch, tangency, missing hit, last-cell integral bounds, cache/model mismatch.
4. [x] Add bounded domain subdivision with normalized widths and complete tree replay. Allow unresolved leaves and report resolved fraction/count without interpreting it as probability or global proof.
5. [x] Run new D study from frozen v3 inputs, measure path preparation versus reused-box evaluation time, compare independent high-accuracy ODE samples as falsification. Include full-domain bounded subdivision and local cases. Stop at explicit budgets and retain evidence.
6. [x] Independent adversarial review, regressions, actual replay and protected hashes. Write Korean report with research relevance, runtime, achieved scope and outstanding limitations; commit code and compact results.

## Intended files

- d3d_optimization/continuous_path_bounds.py: chart vector field, Picard tube, invariant contraction, integrals, trace/replay.
- d3d_optimization/continuous_transport_certificate.py: cache validation, first hit, thermal and area bridge, complete subdivision tree.
- verification/test_continuous_path_bounds.py: analytic and adversarial tests.
- verification/run_continuous_path_study.py: D-only bounded real study, original hashes, timing and source snapshot.
- results/continuous-path-validation-20261102/REPORT.ko.md: measured results and unresolved obligations.

## Scope gates

Only a complete compatible path horizon and guaranteed first hit for every positive-power donor permit continuous-transport inclusion. A partial path can be reused only if the guaranteed hit is within its proven prefix, with explicit prefix scope. No automatic incumbent promotion or change to saved optimized shapes. Final global optimum/gap flags remain false/null unless all separate obligations are actually discharged.

## Adversarial review outcomes

- Reproduced and fixed wrong-physical-branch cache reuse. Cache identity now binds independently inferred donor branch; same-branch guide resampling remains valid.
- Preserved interval incoming axis and variable coordinate-plane stopping times; strict Picard inclusion precedes flux-invariant contraction.
- Coarse first-hit tubes initially resolved only 12/16 IMP donors. Refined possible-hit cells and neighbors, retraced from the exact source, and replayed all final records. No nominal reset or sampled-crossing proof.
- Independent DOP853 audit initially failed by 2.6335e-8 m. Tightening only the ODE solver tolerances reduced the discrepancy below the unchanged comparison threshold; failed attempt retained on D.
- Reviewer independently passed all 17 new tests and found no remaining proof/provenance blocker. Explicit cautions: guide completion is horizon-relative, and reused preparation time is not cold-cache timing.
- Final targeted regression: 117 passed in 5.52 s. Full-domain finite subdivision is an explicit unresolved result, not evidence of infeasibility or global optimality.

## Measured acceptance results

The v2 D study completed in 548.1 s. Four local/wide boxes enclosed all 16 donors each. Cache-reused box evaluations took 3.674-4.131 s. Final cache replay took 123.685 s (IMP) and 127.138 s (OMP). Full-domain subdivision retained three leaves per source; zero were resolved. All path/box/tree replays passed, as did 102,049 independent floating comparisons (falsification only). Verified 1,022 protected hashes, 150 new artifact hashes and 96 source snapshots. Peak bounds remain too loose for candidate ranking or a global gap; no optimized shape or objective change is claimed. See the Korean report and verification.json for evidence and the distinction between reused and cold preparation.
