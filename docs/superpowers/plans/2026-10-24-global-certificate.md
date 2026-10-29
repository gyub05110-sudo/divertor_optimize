# Fixed-model global optimality certificate plan

> Execute inline with executing-plans; test proof primitives first and independently review the certificate before claiming a theorem.

**Goal:** Attempt coefficient-domain global certification for the existing strict (buffer=0) 24-variable model, and explicitly distinguish proof, counterexample, and unresolved optimization.

**Architecture:** A read-only verifier freezes all coefficient-independent geometry/field-path inputs. Outward interval evaluation encloses all allowed B-spline node displacements. Exact rational convex hull/separating-axis predicates exclude intersections for every coefficient vector at once. The production heat/search code and design domain are not changed.

**Scope:** The certificate is about each stated frozen polygonal/path discretization of the current code. It is NOT a certificate for continuous field-line ODE solutions, EFIT uncertainty, or real plasma. Ordinary DOP853/GEOS/SciPy are not formally verified; explicitly separate exact geometric statements from operational evaluator assumptions. Do not set the existing broad `global_optimum_proven` flag from a narrower theorem.

## Literature and approach review

- Rump, *Mathematically Rigorous Global Optimization in Floating-Point Arithmetic*: https://www.tuhh.de/ti3/paper/rump/Ru18.pdf . Global proofs must exclude entire sets, with roundoff included; a successful local/stochastic search is insufficient.
- IntervalArithmetic documentation: https://juliaintervals.github.io/IntervalArithmetic.jl/dev/intro/ . Enclosures require outward rounding. Merely using dense samples is not an enclosure.
- Alternatives: full 24D interval branch-and-bound over spline/ODE/first-hit geometry would require a validated replacement of the evaluation graph; an analytic invariant can close the IMP bound without that implementation. OMP must remain unproven unless a valid lower/upper gap closes.

## Tasks

- [x] Write and run failing tests for outward interval arithmetic, B-spline whole-box displacement enclosures, exact convex separation (including touching/nonseparation), and certificate failure on overlapping envelopes.
- [x] Add `verification/global_certificate_geometry.py`: finite binary64 arithmetic with nextafter outward endpoints; de Boor recurrence over coefficient intervals; exact rational 2D hull/separating witnesses. Use reference.node_parameters, not a substituted uniform parameter grid.
- [x] Add `verification/certify_fixed_model.py`: rebuild the selected model and fingerprint inputs; find the first FIXED wall intersection with the stored separatrix path using exact rational segment predicates; construct prefix hull; bound each moving polyline edge by endpoint displacement enclosures and a stated roundoff box. Verify every separating witness using Fraction arithmetic.
- [x] IMP theorem: if all moving edges exclude the fixed-hit prefix, the witness's selected fixed material, refined root, normal, h0 and q_limit are coefficient-independent. If at least one baseline peak equals this invariant limit and the baseline is feasible, F(c)>=1 and F(0)=1. Report assumptions about GEOS selection and frozen numerical upstream inputs. Do not extend this theorem to the continuous physical problem.
- [x] OMP: report which envelopes overlap (no false certificate). Re-evaluate a saved-candidate perturbation known from earlier numerical challenges, at the same resolutions and geometry constraints. A valid smaller objective refutes optimality of that candidate for that discrete evaluator; it does not certify the new candidate.
- [x] Review exact proof code and outputs adversarially. Run full regressions, preserve legacy results, save certificate JSON, replayable proof data, and a mathematical explanation with explicit failed obligations.

## Pre-implementation adversarial review

1. Bounding node coordinates independently destroys correlation along a normal and may overestimate intersections. Keep the scalar displacement strip plus a small justified roundoff box, then take convex hull.
2. Two fixed zero controls at each end are essential. Enclose the actual de Boor recurrence, including its rounded constant alpha and coefficient-bound tolerance.
3. A prefix must stop at the first fixed intersection, not the historical strike/root. Use exact rational intersection on stored binary64 coordinates. No GEOS-derived prefix will be labelled exact.
4. Convex hull disjointness is sufficient, not necessary. Any nonseparated hull causes a failed proof obligation, never a proof by absence of sampled collisions.
5. The continuous spline between target nodes is not the production polyline. Certify the latter and state that restriction.
6. A stored approximate separatrix can have a tiny gap from the exact strike; this matters. Report that margin and explicitly exclude certification of the true continuous trajectory.
7. A computed positive score difference is not a physical performance gain. Report OMP counterexamples as discrete-model evidence with measured cross-resolution behavior.
8. Do not fabricate an OMP lower bound by dividing normalized tube power by area: current q values and tube power accounting are not the same discretization.

The user has authorized proof work; no geometry extension, new physical assumptions, or additional design authorization is needed for this read-only verification work.

## Execution result

Completed the bounded certification attempt, not the broader continuous-physics proof. Cubic IMP 320 and 2560 whole-box geometry passed; conditional minimum 1.0. Fine minimum exact separation lower bound 5.741024053186218e-12 m. Quintic IMP and cubic OMP lack a fixed-wall intersection, so the invariant proof precondition fails before edge-envelope exclusion. This replaces the planned OMP overlap listing with an explicit failed precondition. OMP incumbent has the same feasible better neighbor at cubic320/2560 and quintic1280. Broad proof flags remain false. See input/D3D/optimized/global-certificate-20261024/PROOF_REPORT.md.

Independent review found no blocker under the stated assumptions. Corrected empty-surface vacuity, added unique operational root check, fixed verifier provenance timing, and retained explicit geometry-only replay scope. Reviewer independently checked fine provenance and ten exact witnesses including the critical imp_sloped first edge; the builder had already completed full replay. Full suite 58 passed; restoration898 passed. A slow early fine run was stopped before producing an artifact; only proof-candidate selection was accelerated, with exact rational acceptance retained, then fine build/full replay were rerun successfully.
