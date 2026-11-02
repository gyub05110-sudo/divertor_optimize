# Geometry-to-transport enclosure implementation plan

> **For agentic workers:** Use executing-plans task-by-task, with independent requesting-code-review checkpoints. Existing user authorization includes implementation after adversarial review.

**Goal:** Enclose every member of an IMP9/OMP16 coefficient box from spline geometry through possible first hits, toroidal connection length, two-point spreading and area-kernel parameters; retain unresolved cases rather than invent a certificate.

**Architecture:** Add a verifier outside the fast production evaluator. Its manifest defines an exact frozen-path model (binary64 anchors, normals, knots, field spline coefficients, source values and cached nodal magnetic quantities become exact constants). Continuous spline walls and a positive continuous spreading integral are enclosed with directed arithmetic. This is not bitwise NumPy execution or a continuous magnetic ODE theorem. Explicit gates prevent promotion to either claim.

**Tech stack:** Python, existing 42-digit directed Decimal intervals, exact rational spline calculations where useful, SciPy only for exporting the already-defined numerical model and independent comparisons.

## Prior delivery (completed before implementation)

- [x] Fetch remote, establish df239e3 as prior publication, audit eight subsequent commits.
- [x] Write docs/WORK_REPORT_2026-11-02.ko.md, run 76 regressions, commit dc01bee and push origin/codex/transport-peak-bounds. Verify remote SHA.

## Adversarial design decisions

Two independent read-only reviewers examined geometry and thermal code before implementation. Critical findings and resulting requirements:

1. **No sampled envelope:** all interval outputs must follow arithmetic inclusion. Samples are counterexample tests only.
2. **First-hit switches:** retain all possible surfaces and path segments before the earliest universally guaranteed crossing. A possible intersection alone never proves arrival for every member.
3. **Normal-graph coordinates:** for each reference v=b-a and frozen normal n use u=cross(p-a,n)/cross(v,n), d=cross(v,p-a)/cross(v,n). This exact invertible coordinate map avoids assuming rounded n is exactly unit/orthogonal. Solve d(path(t))-spline(u(path(t)))=0 by interval exclusion. Uniform opposite signs on connected path portions prove existence even when the intersection crosses a path vertex, but only while the entire connecting path remains inside this same material graph domain; exit/re-entry resets the bracket. Earlier potential roots cannot be discarded based on the midpoint.
4. **All variables:** IMP shared junction angle is an interval affine function of variable 8, with enclosed sin/cos and denominator checks. OMP all16 and IMP all8 offsets remain independent interval inputs. Endpoint constants unrelated to varying junction are frozen explicitly.
5. **Distance identity:** use path edge index/fraction as the primary hit coordinate. Recompute segment lengths and cumulative integrals with directed arithmetic; do not mix GEOS projection with NumPy cumulative lengths.
6. **Toroidal metric:** integrate the declared positive nodal |B|/|Bp| metric piecewise linearly, including the quadratic fractional last segment. Check only prefixes that may contribute to first hits.
7. **Thermal root:** with E=3.5 qu J/kappa0, H=qt/C, y=Tt^(7/4), solve the strictly increasing cubic y³+Ey=H³sqrt(H). Directed bisection brackets roots for all interval inputs. No Brent residual is accepted as proof.
8. **Width:** sigma=sqrt(3 e n chi/kappa0)*(y²+E)^(1/7)*J/(sqrt(y²+E)+y)*average(mapping_over_b). A min/max enclosure of the positive frozen piecewise-linear mapping on the used prefix encloses its continuous transformed-coordinate average. This may be loose but has no quadrature guess.
9. **Area bridge:** enclose positive R, curve speed, full and partial rotation areas. Bound area center divided by total area, intersecting only proven [0,1] domain. Target gradient and magnetic norm come from interval evaluation of the exact exported tensor spline. Zero gradient blocks finite width.
10. **Feasibility vs enclosure:** geometry boxes are supersets; inclusion does not prove every member satisfies curvature/collision constraints. No incumbent promotion or global optimum assertion. Missing arrival, singularity or a numerical budget stop preserves an unresolved record.
11. **Path scope:** frozen RK4 vertices are input constants, not validated continuous trajectories. `continuous_fieldline_enclosure_proven=false`, `production_roundoff_enclosed=false`, `geometry_global_optimum_proven=false` are mandatory. Implementing a validated ODE integrator is a distinct remaining requirement; do not relabel this result as completing it.
12. **Protection:** new D study only; hash old campaign result/snapshot inputs, preserve source and certificate hashes, never resume expired campaign.

## Task 1 — Arithmetic and spline range primitives

Files: create d3d_optimization/enclosure_math.py and verification/test_geometry_transport_bounds.py.

- [x] Write failing tests for exact affine spline values, knot-crossing hulls, sin/cos intervals including turning points, seventh-root enclosure against integer powers, malformed/nonfinite data rejection.
- [x] Implement generic de Boor with interval controls; derivative coefficients use directed arithmetic, and ranges split at every positive knot span. Scalar trigonometric enclosures use Taylor remainder and a Lipschitz radius, never ordinary math.sin as proof.
- [x] Export/test a tensor-spline field adapter; derivative grids must preserve coefficient ordering and correct R/Z axis identity. Test polynomial fields with analytic gradients across knot boundaries.

Concrete API:
```python
root = positive_root(value_interval, 7)
value = spline_range(knots, control_intervals, degree, parameter_interval)
gradient = SplineField(manifest).gradient([R_interval, Z_interval])
```

## Task 2 — Exhaustive first-hit geometry

Files: create d3d_optimization/geometry_transport_bounds.py; extend the above tests.

- [x] RED: analytic plane hit, competing surfaces, earlier fixed wall, knot/path-vertex crossing, interior tangency, no-hit horizon, variable junction, and multiple offset variations.
- [x] Construct interval spline controls from the full box, positive-R surface area bounds, and all path/surface residual candidates. Use interval subdivision to reject only sign-definite residuals. Uniform crossings provide an upper limit on the first hit; retain every earlier possibility.
- [x] Return first-hit union with edge/fraction, curve/u, RZ and existence status. Missing finite-horizon guarantees stay unresolved and cannot enter a normalized kernel.

Concrete output contract:
```python
{"guaranteed_hit": True, "possible_hits": [
 {"edge": 3, "fraction": ["0.2", "0.4"], "surface": 1,
  "u": ["0.1", "0.2"], "point": [["1", "1.1"], ["-1.3", "-1.2"]]}]}
```

## Task 3 — Thermal and area propagation

Files: create d3d_optimization/thermal_box_bounds.py; extend tests.

- [x] RED: constant field analytic L, toroidal scaling, variable final fractions, cubic root residual sign enclosure, small/large positive heat inputs, zero gradient, invalid metric before/after hit.
- [x] Compute exact-model prefix metric/integral enclosures, monotone positive thermal root and mapping average bounds. Propagate every possible first hit; hulls may widen but cannot exclude a member.
- [x] Enclose target field, sigma physical/area and total/partial areas. Emit a complete positive donor kernel only if every positive-power donor has a guaranteed hit and finite positive width. Otherwise emit explicit unresolved reasons without dropping or renormalizing donors.

## Task 4 — Full-domain records and actual studies

Files: create verification/run_geometry_transport_bounds.py; docs/geometry_transport_enclosure.md; results/geometry-transport-enclosure-20261102/REPORT.ko.md.

- [x] Freeze all inputs and bounds, replay deterministic output, reject changed leaves/input/output. Optional geometry bisection retains both children with the identical boundary; unsolved leaves remain in the record.
- [x] Run a bounded real IMP and OMP case using all regional variables, fixed reference, new D output; record interval success/failure per donor and runtime. Compare midpoint/perturbed points only as falsification tests.
- [x] Feed a resolved geometry-derived family to the existing kernel evaluator. Label its lower bound as valid only for the declared frozen-path family; no global gap/candidate feasibility promotion.
- [x] Run targeted regression suite and independent adversarial code review. Repair all important findings and replay saved results. Record actual remaining obligations rather than asserting no possible bugs.

Commands (from repository root, configured Python):
```powershell
python -m pytest verification/test_geometry_transport_bounds.py -q
python -m verification.run_geometry_transport_bounds --output D:/DivertorOptimization/studies/geometry-transport-enclosure-20261102
python -m verification.run_geometry_transport_bounds --replay D:/DivertorOptimization/studies/geometry-transport-enclosure-20261102-v3/IMP_local.json
```

## Acceptance distinction

Success is inclusion for every coefficient vector in the recorded boxes of the explicitly defined model, or a correctly retained unresolved outcome. Useful finite bounds for all full-range boxes are not assumed. Continuous ODE accuracy, continuous source integration, rigorous obstacle/curvature feasibility and a finite global optimization gap remain separate proof obligations. The report must distinguish those from what the new verifier actually establishes.


## Execution record - 2026-11-02

- Implementation: enclosure_math, geometry_transport_bounds, thermal_box_bounds, geometry_transport_certificate, actual model exporter/study runner, independent numerical audit, and regression tests.
- Added validated transformed-coordinate integration (64 panels in the real study), exact straight-segment rotation areas, and a universal wrapped-Gaussian peak bracket when spatial image-series budgets are exceeded.
- Independent review repaired material-domain sign brackets, ambient-context trig negation, unclamped wall endpoints, uncertain knots, and invalid used cache markers. Added targeted regression cases.
- Real v3 study completed in 121.13 s: all 9 IMP / 16 OMP variables active, 16 donors each, local and wide boxes all resolved. Full +/-100 mm boxes remain unresolved and were retained; no feasibility/global-optimality claim.
- Independent numerical falsification: 20 shapes, 2,900 checks passed. Replay verifies IMP/OMP local certificates and both full-domain records. These checks are distinct from the mathematical inclusion argument.
- Original campaign/snapshot: 903 protected file hashes unchanged. Initial image-budget failure, v2 and v3 are separate preserved D studies.
- Final targeted suite: 24 new tests plus 76 existing regressions. Detailed results and limitations: results/geometry-transport-enclosure-20261102/REPORT.ko.md.
- Not completed or claimed: validated continuous field-line ODE tubes, continuous source integration error, full-domain useful finite bounds, geometry-feasibility proof, or certified global optimization gap.
