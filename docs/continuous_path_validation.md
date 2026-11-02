# Continuous-path validation for fast geometry optimization

The research objective is to find divertor shapes that spread deposited power and reduce peak heat flux for a fixed baseline and fixed operating inputs. Continuous validation is a reusable candidate-checking layer. It does not replace the fast evaluator, change the optimization objective to connection length, or modify stored best shapes.

## What is now continuous

The exported C1 tensor spline psi and toroidal constant K are exact defining inputs. In contrast to the earlier polyline verifier, the new path is the continuous solution of the poloidal field-line ODE for that spline. Its toroidal connection length, magnetic thermal integral, and mapping are enclosed along the continuous curve. The nominal polyline selects coordinate planes only and is not assumed to be the solution.

Spline construction/EFIT reconstruction error and continuous source integration remain outside this model. The finite donor source, fixed transport closure, periodic deposition kernel, and prescribed inputs still define a reduced optimization problem. Agreement with measurements is not required.

## Coordinate charts and continuation

Let `(A,q)` be either `(R,Z)` or `(Z,R)`. A chart moves from the incoming axis interval A0 to an exact plane A1 with sign s. Its duration is the interval H=s*(A1-A0), required strictly positive. For each actual incoming state, A0 is fixed and its actual duration belongs to H; no nominal reset is allowed at chart switches.

With `b_pol=(-psi_Z,psi_R)/R`,

\[
A(t)=A_0+st,\qquad q'(t)=s\,b_q/b_A.
\]

The trial rectangle includes every `A0+s*[0,H.hi]`. Strict Picard inclusion of `q0+[0,H.hi]*f(rectangle)` inside the trial q interval proves existence through every actual stopping time. C1 spline regularity and the denominator guard provide local Lipschitz regularity and uniqueness. The checked sign `direction*s*b_A>0` prevents following the opposite physical branch.

The proved endpoint enclosure `q0+H*f(rectangle)` is intersected with an interval Newton enclosure of `psi(A1,q)=psi(source)`. The derivative is evaluated over the entire current endpoint interval. If it includes zero, contraction is skipped. Newton never substitutes for reachability or selects a nominal root. The source, physical direction, field, and source gradient are all bound to the cache.

## Toroidal and thermal integration

Writing `g=|grad psi|` and `B=sqrt(g²+K²)/R`,

\[
{dL\over dt}={\sqrt{g^2+K^2}\over|\psi_q|},\qquad
{dJ\over dt}={B\over B_0}{dL\over dt},\qquad
m={g\over g_{sep}}{B_0\over B}.
\]

Positive interval range integration encloses each full chart cell and every accumulated prefix. This explicitly preserves toroidal winding. A hit anywhere in a cell uses the incoming lower and outgoing upper prefix bounds; it is not assigned a full-cell lower integral.

For spreading, the existing positive two-point/sheath root enclosure is reused with these continuous J and arrival-field bounds. The transformed-coordinate mapping integral enumerates every path cell compatible with each inverse-integral coordinate interval and integrates its mapping range over the entire panel. There is no inferred attenuation or missing-power renormalization.

## First-hit geometry and completeness

For each wall spline, the entire continuous tube is tested against the normal-graph residual. All possible earlier cells remain, including tangency and surface switching. A guaranteed crossing needs uniform opposite endpoint signs while the **un-clipped entire tube** stays in the same material parameter domain. Leaving that domain resets the bracket.

An incomplete late horizon may still certify a first hit already guaranteed inside its proved prefix. A path with no guaranteed first hit cannot supply a normalized donor kernel. Record scope and unresolved reasons remain explicit.

## Reuse and bounded subdivision

`ValidatedPaths` replays path records once, takes private copies, binds source/field/branch, and checks their hashes before reuse. Geometry boxes reuse that cache. The public interface has no `trust=True` bypass.

`cover` retains both children and their common exact split boundary. It prioritizes unresolved leaves; a junction denominator failure can target the IMP junction fraction, otherwise widths are normalized by the original variable domain. This avoids comparing metres with dimensionless fractions. Its finite split budget is not an exhaustive optimization. No child is dropped based on a nominal solution or unproved infeasibility, and no global gap is asserted.

## Evidence and limits

Tests include analytic straight/helical lengths, curved invariant contours, variable chart durations, chart switching, interior C1 knots, near-zero/incorrect direction rejection, tangencies, domain exit/re-entry, missing hits, retained incomplete prefixes, nonconstant mapping, ambient decimal precision, tampered cache, and missing cover children. A separate floating DOP853 audit searches for counterexamples without projecting onto the flux invariant. Its numerical tolerance is not the proof tolerance.

Proof flags concern the exact exported spline and finite donor model. They do not imply every box member meets curvature/collision constraints, a certified feasible incumbent, continuous-source convergence, experimental validity, or global optimality. Runtime and the usefulness of the resulting bounds must be reported before deciding how widely to use the verifier.

```powershell
python -m pytest verification/test_continuous_path_bounds.py -q
python -m verification.run_continuous_path_study --output D:/DivertorOptimization/studies/NEW-UNUSED-DIRECTORY --stride 20 --splits 2
```

All old studies are protected by hashes. New input snapshots, paths, certificates, and numerical audits are stored on D. Repository results contain compact evidence and the Korean research report.
