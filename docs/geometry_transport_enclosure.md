# Geometry-to-transport interval verifier

## Contract and scope

`geometry_transport_certificate.certify_box` encloses every coefficient vector in a supplied box of the **declared frozen-path model**. The exported model contains all IMP 9 or OMP 16 variables. A resolved record includes possible first-hit branches, arrival R/Z, toroidal connection length, magnetic integral, spreading width, rotation area, and periodic-kernel parameters. An unresolved record is a valid refusal to certify, not an infeasibility result.

The model treats exported binary64 path vertices, spline coefficients, normals, source powers, and nodal thermal caches as exact constants. Wall curves are continuous splines; nodal thermal quantities are piecewise linear along the poloidal polyline. The spreading integral is continuous and enclosed, not replaced by a sampled maximum or assumed quadrature error.

This is **not** a validated solution of the continuous magnetic field-line ODE, a bound on NumPy/SciPy execution roundoff, or a certificate that every geometry satisfies curvature and obstacle constraints. It does not certify continuous source integration. The corresponding flags remain false, and `certified_geometry_optimality_gap` remains null. The production evaluator and old campaign are unchanged.

## Inclusion chain

### 1. Wall and first arrival

Each surface is `x(u)=a+u v+n d(u)`, `0<=u<=1`, with `v=b-a`. The frozen normal need not be exactly unit or perpendicular after binary rounding. Coordinates are therefore computed with the invertible determinant formulas

\[
u={\operatorname{cross}(p-a,n)\over\operatorname{cross}(v,n)},\qquad
d={\operatorname{cross}(v,p-a)\over\operatorname{cross}(v,n)}.
\]

The IMP shared junction parameter changes both endpoint slopes through interval sin/cos and a checked positive denominator. It is not frozen at the nominal candidate.

For every path edge and surface, interval de Boor evaluation encloses `d(path(t))-d_spline(u(path(t)))`. Only sign-definite residuals or disjoint boxes are excluded. Remaining subintervals retain tangencies, shared endpoints, and switches between receiving surfaces. A uniform sign change guarantees a hit only along a connected path portion wholly inside the **same** surface's material domain. Exiting and re-entering that domain resets the sign bracket.

The earliest universally guaranteed hit gives a cutoff. Every possible earlier hit remains, even if a nominal ray would hit a different surface. No guaranteed hit before the finite path horizon leaves the donor unresolved. Missing donor power is never removed and renormalized.

### 2. Toroidal length and thermal integral

With poloidal edge length `h`, positive nodal values `v_j,v_{j+1}`, and terminal fraction `t`, the partial integral is

\[
h\left[t(1-t/2)v_j+(t^2/2)v_{j+1}\right].
\]

Directed square roots enclose `h`; directed cumulative integration encloses every prefix. Apply this expression to the frozen metric `|B|/|Bp|` for connection length and to `(B/B0)|B|/|Bp|` for `J`. Thus the length includes toroidal travel in this model. Invalid caches beyond all possible first hits do not invalidate a usable prefix; invalid contributing caches do.

### 3. Thermal root and spreading

Let `qu` be upstream parallel heat flux, `qt=qu B_hit/B0`, and

\[
C={\gamma e n\over2}\sqrt{2e/m_i},\quad H=q_t/C,\quad
E={3.5q_uJ\over\kappa_0},\quad y=T_t^{7/4}.
\]

The closure becomes `y^3+E y=H^(7/2)`. The root is positive, strictly increasing in the right-hand side and decreasing in `E`. Directed bisection at monotone corners, with exact rational residual comparisons when needed, encloses it without trusting a floating root tolerance. Then `Tu=(y²+E)^(2/7)`.

For positive frozen mapping `m`,

\[
\sigma_{flux}=\sqrt{3en\chi_\perp/\kappa_0}
 (y^2+E)^{1/7}{J\over S}\,\overline m,\qquad
S=\sqrt{y^2+E}+y.
\]

The transformed average is

\[
\overline m=\int_0^1m(F^{-1}(p(u)))du,\qquad
p(u)=J(1-u)(1+uE/S^2),
\]

where `F` is the positive cumulative thermal integral. The broad bound is the full prefix mapping hull. Optional interval panels tighten it: for every panel enumerate every compatible inverse segment, bound the inverse fraction by the positive nodal integrand range, enclose the mapping, and integrate that enclosure over the whole panel. Shared exact panel boundaries cover `[0,1]` without gaps. Intersect with the broad bound. Panel counts govern tightness, not validity.

### 4. Field, area, and kernel

The target field uses the exact exported tensor spline (Z rows, R columns), including differentiated coefficients and every intersected knot span. A gradient enclosure containing zero prevents finite width certification. Convert by `sigma_physical=sigma_flux*grad_sep/|grad psi_hit|` and `sigma_area=2*pi*R_hit*sigma_physical`.

The receiving coordinate is the cumulative rotation area divided by the total area. Positive interval integration encloses `2*pi*R*|dx/du|`; straight zero-displacement segments have an analytic area formula. Total area includes the entire declared closed wall, including fixed surfaces. Dependencies are relaxed by interval hulls; this widens bounds but does not remove possibilities.

Each donor uses an exactly normalized infinite periodic Gaussian in area coordinates. Its integral over the period is its prescribed frozen power. For positive area `A`, widths `sigma_i`, and powers `w_i`,

\[
{\sum_iw_i\over A}\le\max_a q(a)
\le\sum_iw_i\left({1\over A}+{1\over\sqrt{2\pi}\sigma_i}\right).
\]

The upper bound follows from positive Fourier coefficients (a single wrapped Gaussian peaks at its center) and bounding the monotone image tails by their integral. It needs no finite image budget. Existing spatial subdivision can tighten these bounds; exceeding its image budget retains the universal bound and records the limitation.

## Global optimization interpretation

A bound on an entire family is not an attained best candidate. A true optimization certificate needs both a lower bound `L` covering **all feasible geometries** and a validated feasible incumbent with objective upper bound `U`. Then `L <= f* <= U` without knowing `f*` or its shape. `U-L` bounds remaining possible improvement; a lower bound need not correspond to any realizable shape.

The current record encloses local boxes and conservatively records unresolved full-domain boxes. `cover` can bisect a box, retaining both closed children and their exact shared boundary. Replaying reconstructs the complete tree; omitted children or changed results are rejected. It does not yet select useful boxes by objective, prove full-domain feasibility, or produce a global optimization gap.

## Verification and reproducibility

The verifier reuses 42-digit directed Decimal interval primitives and rational pi bounds. Rational Taylor remainders enclose trigonometric functions. Replay compares the entire recomputed record, including the embedded input digest and claim flags. This is a reproducible numerical certificate, not a machine-checked formal proof of the Python program.

Tests cover tangencies, first-hit switching, leaving/re-entering a material domain, knot peaks, shared endpoints, all-variable junction movement, toroidal prefixes, extreme thermal root scales, invalid data, zero gradients, tampering, and incomplete cover trees. Independent SciPy scalar roots and quadrature are used only to search for counterexamples; they are not the inclusion oracle.

```powershell
python -m pytest verification/test_geometry_transport_bounds.py -q
python -m verification.run_geometry_transport_bounds --output D:/DivertorOptimization/studies/NEW-UNUSED-DIRECTORY --count 16
python -m verification.run_geometry_transport_bounds --replay D:/DivertorOptimization/studies/geometry-transport-enclosure-20261102-v3/IMP_local.json
python -m verification.audit_geometry_transport_bounds D:/DivertorOptimization/studies/geometry-transport-enclosure-20261102-v3/OMP_wide.json --samples 5
```

The study command refuses existing output directories and the expired campaign directory, snapshots source, and verifies protected inputs before and after execution. CPU directed arithmetic is intentional: the current GPU fast evaluator is not used as a rounding or inclusion oracle.
