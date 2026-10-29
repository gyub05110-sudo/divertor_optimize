# Fixed hypothetical operating-point transport comparison

The user explicitly approved selecting a hypothetical reference operating point and comparing all shapes under the same input/boundary conditions. Proceed through design, independent review, implementation and execution in this session; no additional permission gate. Preserve the historical direct and 2026-10-27 broadening source files and outputs.

## Scope and alternatives

Implement a geometry-dependent, one-way two-point conduction closure feeding the existing conservative endpoint deposition methodology. This is a useful intermediate between an arbitrary sigma(L) law and a coupled 2D SOL solver. A constant-speed advection model requires a prescribed velocity profile and is not selected for the first reference. A coupled conduction/advection/cross-field PDE needs a mesh, particle/momentum closure and substantially broader validation and is outside this increment.

The new model solves each donor's parallel thermal state for its actual first hit, then derives a transverse kernel. It does NOT feed deposited power back into the donor temperature solution, solve particle continuity/momentum, model neutrals/radiation/detachment, or change the equilibrium. State this explicitly in outputs. Thermal endpoint quantities are pre-redistribution tube states, not the final receiving plate plasma temperature.

## Reference point, fixed before any shape search

One independent source leg receives P=1e6 W, lambda_q=0.003 m, upstream electron density nu=1e19 m^-3, chi_perp=0.3 m²/s, deuterium ions with Ti=Te, sheath gamma=7, effective parallel conductivity coefficient kappa0=2000 W m^-1 eV^-7/2. IMP and OMP are independent one-megawatt scenarios, NOT two parts of a prescribed total machine power. These values are hypothetical design inputs, not inferred DIII-D discharge measurements. Source support, normalization, geometry freedoms, field and engineering constraints remain those of the existing case.

No independent Tu or Tt is imposed. All candidate shapes share these input parameters and boundary formulas. Each source's same nominal 3 mm width is the optimization objective. Other operating points are subsequent fixed-candidate sensitivity checks, not extra optimized cases.

## Equations and conventions

Temperatures are eV; q is W/m²; e is elementary charge. Positive B denotes magnitude. Along each source-to-first-hit path, AB=constant and q_parallel(s)=qu B(s)/Bu, with qu=P exp(-x/lambda_q)/source.normalization. Integrated source ray weights are used only for deposition power, not to infer qu from a wall area.

Define J=integral_0^L B/Bu ds. Solve the positive, monotone two-point system:

```
Tu^(7/2) = Tt^(7/2) + (7/2) qu J / kappa0
nt = nu Tu / (2 Tt)
qt = qu Bt/Bu = gamma e nt Tt sqrt(2 e Tt / mi)
T(s)^(7/2) = Tt^(7/2) + (7/2) qu/kappa0 integral_s^L B/Bu ds
```

The bulk static-pressure approximation n(s)=nu Tu/T(s) is distinct from nt: the factor-two transition accounts approximately for the sonic sheath endpoint's dynamic pressure. No continuous Mach/density transition or momentum solution is claimed. The conductivity and total sheath heat coefficient form a conventional effective two-point energy approximation, not a separately resolved electron/ion model.

Use chi_parallel=kappa0 T^(5/2)/(1.5 e n) and the Nille-inspired conductive spreading closure:

```
g(s)=|grad psi| / |partial_R psi_sep|
S_xi = integral g(s) sqrt(chi_perp/chi_parallel(s)) ds
sigma_xi = S_xi / sqrt(2)
```

S is the Eich Gaussian width convention; sigma is its standard deviation. This is an explicitly chosen spatial integral closure, not an exact heat-packet variance solution. In poloidal coordinates g ds=R B/|partial_R psi_sep| ds_p; exploit this cancellation for stable integration near the X-point. This application uses the full midplane-to-target path because the exponential is prescribed at the midplane, not an already broadened measured target profile. Geometry mapping is applied in the width and then the wall Jacobian; do not multiply target flux expansion a second time.

Source: Nille, Sieglin and Eich, https://arxiv.org/html/1610.04148 (conductive broadening integral, explicit omission of convective transport). Two-point pressure/sheath convention: https://www.sciencedirect.com/science/article/am/pii/S2352179123000649 (Variation in the volumetric power and momentum losses in the JET-ILW scrape-off layer).

Independent physics review clarified: chi_perp is an ELECTRON diffusivity with kappa_perp=1.5 e n_e chi_perp. Its numerical value must not be reused with a different heat capacity. The pre-sheath n_t jump must never be smeared over a quadrature interval. The ray-dependent full-path metric is an extension of Nille's separatrix/divertor model, not a reproduction. For constant B and constant mapping a, the analytic check is S=a sqrt(1.5 e nu Tu chi_perp/kappa0) 2L/(Tu^1.75+Tt^1.75). Diagnostics identify cold/weak-conduction/large free-streaming-ratio conditions without asserting detachment or valid fluid flow. No temperature floor, width cap or discarded-power renormalization is allowed.

## Architecture and preservation

- `reference_transport.py`: immutable reference/numerical controls; positive two-point solver; cut-path thermal integration and donor diagnostics.
- `transport_deposition.py`: additive copy of the previously reviewed Gaussian deposition routine accepting explicit donor sigma; import existing wall-mesh/CDF helpers. Preserve original source; record its hash and test equivalence when widths match.
- `transport_evaluator.py`: one-source-width adapter, direct/thermal/frozen-width controls, full receiving-angle and power constraints, absolute MW/m² export. Reuse existing geometry and direct first-hit evidence.
- `compare_reference_transport.py`: immutable run inputs/provenance, matched bounded searches, cross-evaluation and readable artifacts.
- `verification/validate_reference_transport.py`: independent thermal quadrature, full first hits, source power and wall-area audits, resolution comparison and off-reference fixed-candidate evaluations.

Controls: direct mapping; frozen donor widths calculated once on the baseline at the SAME operating point; recomputed thermal widths. Frozen widths isolate the response of the thermal closure to changing geometry without arbitrary calibration. None may alter the fixed input power or renormalize unresolved deposition. Preserve old failed screens, and never overwrite a partial/completed run.

## Validation and acceptance

Analytic constant-B two-point identities, gamma/nu/power dependence, variable-B conservation, finite/positive validation, source normalization, and Gaussian-width conventions require tests before implementation. Reproduce the old kernel for identical sigma. Verify changed geometry changes thermal states but not imposed source power. At least 320/640/1280 donor levels and a 640 quintic field sensitivity are used for final candidates; path steps 0.5/0.25/0.125 mm and 0.25 mm respectively, reflecting the previous study's X-point integration evidence. Kernel density/CDF power difference <=1e-4, independent thermal length and integral differences <=0.5%, absolute peak refinement <=0.5%, peak ratio change <=0.002. All original/refined feasibility screens must pass. An improvement must exceed twice the observed ratio change plus 1e-4. Numerical acceptance never implies physical validation or a global optimum.

Report weighted cold-target fraction, weak-conduction/free-streaming diagnostics, sigma/lambda, unresolved tail/blocked power and two-point energy/sheath residuals rather than hiding model limits. A baseline that fails power or resolution screens is diagnosed and excluded; reference parameters are not tuned after observing optimization outcomes. Off-reference checks halve/double P, nu or chi separately at fixed candidate shape; report validity and ratios without claiming new optima.
