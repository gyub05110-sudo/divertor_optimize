# Fast area-measure boundary-footprint surrogate

## User scope

Fast geometry optimization from an unchanged supplied D3D baseline and arbitrary supported hypothetical operating inputs. Measured agreement is explicitly out of scope. The user authorized diagnose → improve model → optimize. Earlier studies and code remain reproducible; new implementation is additive.

## Diagnosis

`boundary-diagnostic-20261028/diagnostic.json` reconstructs the old OMP320 ledger exactly: branch-tail 0.796088%, blocked 0.131611%, unhit 0%. Every blocked contribution is an anchor-to-cell fan extending outside the concave wall domain; neither LCFS nor solid intersection contributes. Only 0.343082% of input is unresolved on Tt<5 eV donors, so removing cold donors would not explain all missing power. These are exclusions of the old assumed receiver, not measured physical losses or invalid baseline geometry.

## Alternatives and choice

1. Keep the old same-branch kernel, expand/renormalize survivors: rejected, because branch multiplicity and masked chords have no complete arrival law.
2. Full anisotropic 2D transport: useful possible later reference, but not required for the requested fast reduced tool.
3. A closed connected boundary footprint diffusion kernel: selected. First trace the original donor to its first wall hit, then redistribute its footprint along the actual connected wall boundary. This is an explicit engineering closure, not a claim to reconstruct actual plasma trajectories.

Independent review rejected cumulative |dpsi/ds| as receiving distance because an isoflux plateau collapses physical distance and enables optimization shortcuts. Use cumulative toroidal surface area A(s)=integral 2piR ds instead. Positive-length positive-R cells always have positive measure. A periodic Gaussian on a single closed component integrates to one without rescaling surviving weights; wide kernels tend to uniform heat flux, not uniform per wall length.

## Operator and width convention

K(A|Ah,sigmaA)=sum_k exp(-(A-Ah+kAtotal)^2/(2 sigmaA^2))/(sqrt(2pi)*sigmaA).
For each donor, sigma_phys = sigma_flux*abs(partial_R psi_sep)/norm(grad psi_hit), and sigmaA=2piR_hit*sigma_phys. sigma_flux is the preserved one-way thermal path integral. q(A)=sum donor_power*K; each cell's power is the exact CDF interval mass.

This finite-width convention deliberately omits a grazing-incidence 1/sin projection. Geometry/incidence enters donor first hits and powers; sigma_phys is an effective along-wall footprint width. It is a NEW footprint model, not an equivalent implementation of the previous projected Eich kernel. Report this prominently. Do not cap broad widths or floor a zero target gradient: undefined conversions fail explicitly. Record gradient, physical/area widths, broad-kernel power and arrival distances.

## Numerical contracts

- Validate one simple closed positive-R boundary, no disconnected concatenation or passage through solid interiors. LCFS is not an absorbing material obstacle for this boundary model.
- Exact analytic area of each straight cell and of the partial hit segment. Invert area to position analytically.
- Narrow kernels: periodic Gaussian image CDF/density, with omitted probability bounded below double-precision closure tolerance. Broad kernels: controlled Fourier representation; uniform limit handled without subtractive cancellation. No empirical normalization.
- Unhit donor power remains unhit. All hit donor weights are preserved. Nonnegative probabilities and independent integral q dA are tested.
- Keep named-region power and donor-to-region transfers; report power moving over 0.1 m boundary arclength and width/total-area ratios. These expose the surrogate's nonlocal behavior instead of claiming actual paths were solved.
- Arbitrary positive finite P, n, chi, lambda, gamma and kappa are accepted within floating-point support. Numerical overflow/undefined geometry produces a clear error, not a silently adjusted input. Diagnostics about conduction regime do not block optimization.
- Optimization retains original geometry/angle constraints. Plate arclength is additionally integrated on a fixed high-order representation so coarse search cannot exploit the original polygon length underestimate.
- Cache unchanged candidate evaluations; source width is not restricted to the five historical cases. Baseline remains an available candidate even if a bounded search finds no improvement.

## Verification and delivery

Analytic kernel limits; independent numerical quadrature; cyclic seam/reversal; concave and flux-flat boundaries; conservation with unhit donors; unchanged file hashes; fixed-baseline legacy diagnostic; half/double and wider input sweeps; old/new timing; selected OMP shape reevaluated on refined tracing/wall mesh and alternate interpolation. Improvement is claimed only under this same new model and validated reference. No global-optimality or measured-agreement claim. User-facing main command and result summary must identify the model and include runtime.

## Final numerical refinements from implementation review

Continuous maxima are refined between donor/wall samples using normalized brackets. Cyclic seam arithmetic preserves small displacements. Converted widths below floating coordinate resolution and nonfinite output densities fail explicitly. Mean wall-distance estimates disclose max(cell length)/2 uncertainty; the probability beyond 0.1 m is integrated by exact partial-arc CDF intervals. Final acceptance also requires both wall-only refinement cases to pass the full audit, not merely a qmax comparison. The validator is included in frozen study provenance.
