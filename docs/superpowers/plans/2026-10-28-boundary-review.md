# Adversarial review and implementation decisions

User objective: fixed supplied D3D baseline plus arbitrary supported hypothetical inputs → fast geometry optimization. No measurement-agreement gate. Original files are preserved.

## Design review

- Legacy OMP exclusion was reconstructed without altering geometry or thermal inputs. Same-branch Gaussian tails and an anchor-to-cell fan outside a concave wall caused the missing-power ledger; LCFS/solid hits did not cause these exclusions.
- Expanding a branch and renormalizing surviving weights was rejected: no complete arrival law.
- A cumulative absolute-flux boundary metric was rejected: isoflux wall plateaus collapse the coordinate and let a geometry search exploit vanishing distance.
- Selected cumulative toroidal area A=∫2πR ds on one simple connected closed boundary. The periodic CDF supplies a complete normalized kernel without survivor correction. Report this as a new engineering footprint surrogate, not reconstructed plasma trajectories.
- The width convention omits grazing-angle projection and is explicitly different from the prior projected kernel. Thus old/new differences cannot be attributed solely to filling missing mass.

## Implementation review and fixes

1. Narrow Gaussian mixtures can peak between donor-local samples. Added continuous local-maximum refinement in normalized brackets; the regression uses two close donors with a known analytic mode. Independent reviewer swept 150 random heterogeneous mixtures against dense donor-adaptive samples and found no underestimate. Final numerical refinement remains required; no formal global-mode proof is claimed.
2. Arithmetic at the cyclic seam lost small displacements. CDF image endpoints now subtract the center after forming identical periodic endpoints; density uses conditional displacement wrapping without adding and subtracting P/2. Reviewer seam sweep σ/P=1e-12…100: maximum mass error 2.22e-16.
3. A tiny positive sigma could yield infinite q despite a conserved CDF ledger. Reject widths below floating-coordinate resolution; require finite q, cell power, quadrature and total donor power. Never floor/cap widths silently.
4. Cell-midpoint distance estimates were misleading for narrow footprints. Disclose max(cell length)/2 error bound for mean wall distance; evaluate power beyond 0.1 m using exact partial-arc area maps and CDF intervals.
5. Original polygon arclength could allow a candidate to pass a 1.2 length ratio at coarse shape resolution and fail when refined. Integrate the continuous spline with independent 16/32-point Gauss rules per knot span before expensive thermal evaluation.
6. A wall-only refinement initially checked only qmax. Review exposed a potential missed receiving-angle violation. It now runs the full audit for both baseline and candidate, and the final verdict requires both audits and feasibility. Added a regression that rejects a wall-only infeasible candidate.
7. Include the validator in the frozen model provenance so a changed acceptance policy cannot silently attach to the same saved study.

The identified defects had failing regression tests before correction. The final full suite was run again after fixes: 174 passed, 1 skipped. All 169 pre-existing Python files and both prior reference-study provenance records were checked unchanged. Main search and final benchmark provenance were verified after the run.

## Evidence directories

- `input/D3D/optimized/boundary-diagnostic-20261028`: legacy route diagnostic
- `input/D3D/optimized/boundary-optimization-20261028/OMP`: frozen final search and independent refinement
- `input/D3D/optimized/boundary-benchmark-final-20261028`: frozen same-baseline comparison and 10 input variants
- `input/D3D/optimized/boundary-input-smoke-20261028`: independently specified composite inputs, IMP and OMP smoke searches

`boundary-benchmark-20261028` is an earlier exploratory run captured before the final validator fix. Its evidence is retained, but the final report uses `boundary-benchmark-final-20261028`.
