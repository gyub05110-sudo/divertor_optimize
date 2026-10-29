# Adversarial broadening plan review — 2026-10-27

## Disposition

Approved for implementation as an exploratory **connection-length-dependent deposition surrogate**, subject to the implementation and validation gates below. This is not approval of physical transport fidelity, optimized results, or numerical convergence. No results from the new model have been reviewed yet. The parent reports 41 original tests passing and a pre-edit preservation inventory captured.

Reviewed the design and implementation plan dated 2026-10-27, README, and `frozen_field.py`, `frozen_mapping.py`, `independent_sources.py`, `multisurface.py`, and `continuous_boundary.py`.

## Findings already resolved in the written design

1. **Poloidal length is not connection length.** Existing first-hit distance measures poloidal arclength. Resolution: integrate |B|/|Bp| to the actual candidate first hit, including the fractional final segment; preserve full source-to-hit scope in metadata.
2. **An unrestricted wall blur can move heat through material.** Resolution: restrict receiving support to the donor's connected monotone-flux wall branch, with a plasma-domain chord locality check. Explicitly describe this as a safeguard, not a solved diffusion trajectory.
3. **Renormalizing surviving kernels hides missing power.** Resolution: no renormalization; deposited, blocked, branch-tail and unhit components remain separately visible, with an unresolved-power feasibility gate.
4. **A binned heat flux is not the existing analytic direct peak.** Resolution: exact zero-spread dispatch to the original evaluator; nonzero spread uses a separately sampled continuous density with its own convergence requirements.
5. **The separatrix witness is not a finite donor.** Resolution: no assigned finite connection length or power at the exact separatrix; no unsmeared witness q added to broadened peaks. Geometry and incidence evidence remains separate.
6. **A length power law does not determine physical diffusivity.** Resolution: sigma_ref is phenomenological, lambda_q remains the prescribed source width, p=1 is the declared default and p=0.5 an optional sensitivity. Neither exponent is established by the current inputs. The chosen 0.1/0.25/0.5 mm strengths are exploratory assumptions at L_ref=10 m.
7. **Endpoint PFR support is not PFR transport.** Resolution: the model allows adjacent PFR-facing wall deposition without claiming PFR field-line transport, two-leg coupling, plasma response, neutrals, radiation or detachment.
8. **The objective cannot be globally smooth.** Resolution: document deterministic, piecewise-smooth behavior; first-hit changes, branch/visibility changes and a maximum objective retain nonsmooth events.

## Required implementation/validation resolutions

These are planned acceptance gates, not claims that implementation already satisfies them.

### R1 — Independent deposition-density accounting (high)

The component ledger can close by construction while q(s) is wrong. The CDF cell powers and sampled density must use the same donor/cell visibility and branch masks. Independently integrate q(s) 2*pi*R(s) ds and compare with deposited CDF power under quadrature refinement. Do not accept a ledger-only conservation test. Distinguish exact ledger tolerance (1e-10) from the declared, refined quadrature tolerance.

### R2 — Hidden flux extrema and first-hit branch ambiguity (high)

Endpoint flux ordering alone does not certify a monotonic cell. Check interior directional derivatives/roots and subdivide or reject unresolved extrema. Every powered donor must have an unambiguous branch at its actual first hit; ties at branch boundaries need a documented deterministic rule or rejection. Include hidden-extremum and hit-at-extremum tests. No seam wrapping without an explicit topology proof.

### R3 — Visibility and peak masks must remain consistent (high)

Testing only receiver endpoints does not prove every interior chord is inside plasma. Subdivide any uncertain receiver cell or conservatively reject it; at minimum perform an independent finer-cell audit for selected candidates. Do not mask CDF mass one way and density/peak points another way. Report blocked mass and topology/mask changes rather than silently treating them as diffusion losses.

### R4 — Source quadrature and the zero-spread limit (high)

Discrete Gaussian donors can create artificial peak spikes when sigma is small compared with source-cell spacing. The exact sigma=0 dispatch alone does not establish convergence for sigma>0. Refine the source grid and receiver/peak grid separately at the smallest configured spread, including near-separatrix cells where L grows without bound. Show stabilized q/P and shape ranking, not just total power. Never silently cap L.

### R5 — Locality-anchor sensitivity (medium)

The 5 mm poloidal anchor is a numerical/model choice, not a physical transport parameter. Report its value and assess selected candidates with at least two bracketing values (for example 2.5 and 10 mm). Anchor-dependent admissibility or ranking must be reported as surrogate sensitivity. Do not tune the anchor to make a preferred candidate feasible.

### R6 — Powered receiving incidence and one-sided peaks (high)

Apply the explicit shelf incidence floor to every genuinely powered receiving region, not merely old first-hit rays or design segments. Use the same positive-power policy as the original model, and document any numerical tail truncation. Include one-sided samples at branch ends, visibility boundaries and wall corners in peak refinement; the unsmeared SOL witness cannot certify these broadened peaks.

### R7 — Fair comparison and failed baseline handling (high)

Preserve matched budgets and seeds, and report each candidate under every configured model and source width. A baseline failing the surrogate missing-power or visibility gate makes a baseline-relative improvement unqualified; do not renormalize it, relax its gate, or silently omit the failed strength. Report feasibility separately from ratios. Claimed improvement requires the written refinement/variation gate and no worsening at any prescribed lambda.

### R8 — Independent first-hit and length validity (high)

Compare cached candidate hits with an independent complete-wall calculation for selected shapes. Test length quadrature in a varying field in addition to the constant helical analytic case, including fractional final segments. Reject invalid Bp, nonfinite metrics, inadequate tracing, or absent first hits explicitly. Near-X-point path-step refinement is required; a constant-field test cannot expose it.

## Scope limit and practical choice

The endpoint surrogate is the smallest defensible scoped implementation for the current code. A bundle-marching alternative would require new transverse metrics, PFR coverage, material-clipped exchange interfaces and live/absorbed tube bookkeeping. Implementing a simplistic marcher on the current SOL-only ray set would create a stronger but unsupported physical appearance.

The approved deliverable is a reproducible sensitivity comparison with explicit validity limits. If all configured spreads or optimized candidates fail numerical/locality screens, that is a valid finding and must remain visible in the report; it is not permission to claim a successful transport optimization.

## Primary background consulted

- https://arxiv.org/abs/1306.6777 — measured exponential and Gaussian footprint components and an approximate diffusion interpretation.
- https://arxiv.org/abs/1610.04148 — conductive spreading approximation involving perpendicular/parallel diffusivity.

These support distinguishing upstream width from divertor spreading, not the detailed topology or visibility assumptions of this new surrogate.

## Independent specification/physics implementation review

Reviewed `broadening.py`, `broadened_evaluator.py` and their targeted tests after initial implementation. Initial implementation disposition: **changes required before accepting numerical results**. No production code was edited by this reviewer. Targeted test run during concurrent repairs: 30 passed, 1 skipped in 6.24 s; passing tests alone did not expose the counterexamples below.

### I1 — Thin solid between visibility probes (high)

**Confirmed counterexample.** In the linear-field box fixture with 25 mm bottom cells and anchor (2, 0.005), put a thin solid rectangle at R=[1.9529,1.9533], Z=[0.0024,0.0026]. The ray to receiver (1.90625,0) intersects the solid, but the initial endpoint/midpoint checks returned blocked=0 and valid=True; that point received q/P=0.2146338785 and its cell received power 0.0679720984. This violated the intended locality safeguard.

**Applied repair observed in code:** visibility now tests the full triangular fan formed by anchor and cell endpoints against the plasma domain and obstacles, with a line representation for degenerate fans. The fan is the union of all chords to a straight receiving cell, so one common mask can safely apply to CDF mass and every density sample. A thin-obstacle regression was added. This resolves the structural flaw; it remains a geometric locality test, not a physical diffusion solution.

### I2 — Grazing receiving samples excluded by q mask (high)

The initial adapter used q>0 to identify powered incidence samples. At exact grazing, q vanishes through the flux Jacobian even though the receiving cell has positive integrated Gaussian power. A zero-angle endpoint could therefore evade the explicit shelf floor.

**Applied repair observed in code:** incidence now selects samples through their receiving cell's positive integrated power, including zero-q endpoints. A targeted grazing-sample regression was added. This restores the original power-based constraint policy.

### I3 — Paired hidden flux extrema accepted as one branch (high)

**Confirmed counterexample.** Let t=R-1 and psi=t^3/3-0.01*t^2/2+0.000024*t. On a box bottom R=[1,3] with wall_step=0.15, derivative roots R=1.004 and 1.006 both lie inside the first cell [1,1.142857142857]. All initial and final sparse probes miss the reversed interval. `_wall_mesh` returned one bottom branch instead of three. CDF interval power then need not match integration of the absolute Jacobian, and the topology restriction is incorrect.

**Accepted repair in progress, not yet re-reviewed:** split wall intervals at spline R/Z knot crossings and determine all directional-derivative polynomial roots within each spline patch, with independent residual checks and the paired-root regression. Sparse sign probes alone cannot establish monotonicity. This remains the main numerical acceptance blocker at this review checkpoint.

### I4 — Exact-vertex length includes unused downstream metric (medium)

`PathLengths.at` initially selected the next interval at an exact hit vertex, requiring its downstream metric to be finite even when that metric has zero integration weight. This contradicts the stated rule that invalid samples beyond the hit cannot invalidate a valid earlier hit. It causes conservative false rejection rather than incorrect deposited power.

**Accepted repair in progress:** choose the left interval at exact nonzero vertices, handle distance zero explicitly, and add an invalid-after-hit regression. Awaiting re-review.

### I5 — Final qualification gates distinct from coarse feasibility

The initial adapter exposes independent Boole density integration but does not make coarse quadrature error a final qualification claim. Reported coarse density/CDF discrepancy around 2.6e-4 is larger than the proposed final 1e-4 gate, so those coarse evaluations cannot alone establish qualified results.

**Parent-reported validation additions, not independently reviewed at this checkpoint:** max pre-hit flux-coordinate drift/sigma <=0.05, independent Simpson connection-length relative discrepancy <=0.005, density integral/CDF discrepancy <=1e-4, absolute q refinement discrepancy <=0.005, and peak-ratio discrepancy <=0.002. These supplement the design's improvement margin and all-width nonworsening checks. Recheck selected candidates and baselines on refined levels, retaining failure flags.

Reported preflight baseline failures at OMP sigma_ref=0.25 mm and both sources at 0.5 mm are legitimate limits of this surrogate and must remain explicit. A baseline-invalid setting cannot support an accepted baseline-relative optimization improvement, even when its numerical peak ratio is finite.

## Final kernel/adapter re-review — ready to run

**Disposition: ready to run the declared sensitivity study.** No remaining important specification/physics implementation blocker was identified in the repaired kernel and adapter. This is permission to generate and assess numerical evidence, not a finding that selected shapes or this physical closure are validated.

- **I1 resolved:** whole-cell triangular fan checks preserve the same material-aware mask for CDF power and density; the reproduced thin-solid case now has a regression.
- **I2 resolved:** receiving incidence selects positive cell power, retaining zero-q grazing endpoint samples. The adapter regression covers the exclusion failure.
- **I3 resolved for the supported FrozenField interpolants:** wall segments split at spline knot crossings; degree-at-most-5 cubic or degree-at-most-9 quintic directional derivatives are reconstructed in Chebyshev form, tested at disjoint nodes, and rooted when a sign bound does not suffice. Root residuals are checked. Tests cover the paired hidden extrema, actual cubic/quintic FrozenField splines, knot boundaries, and rejection of an unresolved generic nonpolynomial field. These are floating-point numerical checks, not a mathematical proof against arbitrary ill-conditioned roots.
- **I4 resolved:** exact nonzero hit vertices use the preceding interval, and zero distance is explicit. The invalid-downstream-metric regression passes.
- **I5 gates independently read in verifier:** density/CDF power <=1e-4; midpoint-Simpson versus cached-trapezoid connection length <=0.005; path flux drift/sigma <=0.05; absolute peak refinement <=0.005; peak-ratio refinement <=0.002. Complete-wall first hits, independent source normalization, fixed geometry and provenance checks are also present. The final result still depends on actually executing these checks on baseline and selected candidates.

Fresh reviewer test execution: `python -m pytest verification/test_broadening.py verification/test_broadened_evaluator.py -q` — **35 passed, 1 skipped in 5.97 s**.

Two reporting caveats remain intentional rather than blockers:

1. Acceptance currently requires every one of the four refinement-level audits to pass, including the coarse level. A coarse density error above 1e-4 therefore prevents qualification even if later levels converge. Report this as insufficient numerical qualification, separately from physical/surrogate feasibility. Do not retroactively relax the rule to obtain an improvement claim.
2. Anchor sensitivity is a separate robustness flag. A numerically confirmed candidate that changes materially with anchor choice must not be described as robust to the locality closure. The exponent sensitivity is likewise a model-assumption study, not new calibration evidence.

## Adaptive validation resolution review

The parent reported an initial baseline audit at 320 donors and 0.002 m path step with maximum cached-trapezoid versus independent-Simpson length discrepancies of 0.02182 (IMP) and 0.01972 (OMP), exceeding the unchanged 0.005 gate. The reported density/CDF errors were 6.04e-5 and 1.23e-6; flux-drift/sigma was 0.0122 and 0.0025. These figures are parent-reported execution evidence, not independently rerun by this reviewer.

**Approved numerical response:** use the existing `validate_source(..., levels=...)` API to run an explicitly recorded, finer path-step matrix:

```python
levels = [
    ('cubic', 320, 0.0005, 161, 0.002),
    ('cubic', 640, 0.00025, 321, 0.001),
    ('cubic', 1280, 0.000125, 641, 0.0005),
    ('quintic', 640, 0.00025, 321, 0.001),
]
```

This is scientifically legitimate adaptive resolution after a failed error check, subject to applying the complete matrix uniformly to all saved candidates, baselines, models and both sources. Keep the model parameters, thresholds and acceptance logic unchanged. Preserve the original failed diagnostic, rationale, exact executed invocation and recorded actual levels in the execution record. The proposed matrix is not itself evidence that the checks pass.

The parent reports that the main search uses 320 donors and a 0.001 m path step. Describe that search resolution accurately and distinguish it from verification resolution. If finer verification changes candidate ranking, report the change; do not select a convenient partial set of refinement levels. Frozen pre-search code hashes can remain intact because this change uses an existing explicit API parameter rather than a code modification.

## Independent code-quality integration review

Reviewed the additive comparison runner, evaluator adapter, study persistence/provenance helpers, validation, plotting, and their tests against the design and plan. The physics/kernel review above owns the physical assumptions and repaired branch-topology issues. No preserved production code was edited by this reviewer.

### Q1 — Preserve original failed screens through refinement (important, resolved)

The original validator formed its final verdict only from newly evaluated refinement levels. It did not include the saved search/preflight baseline admissibility and cross-evaluated candidate feasibility. Consequently, a failed original screen could be replaced by a positive refined verdict, contrary to the requirement that failed screens remain explicit.

**Repair independently observed:** `gate_original_screen` now gates each model/candidate verdict using the original baseline and candidate evidence. An original failure prevents `numerically_verified_feasible` and `improvement_confirmed`, sets `verification_status=original_screen_failed`, and retains `refined_verification_status` separately. A regression verifies that successful refinement cannot erase a failed original baseline.

### Q2 — Refuse validation artifact overwrites (important, resolved)

A repeated verification previously replaced `validation.json` and refined arrays, allowing an interrupted rerun to destroy a completed verdict.

**Repair independently observed:** `validate_source` now checks for an existing `validation.json` or any model `*_refined.npz` before provenance checks and numerical work, and raises `FileExistsError` for completed or partial validation. A regression confirms that an existing completed validation remains unchanged after the rejected invocation.

### Q3 — Respect requested bounds for historical warm starts (important, resolved)

The OMP historical candidate has a maximum absolute coefficient of approximately 80.0781 mm. The initial runner accepted it without checking the requested coefficient bound, so a legal smaller bound such as 50 mm raised `ValueError` during search or cross-evaluation. The failure was independently reproduced using the actual historical vector and evaluator.

**Repair independently observed:** `select_warm_start` excludes an out-of-bound historical vector and records `outside_requested_coefficient_bound` in manifest metadata. It does not clip or otherwise change the historical candidate. A regression checks both exclusion and unchanged inclusion within the bound.

### Fresh verification and disposition

After the three repairs, this reviewer ran:

- `python -m pytest verification/test_broadening_study.py verification/test_broadening_cli.py -q`: **13 passed**.
- `python -m pytest verification/test_broadening.py verification/test_broadened_evaluator.py -q`: **35 passed, 1 skipped**.

Combined reviewed test scope: **48 passed, 1 skipped**. This is the targeted additive integration/kernel test scope, not a claim about a fresh full-repository suite.

**Final code-quality disposition: approved for bounded exploratory execution and independent numerical validation.** No remaining important or critical issue was found in the reviewed scope. This disposition does not establish numerical convergence, physical validation, a global optimum, or a confirmed improvement. Actual search and refinement results must retain every failed screen and satisfy the declared acceptance gates before any improvement is qualified.

At the time this review was appended, the parent reported a 320-tube search in progress and a uniformly finer trace-step validation planned after a 2 mm trace-step length check failed. Those execution outcomes have not been reviewed here; no pass is implied by this code-quality approval.

## Post-run presentation repair review

After both complete numerical validations, full historical provenance was verified and the original plotting source archived. A report-only repair handles incomplete checkpoints, labels search versus refined unresolved power, and marks anchor-sensitive gains as conditional. Two regressions were observed failing before the repair. The independent code reviewer reran them: **2 passed**, confirmed the archived source hash, and found **no blockers**. The reviewer also read the completed evidence: IMP 5.8600%, IMP 3.7488% conditional on anchor, OMP direct 35.7085%; OMP broadening remains unqualified with refined baseline unresolved power reaching 0.260291%. Historical provenance is unchanged; the explicit presentation amendment records the only later hashed-file edit.
