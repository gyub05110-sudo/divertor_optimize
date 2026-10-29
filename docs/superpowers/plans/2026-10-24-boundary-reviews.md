# Implementation review before code changes

## Plan review

- Coverage: source separation, no transport, user Z cutoff, excluded walls, original optimization recipe, legacy replay and new result paths are explicit.
- Historical kernels remain unchanged. A new adapter owns the correction. Strict geometry correction is already authorized; 5 mm extension is awaiting the user's answer.
- Boundary change is specific to IMP segments 69/70. OMP sampled region boundaries are retained; no claim that every OMP footprint boundary has been reconstructed continuously.
- Tests precede implementation. Numerical acceptance is defined before examining new results.

## Adversarial review

1. **A continuous endpoint still has a fixed tangent.** Fix: include the SOL-side separatrix limit in every search objective. Strict IMP may legitimately give zero improvement; do not manufacture gain.
2. **A historical strike can be shadowed by new geometry.** Fix: use the first intersection of an extended leg with the entire candidate wall before refining its location. Never use an invariant lower bound without checking visibility.
3. **Polyline normals can create false convergence.** Fix: analytic spline derivative at the refined witness, while ordinary finite source rays continue to use the preserved mapper. Independently refine both ray and target resolution.
4. **A numerical witness may be unresolved or on the wrong branch.** Fix: check the refinement is near the actual first intersection and has small ψ residual. An unresolved witness rejects a candidate; it cannot supply zero heat.
5. **A separatrix witness does not bound all interior maxima.** Fix: keep all source rays, test 1280/2560 and quintic, explicitly label measured convergence rather than a global certificate.
6. **Boundary witness could bypass shelf incidence policy.** Fix: apply 1.5 degree floor when the witness hits the shelf, in addition to finite-power-ray checks.
7. **Apparent power conservation may be just normalized bookkeeping.** Fix: full first-hit and fixed-geometry checks plus normalization comparison; never call this experimental validation.
8. **Small improvement could be discretization noise.** Fix: require gain greater than twice observed refinement/interpolation changes plus 0.0001 before setting improvement_confirmed.
9. **User scope could silently expand.** Fix: strict defaults; any PFR support buffer is an explicit run parameter awaiting the answer. The Z cutoff and excluded green wall are unconditional.
10. **Result cannot be reproduced after code changes.** Fix: saved boundary mode/buffer and input hashes, versioned objective, no overwrite of historical runs. Validation reconstructs from saved policy and fails on differences.

Disposition: implement strict correction now. Buffer-dependent search waits for the user's design answer. The review accepts a numerically sound zero-gain IMP result; it does not promise a global optimum or experimental agreement.
