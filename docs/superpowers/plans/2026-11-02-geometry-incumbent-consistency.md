# Geometry consistency and durable qualified candidates

Goal: use the final geometric screen before search heat evaluation and prevent validated candidates from being displaced by coarse rankings or interrupted writes. User authorized both changes on 2026-11-02. Existing D3D, inputs, signed curvature limits and old snapshots remain unchanged. Work in the existing D-drive repository on codex/geometry-incumbent-consistency.

Design: a shared versioned 1e-7 m chord-error screen checks full-wall simplicity, LCFS exclusion and vessel exclusion. Report specific failure reasons. Keep independent dense derivative auditing. A separate immutable, checksummed qualification bundle holds the manifest and selected arrays; the current manifest points to that bundle. Load the best compatible complete bundle on restart. Imported historical qualified geometries are protected revalidation candidates, never automatically qualified under new inputs/policy.

Rejected alternatives: loosening final geometry checks hides real collisions; enlarging the coarse archive alone still permits loss of valid incumbents. Neither is used.

- [x] Reproduce archived OMP20 collision using geometry only; save a small fixture and a failing early-rejection test. Test shared tolerance and positive examples.
- [x] Add GEOMETRY_POLICY and geometry_screen in campaign_geometry.py; use the identical screen in evaluate_candidate and geometric_audit, retain diagnostic reasons in describe.
- [x] Add campaign_incumbents.py: immutable bundles with context/policy/content hashes, atomic publication, corruption refusal, best compatible recovery; test worse proposal, wrong context, interrupted publication and checksum corruption.
- [x] Freeze historical qualified candidate vectors and origin hashes separately from coarse warm seeds. Test same-source/looser-radius handling. Include these and current incumbents independently of the four coarse finalist slots; revalidate before promotion.
- [x] Integrate Qualifier cached-policy checks, protected candidate priority and bundle promotion. Summary must read the manifest-selected bundle rather than potentially stale compatibility exports. Test finalization interruption/resume and legacy import behavior without expensive plasma calculations.
- [x] Run targeted regression and appropriate existing tests; verify old study hashes. Record geometry-only performance and exact reproduction. Explain optimality certificates: fixed domain, feasible upper bound, rigorous global lower bound, epsilon gap; existing fixed-polyline certificate is not a thermal geometry optimum proof.

Validation commands: python -m pytest verification/test_global_campaign.py verification/test_campaign_consistency.py verification/test_tangent_geometry.py verification/test_global_certificate_geometry.py -q. Add integration tests for search early gate and Qualifier selection/promotion. No new campaign launch, no modifications to expired campaign files. Large diagnostics go to D:/DivertorOptimization/studies/geometry-consistency-20261102.
