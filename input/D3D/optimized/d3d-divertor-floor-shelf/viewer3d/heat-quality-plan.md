# Heat mapping coverage and numerical verification

User approved the proposed coverage, extra tracing, independent surface integral,
and peak convergence checks in this side conversation. Work is confined to this
viewer and new derived files. Original optimization records remain the provenance.

1. Test adaptive source-cell selection, explicit connection statuses, independent
   axisymmetric surface integral, and peak/location/improvement convergence.
2. Add an isolated verifier importing preserved case math. Compare uniform source
   resolutions, then refine cells whose target arrivals are sparse or switch
   components, and cells near peaks. Inspect both baseline and retained candidate.
3. Independently invert target psi to a source coordinate and trace that coordinate
   to certify first-hit connectivity at fine surface cell midpoints. Integrate the
   resulting q over real axisymmetric cell areas; compare with assigned source
   power. Keep blocked/outside/unresolved states distinct; never infer zero from
   missing samples. Refine surface cells near connectivity transitions and peaks.
4. Publish separate quality artifacts with exact shape fingerprints. Show audited
   target patches, explicit status colors and tooltips; expose convergence and
   integral residuals. Historical candidates retain a clear un-audited label.
5. Run unit tests, execute verification, and exercise the local viewer in a browser.
   Report any failed numerical criteria rather than claiming physical validation.

Criteria: peak and candidate/baseline ratio changes <=1%, peak location <=2 mm,
target integral residual <=1% of that component's assigned power (absolute floor
1e-5 of unit source power), unmapped power <=0.1%. Shelf-only 1.5 deg policy
unchanged. These are numerical acceptance criteria, not experimental uncertainty.
