# Additional heat mapping verification

This adds numerical coverage checks for the registered baseline and retained final
shape. Search logs, original manifests, coefficients, and archived math are intact.
Only viewer code and new derived quality artifacts are affected.

The source remains unit outer-branch power at the fixed 1400 ms equilibrium,
with widths 1, 2, 3, 5, 10 mm. The shelf-only 1.5 degree constraint is unchanged.
Uniform source resolutions 320/640/1280 are followed by adaptive refinement of
large arrival gaps, component transitions, and high-q cells. The final shared
source partition contains 2502 cells.

An independent target-space check evaluates surface cells initially <=0.5 mm
long, then subdivides near peaks/connectivity boundaries and inadequate integrals.
For each midpoint it inverts psi at the upstream plane and traces that exact
source coordinate. A point receives q only when this trace reaches the point
and the point coincides with the first material hit within 0.1 mm. Integration
uses q times the actual axisymmetric surface cell area, independent of the
source tube's assigned weight. The smallest refined cells are about 0.1 mm.

Status colors identify direct hits, another surface hit first, outside the
prescribed source range, and unresolved trajectory connectivity. Non-direct
patches retain null q; they are never silently converted into known zero heat.
The direct-heat integral includes only certified direct patches. Any unresolved
target length prevents an overall pass. Background wall gray remains contextual;
other wall heat retains finer forward samples rather than a target-space audit.

Checks passed for this case:

- Both shapes: source accounting, peak/position/ratio convergence, independent
  target surface integrals, surface integral convergence, and shelf angle scope.
- No unresolved target patches under this finite outer-branch model.
- Maximum target integral residual: baseline 0.2632%, retained final 0.2618%.
  Small components use an absolute denominator floor of 1e-5 of input power.
- Fine q/P peak ratios: 0.60374445 for widths 1–5 mm, 0.64532018 for 10 mm.
  Reductions are approximately 39.6% and 35.5%, respectively.
- Added tests and Chrome coverage/status/history switching checks pass.

Numerical thresholds are 1% for peak and improvement-ratio changes, 2 mm for
peak position, 1% for target integrals and integral refinement changes, and 0.1%
for unmapped source power. Every criterion must pass before `verified` is true.
The independent target probe's peak and location also agree with fine source
mapping within the same thresholds. Shape fingerprints prevent stale audits
from attaching to different history snapshots. Historical intermediate candidates
remain explicitly marked as not additionally audited.

These checks concern the stored polygon geometry and prescribed direct heat
mapping. They do not refine the underlying geometric spline or certify experimental
heat flux, radiation, cross-field transport, or campaign hardware.

## Reproduce

From this directory, run `python verify_heat_quality.py` in an environment with
numpy/scipy/shapely/pyarrow. On this host, pyarrow is in the optional inspection
directory `%TEMP%/codex-d3d-inspect` and may need that directory on PYTHONPATH.
The verifier imports the preserved case math under an isolated namespace.

Generated artifacts are `quality/<case-id>.json` and `quality_catalog.json`.
`python quality_result.py` rechecks stored results and shape compatibility without
re-running trajectories. Tests: `python -m unittest test_heat_quality.py
test_quality_result.py test_history_data.py`; `node --test *.test.mjs`.
Chrome UI validation: `node quality_browser_check.cjs` against the existing server.

The original history score is labelled as the score recorded during search.
The additional verification score is displayed separately; primary heat metrics
for audited shapes use the finer result.

## Source range extension audit

`audit_source_domain.py` checks the public wall at the upstream plane against the
EFIT grid. The grid permits 288.68885 mm beyond the separatrix, but the continuous
plasma-side interval ends at the public wall after 98.77583 mm. Proposed 150/200 mm
sources start in material/exterior and are not physically valid for this case.
No extended outside-wall trajectories were launched. The legacy 100 mm source
already includes the complete contiguous plasma-side upstream interval. Its small
outside-wall weight is retained in the original explicitly unmapped accounting.

All target patches previously flagged outside the source range have psi on the
separatrix side of the upstream source start. None correspond to flux beyond the
100 mm upper cutoff. Widening the source therefore cannot resolve those patches.
The viewer now explains this distinction and includes the exact physical limit.
Patches retain null q; separatrix-side flux alone does not establish zero actual
heat or uniquely identify private-flux versus other magnetic regions.

An optional `--source-max-mm` is now accepted by the verifier. Explicit ranges
must fit inside the public wall and EFIT grid; 150 mm is rejected before tracing
or changing result files. Omitting the option retains the original source for
reproducible comparison.

The audit also integrates a hypothetical exponential source out to 200 mm using
4096 fine intervals (4-point Gaussian quadrature per cell), purely to quantify
tail weights. This is not an expanded physical plasma calculation. At width
10 mm, the 98.77583–100 mm interval contains about 5.64e-6 of that hypothetical
source power; 100–200 mm about 4.30e-5. Both are outside the public material wall.
