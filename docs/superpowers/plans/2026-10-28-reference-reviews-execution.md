# Reference transport: review and execution evidence

## Authorization and preserved state

The user approved choosing a hypothetical fixed operating point and completing the comparison. All 46 pre-existing Python files are inventoried in `2026-10-28-reference-preservation.json` and remain byte-identical. No Git checkout exists. New source, tests and outputs are additive. Neither 2026-10-27 numerical results nor their provenance is overwritten.

## Physics review

Independent reviewer approved a one-way effective one-temperature two-point conduction closure for a hypothetical point, with explicit restrictions. Source q=P exp(-x/lambda)/normalization, not donor power divided by target area. The electron heat capacity is 1.5 e n; chi_perp=.3 uses this convention. Bulk density at target is twice the separate sonic sheath density; never smear the pressure jump across a cell. Temperature diagnostics are pre-redistribution; final Gaussian wall loads do not feed back to the sheath solution. Full-path/ray-dependent mapping extends Nille's model rather than replicating it. Cold/weak-conduction diagnostics cannot establish detachment or flow validity.

Reviewer supplied exact constant-B spreading regression. The initial poloidal Gauss quadrature under-resolved a very cold target by 0.34% in an added analytic regression. Transforming to y=T^(7/4) removes the integrable endpoint singularity without floors or caps:

```
S = 2 C J / (Tu^(7/4)+Tt^(7/4)) * mean_y[g Bu/B]
C = sqrt(1.5 e nu Tu chi_perp/kappa0)
```

The inverse cumulative magnetic integral locates quadrature points. The cold regression failed before the change and passed after it.

## Code review and repairs

Independent code review identified two important issues, both reproduced with failing regressions and fixed before the main run:

1. Interior quadrature missed the bulk target maximum of the free-streaming diagnostic for a cold tail. Exact upstream and pre-sheath bulk endpoints are now included. The analytic bulk target ratio is gamma/2 sqrt(me/mi).
2. The nominal 3 mm adapter initially inherited the original five-width missing-power/feasibility aggregate. It now reconstructs nominal-width power, angle, boundary and illumination constraints. A failure at an unused source width cannot reject the nominal scenario.

Actual deposited-width statistics are distinguished from recomputed thermal-width diagnostics in frozen mode. The new explicit-width kernel was compared against the original reviewed kernel: q, CDF mass, geometry and all ledger arrays match exactly for equal widths.

The integration reviewer found the scalar acceptance, original failed-screen gate, before-search provenance, no-overwrite guard and fixed-candidate sensitivity implementation sound. A separate constant-field fractional-hit independent integration check agreed with production sigma, L and J. Final disposition: approved for freeze and bounded execution, no important/critical blocker. Report now includes free-streaming diagnostic and candidate-by-model cross-evaluation.

## Initial numerical evidence

Smoke preflight: `input/D3D/optimized/reference-preflight-20261028`, 48 paths. This is wiring evidence only.

Independent 320-path audit at 0.5 mm tracing and 1 mm receiving spacing is saved in `input/D3D/optimized/reference-initial-audit-20261028`:

- IMP: feasible, audit passed; unresolved power 5.37863e-6; sigma relative discrepancy 3.68473e-5; L discrepancy 0.00126325; magnetic heat integral discrepancy 0.00124036; density/CDF power difference 1.82523e-9.
- OMP: baseline infeasible, unresolved power 0.00927699, exceeding the unchanged 0.001 gate. Sigma/L/J integration discrepancies were 5.97063e-5 / 0.000871751 / 0.000930261. Numerical integration accuracy does not rescue deposition inadmissibility. Thermal/frozen optimization is excluded; the point is not tuned to make OMP pass.

Before main run: full `python -m pytest verification -q` **132 passed, 1 skipped in 21.93 s**. The final presentation-only additions were separately checked by the report regression before freeze.

## Main run, fixed before candidate search

```
python compare_reference_transport.py --mode optimize --source both --tubes 320 --nodes 161 --step-mm 0.5 --wall-step-mm 1 --quadrature 128 --sobol-power 7 --population 24 --generations 8 --output input/D3D/optimized/reference-comparison-20261028
```

P=1MW per independent source, nu=1e19 m^-3, electron chi_perp=.3 m²/s, lambda=3mm, gamma=7, kappa0=2000. Equal seed/proposal budget and historical warm starts for each model. Actual evaluation counts and early termination are saved, not assumed equal. Selected shapes remain bounded exploratory candidates, not proven global optima.

Verification uses the declared `LEVELS` matrix in `verification/validate_reference_transport.py`, with independent full-wall first hits, Simpson-sampled thermal integrals, area-power integration and fixed geometry. Sensitivities halve/double P, n or chi one at a time and vary the local deposition anchor at fixed shape. Final outcomes are recorded after completion.

## Completed search and sensitivity evidence

Output: `input/D3D/optimized/reference-comparison-20261028`. Main search exited successfully. IMP elapsed 1321.9 s; OMP 58.9 s. IMP direct stopped after 153 distinct evaluations with the baseline; frozen and thermal each evaluated 321 candidates and selected exactly equal coefficient vectors. Raw frozen/thermal ratios were 0.871194423/0.871216185; maximum plate displacement was 7.343126 mm. These are search values, superseded by refined results below when available. OMP frozen and thermal correctly skipped search because their original baseline was inadmissible.

OMP's new direct candidate fails the refined geometry gate. An additive diagnostic `OMP/selected_geometry_diagnostic.json` identifies the floor arc-length ratio: 1.199998044 at 161 nodes, 1.200048004 at 321 and 1.200060500 at 641. The pre-existing 1.2 limit was not relaxed. The historical direct candidate does pass all four validation levels, with a nominal-3-mm reduction of 40.617395%; this is a previously available candidate under the new nominal objective, not a new transport result.

IMP fixed-shape sensitivities at cubic640: all reference baselines are feasible. The selected thermal candidate is feasible at double power (ratio 0.914996), half density (0.938820), and half chi (0.897067). It fails the unchanged 0.1% unresolved-power gate at half power (0.199640% unresolved), double density (0.361851%), and double chi (0.210980%). A 2.5 mm deposition anchor also fails (0.400590%), while 10 mm passes (0.021184%). These are single-resolution diagnostics, not reoptimized or fully validated operating ranges. A lower peak in a failed variant is not accepted evidence.

After the final presentation changes and code freeze, full `python -m pytest verification -q`: **132 passed, 1 skipped in 27.67 s**. All 46 original Python hashes and both before-search code/input inventories matched. Numerical implementation and report generator remain byte-identical to the main run's captured versions; later changes are documentation only.

## Independent interpretation review

The physics reviewer inspected the saved search, sensitivity and available validation JSONs. Frozen-vs-thermal ratio differences reverse sign between 320 (+2.1762e-5 thermal minus frozen) and 640 (-2.6590e-5), while the coarse-to-640 thermal ratio changes by about 0.00307. Thus the extra effect of recalculating transport is not resolved; identical selected coefficients should be reported without claiming a benefit or harm from thermal-state recalculation. Inclusion of broadening changes the preferred shape in this bounded search, but this is a different claim from thermal-response optimization.

The reviewer also emphasized that the 2.5/10 mm anchor sensitivity changes feasibility while the peak ratio remains identical to displayed precision: peak insensitivity does not establish design robustness. The OMP failed length constraint concerns material plate arc length, not magnetic connection length. With conduction dominance unverified, the study remains a conditional numerical comparison rather than validated plasma transport or global shape optimality.

## Final verification and delivery

Both independent validation processes exited 0 and saved complete four-level records plus all eight sensitivities. IMP selected frozen and thermal candidates pass every geometry/power/independent-integration screen and the predeclared resolution criteria. At cubic1280:

- Thermal baseline peak 2.541298766 MW/m²; selected 2.203046640 MW/m²; reduction **13.310207%**; unresolved selected power **0.094308332%**.
- Frozen selected peak 2.203114020 MW/m²; reduction **13.307556%**. Selected coefficients are identical to thermal; maximum plate displacement 7.343126 mm.
- Thermal maximum ratio change 0.001247784 <0.002; maximum absolute-peak relative change 0.001432627 <0.005. All four rows are feasible and independently audited.
- Thermal power-weighted L: 27.903091→27.988542 m; sigma: 0.251482829→0.251463785 mm. This modest state/width response does not support a resolved additional optimization advantage over frozen widths.
- IMP direct selected the baseline; no improvement claim. OMP new direct candidate remains rejected and thermal/frozen remain originally inadmissible; the historical direct comparison is recorded separately above.

The summed four-level evaluation times were 1257.19 s (IMP) and 682.27 s (OMP); these sums exclude the separately executed sensitivity block. Source validation records, combined `validation_summary.json`, regenerated `REPORT.md` and figures, plus `INTERPRETATION.md` preserve the final findings. Source/input provenance was reverified after both validations, and all 46 original Python files remain unchanged. No main-run numerical or plotting source was changed after capture. The full post-freeze suite result remains 132 passed / 1 skipped; subsequent edits are documentation only.
