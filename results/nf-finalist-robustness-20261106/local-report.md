# OMP same-640 discrepancy on the frozen wall

At χ = 0.3 m²/s, a roughly 0.1 mm shift in the landing points of near-separatrix donors produces the 2.6930% local discrepancy on the candidate floor. The 1% v3 limit fails. The maximum occurs on the falling side of a local peak, at wall arclength 4.948040 m, (R, Z) = (1.225701, −1.367985) m. It falls between the two interpolants’ floor intersections with their respective separatrix flux levels, at 4.947967 and 4.948056 m. Wall arclength starts at the first stored vertex and follows the archived wall orientation.

## Where the local change comes from

At that location q rises from 0.752453 to 0.791321 MW/m². Each method supplies the saved first-hit centers, physical standard deviations σ and independently reconstructed source weights. We evaluated all eight cubic/quintic combinations of these three inputs. For each donor the table averages the signed effect over the six possible orders of substitution. The center factor includes the hit radius in the production area width 2πRσ. The width factor changes the physical σ. All substitutions use the same physical coordinates.

| Quintic minus cubic (MW/m²) | At the v3 maximum | Native cell average |
|---|---:|---:|
| First-hit centers | +0.03891214 | +0.03872960 |
| Physical widths | −0.00005560 | −0.00014169 |
| Source weights | +0.00001135 | +0.00001221 |
| Full kernel difference | +0.03886789 | +0.03860012 |

Native cell 10709 covers s = 4.947724–4.948068 m on floor segment 1274. Its saved cell averages are 0.811258 and 0.849858 MW/m², the largest native difference. The pointwise kernel difference agrees with the piecewise-linear v3 discrepancy within 5 × 10⁻⁸ MW/m² there. Source weights were rebuilt by integrating the production source on the equilibrium independently for each interpolant. They match the archived donor-region sums within 5.3 × 10⁻¹⁸ of assigned power. The wrapped Gaussian interval masses reproduce every exported native cell power exactly in this run, across all eight exports.

* Zero-based donor indices 8–25 span launch offsets 0.01748–0.15687 mm from the OMP separatrix and contribute +0.037075 MW/m², 95.4% of the signed point difference. The ten largest contributors are indices 12–21. They contribute +0.028188 MW/m² and their landing positions move 93.0–99.9 μm in the same arclength direction.
* Donor 16 contributes +0.003402 MW/m², the largest individual change. Its hit moves +96.562 μm along floor segment 1272. Its physical σ changes from 1.495478 to 1.497497 mm. Across all 640 donors, the median absolute displacement is 2.65 μm and the maximum is 141.34 μm. Six donors cross to an adjacent discretized segment. Both interpolants place the before/after hits of every donor within the same monotone wall-flux arc when judged on that interpolant’s own wall map.
* An independent scalar check solved ψ(wall) = ψ(launch) on donor 16’s observed arc. It predicts +96.568 μm and each root lies within 0.024 μm of the traced hit. Switching the field interpolant at fixed launch coordinates contributes +97.357 μm to that root shift, averaged over both launch choices. Switching only the launch coordinates contributes −0.789 μm. The same donor on the baseline moves approximately +28.1 μm. The cubic wall-tangential flux slope for this donor has magnitude 0.03198 on the candidate versus 0.06055 on the baseline, in the equilibrium’s flux units per metre. These comparisons establish sensitivity to the interpolated equilibrium at this floor location. They do not select a physically correct interpolant.

## Global peak and area

All relative peak and A₂ changes below use cubic as the denominator. The v3 local error uses the larger peak. A₂ uses the saved cell powers and physical areas. Both candidates fail only the v3 local limit. Both baselines pass all v3 distribution checks.

| χ (m²/s) | Wall | v3 local | v3 L¹ | Cubic → quintic qpeak (MW/m²) | qpeak change | A₂ change | Peak arclength shift |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0.3 | Candidate | 2.6930% | 0.2729% | 1.440223 → 1.443282 | +0.2124% | +0.03155% | +71.268 mm |
| 0.3 | Baseline | 0.7082% | 0.1147% | 2.730737 → 2.729744 | −0.03636% | +0.02702% | +0.02290 mm |
| 0.6 | Candidate | 1.8916% | 0.2707% | 1.449940 → 1.450742 | +0.05530% | +0.03630% | +0.03562 mm |
| 0.6 | Baseline | 0.5044% | 0.1110% | 2.669777 → 2.668786 | −0.03713% | +0.02864% | +0.02114 mm |

On the χ = 0.3 candidate, two near-tied local maxima exchange order. Near s = 4.874 m, the local maximum rises from 1.440223 to 1.441051 MW/m². Near s = 4.945 m, it rises from 1.439571 to 1.443282 MW/m². The cubic global maximum is at (R, Z) = (1.278733, −1.419590) m, and the quintic global maximum is at (1.228175, −1.369951) m. The individual local maxima move +0.063 and +0.097 mm. A₂ rises from 0.934794 to 0.935088 m². At χ = 0.6 the locations of all first hits are the same as in the χ = 0.3 case, and the larger footprint widths give a smaller discrepancy at the same floor edge. The fixed-point contributions there are +0.027510 MW/m² from hit centers, −0.000078 from physical widths and +0.000011 from weights.

## Candidate identity and scope

Both finite-grid OMP searches selected index 0, ID 437d6200…bb37fb. Historical OMP20 stage 18 has ID 0eb363d2…755cc1. The coefficient vectors differ by at most 6.94 × 10⁻¹⁸. At matched χ = 0.3 the maximum physical wall-vertex difference is 2.22 × 10⁻¹⁶ m. For each method, every historical/current first-hit segment matches and the maximum difference in saved q samples is less than 3.8 × 10⁻¹³ MW/m². The measurements here use the finite-grid exports. The two archived identities remain separately identified in the JSON evidence.

The v3 metrics describe reconstructed saved profiles on the same physical support. The counterfactuals identify the inputs responsible for the observed difference at the measured location. Higher donor counts, changes to the transport model and experimental validation are separate questions. No geometry, operating point, normalization or acceptance threshold was changed.

## Replay

From the repository root, with the project's Python dependencies including PyArrow available:

```powershell
python -m verification.nf_local_attribution --study D:/DivertorOptimization/studies/nf-method-20261105 --historical D:/DivertorOptimization/studies/omp-stage-qualification-20261104/R020/validation/0eb363d2b21c1d391948f380244813840964a06edab792fffca167be74755cc1 --output D:/DivertorOptimization/tmp/nf-local-attribution-20261106
python -m pytest verification/test_nf_local_attribution.py -q
```

The JSON includes all factorial corners, native replay residuals, source and archive hashes, physical peak locations, and independent flux-root checks. The CSV contains all 2,560 paired donor rows, including both hits, widths, gradients, weights, connection lengths and signed point/cell contributions.
