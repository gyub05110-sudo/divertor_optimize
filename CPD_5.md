# CPD_5: fixed FreeGS equilibrium and three-hour divertor search

The supplied device geometry replaces D3D geometry. Existing D3D studies and
source files are preserved. Raw YAML, a losslessly zipped STEP and their hashes
are in [input/CPD_5](input/CPD_5). The node contours define the axisymmetric
computational section; the STEP remains the original CAD reference.

## Fixed operating point

The user approved deriving CPD PF/OH currents in FreeGS while preserving the
available D3D operating inputs, with a simple fitted profile and a
geometry-compatible lower-single-null plasma boundary.

| Quantity | Definition / value |
|---|---|
| D3D reference time | 1400 ms, exact EFIT frame 60 |
| Plasma current | 1.388124068 MA |
| Toroidal field | 2.077512221 T: D3D vacuum field at D3D magnetic-axis R, applied at CPD R0=4.007 m |
| Normalized beta target | 0.616604388; documented vacuum-reference convention |
| FreeGS result betaN | 0.616585335 (0.00309% relative residual) |
| Axis pressure | 46,500.29 Pa, simple alpha_m=1 / alpha_n=2 profile |
| Lower X-point | R=3.550005 m, Z=-1.799997 m |
| Frozen field | 129×129 psi grid; SHA256 `70b1c6caf974a4628afcf0cd6584ac3ac127c7bfd7c67430670199e2c5b3501d` |

Full coil **ampere-turns**, shape controls, source-value extraction and numerical
checks are saved in [input/CPD_5/equilibrium](input/CPD_5/equilibrium). Electrical
winding counts and engineering current-density limits were not supplied.
The exact locally modified FreeGS sources and their license are preserved in
`input/CPD_5/freegs-source.zip`; this is more precise than a package version alone.

The 65→129 comparison gives a 0.338 mm LCFS boundary difference and a 0.3175%
95th-percentile relative poloidal-field difference near the separatrix. These
are numerical checks, not experimental or engineering qualification. The EFIT
betaN convention is not fully available and the adopted convention is explicit
in the equilibrium report. The source thermal closure and equilibrium pressure
profile are distinct prescribed reduced models, as in the original study.

## Geometry and heat transport

Blanket and divertor are independent solids. Actual gaps are preserved, including
19.43 mm blanket/divertor and 4.76 mm minimum divertor/vessel clearance. The
axis-connected vacuum has three boundary loops. First hits compete across all
loops. The existing conservative area-measure Gaussian footprint is applied
separately on each loop with unchanged donor weights and a global power ledger.

This is a phenomenological **surface-footprint model**, not a perpendicular
plasma transport solver. It permits redistribution around a component onto
rear surfaces but not between disconnected boundary loops. Front/rear deposited
power is reported explicitly. No missing-power renormalization or width cap is
applied. Toroidal winding is retained in the connection and conduction integrals.

IMP and OMP are independent studies, each at 1 MW, not a combined 2 MW loading:
density 1e19 m^-3, chi_perp 0.3 m²/s, lambda_q 3 mm, gamma 7 and kappa0 2000.

Ten original divertor faces (8→9 through 17→18, including the dome) each have
four cubic-spline normal-displacement coefficients. Original node positions and
face endpoint tangents remain fixed. Other faces, blanket, shield and vessel
are fixed. Bounds are ±100 mm per coefficient, per-face length ratio ≤1.2 and
smooth-face curvature radius ≥20 mm. Existing sharp node corners are preserved
and exempt; this is not a globally smooth or globally optimized shape family.
Invalid polygons, collisions, vessel escape and LCFS intrusion are rejected.

## Baseline and execution

At 320 donors and 1 mm path/receiver/arc settings, the untouched baseline gives
IMP **3.45294 MW/m²** and OMP **2.25648 MW/m²**. Relative changes from the faster
160-donor settings are 0.0111% and 0.0143%. Both resolve all prescribed power.
Further refinement to 640 donors, 0.5 mm paths/receiver/arcs and quadrature 256
gives IMP 3.45310 and OMP 2.25657 MW/m² (0.00450% and 0.00355% changes from
320 donors). These are reference values, not optimization improvements. Full evidence is in
[results/cpd5-20261103](results/cpd5-20261103).

```powershell
python run_cpd_baseline.py --equilibrium input/CPD_5/equilibrium/equilibrium.npz --output D:/DivertorOptimization/studies/cpd5-baseline --count 320 --path-step .001 --wall-step .001 --arc-step .001 --quadrature 128
python run_cpd_campaign.py --equilibrium input/CPD_5/equilibrium/equilibrium.npz --output D:/DivertorOptimization/studies/cpd5-campaign --seconds 10800 --reserve 900
```

The three-hour clock starts after baseline preparation and includes final
candidate re-evaluation. Checkpoints are written during search. Candidate ranking
uses global peak at fixed source power, with footprint and unresolved-power
guards. Final baseline/candidate comparisons use the same frozen field and
refined numerical settings. An unqualified candidate is not reported as a
validated improvement. No global optimum certification is claimed.

To reproduce the solver exactly, extract `freegs-source.zip` and put its parent
directory on `PYTHONPATH`. The Python numerical dependencies are documented in
`requirements-repro.txt`; CPD additionally uses PyYAML and the archived FreeGS.
