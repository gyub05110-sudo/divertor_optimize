---
pretty_name: Fusion Equilibrium Challenge
license: cc-by-4.0
language:
  - en
tags:
  - fusion
  - plasma-physics
  - tokamak
  - scientific-ml
  - equilibrium-reconstruction
  - image-regression
size_categories:
  - 1K<n<10K
configs:
  - config_name: diii_d_train
    data_files:
      - split: train
        path: data/diii_d_train/*.parquet
  - config_name: diii_d_public_test
    data_files:
      - split: public_test
        path: data/diii_d_public_test/*.parquet
  - config_name: mast_public_test
    data_files:
      - split: public_test
        path: data/mast_public_test/*.parquet
---

# ⚛️ The Fusion Equilibrium Challenge

Predict the **shape of a fusion plasma** (the magnetic equilibrium, $\psi$) from
control inputs and diagnostics alone — *without* magnetic sensors. Data comes from
two tokamaks:

- **DIII-D** — General Atomics tokamak (San Diego, USA)
- **MAST** — Mega Ampere Spherical Tokamak (Culham, UK)

In data-science terms this is an **image-regression / control problem**: predict a
2-D poloidal flux map (`efit_psirz`) at each EFIT timestep from coil currents
(the *actuators*), Thomson-scattering temperature/density profiles (the *sensors*),
and the machine's fixed **coil and diagnostic geometry**.

> One row per shot. All time-series, profiles, and flux maps are stored as nested
> arrays within that row. **Every time array is in milliseconds** (both machines).

## ⚠️ Diverted plasmas only

**Every timeslice in this dataset is a diverted plasma.** Limited frames — where the plasma's
boundary is set by physical contact with a material surface rather than by a magnetic X-point —
have been removed from `efit_times` and from every column that lives on it, on **both machines
and all splits**.

This makes the problem well-posed. On a limited plasma the last closed flux surface is determined
by the wall, and wall geometry is not contained in ψ(R,Z) — so a boundary extracted from the flux
map alone systematically overshoots (91% of DIII-D limited frames overshot plasma volume by more
than 0.5 m³). Scoring a boundary you cannot infer is not a useful test, so those frames are gone
rather than silently down-weighted.

How the frames were classified — the same physical criterion, each machine's best evidence:

| | criterion | kept |
|---|---|---|
| **DIII-D** | `magnetics_dsep > 0`. That column is EFIT's a-file `DSEP`, a separatrix↔limiter clearance whose **sign encodes the configuration** (`>0` diverted gap, `<0` limited). 0% NaN over the corpus. | 88.0% of frames |
| **MAST** | EFIT++ reports X-point coordinates (`efm/xpoint{1,2}`) only on frames where it found one. Cross-checked against the flux itself: `\|ψ(X-point) − ψ_boundary\| / \|ψ_boundary − ψ_axis\|` agrees on **97.3%** of 184,396 frames, median offset 0.0006. | 66.4% of frames |

The full-rate input time series (`magnetics_*`, `thomson_*`) are **not** subset — you
still get the complete discharge and resample onto `efit_times` as before. `efit_times` is simply
no longer contiguous.

> **A caution on MAST's `magnetics_dsep`.** It is *not* the same quantity as DIII-D's. MAST ships
> `esm/dr_sep_out`, the upper/lower divertor **balance** (δR_sep), which straddles zero for
> perfectly ordinary diverted plasmas. **Its sign is not a diverted/limited indicator** — do not
> use it as one. It remains in `train` as context only.

## 🏆 The competition

This dataset backs the **Fusion Equilibrium Challenge**, a NeurIPS 2026 competition. You do not
have to enter to use the data — but if you want to:

- **Leaderboard & submission:** [Codabench competition #17456](https://www.codabench.org/competitions/17456/)
- **Starter kit** (baselines, submission skeleton, local validator): [Sophelio/fusion-equilibrium-challenge-starter](https://github.com/Sophelio/fusion-equilibrium-challenge-starter)
- **Website** (full rules, scoring, glossary): <https://fusion-equilibrium-challenge.sophelio.io>

There are **two challenges**, and each grants exactly one award — the mapping is one to one. Both
share this dataset, one submission format, and one leaderboard (the challenges are two sortable
columns on it).

| | Challenge 1 — Intra-machine | Challenge 2 — Cross-machine |
|---|---|---|
| **What you do** | Reconstruct the equilibrium on DIII-D, the machine you trained on. | Take that model zero-shot to MAST — which has **no training data at all**. |
| **Scored by** | `S_model` on the hidden DIII-D test set | `G_ratio = S_MAST / S_DIII-D` |
| **Award** | **Award #1** — $500 | **Award #2** — $500 |

`S_model = 0.55·R²ψ + 0.15·R²{q95,βN} + 0.10·(1 − D_LCFS) + 0.20·Consistency`, where
`Consistency` scores the seven ψ-derived scalars — magnetic axis `R_axis`/`Z_axis`, elongation
`κ`, triangularities `δ_top`/`δ_bot`, plasma `volume`, and internal inductance `li` — computed
**from your submitted flux map** against the same computation on the true flux. You submit only
`psirz`, `q95`, and `betaN`. Challenge 2 needs predictions for **both** machines: a DIII-D-only
entry scores `G_ratio = 0`. A **compute-light badge** (trains in under 2 h on one GPU/CPU)
applies across both challenges.

| Phase | Window |
|---|---|
| **Development** — public leaderboard, iterate freely | **Jul 27 – Oct 18, 2026** |
| **Final** — blind private fold, decides the winners | **Oct 19 – Oct 26, 2026** |

Note that `mast_train` does not exist, and that is deliberate: **MAST has zero training shots.**
That is the whole point of Challenge 2.

## Splits & configs

Each (machine, role) is its **own config**, so every config has a clean schema —
the two machines have different coil layouts, and the test configs genuinely omit
the target column rather than null-filling it:

| Config               | Split         | Shots | Contents |
|----------------------|---------------|-------|----------|
| `diii_d_train`       | `train`       | 7,041 | **Full** — includes all targets (`efit_psirz`, EFIT scalars) |
| `diii_d_public_test` | `public_test` | 874   | **Inputs only** — all EFIT-derived columns withheld |
| `mast_public_test`   | `public_test` | 1,206 | **Inputs only** — all EFIT-derived columns withheld |

Shot counts are lower than the previous release because the diverted-only re-cut removed shots
left with no diverted frames at all: DIII-D 9,113 → 8,794 shots and 2,216,380 → 1,949,000
timeslices; MAST 2,416 → 2,414 shots and 207,066 → 137,563 timeslices.

```python
from datasets import load_dataset

REPO = "Sophelio/fusion-equilibrium-challenge"

# DIII-D training data (with targets)
train = load_dataset(REPO, "diii_d_train", split="train")

# Public test (inputs only — predict efit_psirz at each efit_times step)
d3d_test  = load_dataset(REPO, "diii_d_public_test", split="public_test")
mast_test = load_dataset(REPO, "mast_public_test",  split="public_test")

import numpy as np
shot = train[0]
psirz = np.array(shot["efit_psirz"])   # (T, 65, 65) for DIII-D
times = np.array(shot["efit_times"])   # (T,) ms — predict a flux map at each
```

**Notes on the design:**
- **No MAST training data.** Cross-machine generalization (DIII-D → MAST) is
  zero-shot by construction — the `mast` config has only a test split.
- **All EFIT-derived quantities are withheld on the test splits.** Anything computed by
  the EFIT reconstruction leaks the equilibrium you must predict, so every one is removed
  from `public_test`: the flux map `efit_psirz`, the EFIT scalars (`efit_beta_n`,
  `efit_li`, `efit_q95`, `efit_r_axis`, `efit_z_axis`, plus the auxiliary/boundary
  scalars), and `magnetics_dsep` (EFIT-derived despite its `magnetics_` prefix — and on DIII-D
  its sign is exactly the classification that defined this corpus). Only the coordinate
  scaffolding is kept — `efit_times` and
  `efit_grid_R/Z` (both machines) — so you know exactly which timestamps and grid to predict on.
- A **private test split** is held back entirely for final scoring and is not part of
  this dataset.

## The targets — what you predict

The **primary target** is the 2-D flux map. You *submit* only the flux map plus `q95` and `betaN` —
the two scalars a flux map cannot contain (they need the toroidal field function `F(ψ)` and the
pressure profile `p(ψ)`, neither of which is in the flux geometry). The dataset still ships all
**EFIT scalar labels** (one value per `efit_times` step) in `train` as supervision — withheld on
the test splits, same as the flux map — but at scoring time the axis, shape and `li` are
**derived from your submitted flux map** and scored for consistency against the same derivation on
the true flux. A good ψ is what earns them.

| Key | Shape | Description |
|-----|-------|-------------|
| `efit_psirz`  | `(T, 65, 65)` (both machines) | **Primary target** — poloidal flux map, a 2-D image at each timestep (V·s/rad). *Withheld on test.* |
| `efit_q95`    | `(T,)` | Safety factor at 95% flux surface — **submitted scalar**. *Withheld on test.* |
| `efit_beta_n` | `(T,)` | Normalised beta — **submitted scalar**. *Withheld on test.* |
| `efit_li`     | `(T,)` | Internal inductance ℓi (supervision; derived from your ψ at scoring). *Withheld on test.* |
| `efit_r_axis` / `efit_z_axis` | `(T,)` | Magnetic-axis R/Z (m) (supervision; derived from your ψ at scoring). *Withheld on test.* |
| `magnetics_dsep` | `(T,)` | Context only, **not scored**. On DIII-D it is EFIT's a-file `DSEP` (separatrix↔limiter clearance, now `>0` on every shipped frame — it is what defined the diverted filter). On MAST it is a different quantity entirely: `esm/dr_sep_out`, divertor balance. *Withheld on test.* |
| *(boundary — both machines)* `efit_lcfs_n`, `efit_lcfs_r`, `efit_lcfs_z` | `(T,)` / `(T, N)` | LCFS boundary contour + valid-point count — provided in `train` as bonus context. *Withheld on test.* |
| *(aux — MAST only)* `efit_beta_p`, `efit_elongation`, `efit_psi_axis`, `efit_psi_boundary`, `efit_geom_axis_r/z` | `(T,)` | Extra EFIT scalars, MAST only (the DIII-D fetch did not include them). ⚠️ **Not in any released config** — they are EFIT-derived so they are withheld on test, and there is no `mast_train`. You will only see them in the starter kit's `parquet_data/` demo shots. |

> **`dsep` is no longer scored at all.** It was one of eight consistency scalars in the previous
> metric and has been dropped, for a reason worth stating plainly: the two machines' `dsep`
> columns are **not the same physical quantity**. DIII-D's is a clearance whose sign encodes the
> configuration; MAST's is divertor balance (δR_sep). Scoring one functional against both labels
> was measuring different things per machine. With the corpus now diverted-only, the
> classification job it also did is gone. Seven scalars remain in `Consistency`.

**Scaffolding (kept on every split):**

| Key | Shape | Description |
|-----|-------|-------------|
| `efit_times`  | `(T,)` | Timestamps (ms) for the target images. Align all inputs to these times. |
| `efit_grid_R` / `efit_grid_Z` | `(65,)` | Physical R/Z (m) labelling the flux-map columns/rows (both machines). |

Both machines ship a clean `65×65` flux map with **no NaNs**. (MAST's upstream EFIT
stores `psirz` on a doubled 129-column R grid — 65 real columns interleaved with 64
empty ones; this is *not* a central-column hole but a grid artifact. We drop the
empty columns to recover the dense 65×65 grid the data is actually defined on:
MAST R ∈ [0.06, 2.0] m, Z ∈ [−2.0, 2.0] m.)

> **The two machines' flux maps use opposite sign conventions.** On MAST the magnetic axis is the
> **maximum** of ψ (`efit_psi_axis > efit_psi_boundary` on 100% of frames); on DIII-D it is the
> **minimum** (99.98% of 1,559,340 train frames, with no frame the other way). Amplitudes differ
> too — a DIII-D flux map has a median frame peak-to-peak of 0.61 Wb/rad against MAST's 0.21, so
> DIII-D ψ swings roughly **3×** wider. This is provenance, not physics: EFIT++/FAIR-MAST stores ψ
> with the opposite sign to DIII-D's EFIT, and both machines run positive plasma current on
> essentially every shot. **The scorer normalizes the global sign of your submitted flux map**, so
> a correctly-shaped prediction in either convention scores the same. The ~3× amplitude difference
> is real and is yours to handle.

## Inputs

**DIII-D actuators** — 18 shaping coils `magnetics_F{1-9}{A,B}` (kA·turn; ~77% of in-window peaks
stay within ±200, the extremes reach ~±580),
`magnetics_ECOILA` (ohmic/central solenoid, kA), `magnetics_bcoil` (toroidal field, kA),
`magnetics_plasma_current` (Ip, kA).
**MAST actuators** — 10 poloidal coils `magnetics_p{2-6}{l,u}_current` (kA·turn; the set is far
more spread than DIII-D's — per-coil median in-window peaks run from ~3 (`p3`, `p6`) to ~185
(`p2`), extremes ~±600), `magnetics_sol_current` (central solenoid, kA),
`magnetics_tf_current` (toroidal field, kA), `magnetics_efps_current` (error-field correction, kA),
`magnetics_plasma_current` (Ip, kA).

**Both machines** — **Thomson scattering** core + edge profiles
(`thomson_core_*`, `thomson_edge_*`: electron temperature `Te` in eV and density
`ne` in m⁻³, with time bases). Each system carries a single spatial-coordinate array,
and the axis differs by machine: `thomson_core_R` is radius (DIII-D constant ≈ 1.94 m
vertical chord; MAST per-channel R), and `thomson_edge_spatial` is **Z** on DIII-D
(−0.09 to −0.04 m across the corpus, near-constant within a shot) but **R** on MAST
(≈ 1.3–1.5 m). There is no separate `*_z`/`tan_z` column.

> **Note:** `magnetics_dsep` used to be listed here as an input. It is EFIT-derived, so it is
> present in `train` as context and withheld on test.

Time bases: each machine exposes a shared `magnetics_time`; on DIII-D, Ip sits on
its own ADC (`magnetics_plasma_current_times`). In `train`, `magnetics_dsep_times`
equals `efit_times` on every shot (it and `magnetics_dsep` are withheld on test).

> **DIII-D's `magnetics_time` has two sampling rates.** 4,932 of 7,041 train shots (70.0%) carry a
> ~20 kHz base — 480,256 samples at dt = 0.049927 ms — and the remaining 2,109 (30.0%) a 2 kHz base
> of 49,152 or 50,176 samples at dt = 0.5 ms. The public-test mix is the same within a few points
> (586 shots, 67.0%, at ~20 kHz). Six distinct array lengths occur in all, so do **not** hard-code a
> sample count or a rate — read `magnetics_time` per shot. The record spans roughly −9.4 s to
> +20.1 s (median duration 24.0 s), so only about 20% of its samples fall inside the `efit_times`
> window.

> **MAST's `magnetics_time` has two populations.** 1,092 of 1,206 test shots carry `(30000,)`, a
> uniform 0.2 ms / 5 kHz grid. The 114 early-campaign shots carry `(15482,)`, which is the
> **union of two acquisition grids** — the poloidal set at 5 kHz over a shorter window, the
> toroidal field coil at 1 kHz over the full record — so on those shots each column is `NaN` on
> the samples belonging to the other grid. No data is missing; both native grids fully cover the
> plasma window on every shot. This is the layout FAIR-MAST's own `amc` group ships.

### Missing values: always `NaN`, never a parquet null

Every floating-point column in this release is **null-free at every nesting depth**, and there
are no `±inf` values anywhere — wherever a value is absent, it is a plain `NaN`. So every column
converts cleanly to a dense float array, on both machines, with no `object`-dtype surprises:

```python
v = np.asarray(shot["magnetics_tf_current"], dtype=float)
```

Which columns carry `NaN`, how much, and what to do about it are yours to find out and decide —
that is part of the problem, not a defect in the data.

### 🆕 Machine geometry — where the coils and chords actually are

A coil current only tells you something once you know **where the coil is**, and a Thomson
measurement only localises a profile once you know **where the chord is**. Both are now shipped,
on every row of every split (they are inputs, so nothing is withheld):

| Key | Shape | Description |
|-----|-------|-------------|
| `coil_name` | `(C,)` str | Coil (or conductor element) identifier |
| `coil_input_column` | `(C,)` str | **The join key** — the current column this geometry belongs to, e.g. `magnetics_F1A` |
| `coil_R` / `coil_Z` | `(C,)` | Centre of the conductor rectangle (m) |
| `coil_width` / `coil_height` | `(C,)` | Its radial / vertical extent (m) |
| `coil_angle1` / `coil_angle2` | `(C,)` | Parallelogram **skew angles** (degrees) — see below |
| `thomson_chord_name` | `(N,)` str | Chord identifier, IMAS-style (`TS_core_*`, `TS_tangential_*`, `TS_divertor_*`) |
| `thomson_chord_R` / `thomson_chord_Z` | `(N,)` | Chord position (m) |

**The two machines describe their coil sets at different granularity, and that is not an
inconsistency.** DIII-D ships **C = 19 lumped rectangles** (18 F-coils + the `ECOILA` solenoid) —
the representation EFIT's own `mhdin.dat` uses. MAST ships **C = 812 individual conductor
elements** from the FAIR level-2 `pf_active` IDS — one row per conductor turn, the solenoid alone
accounting for 656 — so a MAST coil's turn count is how many rows share its `coil_input_column`
(P2 = 20, P3 = 8, P4 = 23, P5 = 23, P6 = 4). Use `coil_input_column` to join either one back to
the currents; do not read the row counts as a discrepancy.

**There is no `coil_turns` column, and you do not need one.** Coil currents are already in
**ampere-turns**: each DIII-D F-coil's 58 or 55 turns is folded into `magnetics_F*`, and MAST's
`p{2-6}{l,u}_current` are natively `kA·turn`. That is the quantity a Green's-function calculation
wants — total ampere-turns per rectangle — so `F4A` and MAST `p4u` are directly comparable: median
in-window peak 137 kA·turn and 171 kA·turn respectively, over every shot in each config. Plasma
current is in **kA** on both machines. Four channels stay in plain
**kA** rather than kA·turn because no trustworthy turn count exists for them: DIII-D `ECOILA` and
`bcoil`, MAST `sol_current` and `tf_current` (and `efps_current`, which is not a coil).

**Six DIII-D F-coils are parallelograms, not rectangles.** `coil_angle1` / `coil_angle2` are
EFIT's `AF` / `AF2` shear angles in degrees, and they are non-zero on `F5A`/`F5B` (±45°),
`F6A`/`F6B` (±92.4°) and `F7A`/`F7B` (±108.06°). For those six, `(R, Z, width, height)` alone
under-describes the conductor.

> **`0.0` means "no skew — a plain axis-aligned rectangle."** It does *not* mean the sides lie
> flat. This is EFIT's own convention: its geometry code branches on
> `if angle1 == 0 and angle2 == 0` to emit a rectangle rather than a polygon, and on write-out it
> normalises `90 → 0`, i.e. **90° and 0° denote the same unskewed coil**. We ship the canonical
> `0`. So the 13 DIII-D coils at `0.0` and all 812 MAST elements at `0.0` are making the identical
> statement — *this conductor is a rectangle* — and you can treat the two machines' geometry with
> one code path. MAST's zeros are not a missing-value sentinel: its elements come from IMAS
> `element[:].geometry.rectangle`, which has no skew degree of freedom, so there is no angle
> upstream that could have been lost.

To build the conductor outline: if both angles are 0, use the rectangle
`(R ± width/2, Z ± height/2)`. Otherwise shear that rectangle — `angle1` tilts one pair of sides,
`angle2` the other — which is what EFIT does when it converts these six coils to polygons.

*Not covered:* `magnetics_bcoil` (DIII-D toroidal field) and MAST's `magnetics_tf_current` /
`magnetics_efps_current` have no poloidal-plane rectangle, so they appear in no geometry row —
19 of 21 DIII-D current columns and 11 of 14 MAST ones are covered.

*One caveat on `ECOILA`:* its row is a single lumped rectangle for the ohmic solenoid, and its
turn convention is genuinely ambiguous — EFIT's machine file models the same group as 48
single-turn elements spanning R 0.6875–0.7625 m, Z ±1.7255 m (an envelope the shipped rectangle
reproduces to ~1.6 cm), and `ECOILB` is a second co-located group not shipped here. That is why
`magnetics_ECOILA` is left in **kA** rather than kA·turn: do not apply an ampere-turn multiplier
to it. The 18 F-coil rows carry no such ambiguity (see below).

**Chord positions are per shot on DIII-D and they genuinely vary** — 22 distinct channel-name
layouts in train (19 in public test) across 6 distinct subsystem
layouts and channel counts from 59 to 138 across the corpus, mostly because the divertor system
was reconfigured between campaigns. Do not cache one shot's chord array and reuse it.

This closes a real gap in the previous release. DIII-D's core Thomson is a **vertical** laser, so
`thomson_core_R` is a constant ≈ 1.94 m and the informative coordinate — Z — was never shipped;
the tangential system is the mirror image (`thomson_edge_spatial` is its nearly-constant Z, and
its informative R was missing); and the **divertor** subsystem had no shipped counterpart at all.
`thomson_chord_R`/`_Z` carry all three subsystems with both coordinates.

MAST's Thomson is a **horizontal midplane** laser: per-channel R (already in `thomson_core_R` and
`thomson_edge_spatial`) and `thomson_chord_Z = 0` by construction. There is no Z upstream to
ship — FAIR level-2 exposes only `thomson_scattering.channel[:].position.r`.

```python
import numpy as np
shot = train[0]
# vacuum field contribution of each coil, using the shipped geometry
geom = {c: np.array(shot[c]) for c in ("coil_R", "coil_Z")}
for name, col in zip(shot["coil_name"], shot["coil_input_column"]):
    # already ampere-turns for the F-coils / MAST P-coils -- no turn multiplier needed
    amp_turns = np.array(shot[col])        # this coil's kA*turn time series
```

👉 **Full signal dictionary, machine differences, and a beginner-friendly guide** are in the
[starter kit](https://github.com/Sophelio/fusion-equilibrium-challenge-starter) — its `README.md`,
`MODELING_GUIDE.md`, and the `fusion_data_provider.py` visualizer plugin.

## License & acknowledgement

**Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)** — use it for the
competition or for any research, with attribution. The private (blind) test fold and its ground
truth are never released.

Work supported by the U.S. Department of Energy, Office of Science, Office of Fusion
Energy Sciences, using the DIII-D National Fusion Facility, a DOE Office of Science
user facility, under Award No. DE-FC02-04ER54698, and Awards DE-SC0024426,
DE-SC0024499, DE-SC0024409, and DE-SC0024571. See the full challenge README for the
complete DOE disclaimer.

## Citation

```bibtex
@dataset{fusion_equilibrium_challenge,
  title  = {The Fusion Equilibrium Challenge: Predicting Plasma Shape from Control Inputs},
  author = {Michoski, Craig and Waller, Matthew and Sammuli, Brian and Boyes, William
            and Clark, Mitchell and Smith, Sterling and Nakkina, Tapan Ganatma
            and Hatch, David and Nazikian, Raffi},
  year   = {2026},
  note   = {Sophelio and General Atomics}
}
```
