# Research archive · 2026-11-05

The [GitHub Release](https://github.com/gyub05110-sudo/divertor_optimize/releases/tag/research-archive-20261105) stores the original files behind the NF method, first-wall routing and FIREFLY studies. Its tag points to [the source snapshot](https://github.com/gyub05110-sudo/divertor_optimize/commit/7a889e86854ae2285c92df262aa565e25c94c0f4). The date refers to the study series. Original reports and result bytes retain their historical values. Refer to the [current NF method report](../nf-method-20261105/REPORT.ko.md) for subsequent reporting corrections.

## Read the FIREFLY results

- [Baseline and midplane-source comparison](firefly/REPORT.ko.md)
- [Baseline and two OMP geometries on common meshes](firefly/common_rank/REPORT.ko.md)
- [Blinded screen of 80 candidate geometries](firefly/blind_screen/REPORT.ko.md) and [sampling protocol](firefly/blind_screen/PROTOCOL.md)

## Download and verify

Download all 13 ZIP assets and `inventory.json` from the release. Each ZIP restores paths relative to a common `DivertorOptimization` directory. The [member inventory](inventory.json) lists each ZIP's size and SHA-256 and every original file's relative path, size and SHA-256. The archive contains 2,139 files totaling 5,216,175,126 bytes before compression.

Use Python 3.10 or later. The first command validates every ZIP and member without writing extracted data. The second verifies and restores each file into a new destination. The destination must not already exist.

```console
python verify_archive.py /path/to/downloaded/assets
python verify_archive.py /path/to/downloaded/assets --restore /path/to/new/DivertorOptimization
```

The restoration tool reads the inventory beside itself by default. For downloads without a repository checkout, place the release's `inventory.json` beside the release copy of `verify_archive.py`. The archives require approximately 6 GB of free space after download. The archived scripts refer to `D:/DivertorOptimization` and its `worktrees/d3d-report-correction-20261104` checkout. Use that layout or update the path arguments in a working copy. Restore into an empty location first if studies already exist there.

| Asset | Files | ZIP bytes |
| --- | ---: | ---: |
| `analysis-and-donor-inputs.zip` | 39 | 40,978,899 |
| `firefly-baseline.zip` | 297 | 253,120,881 |
| `firefly-blind-screen.zip` | 945 | 232,416,475 |
| `firefly-common-rank.zip` | 293 | 363,159,131 |
| `nf-method.zip` | 316 | 84,118,497 |
| `preimpact-firstwall.zip` | 106 | 11,746,819 |
| `screen-campaign-inputs.zip` | 137 | 14,858,493 |
| `screen-ledger-imp-r020.zip` | 1 | 69,463,315 |
| `screen-ledger-imp-r050.zip` | 1 | 179,154,679 |
| `screen-ledger-imp-r100.zip` | 1 | 789,426,863 |
| `screen-ledger-omp-r020.zip` | 1 | 28,320,087 |
| `screen-ledger-omp-r050.zip` | 1 | 37,011,684 |
| `screen-ledger-omp-r100.zip` | 1 | 69,730,351 |

The three FIREFLY archives contain native NetCDF meshes and heat results, original configurations and inputs, analysis code, diagnostic run records, tables and reports. The baseline ZIP contains the exact locally patched shared library and Python extension from its protected checksum list. The analysis archive includes the three audited OMP input arrays with their original paths, two baseline donor models, the two historical IMP donor models, and the prior candidate-validation artifacts inspected by the submission audit. The NF and preimpact archives retain original reports, numerical arrays, routing frames and campaign evidence. The campaign archive contains the scientific input and source snapshot for score replay. Separate ledger archives preserve the six source/radius histories, including the IMP 100 mm stratum with zero heat evaluations.

## Software provenance

Use the repository requirements for Python packages. The frozen campaign inventories record dependency file hashes. Package caches, third-party installs, registry layers, operating-system images, Git metadata and temporary test directories are excluded. The original directories used a generic local account and Windows/WSL paths.

The [FIREFLY build instructions](firefly/REPRODUCE.md) identify the source revision, image digest, adapter patch and f2py compatibility change. Recreate external components before rerunning native solvers. The original instructions describe a workstation on which the Linux image and upstream checkouts already existed. For historical reduced-model code comparisons in `analyses/transport-independence-20261105/audit.py`, the source revision recorded in the original studies is `dc8ab43a9a930ca155813cb53a448bba39f6291f`.

The FIREFLY checksum manifests cover 161 baseline files, 187 common-rank files and 682 blind-screen files. One protected equilibrium parquet and five first-wall verification scripts are supplied by the tagged source checkout. Their committed bytes match the hashes recorded in the original studies.

| Component | Original reference |
| --- | --- |
| FIREFLY upstream | [85a1e6c006e8af0c6a6cb6775ff78efee61af47e](https://gitlab.com/hfrerichs/firefly/-/tree/85a1e6c006e8af0c6a6cb6775ff78efee61af47e), then apply `studies/firefly-baseline-20261105/firefly-midplane.patch` and use LF line endings for `autowrap` |
| FLARE runtime | `hfrerichs/flare:v1.4.7`, registry manifest `sha256:cc26d9250125d515eabc6d1a43eb5de98de11f1463770f36f66c3d2f5432b41c` |
| FLARE source checkout used for reference | [80eec5423fb74ae74870d04c919dd83ca9c07a42](https://gitlab.com/hfrerichs/flare/-/tree/80eec5423fb74ae74870d04c919dd83ca9c07a42) |
| Public example-shot comparison repository | [823f6eef1f527e1b5566a8fc1e952da6e2fadf43](https://github.com/Sophelio/fusion-equilibrium-challenge-starter/tree/823f6eef1f527e1b5566a8fc1e952da6e2fadf43). Clone its history into `common_rank/public-starter` to rerun the Git LFS pointer check. The five example parquet files and their recorded fingerprints are in the common-rank archive. |
