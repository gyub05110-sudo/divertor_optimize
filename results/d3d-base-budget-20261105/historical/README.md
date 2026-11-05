# Historical comparisons for the DIII-D base-budget trial

Run the read-only audit from the repository root, replacing `TRIAL` with the directory containing the new run's `protocol.json`, `OMP/`, `IMP/` and `audit/`. Python 3.13 and its standard library suffice.

```text
python -m verification.nf_base_budget_audit --study TRIAL --studies results/d3d-base-budget-20261105/historical
```

`catalog.json` records the source-file SHA-256 hashes, category memberships, the earlier cubic-1280 baselines and their operating inputs. The audit rebuilds the old 56 grid vectors from the repository anchors and reads published winners from their versioned manifests. The catalog adds the remaining 51 vectors. Identifiers include the source, radius and, for IMP, junction fraction. The reported distances use only spline displacement coefficients, in metres.

The original archive and a clone assign different filesystem timestamps. The source-run audit recorded when the six panels and earliest thermal pair were saved and bound those times to file hashes. A clone uses that record to verify the original order. Per-pair hashes in the audit bind the cost replay to the saved evaluations.

To reproduce the compact reference from the full study archive, choose a new output directory.

```text
python -m verification.nf_base_budget_audit --studies ARCHIVED_STUDIES --export-history NEW_DIRECTORY
```
