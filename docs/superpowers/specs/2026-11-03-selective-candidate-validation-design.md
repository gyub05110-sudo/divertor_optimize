# Selective candidate validation

The user approved the design in this chat: keep fast search, spend additional validation only on finalists, measure total time to a reliable improvement, and defer global optimality certification.

## Scope and acceptance

1. Revalidate baseline and saved IMP20/50, campaign OMP50 and historical OMP50 candidates with the current common geometry policy and the existing four-resolution independent audit. Preserve original studies and report exclusions.
2. Diagnose the 16-donor continuous model at a point before changing geometry domains. Record per-donor arrival, L/J, mapping/width bounds, analytical versus spatial peak bounds, preparation/replay/evaluation times.
3. Add a bounded spatial peak checker that uses donor centers, supports broad intervals without dropping power, and returns an honest unresolved result on budget exhaustion. Preserve existing certificate replay.
4. Add selective path refinement chosen from contributions to peak uncertainty, retaining original source and all hit branches. Integrate the resulting kernel with the spatial checker. A finite budget must stop even when intervals still overlap.
5. Run measured studies and relevant regressions. Establish numerical resolution qualification at 320/640/1280 and quintic for real candidates; continuous certification remains explicitly a finite-donor study unless measured scalability supports expansion.

## Architecture

The production search and old proof APIs remain unchanged. New modules provide bounded comparison and refinement; new verification scripts run read-only input studies under a new D: directory. Existing path caches are replayed before reuse. All old certificates/results remain intact.

Peak brackets from different source resolutions or model definitions must never be intersected. Frozen-kernel certificates are not continuous-transport certificates. Certified candidate ranking concerns the declared transport model, not geometry feasibility or global optimality.

## Budget policy

Check wall-clock budgets at bounded units of work; report any single-operation overrun and the smallest noninterruptible unit. Keep best valid evidence on exhaustion. A result may be improved, worse, indistinguishable at requested tolerance, or unresolved. No guessed lower bound, power renormalization, width clipping, or dropping ambiguous branches.

## Validation

Analytic narrow/broad periodic kernels, separate versus colocated heat sources, cache/input identity, near ties, incomplete hits, retained branches, deterministic replay, and real DIII-D candidate audits. Compare equal accuracy runtime and bound widths, not only one-call timing. New code must not add calls to production search loops.
