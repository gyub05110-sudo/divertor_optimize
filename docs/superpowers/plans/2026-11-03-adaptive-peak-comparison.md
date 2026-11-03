# Adaptive Peak Comparison Implementation Plan

> **For agentic workers:** Use subagent-driven-development for independent review and execute the approved work without another approval gate.

**Goal:** Resolve finalist rankings with less redundant interval work and replay.

**Architecture:** Add an owned incremental builder and adaptive scheduler in `d3d_optimization/adaptive_peak.py`; reuse unchanged selective peak math and replay. Add a read-only study runner with the existing qualification acceptance gate, immutable code snapshot and input manifest checks.

**Tech Stack:** Python, existing directed Decimal bounds, pytest, recorded D: studies.

## Task 1 — Incremental evidence

- [x] Write `verification/test_adaptive_peak.py` tests calling `PeakBuilder(kernel(), 'same')`, advancing one operation at a time, then asserting `replay(builder.certificate())`, monotonic bounds and caller-model isolation. Verify missing API fails.
- [x] Implement `PeakBuilder` with private copied model, `_initial`, lazily constructed `SpatialBounds`, deterministic probe queue, then `_operation` worst-upper splits. Export existing v1 replayable evidence. Stop at coordinate precision and bounded operations without losing cells.
- [x] Run `python -m pytest verification/test_adaptive_peak.py verification/test_selective_peak.py -q` and require all pass.

## Task 2 — Shared adaptive comparison

- [x] Add tests for separated peaks resolving early, identical kernels remaining unresolved at a cap, zero time, invalid identities, and a shared baseline replayed only once per final record.
- [x] Implement `rank(models, problem_ids, pairs, seconds=136, max_operations=4096)`. Validate pairs/identities before numerical work. Refine only unresolved participants; select maximum `gap * unresolved_degree`, with deterministic ties. Save operation scheduling and decision timings. Independently replay each final record once, then return decisions from those verified records.
- [x] Test near ties, tampered comparison evidence, mutable caller inputs, operation budget and final replay failure propagation.

## Task 3 — Real measurement and delivery

- [x] Add `verification/run_adaptive_peak_comparison.py`: verify recorded input hashes and explicit source-scope audit when needed; load exact prior candidate kernels; snapshot source; execute one shared session; preserve all certificates/timings/manifests in a new D: directory. Keep old runner unchanged.
- [ ] Run at cubic_320 then cubic_1280 with 136 seconds generation and the same five comparisons; verify every input and result hash. Include all mandatory replay and preparation in end-to-end timing.
- [ ] Obtain independent specification/correctness review and fix blockers. Run existing 172-test selection plus new tests, and `git diff --check`.
- [ ] Report whether equal decisions were reached faster, unresolved cases, and limits. Commit and push authorized branch updates after verification; no global optimum or new optimized geometry claim.
