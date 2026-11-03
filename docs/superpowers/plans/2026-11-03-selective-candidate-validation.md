# Selective Candidate Validation Implementation Plan

> **For agentic workers:** Use superpowers:subagent-driven-development or superpowers:executing-plans to execute and review each component. The user has authorized execution; no additional design approval is required.

**Goal:** Preserve fast search while making finalist comparisons reliable within explicit compute budgets.

**Architecture:** Add bounded kernel and continuous-refinement modules without changing historical replay or production search. Run independent saved-candidate requalification in a new D: study. Protect all consumed historical files with hashes.

**Tech Stack:** Python, existing Decimal interval arithmetic, NumPy/SciPy, pytest; CPU validation and existing optional GPU evaluator.

## Task 1 — Saved-candidate requalification
- [x] Add verification/run_candidate_requalification.py, using campaign_geometry.evaluate_candidate and campaign_validation.Qualifier or the same audit/verdict pipeline.
- [x] Baseline and candidates share originals per source/resolution; historical OMP50 is freshly checked. Reject invalid candidates before thermal work.
- [x] Record each preparation/evaluation/audit time, input identity and reference, geometry policy, numerical qualification, full power ledger, source hashes. Save new candidate arrays and frozen kernels for later comparisons.
- [x] Test candidate discovery/identity and preservation. Run full real qualification in a fresh D: directory; do not overwrite/promote old results.

## Task 2 — Bounded spatial peak and comparison
- [x] Write analytic behavior tests in verification/test_selective_peak.py. Verify they fail before implementation.
- [x] Add d3d_optimization/selective_peak.py. Use safe image truncation or broad-kernel uniform/Fourier bounds; combine independent donor bounds with spatial subdivisions. Carry analytic universal lower/upper bounds and retain all spatial cells.
- [x] Budget by elapsed seconds and splits. Save actual completed split count for deterministic replay. Compare brackets only for matched problem identities and distinct candidate identities.
- [x] Test known Gaussian peak, separated vs coincident sources, huge sigma range, zero time, near tie, tampering, and input mismatch. Run original peak tests too.

## Task 3 — Continuous diagnostics/refinement
- [x] Read prior local/wide bounds, then record new point bounds and per-donor uncertainty from saved replayed paths. Do not treat nominal diagnostics as bounds.
- [x] Add d3d_optimization/selective_transport.py and verification/run_selective_transport_study.py. Rank donors/cells by bound contribution, refine hit neighborhoods and dominant J increments, re-trace from exact source, replay changed paths, and preserve earlier valid kernel certificates.
- [x] Retain geometry box and every donor. Never reset to nominal state or mix certificates across models. Record skipped/failed attempts and budget overruns.
- [x] Test analytic contraction, retained source/branches, budget stop, missing-hit handling, cache tampering; then run real IMP/OMP 16-donor point experiments with identical budgets.

## Task 4 — Integration, evidence and review
- [x] Feed accepted requalification kernels into bounded candidate comparison; measure spatial-comparison costs at 320 and 1280; cover 640 in numerical requalification. Clearly separate frozen distribution accuracy from continuous path accuracy.
- [x] Run original 117-test command plus new tests and appropriate campaign regressions.
- [x] Independent code review; fix substantive findings. Recheck protected source/results hashes.
- [x] Write Korean results with actual candidate verdicts, runtime breakdown, bound contractions, unresolved cases, and whether the added validation should be enabled by default. No global optimum claim.

## Execution outcome

Completed first measured cycle; see results/selective-candidate-validation-20261103/REPORT.ko.md. No production search-loop integration, automatic candidate promotion, new optimization campaign, or global proof. Original broad source-guard failure is preserved with a separately verified unused-code exception. Continuous point bounds remain too wide for ranking; therefore validation stays finalist-only and optional. Same-horizon selected-cell refinement was separately replayed after fixing exact terminal-plane rounding. Wall-time budgets stop between atomic operations and report overruns.
