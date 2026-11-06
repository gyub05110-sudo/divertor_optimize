# Independent code/protocol review

Reviewed 2026-11-06 before thermal campaign launch. Scope: `nf_unused_inputs.py`, `nf_unused_live.py`, their tests, frozen `nf_base_live`/`nf_base_budget` loop, evaluator/input/GPU call paths, protocol, selection, preflight, downloads, and source audit. Source files were not changed by this reviewer; only synthetic tests and read-only hash checks were run.

Verdict: no remaining issue that blocks freezing and running the registered campaign.

## Finding resolved before freeze

- **P2 — enforce the preflight-to-runtime provenance chain.** Initial `load_cases` freshly hashed staged inputs but did not compare those bytes with preflight's frozen `inputs`; a consistently replaced parquet plus manifest could therefore be accepted before first launch. The revised implementation compares normalized input hashes, the exact download entry, raw byte count and SHA-256, and raw/staged parquet equality. `prepare_campaign` now also verifies preflight's protocol, selection, download, original-input and two preparation-driver hashes. Source re-review confirms the fix. New regression tests cover consistent parquet/manifest replacement and a changed preflight selection hash.

## Evidence and design checks

- Independently checked the actual preflight protocol/selection/download hash chain, both preparation-driver hashes, all four original input hashes, all four input hashes for each of the three selected roots, and raw/staged parquet identity: all matched.
- Before the provenance fix, this reviewer ran `python -m pytest verification/test_nf_base_live.py verification/test_nf_unused_live.py verification/test_nf_unused_inputs.py -q -p no:cacheprovider`: **23 passed in 6.21 s**, using a D: temporary root. No thermal/GPU evaluation ran. After the fix, the parent reported **15 live/input tests passed in 4.45 s**; those updated checks were not rerun by this reviewer.
- Selection is source/geometry based. The original wall, supplied geometry bytes, OMP endpoints, engineering policy, selected time and operating reference are preserved. The source audit documents the anonymous-record limitation and array-based novelty checks; the study must not claim established physical shot numbers.
- Launch time originates in the parent immediately before subprocess creation. Imports, setup, baseline calculation, generation/rejections, coarse/fine scores and first GPU initialization are charged. Completed results after 200 s are discarded. This is an information cutoff with measured in-flight overrun, not a hard process kill. A baseline fallback is admitted only after its on-time measurement.
- Each arm is a fresh child with its own initially empty CuPy/CUDA disk caches on D:. GPU scoring returns through `cp.asnumpy`, synchronizing scored work. No previously measured candidate heat values or saved winners enter the policy. Fixed-path reuse occurs within an arm.
- Both policies use the same RNG seed and geometry proposal/filter stream. Prefix identity and attempt counts are audited; empty prefixes are explicitly uninformative. Arm order alternates. Failures or missing worker records do not suppress the second arm. Existing pair directories and exclusive outputs prevent silent overwrite or rerun.
- On-time eligible fine events reconstruct incumbent, best new candidate, first 5% improvement and counts. Process exit/record status and serial ordering are checked. Baseline equality is checked when both baselines exist. Runtime package metadata is separately recorded in `environment.json`, as confirmed by the parent.

## Nonblocking robustness note

`audit_pair` currently parses worker JSON outside its guarded validation. A worker interrupted during exclusive JSON serialization could leave a truncated file, causing audit to raise and stop later cases instead of producing an explicit failed-pair summary. Normal numerical failures and missing JSON are handled. Harden malformed-output handling in a future additive revision (or before freeze if already underway); preserve the partial file and do not rerun a failed thermal arm to conceal it. This is not evidence of a failure in the current, not-yet-started campaign.

Review limitations: this review establishes implementation consistency and input integrity, not physical validation, a universal optimizer-speed advantage, or an independent-discharge statistical sample larger than the three selected records. No thermal reruns or outcome-dependent selection were performed by the reviewer.
