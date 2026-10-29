# Review and verification record

Independent reviewer found no material blockers before the final run. Reviewed gpu_footprint.py, accelerated_transport.py, regional_transport.py, optimize_regional_transport.py and the new tests.

- Three new tests observed RED before implementation, then all passed.
- Full suite: 177 passed, 1 skipped (28.17 seconds).
- Reviewer independently repeated regional tests and real IMP/OMP perturbations: CPU/GPU peak difference zero, feasibility identical, opposite-region movement zero.
- Reviewer checked every named-surface export index against actual surface_segments and target coordinates for both sources; all slices exactly matched.
- Accelerated implementation preserves the original CDF power/geometry/thermal/continuous-peak logic. Only narrow periodic Gaussian density summation runs on GPU. Broad kernels remain on CPU.
- Actual OMP320 preflight evaluates three small perturbations: maximum profile error relative to peak <5e-16, peak difference zero; median speedup 1.459x for evaluation, not an end-to-end search claim.
- GPU installation used an isolated D target. Existing Torch is CPU-only and had conflicting OpenMP libraries; it is not used. No unsupported duplicate-OpenMP override was enabled for production.
- Regional constraints: only IMP's two surfaces or OMP's four surfaces vary; global wall peak is minimized to discourage moving heat onto an uncounted receiver.
- Final validation uses the preserved CPU implementation, independent first-hit and thermal integration, four standard mesh/interpolation levels and wall-only refinement. Acceptance thresholds are unchanged. Regional invariance is additionally verified from the finest full-wall coordinates.

See each source/validation.json and source/regional_validation.json for actual final verdicts. Neither measured-data agreement nor global optimality is claimed.

Qualification review: existing target and screen destinations are now checked before any screen write. A rerun against qualified/IMP raised FileExistsError and left the screen SHA256 unchanged. The exact IMP qualifier before this nonnumerical guard fix is retained with its original provenance; the current guard-protected script is separately hashed. Original lowest-peak candidates and all their validation evidence remain intact.

Final delivery: qualified/IMP case 250 and qualified/OMP case 236 passed all four CPU validation levels and unchanged thresholds. Both fixed-region coordinate errors are exactly zero; delivered power is 1 MW each with no unhit or unresolved power. verify_delivery.py passed exact CSV/array comparisons, source and snapshot hashes, screening evidence hashes and all 180 pre-existing Python hashes. PNGs were visually inspected; the transport plot axis was restricted to the visible source-offset range to remove unused vertical space.
