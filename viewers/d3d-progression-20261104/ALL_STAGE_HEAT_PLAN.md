# All-stage heat visualization

Goal: every archived baseline, running-record stage and refined result has its own color and q/R/Z probe in both views. Geometry, inputs, original milestones and qualification labels stay unchanged. This is display-only evaluation, not a new optimization or qualification campaign.

Design: add a hash-bound stage profile index beside the immutable data/leg files. Reuse fixed 320-path preparation per source, evaluate each unique archived candidate once, retain raw NPZ and provenance on D:, and publish a profile only after geometry, original peak and power checks pass. Existing refined profiles remain preferred. Use one fixed source color scale across stages.

Alternatives considered: borrowing the final profile is physically wrong; interpolating distributions between records does not evaluate changed arrival locations. Re-evaluating the exact archived candidate at its original resolution is the selected approach.

- [x] Add failing tests for stage profile lookup, identity/resolution isolation and conservation/shape/peak rejection.
- [x] Implement resumable export_stage_profiles.py and data checks. Verify original fixed input hashes, compare CPU/GPU q on one case per source, and reuse prepared paths. Save every candidate NPZ and JSON before advancing checkpoint.
- [x] Evaluate all 237 milestones, deduplicating only identical candidate IDs under the same fixed input/model/grid. Compare computed and logged peak at relative 1e-7; stop on mismatch rather than label a different model as the original result.
- [x] Wire the independent index into heat-probe.mjs/app.mjs; distinguish 320-path display reevaluation from 1280-path refined results. Retain latest-request protection against stale async profiles. Limit in-memory profile cache for long playback.
- [x] Verify every milestone profile exists and hashes match; independently recompute display-bin power, confirm wall coordinates and each record's original peak. Verify data/legs/original final artifacts were not changed.
- [x] Browser QA: intermediate IMP/OMP, all runs/radii, baseline/overlay/candidate modes, final-to-stage switching, both 3D and R-Z probes, color toggle and rapid stepping. Open actual screenshots and check console.
- [x] Update README and write machine-readable verification evidence. No GitHub push.

Adversarial checks: no final q reused on coarse geometry; no source swap; no null-to-zero conversion; no stale final profile after step changes; no 1280 label on 320 calculations; no qualification promotion; missing/failing cache must fail closed; keep baseline and final selection metrics untouched. Full 320/640/1280 convergence is not claimed for intermediate visualization.
