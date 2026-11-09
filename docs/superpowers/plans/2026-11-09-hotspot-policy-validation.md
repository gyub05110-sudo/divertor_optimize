# Hotspot policy validation implementation plan

The approved design is frozen in `results/nf-stage4-20261109/protocol.ko.md`. Use independent tasks for cohort registration, geometry, the public search, and the timed benchmark.

- [ ] Add failing fixed-wall, endpoint and legacy-support tests. Implement the geometric early-exit witness in `d3d_optimization/user_precheck.py` and check the Stage 3 segments 4/5 case and lower-X failures using development data only.
- [ ] Add failing family-mask, tie, abstention, seed and deadline tests. Extend `user_search.py` and `optimize_divertor.py`; integrate identical preflight and pre-search prediction in all screened arms. Exercise known inputs.
- [ ] Implement and test prospective exclusion, deduplication, freeze and replay in `verification/stage4_inputs.py`. Freeze and verify 24 sources. Do not replace missing frames.
- [ ] Test a 240-slot balanced schedule and independent forced controls before implementing `verification/stage4_campaign.py`. Bind raw inputs, code, interpreter and order. Preserve process receipts, shocks to elapsed time and resumability. Run full tests and publish the sealed protocol before timed heat.
- [ ] Run the schedule through the public entry point, verify event clocks, fixed masks and physical exports, then compute three-seed product comparisons and one-seed multiclass/always-floor controls in `verification/stage4_report.py`. Re-evaluate preset sensitivity triggers and distinguish eligible decisions and failed inputs. Publish Korean results, raw evidence inventory and rerun instructions.
