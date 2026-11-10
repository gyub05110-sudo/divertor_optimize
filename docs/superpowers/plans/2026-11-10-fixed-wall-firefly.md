# Fixed-wall FIREFLY comparison implementation plan

The approved questions, population, compute budget and precision rule are frozen in `results/nf-stage5-firefly-20261110/protocol.ko.md`. Use the supplied Stage 4 worktree and independent input, native, and metrics tasks.

- [x] Register 51 roles against the original thermal and geometry receipts, content hashes and protocol before any prospective FIREFLY heat. Verify exact walls, cases and optimizer seeds.
- [x] Export the supplied equilibria, individual donor locations/weights and fixed walls in `verification/stage5_inputs.py`. Validate Stage 4 source, field conventions, LCFS, original native peak, and exact wall coordinates.
- [x] Test and validate the FIREFLY/FLARE bridge in `verification/stage5_native*` on the nominated shelf and floor examples. Make a shared mesh per input and preserve field, source, wall, lower-leg and power accounting checks.
- [x] Run the sealed fixed-wall matrix with six paired seeds, native MPI and original conditions; store atomic raw output, audit and resumable receipts on D.
- [x] Reconstruct comparable ring powers, peaks, areas, and paired uncertainty using `verification/stage5_metrics.py`. Check the prediction and actual returned walls separately, including the 30757 reversal. Publish verified outputs, Korean results, reproducer and original external evidence location.
