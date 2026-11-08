# DIII-D prospective benchmark implementation plan

Use the approved three-arm protocol in results/nf-stage3-20261108/protocol.ko.md. The deliverables are a prospective registration, timed real runs, audited checkpoint/5% outcomes, fixed-design ablations, and an inspectable report.

- [x] Input registration in verification/stage3_inputs.py. Test deterministic ranking, exclusions, exact time and byte/array duplicate detection before retrieving selected data. Record official pinned sizes and SHA-256. Preserve all 12 original slots.
- [x] Commission the direct fine strategy on previously seen inputs in isolation. Test that no coarse evaluator is constructed, late fine observations never win, and existing screened runs retain their original entry path.
- [x] Implement verification/stage3_campaign.py with an immutable job plan, input/source hashes, balanced serial order, separate process receipt, non-overwriting recovery and actual GPU commands. Test three methods leading once per shot, excluded registrations, immutable case binding, resume and replay of sealed evidence. Run the full test suite and a known-input direct smoke before starting the held-out clock.
- [x] Seal registration and source hashes; execute all scheduled jobs on D: and retain process failures separately. Review partial status, preserve interrupted attempts and diagnose technical interruptions without switching files or retuning.
- [x] Reconstruct 200/400/600 s decisions and time to first ≥5% reduction from the raw journals, retaining taxonomy and denominators. Re-evaluate registered full-winner ablations after the timed phase. Confirm hashes and report outcome scope in Korean.
