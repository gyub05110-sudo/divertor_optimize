# Final outcome audit

Independent saved-file audit passed: 2,147 checks, no discrepancies. All 18 raw arm records and nine paired summary rows were reconstructed with standard-library arithmetic; no production audit/summary imports or thermal evaluations were used. See `final-outcome-audit.json` and `audit-temp/final_outcome_audit.py` for evidence.

All registered source/code snapshot hashes, original material-wall inputs, staged inputs, raw downloads, and protocol hashes remained unchanged. The final supervisor hash is `61dc65d0e3ad39432f8a248948abdae80f82496e4d108358e4833d3e019d1bcd`; its pre-freeze malformed-JSON handling fix supersedes the earlier nonblocking note in `code-review.md`.

Screening produced four strictly lower final peaks, five exact ties, and no losses. The median paired peak difference is 0%. The five ties are also the only pairs within 0.1% or 1%; there are no additional approximate ties.

| Record | Seed | Screen vs direct peak (%) | Direct peak reduction (%) | Screen peak reduction (%) | Direct A2 change (%) | Screen A2 change (%) |
|---|---:|---:|---:|---:|---:|---:|
| deb842cf1e | 20611 | 0.0000 | 42.0033 | 42.0033 | 66.6661 | 66.6661 |
| deb842cf1e | 20627 | -51.1423 | 0.0000 | 51.1423 | 0.0000 | 92.4173 |
| deb842cf1e | 20639 | -14.9989 | 0.0000 | 14.9989 | 0.0000 | 16.9374 |
| fad1624516 | 20611 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| fad1624516 | 20627 | -9.5048 | 0.0000 | 9.5048 | 0.0000 | 0.6492 |
| fad1624516 | 20639 | -7.5211 | 2.5337 | 9.8642 | -6.2826 | 0.0762 |
| 4d14a35fcb | 20611 | 0.0000 | 20.7396 | 20.7396 | 1.5586 | 1.5586 |
| 4d14a35fcb | 20627 | 0.0000 | 13.3378 | 13.3378 | 5.8686 | 5.8686 |
| 4d14a35fcb | 20639 | 0.0000 | 32.8082 | 32.8082 | 14.2929 | 14.2929 |

Peak reduction and A2 change columns use each record's original wall as denominator. Screen-vs-direct peak uses the direct arm's accepted final peak; a negative value favors screening.

Screening reduced peak and increased A2 in eight of nine arms; one retained the original wall. Direct reduced peak in five arms, with increased A2 in four. Its `fad1624516 / 20639` accepted design reduced peak by 2.5337% but reduced A2 by 6.2826%. Screening's A2 gains for `fad1624516 / 20627` and `/20639` are only 0.6492% and 0.0762%; numerical sensitivity must remain visible.

Screening reached a 5% reduction in eight arms, direct in four. Four pairs reached the target only with screening; one reached it with neither. Among the four pairs where both reached it, direct was earlier in three and screening in one. The result supports a pilot endpoint advantage, not a universal first-improvement speed advantage.

The 200-second rule was applied to completed measurements. Ten late fine evaluations were discarded without measurement leakage. Maximum recorded search duration was 214.6681179 s and maximum process duration 215.1404563 s, because in-flight work may finish after the cutoff. Phase-duration sums include such discarded work. Screening completed 108 coarse and 26 fine candidate scores; direct completed 42 fine scores.

This is a pilot on three geometry-admitted unused anonymous public equilibrium records, repeated with three paired seeds each. Nine seed pairs are not nine independent discharges, and numeric physical shot numbers are not established. It does not demonstrate a global optimum, general operating-domain coverage, physical q(s) accuracy, or FIREFLY validation. Posthoc quintic640 sensitivity is separate from the frozen primary comparison and may qualify the small A2 gains without replacing primary outcomes.
