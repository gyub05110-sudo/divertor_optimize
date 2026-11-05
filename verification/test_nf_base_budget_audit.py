from verification.nf_base_budget_audit import replay_scores


def test_one_batch_produces_one_refined_decision_only_when_it_finishes():
    baseline = {"peak": 10., "A2": 2.}
    rows = []
    for i, coarse_peak in enumerate((7., 6., 7., 6.)):
        rows.append(dict(i=i, generation_s=2.,
                         coarse=dict(seconds=1., eligible=True, peak=coarse_peak),
                         fine=dict(seconds=14., eligible=i != 0, peak=9.-i, A2=3.)))
    setup = {"low_s": 4., "high_s": 10., "baseline": baseline}
    def run(policy, limit):
        return replay_scores(rows, setup, policy, batch=4, budget=limit)[0]
    assert run("screen", 39.99)["selected"] is None
    screen = run("screen", 40.)
    assert screen["fine_ids"] == [1] and screen["first_5pct_s"] == 40.
    assert screen["seen_coarse"] == [0, 1, 2, 3]
    direct = run("direct", 42.)
    assert direct["fine_ids"] == [0, 1] and direct["eligible_fine"] == 1
    assert direct["first_5pct_s"] == 42.


def test_missing_high_setup_prevents_even_baseline_selection():
    setup = {"low_s": 1., "high_s": 201., "baseline": {"peak": 9., "A2": 2.}}
    result, _ = replay_scores([], setup, "direct", batch=4, budget=200.)
    assert result["peak"] is None and result["A2"] is None
    assert result["elapsed_last_completed_s"] == 0.
