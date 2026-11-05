import pytest
from verification.nf_base_budget import replay


def sample(i, c, f, low=1, high=3, generation=1):
    return dict(i=i, generation_s=generation,
                coarse=dict(seconds=low, eligible=c is not None, peak=c),
                fine=dict(seconds=high, eligible=f is not None, peak=f, A2=2.0))


def test_shortlist_uses_only_seen_coarse_results_and_charges_full_batch():
    rows = [sample(0, 1.3, 1.1), sample(1, 1.2, 1.4),
            sample(2, None, .8), sample(3, 1.15, 1.0)]
    fast = replay(rows, "screen", 2, 3, 1.5, budget=16)
    assert fast["fine_ids"] == [3] and fast["peak"] == 1.0
    assert fast["first_5pct_s"] == 16
    late = replay(rows, "screen", 2, 3, 1.5, budget=15.99)
    assert late["fine_ids"] == [] and late["peak"] == 1.5
    assert late["A2"] is None


def test_direct_search_counts_only_completed_fine_improvements():
    rows = [sample(0, 1.1, None), sample(1, 1.1, 1.2), sample(2, 1.0, .8)]
    direct = replay(rows, "direct", 2, 3, 1.5, budget=11)
    assert direct["fine_ids"] == [0, 1] and direct["selected"] == 1
    assert direct["first_5pct_s"] == 11 and direct["eligible_fine"] == 1
    pending = replay(rows, "direct", 2, 3, 1.5, budget=1)
    assert pending["peak"] is None and pending["first_5pct_s"] is None
