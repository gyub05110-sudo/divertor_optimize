import pytest

from verification.nf_finalist_decisions import compare_peaks, hot_intervals


def test_paired_improvement_can_have_overlapping_numerical_ranges():
    a = {"cubic": 1.0, "quintic": 1.03}
    b = {"cubic": 1.01, "quintic": 1.04}
    result = compare_peaks(a, b)
    assert result["paired_order"] == "first_lower"
    assert result["interval_order"] == "overlap"
    assert result["paired_delta_fraction"][1] < 0


def test_reversal_is_distinguished_from_separation():
    assert compare_peaks({"a": 2., "b": 1.}, {"a": 1.9, "b": 1.1})["paired_order"] == "changes_or_ties"
    assert compare_peaks({"a": 3., "b": 4.}, {"a": 1.9, "b": 2.})["interval_order"] == "second_lower"
    with pytest.raises(ValueError):
        compare_peaks({"a": 3.}, {"b": 2.})


def test_hot_bands_preserve_two_separated_peaks_and_interpolate_edges():
    import numpy as np
    u = np.arange(5.)/4
    q = np.array([0., 2., 0., 2., 0.])
    np.testing.assert_allclose(hot_intervals(u, q, .5), [[.125,.375],[.625,.875]])
    np.testing.assert_allclose(hot_intervals(u, np.ones(5), .99), [[0.,1.]])
