import pytest
from verification.nf_unused_inputs import nearest_time


def test_equal_distance_uses_earlier_stored_time_without_interpolation():
    assert nearest_time([1420., 1380., 1500.]) == (1, 1380.)
    assert nearest_time([1360., 1400., 1440.]) == (1, 1400.)


@pytest.mark.parametrize("times", ([1300., 1500.], [1400., 1400.], [], [float("nan")]))
def test_invalid_time_does_not_silently_choose_another_frame(times):
    with pytest.raises(ValueError):
        nearest_time(times)

