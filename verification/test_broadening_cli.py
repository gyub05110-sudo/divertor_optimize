import pytest
from compare_broadening import parse_args, model_key


def test_cli_requires_output_and_finite_widths():
    with pytest.raises(SystemExit):
        parse_args([])
    for value in ['nan', '-1', 'inf']:
        with pytest.raises(SystemExit):
            parse_args(['--output','unused','--sigma-mm',value])


def test_cli_preserves_separate_sources_and_explicit_zero_baseline():
    args = parse_args(['--output','unused','--sigma-mm','.1','.25'])
    assert args.source == 'both'
    assert args.sigma_mm == [0., .1, .25]
    assert model_key(0) == 'direct'
    assert model_key(.25) != model_key(.1)


def test_cli_rejects_invalid_numerical_controls():
    for option, value in [('--tubes','0'),('--generations','-1'),('--population','3'),
                           ('--wall-step-mm','nan'),('--anchor-mm','0')]:
        with pytest.raises(SystemExit):
            parse_args(['--output','unused',option,value])


def test_out_of_bound_historical_candidate_is_skipped_not_clipped():
    from compare_broadening import select_warm_start
    import numpy as np
    c,reason=select_warm_start({'verified_coefficients_m':[.05]*24},24,.01)
    assert c is None and reason=='outside_requested_coefficient_bound'
    c,reason=select_warm_start({'verified_coefficients_m':[.005]*24},24,.01)
    np.testing.assert_array_equal(c,np.full(24,.005))
