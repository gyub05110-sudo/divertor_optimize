from pathlib import Path
import pytest
from d3d_optimization.reference_study import acceptance, verify_preserved
from compare_reference_transport import parse_args


def row(ratio=.9,peak=10.):
    return dict(feasible=True,audit_passed=True,peak_q=[ratio*peak],baseline_peak_q=[peak],peak_ratios=[ratio])


def test_single_width_acceptance_preserves_original_failures_and_absolute_checks():
    good=[row() for _ in range(4)]
    assert acceptance(good,True)['improvement_confirmed']
    assert not acceptance(good,False)['improvement_confirmed']
    good[2]=row(peak=20)
    assert not acceptance(good,True)['improvement_confirmed']
    good[2]=row();good[0]['audit_passed']=False
    assert not acceptance(good,True)['improvement_confirmed']


def test_preexisting_code_is_unchanged():
    assert verify_preserved(Path(__file__).resolve().parents[1])['unchanged']


def test_cli_reference_inputs_and_nonfinite_rejection():
    args=parse_args(['--output','unused'])
    assert args.power_mw==1 and args.density==1e19 and args.chi==.3
    with pytest.raises(SystemExit):parse_args(['--output','unused','--power-mw','nan'])
    with pytest.raises(SystemExit):parse_args(['--output','unused','--population','2'])
