from d3d_optimization.broadening_study import save
from verification.plot_reference_transport import write_report


def test_report_accepts_partial_validation_and_never_calls_it_complete(tmp_path):
    save(tmp_path/'IMP'/'manifest.json',dict(source='IMP',reference={'power_W':1e6,'density_m3':1e19,'chi_perp_m2_s':.3,'lambda_q_m':.003},
        mode='preflight',smoke_test=False,metadata={'tubes':48},models={'thermal':{'baseline':{'feasible':False,'missing_power':.1}}},
        candidates={'baseline':[]},cross_evaluation={}))
    save(tmp_path/'IMP'/'validation.json',{'complete':False,'comparisons':{},'levels':[]})
    write_report(tmp_path)
    s=(tmp_path/'REPORT.md').read_text(encoding='utf-8')
    assert '검증 진행 중' in s
    assert '0.1' in s or '10%' in s
    assert 'one-way' in s
