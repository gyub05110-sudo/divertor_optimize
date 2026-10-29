"""Reporting must tolerate validation checkpoints without implying a pass."""
from d3d_optimization.broadening_study import save
from verification import plot_broadening as plot


def test_report_handles_partial_validation(tmp_path, monkeypatch):
    model = 'sigma_0.1mm'
    row = dict(baseline_admissible=True, missing_power=0.,
               feasible=True, worst_peak_ratio=.9)
    manifest = dict(source='IMP', mode='optimize', smoke_test=False,
                    metadata={'tubes': 320}, candidates={model+'_best': []},
                    models={model: {'config': {'sigma_ref_m': .0001}, 'baseline': row}},
                    cross_evaluation={model: {model+'_best': row}})
    save(tmp_path/'IMP'/'manifest.json', manifest)
    save(tmp_path/'IMP'/'validation.json', {'complete': False, 'comparisons': {}})
    monkeypatch.setattr(plot, 'source_figure', lambda *args: None)
    plot.write_report(tmp_path)
    report = (tmp_path/'REPORT.md').read_text(encoding='utf-8')
    assert '검증 전' in report
    assert '| 확인 |' not in report


def test_report_labels_search_power_and_reports_failed_refined_baseline(tmp_path, monkeypatch):
    model = 'sigma_0.1mm'
    row = dict(baseline_admissible=True, missing_power=.0007,
               feasible=True, worst_peak_ratio=.85)
    manifest = dict(source='OMP', mode='optimize', smoke_test=False,
                    metadata={'tubes': 320}, candidates={model+'_best': []},
                    models={model: {'config': {'sigma_ref_m': .0001}, 'baseline': row}},
                    cross_evaluation={model: {model+'_best': row}})
    save(tmp_path/'OMP'/'manifest.json', manifest)
    save(tmp_path/'OMP'/'validation.json', {
        'complete': True,
        'levels': [{'models': {model: {'baseline': dict(row, feasible=False, missing_power=.0026)}}}],
        'comparisons': {model: {model+'_best': {
            'verification_status': 'incomplete_or_invalid_evidence', 'improvement_confirmed': False}}}})
    monkeypatch.setattr(plot, 'source_figure', lambda *args: None)
    plot.write_report(tmp_path)
    report = (tmp_path/'REPORT.md').read_text(encoding='utf-8')
    assert '탐색 미해결 전력' in report
    assert '검증 기준 형상 적합=False' in report
    assert '0.26%' in report
