from pathlib import Path
import numpy as np
import pytest
from d3d_optimization.broadening_study import acceptance, preserve_check, save, load_json

ROOT = Path(__file__).resolve().parents[1]


def level(ratio=.8, audit=True, absolute=10.):
    return dict(feasible=True, audit_passed=audit, peak_ratios=[ratio]*5,
                peak_q=[absolute*ratio]*5, baseline_peak_q=[absolute]*5)


def test_resolution_screen_rejects_unresolved_or_grid_dependent_gains():
    good = [level(.801), level(.8002), level(.8001), level(.8003)]
    assert acceptance(good)['improvement_confirmed']
    bad = [level(.98), level(.99), level(1.), level(1.)]
    assert not acceptance(bad)['improvement_confirmed']
    good[-1]['audit_passed'] = False
    assert not acceptance(good)['numerically_verified_feasible']


def test_ratio_cancellation_cannot_hide_absolute_peak_nonconvergence():
    rows = [level(absolute=8),level(absolute=9),level(absolute=10),level(absolute=10)]
    assert not acceptance(rows)['resolution_screen_passed']


def test_nonfinite_or_incomplete_evidence_fails_closed():
    assert not acceptance([])['improvement_confirmed']
    rows = [level() for _ in range(4)]
    rows[2]['peak_ratios'][0] = float('nan')
    assert not acceptance(rows)['improvement_confirmed']


def test_preserved_original_production_files():
    assert preserve_check(ROOT)['unchanged']


def test_json_no_nonfinite_and_no_overwrite(tmp_path):
    p = tmp_path/'data.json'
    save(p, {'value': 1}, exclusive=True)
    with pytest.raises(FileExistsError):
        save(p, {'value':2}, exclusive=True)
    assert load_json(p)['value'] == 1
    with pytest.raises(ValueError):
        save(tmp_path/'nan.json', {'value': float('nan')})


def test_provenance_detects_model_input_and_missing_inventory(tmp_path, monkeypatch):
    import d3d_optimization.broadening_study as study
    root, output = tmp_path/'root', tmp_path/'run'
    root.mkdir()
    (root/'model.py').write_text('original')
    (root/'input.dat').write_text('input')
    snapshot = 'preserve.json'
    save(root/snapshot, {'model.py': study.digest(root/'model.py')})
    monkeypatch.setattr(study,'PRESERVATION',snapshot)
    monkeypatch.setattr(study,'REQUIRED',('input.dat',))
    monkeypatch.setattr(study,'STUDY_FILES',('model.py',snapshot))
    study.capture_provenance(root,output)
    study.verify_provenance(root,output)
    (root/'input.dat').write_text('changed')
    with pytest.raises(ValueError,match='changed'):
        study.verify_provenance(root,output)
    (root/'input.dat').write_text('input')
    record=load_json(output/'provenance.json')
    record['files'].pop('input.dat')
    save(output/'provenance.json',record)
    with pytest.raises(ValueError,match='incomplete'):
        study.verify_provenance(root,output)


def test_independent_length_integration_matches_variable_metric():
    from verification.validate_broadening import independent_lengths
    from shapely.geometry import LineString
    from d3d_optimization.frozen_mapping import Hits
    class Field:
        def magnetic(self,p):
            z=np.asarray(p)[:,1]
            metric=2+z*z
            return np.column_stack([np.ones(len(p)),np.zeros(len(p)),np.sqrt(metric*metric-1)])
    paths=[LineString([[2,0],[2,.5],[2,1]])]
    # Integrate 2+z^2 to .7, including a fractional final segment.
    value=independent_lengths(Field(),paths,Hits(np.array([0]),np.array([[2,.7]]),np.array([.7])))
    assert value[0]==pytest.approx(1.4+.7**3/3,rel=1e-13)


def test_original_failed_screen_is_not_erased_by_refinement():
    from d3d_optimization.broadening_study import gate_original_screen
    good=acceptance([level() for _ in range(4)])
    verdict=gate_original_screen(good, {'baseline_admissible':False}, {'feasible':True})
    assert not verdict['improvement_confirmed']
    assert verdict['verification_status']=='original_screen_failed'
    assert verdict['refined_verification_status']=='passed'


def test_verifier_refuses_to_overwrite_completed_or_partial_validation(tmp_path):
    from verification.validate_broadening import validate_source
    (tmp_path/'validation.json').write_text('{"complete":true}')
    with pytest.raises(FileExistsError):
        validate_source(tmp_path)
    assert load_json(tmp_path/'validation.json')['complete']
