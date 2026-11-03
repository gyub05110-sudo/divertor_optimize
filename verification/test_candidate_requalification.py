import importlib
import json
from pathlib import Path
import pytest
from d3d_optimization.campaign_geometry import Candidate


def api():
    assert importlib.util.find_spec('verification.run_candidate_requalification'), 'candidate requalification runner missing'
    return importlib.import_module('verification.run_candidate_requalification')


def manifest(tmp_path, data):
    path=tmp_path/'manifest.json';path.write_text(json.dumps(data));return path


def test_campaign_loader_preserves_identity_reference_and_provenance(tmp_path):
    c=Candidate('IMP',20,tuple([0.]*8+[.4]));reference={'power_W':1e6}
    p=manifest(tmp_path,dict(candidate=c.record(),reference=reference,selected={'peak_MW_m2':2.3}))
    row=api().load_candidate(p,'IMP',20,reference)
    assert row['candidate']==c.record() and row['identity']==c.identity
    assert row['origin_sha256']==api().digest(p)
    assert row['historical_peak_MW_m2']==2.3
    assert row['requires_current_validation'] is True


def test_historical_omp_loader_keeps_exact_sixteen_coefficients(tmp_path):
    values=[i*.001 for i in range(16)];reference={'power_W':1e6}
    p=manifest(tmp_path,dict(radius_mm=50,reference=reference,coefficients_m=values+[0.]*8))
    row=api().load_candidate(p,'OMP',50,reference)
    assert row['candidate']==Candidate('OMP',50,tuple(values)).record()
    assert row['identity']==Candidate('OMP',50,tuple(values)).identity


@pytest.mark.parametrize('change',['reference','source','radius','opposite'])
def test_loader_rejects_mismatched_inputs_and_geometry(tmp_path,change):
    reference={'power_W':1e6};c=Candidate('OMP',50,tuple([0.]*16))
    data=dict(reference=reference,candidate=c.record())
    if change=='reference':data['reference']={'power_W':2e6}
    elif change=='source':data['candidate']=Candidate('IMP',50,tuple([0.]*8+[.4])).record()
    elif change=='radius':data['candidate']['radius_mm']=20
    else:data.pop('candidate');data.update(radius_mm=50,coefficients_m=[0.]*23+[.01])
    with pytest.raises(ValueError):api().load_candidate(manifest(tmp_path,data),'OMP',50,reference)


def test_new_output_refuses_existing_and_protected_tree(tmp_path):
    protected=tmp_path/'old';protected.mkdir();(protected/'keep').write_text('unchanged')
    with pytest.raises(ValueError):api().new_output(protected/'new',[protected])
    with pytest.raises(FileExistsError):api().new_output(protected,[])
    assert (protected/'keep').read_text()=='unchanged'
    created=api().new_output(tmp_path/'new',[protected]);assert created.is_dir()


def test_preservation_reports_mutation_without_overwriting(tmp_path):
    p=tmp_path/'input';p.write_text('before');hashes={str(p):api().digest(p)}
    assert api().check_preservation(hashes)['unchanged']
    p.write_text('after');report=api().check_preservation(hashes)
    assert not report['unchanged'] and report['changed']==[str(p)]
    assert p.read_text()=='after'


def test_source_runs_all_four_levels_and_wall_refinement_even_after_failed_audit(tmp_path,monkeypatch):
    """Cheap orchestration fixture: expensive numerical operations are substituted."""
    from types import SimpleNamespace
    m=api();reference={'power_W':1e6,'density_m3':1e19,'chi_perp_m2_s':.3,'lambda_q_m':.003,'gamma':7.,'kappa0':2000.}
    original=SimpleNamespace(metadata={'source':'IMP'},baseline={'elapsed_seconds':.1},direct=None,
                             wall_step=.001,thermal_paths=None)
    prepared=[]
    monkeypatch.setattr(m,'prepare_boundary',lambda *a,**kw:original)
    class SharedQualifier:
        def __init__(self,*args):self.reference=reference
        def original(self,index):prepared.append(index);return original
    monkeypatch.setattr(m,'Qualifier',SharedQualifier)
    monkeypatch.setattr(m,'BoundaryEvaluator',lambda *a,**kw:original)
    monkeypatch.setattr(m,'evaluator',lambda *a:original)
    monkeypatch.setattr(m,'evaluate_candidate',lambda *a:{'geometry_valid':True,'geometry_screen':{'passed':True}})
    monkeypatch.setattr(m,'audit',lambda *a:{'audit_passed':False,'feasible':True,'peak_q':2.,'peak_MW_m2':2.})
    monkeypatch.setattr(m,'geometric_audit',lambda *a:{'passed_independent':True})
    # This fixture checks orchestration, not repeated atomic disk replacement.
    monkeypatch.setattr(m,'save',lambda *a:None)
    exported=[];monkeypatch.setattr(m,'export_evidence',lambda folder,*a:exported.append(folder))
    c=Candidate('IMP',20,tuple([0.]*8+[.4]));rows={name:dict(candidate=c.record(),identity=c.identity) for name in ['a','b']}
    report=dict(timings=[],candidates={},baselines={},baseline_wall_refinement={})
    m.run_source(tmp_path,'IMP',rows,tmp_path,dict(reference=reference,backend='cpu'),report)
    assert prepared==[1,2,3]
    assert len(report['baselines']['IMP'])==4
    for record in report['candidates'].values():
        assert len(record['levels'])==4 and record['complete']
        assert 'wall_only_refinement' in record
        assert not record['verdict']['resolution_screen_passed']
    assert len(exported)==15  # Five levels for two candidates and shared baseline.


def test_code_guard_tracks_consumed_modules_and_ignores_unimported_neighbors(tmp_path):
    from types import SimpleNamespace
    loaded=tmp_path/'loaded.py';loaded.write_text('pass')
    neighbor=tmp_path/'unimported.py';neighbor.write_text('pass')
    external=tmp_path.parent/'external.py'
    modules={'loaded':SimpleNamespace(__file__=str(loaded)),
             'external':SimpleNamespace(__file__=str(external)),
             'builtin':SimpleNamespace()}
    assert api().consumed_code_paths(tmp_path,modules)==[loaded.resolve()]
