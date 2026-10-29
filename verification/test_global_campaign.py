import importlib.util
import json
import numpy as np
import pytest


def api():
    assert importlib.util.find_spec('d3d_optimization.campaign_state'), 'durable campaign state missing'
    from d3d_optimization import campaign_state
    return campaign_state


def test_atomic_roundtrip_corrupt_current_and_previous_recovery(tmp_path):
    m=api(); path=tmp_path/'state.json'
    m.atomic_save(path,{'value':1}); m.atomic_save(path,{'value':2})
    assert m.read_state(path)=={'value':2}
    path.write_text('{broken')
    assert m.read_state(path)=={'value':1}
    with pytest.raises(ValueError): m.atomic_save(path,{'value':float('nan')})


def test_mid_generation_resume_replays_identical_rng_and_proposals(tmp_path):
    m=api(); bounds=np.array([[-1,1],[-2,2]])
    a=m.new_population([[.1,.2],[.2,.3],[-.1,.4],[.5,-.2],[.6,.1]],77)
    for _ in range(13):
        x=m.propose(a,bounds); m.accept(a,x,float(x@x))
    path=tmp_path/'state.json'; m.atomic_save(path,a); b=m.read_state(path)
    for _ in range(17):
        xa=m.propose(a,bounds); xb=m.propose(b,bounds)
        np.testing.assert_array_equal(xa,xb)
        m.accept(a,xa,float(xa@xa)); m.accept(b,xb,float(xb@xb))
    assert a==b


def test_lock_excludes_second_writer_and_releases(tmp_path):
    m=api(); path=tmp_path/'lock'
    with m.process_lock(path):
        with pytest.raises(RuntimeError):
            with m.process_lock(path): pass
    with m.process_lock(path): pass


def test_absolute_deadline_not_extended_by_resume():
    m=api(); config={'deadline_epoch':1000,'search_deadline_epoch':900}
    assert m.should_stop(config,901,search=True)
    assert not m.should_stop(config,901,search=False)
    assert m.should_stop(config,1000,search=False)


def test_transient_windows_replace_lock_is_retried(tmp_path,monkeypatch):
    m=api(); path=tmp_path/'state.json'; real=m.os.replace; attempts=[]
    def transient(a,b):
        attempts.append(1)
        if len(attempts)<3:raise PermissionError('transient scanner lock')
        return real(a,b)
    monkeypatch.setattr(m.os,'replace',transient)
    m.atomic_save(path,{'value':42})
    assert m.read_state(path)=={'value':42}


def test_checkpoint_rotation_gap_recovers_previous(tmp_path,monkeypatch):
    m=api();path=tmp_path/'state.json';m.atomic_save(path,{'value':1});real=m.os.replace
    def interrupted(a,b):
        if str(a).endswith('.tmp'):raise OSError('simulated process interruption after rotation')
        return real(a,b)
    with monkeypatch.context() as patch:
        patch.setattr(m.os,'replace',interrupted)
        with pytest.raises(OSError):m.atomic_save(path,{'value':2})
    assert not path.exists()
    assert m.read_state(path)=={'value':1}
    m.atomic_save(path,{'value':3})
    assert m.read_state(path)=={'value':3}


def test_candidate_identity_includes_junction_and_restores_angle():
    assert importlib.util.find_spec('d3d_optimization.campaign_geometry'), 'candidate wrapper missing'
    from d3d_optimization.campaign_geometry import Candidate
    a=Candidate('IMP',50,tuple([0.]*8+[.3])); b=Candidate('IMP',50,tuple([0.]*8+[.7]))
    assert a.identity!=b.identity
    assert Candidate.from_record(a.record())==a
    np.testing.assert_array_equal(a.coefficients[:16],np.zeros(16))
    assert a.junction_fraction==.3
    with pytest.raises(ValueError): Candidate('IMP',50,tuple([0.]*8+[1.1]))


def test_qualification_shortlist_contains_distinct_basins():
    from d3d_optimization import campaign_validation as module
    assert hasattr(module,'diverse_shortlist'), 'diversity must reach the actual qualification shortlist'
    rows=[dict(objective=1+i*.001,candidate=dict(vector=[i*1e-5,0])) for i in range(20)]
    rows.append(dict(objective=1.1,candidate=dict(vector=[.05,.05])))
    result=module.diverse_shortlist(rows,6,np.array([.01,.01]))
    assert result[0]==rows[0]
    assert rows[-1] in result


def test_finalization_resume_keeps_qualified_incumbent_at_deadline(tmp_path):
    from types import SimpleNamespace
    from d3d_optimization.campaign_validation import Qualifier
    from d3d_optimization.campaign_state import public_save
    incumbent=dict(fingerprint='test',source='IMP',radius_mm=20,status='qualified',reference={},
                   qualification_attempts=[],selected={'peak_MW_m2':2.},
                   validation={'complete':True,'verdict':{'resolution_screen_passed':True}})
    public_save(tmp_path/'manifest.json',incumbent)
    q=object.__new__(Qualifier);q.source='IMP';q.config={'fingerprint':'test','reference':{},'deadline_epoch':0}
    case=SimpleNamespace(folder=tmp_path,mm=20,state={'best':None},e=SimpleNamespace(metadata={'source':'IMP'}))
    q.finalize([case])
    assert json.loads((tmp_path/'manifest.json').read_text())==incumbent


def test_long_run_sleep_inhibition_is_scoped():
    m=api();assert hasattr(m,'keep_awake'), 'worker must inhibit automatic sleep while computing'
    with m.keep_awake():pass


def test_pilot_warm_candidates_only_flow_to_same_source_and_looser_radius(tmp_path):
    from d3d_optimization.campaign_search import warm_vectors
    (tmp_path/'campaign_warm_candidates.json').write_text(json.dumps({'candidates':[
        {'source':'IMP','radius_mm':50,'vector':[0.]*8+[.4]},
        {'source':'OMP','radius_mm':20,'vector':[0.]*16}]}))
    assert len(list(warm_vectors(tmp_path,'IMP',20)))==1
    assert list(warm_vectors(tmp_path,'IMP',100))==[]
