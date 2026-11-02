import copy
import importlib
import json
from types import SimpleNamespace
import numpy as np
import pytest


def api():
    assert importlib.util.find_spec('d3d_optimization.campaign_incumbents'), 'qualified store missing'
    return importlib.import_module('d3d_optimization.campaign_incumbents')


def example(peak=2.,fraction=.5):
    from d3d_optimization.campaign_geometry import Candidate,GEOMETRY_POLICY
    c=Candidate('IMP',20,tuple([0.]*8+[fraction]))
    ctx=dict(fingerprint='same-inputs',source='IMP',radius_mm=20,reference={'power_W':1e6},
             qualification_policy=api().QUALIFICATION_POLICY)
    audit={'audit_passed':True}
    levels=[dict(method=m,tubes=n,rows={'selected':audit,'baseline':audit},
                 geometry_audit={'passed_independent':True,'geometry_screen':{'policy':GEOMETRY_POLICY,'passed':True}})
            for m,n in [('cubic',320),('cubic',640),('cubic',1280),('quintic',640)]]
    for level in levels:
        level['rows']={'selected':dict(audit,peak_MW_m2=peak),'baseline':dict(audit,peak_MW_m2=3.)}
    manifest=dict(**ctx,status='qualified',candidate=c.record(),selected={'peak_MW_m2':peak},
                  baseline={'peak_MW_m2':3.},cpu_gpu_peak_relative_difference=0.,
                  validation=dict(complete=True,fingerprint=ctx['fingerprint'],qualification_policy=ctx['qualification_policy'],
                                  candidate=c.record(),identity=c.identity,levels=levels,
                                  wall_only_refinement={'rows':{'selected':audit,'baseline':audit},'max_peak_relative_change':0.},
                                  verdict={'resolution_screen_passed':True}))
    return ctx,manifest


def write_artifacts(folder):
    np.savez(folder/'selected_refined.npz',q=np.array([1.,2.]))


def test_worse_candidate_and_restart_do_not_lose_best(tmp_path):
    ctx,m=example();store=api().QualifiedStore(tmp_path,ctx)
    first=store.publish(m,write_artifacts)
    _,worse=example(3.,.6);store.publish(worse,write_artifacts)
    assert api().QualifiedStore(tmp_path,ctx).best()['selected']['peak_MW_m2']==2.
    assert (tmp_path/first['artifact_bundle']/'selected_refined.npz').is_file()
    assert len(store.records())==2


def test_interrupted_publish_keeps_previous_candidate_and_arrays(tmp_path):
    ctx,m=example();store=api().QualifiedStore(tmp_path,ctx);first=store.publish(m,write_artifacts)
    def interrupted(folder):
        write_artifacts(folder)
        raise OSError('simulated crash before publication')
    _,better=example(1.,.7)
    with pytest.raises(OSError):store.publish(better,interrupted)
    assert store.best()==first
    assert len(store.records())==1


def test_corrupt_bundle_is_reported_not_silently_replaced(tmp_path):
    ctx,m=example();store=api().QualifiedStore(tmp_path,ctx);record=store.publish(m,write_artifacts)
    (tmp_path/record['artifact_bundle']/'selected_refined.npz').write_bytes(b'corrupt')
    with pytest.raises(ValueError,match='checksum'):store.best()


@pytest.mark.parametrize('change',['context','policy','incomplete','identity','geometry','parity'])
def test_unverified_or_incompatible_manifest_cannot_be_promoted(tmp_path,change):
    ctx,m=example();store=api().QualifiedStore(tmp_path,ctx)
    if change=='context':m['reference']={'power_W':2e6}
    elif change=='policy':m['validation']['qualification_policy']='legacy'
    elif change=='incomplete':m['validation']['levels'].pop()
    elif change=='identity':m['validation']['identity']='not-the-candidate'
    elif change=='geometry':m['validation']['levels'][0]['geometry_audit']['passed_independent']=False
    elif change=='parity':m['cpu_gpu_peak_relative_difference']=.1
    with pytest.raises(ValueError):store.publish(m,write_artifacts)
    assert store.best() is None


def test_protected_historical_candidates_survive_shortlist_limit(tmp_path):
    from d3d_optimization.campaign_geometry import Candidate
    m=api(); candidates=[]
    for j in range(6):
        c=Candidate('IMP',50,tuple([0.]*8+[.2+j*.1]))
        candidates.append({'candidate':c.record(),'origin':'historical','origin_sha256':str(j),'historical_peak':2+j*.1})
    (tmp_path/'campaign_protected_candidates.json').write_text(json.dumps({'candidates':candidates}))
    protected=m.protected_candidates(tmp_path,'IMP',20)
    assert len(protected)==6
    assert all(c['candidate']['radius_mm']==20 for c in protected)
    assert m.protected_candidates(tmp_path,'IMP',100)==[]
    assert m.protected_candidates(tmp_path,'OMP',20)==[]
    assert not any(c.get('status')=='qualified' for c in protected)


def test_finalizer_tries_protected_candidates_even_with_empty_search_archive(tmp_path,monkeypatch):
    from d3d_optimization.campaign_validation import Qualifier
    from d3d_optimization.campaign_geometry import Candidate
    root=tmp_path/'snapshot';root.mkdir();folder=tmp_path/'IMP/R020';folder.mkdir(parents=True)
    candidates=[Candidate('IMP',20,tuple([0.]*8+[.2+i*.1])) for i in range(6)]
    (root/'campaign_protected_candidates.json').write_text(json.dumps({'candidates':[
        dict(candidate=c.record(),origin='old',origin_sha256=str(i),historical_peak=2+i*.1) for i,c in enumerate(candidates)]}))
    q=object.__new__(Qualifier);q.root=root;q.source='IMP'
    q.config=dict(fingerprint='test',reference={},deadline_epoch=float('inf'))
    novel=Candidate('IMP',20,tuple([0.]*8+[.95]))
    q.shortlist=lambda *a,**kw:[dict(candidate=c.record()) for c in candidates+[novel]]
    attempted=[]
    def reject(case,candidate):
        attempted.append(candidate.identity)
        return {'complete':True,'verdict':{'resolution_screen_passed':False}}
    q.validate=reject
    case=SimpleNamespace(folder=folder,mm=20,state={'best':None},e=SimpleNamespace(metadata={'source':'IMP'}))
    q.finalize([case])
    assert attempted==[c.identity for c in candidates+[novel]]
    report=json.loads((folder/'manifest.json').read_text())
    assert len(report['qualification_attempts'])==7


def test_collect_protected_finds_bundle_after_pointer_crash(tmp_path):
    ctx,m=example();folder=tmp_path/'results/old/IMP/R020'
    api().QualifiedStore(folder,ctx).publish(m,write_artifacts)
    assert not (folder/'manifest.json').exists()
    rows=api().collect_protected(tmp_path)
    assert len(rows)==1 and rows[0]['candidate']==m['candidate']


def test_finalizer_recovers_published_bundle_when_pointer_write_was_interrupted(tmp_path):
    from d3d_optimization.campaign_validation import Qualifier
    ctx,m=example();best=api().QualifiedStore(tmp_path,ctx).publish(m,write_artifacts)
    q=object.__new__(Qualifier);q.source='IMP';q.root=tmp_path
    q.config=dict(fingerprint=ctx['fingerprint'],reference=ctx['reference'],deadline_epoch=0)
    case=SimpleNamespace(folder=tmp_path,mm=20,state={'best':None},e=SimpleNamespace(metadata={'source':'IMP'}))
    q.finalize([case])
    recovered=json.loads((tmp_path/'manifest.json').read_text())
    assert recovered['selected']==best['selected']
    assert recovered['artifact_bundle']==best['artifact_bundle']


def test_summary_uses_published_bundle_not_stale_pointer_or_arrays(tmp_path):
    from summarize_global_campaign import selection
    ctx,m=example();store=api().QualifiedStore(tmp_path,ctx)
    first=store.publish(m,write_artifacts)
    (tmp_path/'manifest.json').write_text(json.dumps(first))
    (tmp_path/'selected_refined.npz').write_bytes(b'outdated convenience export')
    _,better=example(1.,.6);second=store.publish(better,write_artifacts)
    selected,path=selection(tmp_path,ctx)
    assert selected['selected']['peak_MW_m2']==1.
    assert path==tmp_path/second['artifact_bundle']/'selected_refined.npz'
    assert np.load(path)['q'].tolist()==[1.,2.]


def test_collect_protected_imports_manifests_not_only_coarse_progress(tmp_path):
    ctx,m=example();warm=tmp_path/'old';folder=warm/'IMP/R020';folder.mkdir(parents=True)
    (folder/'manifest.json').write_text(json.dumps(m))
    (folder/'progress.json').write_text(json.dumps({'best':None}))
    rows=api().collect_protected(tmp_path,warm)
    assert rows[0]['candidate']==m['candidate']
    assert rows[0]['origin_sha256']==api().digest_file(folder/'manifest.json')


def test_store_rejects_peak_not_matching_validated_candidate(tmp_path):
    ctx,m=example();m['selected']['peak_MW_m2']=.1
    with pytest.raises(ValueError):api().QualifiedStore(tmp_path,ctx).publish(m,write_artifacts)


def test_finalizer_promotes_new_stricter_winner_and_revalidates_looser_cases(tmp_path,monkeypatch):
    from d3d_optimization import campaign_validation as module
    from d3d_optimization.campaign_geometry import Candidate
    ctx,_=example();q=object.__new__(module.Qualifier);q.source='IMP';q.root=tmp_path
    q.config=dict(fingerprint=ctx['fingerprint'],reference=ctx['reference'],deadline_epoch=float('inf'),backend='cpu')
    cases=[]
    for mm in [100,50,20]:
        folder=tmp_path/f'R{mm:03d}';folder.mkdir()
        cases.append(SimpleNamespace(folder=folder,mm=mm,state={'best':None},e=SimpleNamespace(metadata={'source':'IMP'})))
    winner=Candidate('IMP',100,tuple([0.]*8+[.6]))
    q.shortlist=lambda case,*a,**kw:[{'candidate':winner.record()}] if case.mm==100 else []
    attempted=[]
    def validate(case,candidate):
        attempted.append(case.mm)
        _,m=example(2.,.6);v=m['validation'];v['candidate']=candidate.record();v['identity']=candidate.identity
        return v
    q.validate=validate
    original=SimpleNamespace(baseline={'peak_q':3.},describe=lambda r:{'peak_MW_m2':r['peak_q']})
    q.original=lambda i:original
    monkeypatch.setattr(module,'evaluator',lambda *a:original)
    monkeypatch.setattr(module,'evaluate_candidate',lambda *a:{'feasible':True,'peak_q':2.,'deposition':{'q':np.array([2.])}})
    monkeypatch.setattr(module,'write_arrays',lambda folder,*a:write_artifacts(folder))
    q.finalize(cases)
    assert attempted==[100,50,20]
    for case in cases:
        assert api().QualifiedStore(case.folder,api().context(q.config,'IMP',case.mm)).best()['selected']['peak_MW_m2']==2.


def test_old_validation_cache_cannot_skip_current_policy(tmp_path):
    from d3d_optimization.campaign_validation import Qualifier
    ctx,m=example();c=m['candidate'];identity=m['validation']['identity']
    from d3d_optimization.campaign_geometry import Candidate
    candidate=Candidate.from_record(c)
    folder=tmp_path/'validation'/identity;folder.mkdir(parents=True)
    old=copy.deepcopy(m['validation']);old['qualification_policy']='legacy'
    (folder/'validation.json').write_text(json.dumps(old))
    q=object.__new__(Qualifier);q.config={'fingerprint':ctx['fingerprint'],'deadline_epoch':float('inf')}
    def recompute(i):raise RuntimeError('current qualification requested')
    q.original=recompute
    with pytest.raises(RuntimeError,match='current qualification requested'):
        q.validate(SimpleNamespace(folder=tmp_path,mm=20),candidate)


def test_snapshot_freezes_protected_geometry_and_its_provenance(tmp_path,monkeypatch):
    import run_global_campaign as runner
    root=tmp_path/'repo';root.mkdir();(root/'input.txt').write_text('fixed reference')
    _,m=example();warm=tmp_path/'warm';case=warm/'IMP/R020';case.mkdir(parents=True)
    (case/'manifest.json').write_text(json.dumps(m))
    (warm/'campaign.json').write_text('{}')
    monkeypatch.setattr(runner,'ROOT',root)
    def capture(output,extra):
        output.mkdir();(output/'provenance.json').write_text(json.dumps({'files':{'input.txt':api().digest_file(root/'input.txt')}}))
    monkeypatch.setattr(runner,'capture',capture)
    monkeypatch.setattr(runner,'runtime_identity',lambda backend:{'backend':backend})
    output=tmp_path/'new'
    config=runner.snapshot(output,dict(warm_study=str(warm),hours=1,backend='cpu',validation_reserve_seconds=60,search_seconds=None))
    path=output/'snapshot/campaign_protected_candidates.json'
    rows=json.loads(path.read_text())['candidates']
    assert rows[0]['candidate']==m['candidate']
    frozen=json.loads((output/'frozen_files.json').read_text())
    assert frozen[path.name]==api().digest_file(path)
    assert config['qualification_policy']==api().QUALIFICATION_POLICY
    assert rows[0]['origin_sha256']==api().digest_file(case/'manifest.json')
