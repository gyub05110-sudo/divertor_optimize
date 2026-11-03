import copy
import json
from pathlib import Path
import pytest


@pytest.mark.parametrize('update',[
    {'complete':False}, {'failure':{'message':'failed'}},
    {'preservation':{'unchanged':False}}, {'code_preservation':{'unchanged':False}},
])
def test_rejects_incomplete_or_unpreserved_requalification(update):
    from verification.run_selective_peak_comparison import require_accepted_report
    report={'complete':True,'preservation':{'unchanged':True},'code_preservation':{'unchanged':True}}
    report.update(update)
    with pytest.raises(ValueError):require_accepted_report(report)


def test_accepts_completed_preserved_requalification():
    from verification.run_selective_peak_comparison import require_accepted_report
    require_accepted_report({'complete':True,'preservation':{'unchanged':True},'code_preservation':{'unchanged':True}})


def test_problem_identity_separates_source_resolution_and_reference():
    from verification.run_selective_peak_comparison import problem_identity
    r={'reference':{'power_W':1e6},'baselines':{'IMP':[
        {'level':'cubic_320','metadata':{'shot_sha256':'abc'}},
        {'level':'cubic_1280','metadata':{'shot_sha256':'abc'}}],
        'OMP':[{'level':'cubic_320','metadata':{'shot_sha256':'abc'}}]}}
    a=problem_identity(r,'IMP','cubic_320')
    assert a!=problem_identity(r,'IMP','cubic_1280')
    assert a!=problem_identity(r,'OMP','cubic_320')
    b=copy.deepcopy(r);b['reference']['power_W']=2e6
    assert a!=problem_identity(b,'IMP','cubic_320')
    with pytest.raises(ValueError):problem_identity(r,'IMP','absent')


@pytest.fixture
def audited_study(tmp_path):
    """An immutable study with one guarded import and one historical supplement."""
    from verification.run_optimality_bounds import digest
    study=tmp_path/'study';source=study/'source';checkout=tmp_path/'checkout'
    allowed=['d3d_optimization/continuous_path_bounds.py',
             'd3d_optimization/continuous_transport_certificate.py']
    consumed=['verification/run_candidate_requalification.py','d3d_optimization/transport.py']
    hashes={}
    for name in allowed+consumed:
        for root in (source,checkout):
            path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text('# preserved '+name+'\n')
        hashes[name]=digest(source/name)
    origin=tmp_path/'historical'/'transport.py';origin.parent.mkdir()
    origin.write_bytes((source/consumed[1]).read_bytes())
    guarded={str(checkout/name):hashes[name] for name in allowed+[consumed[0]]}
    after=dict(guarded)
    for name in allowed:
        (checkout/name).write_text('# concurrently edited unused module\n')
        after[str(checkout/name)]=digest(checkout/name)
    report={'complete':True,'failure':None,'code_sha256':guarded,
            'code_preservation':{'unchanged':False,'changed':[str(checkout/name) for name in allowed],
                                 'after_sha256':after},
            'protected_sha256':{str(origin):hashes[consumed[1]]},
            'preservation':{'unchanged':True,'after_sha256':{str(origin):hashes[consumed[1]]}}}
    supplement={'dependencies':[{'relative_path':consumed[1],'protected_origin':str(origin),
                                  'working_copy_sha256':hashes[consumed[1]],
                                  'snapshot_sha256':hashes[consumed[1]]}]}
    (study/'supplemental_snapshot_dependencies.json').write_text(json.dumps(supplement))
    evidence=[{'path':consumed[0],'kind':'original_code_guard',
               'before_path':str(checkout/consumed[0]),'after_path':str(checkout/consumed[0])},
              {'path':consumed[1],'kind':'guarded_campaign_snapshot_supplement',
               'before_path':str(origin),'after_path':str(origin),
               'supplemental_snapshot_sha256':hashes[consumed[1]]}]
    audit={'schema':'candidate-requalification-code-scope-audit.v1','accepted':True,
           'historical_preservation_unchanged':True,'static_changed_module_references':[],
           'changed_code_paths':report['code_preservation']['changed'],
           'allowed_unused_changed_names':allowed,'snapshot_code_sha256':hashes,
           'consumed_snapshot_sha256':{name:hashes[name] for name in consumed},
           'consumed_code_sha256':{str(checkout/name):hashes[name] for name in consumed},
           'consumed_code_after_sha256':{str(checkout/name):hashes[name] for name in consumed},
           'source_evidence':evidence,
           'supplemental_dependencies_sha256':digest(study/'supplemental_snapshot_dependencies.json'),
           'import_probe':{'exit_code':0,'cwd':str(source),
                           'imported_local_paths':[str(source/name) for name in consumed]}}
    audit_path=study/'code_scope_audit.json'
    def persist():
        (study/'summary.json').write_text(json.dumps(report))
        (study/'artifact_sha256.json').write_text(json.dumps({'summary.json':digest(study/'summary.json')}))
        audit['summary_sha256']=digest(study/'summary.json')
        audit['artifact_index_sha256']=digest(study/'artifact_sha256.json')
        audit_path.write_text(json.dumps(audit))
    persist()
    return study,checkout,report,audit,audit_path,persist


def test_explicit_unused_code_exception_accepts_immutable_evidence_after_checkout_drift(audited_study):
    from verification.run_selective_peak_comparison import require_accepted_report,verify_code_scope
    study,checkout,report,audit,audit_path,_=audited_study
    with pytest.raises(ValueError):require_accepted_report(report)
    # Later edits to the mutable checkout do not invalidate the completed run.
    (checkout/'verification/run_candidate_requalification.py').write_text('# later work\n')
    checked=verify_code_scope(study,report,audit_path)
    assert str(audit_path.resolve()) in checked
    assert str(study/'source/verification/run_candidate_requalification.py') in checked
    assert str(checkout/'verification/run_candidate_requalification.py') not in checked


@pytest.mark.parametrize('filename',[
    'summary.json','artifact_sha256.json','supplemental_snapshot_dependencies.json',
    'source/verification/run_candidate_requalification.py',
    'source/d3d_optimization/continuous_path_bounds.py',
])
def test_audited_exception_rejects_modified_immutable_files(audited_study,filename):
    from verification.run_selective_peak_comparison import verify_code_scope
    study,_,report,_,audit_path,_=audited_study
    with (study/filename).open('a') as stream:stream.write(' ')
    with pytest.raises(ValueError,match='hash mismatch'):verify_code_scope(study,report,audit_path)


@pytest.mark.parametrize('condition',[
    'failed_report','incomplete_report','changed_historical_input','historical_file_tamper',
    'audit_denied','wrong_schema','failed_import','wrong_import_root','import_outside_snapshot','static_reference',
    'unexpected_changed_module','consumed_changed_module','consumed_before_after_mismatch',
    'guarded_before_after_mismatch','missing_source_evidence','snapshot_path_escape',
])
def test_audited_exception_rejects_invalid_provenance(audited_study,condition):
    from verification.run_selective_peak_comparison import verify_code_scope
    study,checkout,report,audit,audit_path,persist=audited_study
    imported='verification/run_candidate_requalification.py'
    if condition=='failed_report':report['failure']={'message':'numerical run failed'}
    elif condition=='incomplete_report':report['complete']=False
    elif condition=='changed_historical_input':report['preservation']['unchanged']=False
    elif condition=='historical_file_tamper':Path(next(iter(report['protected_sha256']))).write_text('changed')
    elif condition=='audit_denied':audit['accepted']=False
    elif condition=='wrong_schema':audit['schema']='other'
    elif condition=='failed_import':audit['import_probe']['exit_code']=1
    elif condition=='wrong_import_root':audit['import_probe']['cwd']=str(checkout)
    elif condition=='import_outside_snapshot':audit['import_probe']['imported_local_paths'][0]=str(checkout/imported)
    elif condition=='static_reference':audit['static_changed_module_references']=[{'reference':'continuous_path_bounds'}]
    elif condition=='unexpected_changed_module':report['code_preservation']['changed'][0]=str(checkout/'other.py')
    elif condition=='consumed_changed_module':
        name=audit['allowed_unused_changed_names'][0]
        audit['consumed_snapshot_sha256'][name]=audit['snapshot_code_sha256'][name]
    elif condition=='consumed_before_after_mismatch':audit['consumed_code_after_sha256'][str(checkout/imported)]='bad'
    elif condition=='guarded_before_after_mismatch':report['code_preservation']['after_sha256'][str(checkout/imported)]='bad'
    elif condition=='missing_source_evidence':audit['source_evidence'].pop()
    elif condition=='snapshot_path_escape':audit['snapshot_code_sha256']['../outside.py']='bad'
    persist()
    with pytest.raises(ValueError):verify_code_scope(study,report,audit_path)


def test_audit_cannot_omit_guarded_snapshot_to_hide_tampering(audited_study):
    from verification.run_selective_peak_comparison import verify_code_scope
    study,_,report,audit,audit_path,persist=audited_study
    name='d3d_optimization/continuous_path_bounds.py'
    (study/'source'/name).write_text('# altered immutable source\n')
    del audit['snapshot_code_sha256'][name]
    persist()
    with pytest.raises(ValueError):verify_code_scope(study,report,audit_path)


def test_audit_import_inventory_cannot_hide_consumption_of_changed_module(audited_study):
    from verification.run_selective_peak_comparison import verify_code_scope
    study,_,report,audit,audit_path,persist=audited_study
    audit['import_probe']['imported_local_paths'].append(str(study/'source/d3d_optimization/continuous_path_bounds.py'))
    persist()
    with pytest.raises(ValueError):verify_code_scope(study,report,audit_path)
