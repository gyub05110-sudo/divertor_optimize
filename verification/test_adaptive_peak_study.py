import json
import pytest
from verification.run_optimality_bounds import digest
from verification.test_selective_peak import kernel


def frozen_study(tmp_path):
    report={'complete':True,'preservation':{'unchanged':True},
            'code_preservation':{'unchanged':True},'reference':{'power_W':1e6},
            'baselines':{s:[{'level':'cubic_320','metadata':{'source':s}}] for s in ('IMP','OMP')},
            'candidates':{name:{'candidate':{'source':source},'verdict':{'resolution_screen_passed':ok}}
                for name,source,ok in [('candidate_IMP','IMP',True),('historical_OMP50','OMP',True),
                                       ('campaign_OMP50','OMP',True),('rejected','IMP',False)]}}
    (tmp_path/'summary.json').write_text(json.dumps(report))
    for name in ('baseline_IMP','baseline_OMP','candidate_IMP','historical_OMP50','campaign_OMP50'):
        p=tmp_path/name/'cubic_320'/'kernel.json';p.parent.mkdir(parents=True)
        p.write_text(json.dumps(kernel()))
    paths=list(tmp_path.rglob('*.json'))
    (tmp_path/'artifact_sha256.json').write_text(json.dumps({p.relative_to(tmp_path).as_posix():digest(p) for p in paths}))
    return tmp_path


def test_loads_only_qualified_same_resolution_kernels_and_shared_pairs(tmp_path):
    from verification.run_adaptive_peak_comparison import load_inputs
    inputs,protected=load_inputs(frozen_study(tmp_path),'cubic_320')
    assert set(inputs['models'])=={'baseline_IMP','baseline_OMP','candidate_IMP','historical_OMP50','campaign_OMP50'}
    assert inputs['problem_ids']['candidate_IMP']==inputs['problem_ids']['baseline_IMP']
    assert inputs['problem_ids']['baseline_OMP']!=inputs['problem_ids']['baseline_IMP']
    assert len(inputs['pairs'])==4
    assert ('historical_OMP50','campaign_OMP50') in inputs['pairs']
    assert str(tmp_path/'artifact_sha256.json') in protected


def test_rejects_tampered_kernel_before_any_work(tmp_path):
    from verification.run_adaptive_peak_comparison import load_inputs
    frozen_study(tmp_path)
    (tmp_path/'candidate_IMP/cubic_320/kernel.json').write_text('{}')
    with pytest.raises(ValueError,match='hash'):load_inputs(tmp_path,'cubic_320')


def test_rejects_failed_qualification_even_with_fresh_manifest(tmp_path):
    from verification.run_adaptive_peak_comparison import load_inputs
    frozen_study(tmp_path);p=tmp_path/'summary.json'
    report=json.loads(p.read_text());report['failure']={'message':'rejected'}
    p.write_text(json.dumps(report));index=tmp_path/'artifact_sha256.json'
    recorded=json.loads(index.read_text());recorded['summary.json']=digest(p)
    index.write_text(json.dumps(recorded))
    with pytest.raises(ValueError,match='requalification'):load_inputs(tmp_path,'cubic_320')
