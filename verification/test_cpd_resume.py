import json
from types import SimpleNamespace
import numpy as np
import pytest
from run_cpd_campaign import PatternProposals
from run_cpd_resume import remaining_budget, replay_proposals


def evaluator():
    return SimpleNamespace(deformation=SimpleNamespace(dimension=8,bound=.1,
        slices={'a':slice(0,4),'b':slice(4,8)},
        curves={'a':SimpleNamespace(reference_length_m=1.),
                'b':SimpleNamespace(reference_length_m=.2)}))


def test_budget_is_original_minus_consumed_not_new_three_hours():
    assert remaining_budget(dict(requested_seconds=10800,campaign_elapsed_seconds=5674.955197399948)) == pytest.approx(5125.044802600052)
    for spent in (-1,10800,float('nan')):
        with pytest.raises(ValueError):
            remaining_budget(dict(requested_seconds=10800,campaign_elapsed_seconds=spent))


def test_log_replay_restores_rng_adaptation_and_next_proposal():
    e=evaluator();p=PatternProposals(e,{},71);best=np.zeros(8);rows=[]
    for i in range(200):
        candidate=p.next(best);accepted=(i%7==0)
        rows.append(dict(source='IMP',count=160,evaluation=i+1,
                         coefficients=candidate.tolist(),accepted=accepted))
        if accepted:best=candidate
        p.feedback(accepted)
    restored=replay_proposals(e,{},71,'IMP',rows,200,best)
    assert restored.iteration == p.iteration
    assert restored.stagnation == p.stagnation
    assert restored.scale == p.scale
    np.testing.assert_array_equal(restored.next(best),p.next(best))
    rows[22]['coefficients'][0]+=.01
    with pytest.raises(ValueError,match='proposal'):
        replay_proposals(e,{},71,'IMP',rows,200,best)


def test_missing_log_tail_cannot_resume_from_unmatched_checkpoint():
    with pytest.raises(ValueError,match='count'):
        replay_proposals(evaluator(),{},1,'IMP',[],1,np.zeros(8))


def test_atomic_checkpoint_recovers_from_transient_windows_file_lock(tmp_path,monkeypatch):
    import run_cpd_campaign as campaign
    original=campaign.os.replace;attempts=[]
    def locked(a,b):
        attempts.append(1)
        if len(attempts)<3:raise PermissionError('temporary sharing violation')
        return original(a,b)
    monkeypatch.setattr(campaign.os,'replace',locked)
    monkeypatch.setattr(campaign,'sleep',lambda seconds:None,raising=False)
    campaign.atomic_json(tmp_path/'state.json',{'same_physics':True})
    assert json.loads((tmp_path/'state.json').read_text())=={'same_physics':True}
    assert len(attempts)==3


def test_persistent_file_lock_is_not_silently_ignored(tmp_path,monkeypatch):
    import run_cpd_campaign as campaign
    attempts=[]
    def locked(*args):
        attempts.append(1);raise PermissionError('persistent')
    monkeypatch.setattr(campaign.os,'replace',locked)
    monkeypatch.setattr(campaign,'sleep',lambda seconds:None)
    with pytest.raises(PermissionError):campaign.atomic_json(tmp_path/'state.json',{})
    assert len(attempts)==6
