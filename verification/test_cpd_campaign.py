"""Deadline and scientific-acceptance contracts without launching a campaign."""
import importlib
import json
from types import SimpleNamespace

import numpy as np
import pytest


def campaign():
    return importlib.import_module('run_cpd_campaign')


def description(peak=10., area=2., unhit=0., feasible=True, hit=.999):
    return dict(feasible=feasible, errors=[], peak_MW_m2=peak, area_90_m2=area,
                divertor_first_hit_power_fraction=hit,
                ledger=dict(prescribed=1., deposited=1-unhit, unhit=unhit,
                            unresolved=unhit, balance_error=0.))


def test_candidate_gates_cannot_reward_missing_power_or_narrower_footprints():
    baseline=description()
    assert campaign().candidate_reasons(description(peak=9.),baseline) == []
    assert 'unhit_power_increased' in campaign().candidate_reasons(description(peak=1.,unhit=1e-6),baseline)
    assert 'footprint_area_decreased' in campaign().candidate_reasons(description(peak=1.,area=1.),baseline)
    assert 'evaluation_infeasible' in campaign().candidate_reasons(description(peak=1.,feasible=False),baseline)
    bad=description(peak=1.);bad['ledger']['deposited']=.9
    assert 'power_balance_failed' in campaign().candidate_reasons(bad,baseline)


def test_qualification_requires_gain_larger_than_resolution_uncertainty():
    rows={320:dict(baseline=description(10.),candidate=description(9.5)),
          640:dict(baseline=description(10.),candidate=description(9.99))}
    verdict=campaign().qualification_verdict(rows)
    assert not verdict['qualified_improvement']
    assert 'gain_not_above_numerical_uncertainty' in verdict['reasons']
    rows[640]['candidate']=description(9.49)
    assert campaign().qualification_verdict(rows)['qualified_improvement']
    assert not campaign().qualification_verdict({320:rows[320]})['qualified_improvement']


class FakeClock:
    def __init__(self): self.now=0.
    def __call__(self): return self.now
    def advance(self,seconds): self.now+=seconds


class FakeFactory:
    def __init__(self,clock,*,bad_baseline=False,bad_candidates=False):
        self.clock=clock;self.calls=[];self.bad_baseline=bad_baseline
        self.bad_candidates=bad_candidates
    def __call__(self,source,count):
        self.calls.append((source,count,self.clock()))
        self.clock.advance(.02)
        return FakeEvaluator(self,source,count)


class FakeEvaluator:
    def __init__(self,factory,source,count):
        self.factory=factory;self.count=count;self.source_name=source
        self.metadata=dict(equilibrium_sha256='fixed-field',geometry_sha256='fixed-geometry',
                           source=source,count=count,preparation_seconds=.02)
        self.deformation=SimpleNamespace(dimension=8,bound=.1,
            slices={'face_a':slice(0,4),'face_b':slice(4,8)},
            curves={'face_a':SimpleNamespace(reference_length_m=1.),
                    'face_b':SimpleNamespace(reference_length_m=.5)})
    def evaluate(self,c=None):
        self.factory.clock.advance(.04)
        baseline=c is None
        peak=10. if baseline else 10.-10*float(np.sum(c))
        desc=description(peak=peak,unhit=1e-4 if baseline and self.factory.bad_baseline else 0.,
                         feasible=not (not baseline and self.factory.bad_candidates))
        return dict(deposition={},description=desc,coefficients=c)
    def describe(self,result): return result['description'].copy()
    def export(self,path,result): np.savez_compressed(path,coefficients=np.asarray([] if result['coefficients'] is None else result['coefficients']))


def test_bad_baseline_never_starts_campaign_clock_or_search(tmp_path):
    clock=FakeClock();factory=FakeFactory(clock,bad_baseline=True)
    report=campaign().CampaignRunner(factory,tmp_path,seconds=2.,reserve=.8,clock=clock).run()
    assert report['status'] == 'baseline_rejected'
    assert report['campaign_started'] is False
    assert all(count==160 for _,count,_ in factory.calls)
    assert json.loads((tmp_path/'study.json').read_text())['status']=='baseline_rejected'


def test_deadline_reserves_qualification_and_records_fixed_field(tmp_path):
    clock=FakeClock();factory=FakeFactory(clock)
    report=campaign().CampaignRunner(factory,tmp_path,seconds=2.,reserve=.8,clock=clock,seed=4).run()
    assert report['campaign_started'] is True
    assert report['campaign_elapsed_seconds'] <= 2.
    assert all(s in report['sources'] for s in ('IMP','OMP'))
    assert any(count==320 for _,count,_ in factory.calls)
    assert any(count==640 for _,count,_ in factory.calls)
    assert report['equilibrium_sha256']=='fixed-field'
    assert report['search_evaluations'] > 0
    assert report['status'] == 'completed'
    for source in report['sources'].values():
        assert source['baseline']['peak_MW_m2'] == 10.
        assert source['best_search']['peak_MW_m2'] <= 10.
        assert len(source['best_coefficients'])==8
        assert source['qualification']['completed']
    assert json.loads((tmp_path/'status.json').read_text())['status']=='completed'


def test_invalid_candidates_never_replace_best_and_short_budget_reports_unqualified(tmp_path):
    clock=FakeClock();factory=FakeFactory(clock,bad_candidates=True)
    report=campaign().CampaignRunner(factory,tmp_path,seconds=.15,reserve=.04,clock=clock).run()
    assert report['status']=='time_expired_without_qualification'
    assert report['campaign_elapsed_seconds'] <= .15
    for source in report['sources'].values():
        assert source['best_search']['peak_MW_m2']==10.
        assert np.count_nonzero(source['best_coefficients'])==0
        assert not source['qualification']['qualified_improvement']


def test_seed_reproduces_search_candidates_with_injected_clock(tmp_path):
    reports=[]
    for i in range(2):
        clock=FakeClock();factory=FakeFactory(clock)
        reports.append(campaign().CampaignRunner(factory,tmp_path/str(i),seconds=2.,reserve=.8,clock=clock,seed=19).run())
    for source in ('IMP','OMP'):
        assert reports[0]['sources'][source]['best_coefficients']==reports[1]['sources'][source]['best_coefficients']
        assert reports[0]['sources'][source]['evaluations']==reports[1]['sources'][source]['evaluations']


def test_failed_candidate_export_does_not_record_an_accepted_update(tmp_path,monkeypatch):
    original=FakeEvaluator.export
    def fail_candidate(self,path,result):
        if result['coefficients'] is not None and np.any(result['coefficients']):
            raise OSError('injected checkpoint write failure')
        return original(self,path,result)
    monkeypatch.setattr(FakeEvaluator,'export',fail_candidate)
    clock=FakeClock();factory=FakeFactory(clock)
    report=campaign().CampaignRunner(factory,tmp_path,seconds=.3,reserve=.08,clock=clock).run()
    rows=[json.loads(line) for line in (tmp_path/'candidates.jsonl').read_text().splitlines()]
    assert rows
    assert all(not row.get('accepted',False) for row in rows)
    assert all(data['accepted_updates']==0 for data in report['sources'].values())


def test_unexpected_nonpreemptible_overrun_is_explicit_and_checkpoint_survives(tmp_path,monkeypatch):
    original=FakeEvaluator.evaluate
    def slower_candidate(self,c=None):
        if c is not None: self.factory.clock.advance(.5)
        return original(self,c)
    monkeypatch.setattr(FakeEvaluator,'evaluate',slower_candidate)
    clock=FakeClock();factory=FakeFactory(clock)
    report=campaign().CampaignRunner(factory,tmp_path,seconds=.3,reserve=.08,clock=clock).run()
    assert report['status']=='time_expired_without_qualification'
    assert report['deadline_overrun_seconds'] > 0
    assert json.loads((tmp_path/'study.json').read_text())['deadline_overrun_seconds'] > 0
    assert all(not data['qualification']['qualified_improvement'] for data in report['sources'].values())


def test_accepted_generations_keep_matching_immutable_npz_files(tmp_path):
    clock=FakeClock();factory=FakeFactory(clock)
    report=campaign().CampaignRunner(factory,tmp_path,seconds=2.,reserve=.8,clock=clock).run()
    rows=[json.loads(line) for line in (tmp_path/'candidates.jsonl').read_text().splitlines()]
    accepted=[row for row in rows if row.get('accepted')]
    assert len(accepted)>=2
    assert len({row['artifact'] for row in accepted}) == len(accepted)
    for row in accepted:
        with np.load(row['artifact']) as arrays:
            np.testing.assert_array_equal(arrays['coefficients'],row['coefficients'])
    for source,data in report['sources'].items():
        with np.load(data['best_search_artifact']) as arrays:
            np.testing.assert_array_equal(arrays['coefficients'],data['best_coefficients'])


@pytest.mark.parametrize('key,replacement',[
    ('reference',{'power_W':2e6}),('source_RZ_m',[3.5,0.]),('source_support_m',.2)])
def test_fixed_physics_and_source_metadata_reject_resolution_drift(tmp_path,key,replacement):
    runner=campaign().CampaignRunner(lambda *args:None,tmp_path)
    metadata=dict(equilibrium_sha256='field',geometry_sha256='geometry',source='IMP',
                  reference={'power_W':1e6},source_RZ_m=[3.,0.],source_support_m=.1)
    runner._check_metadata(SimpleNamespace(metadata=metadata))
    changed=dict(metadata);changed[key]=replacement
    with pytest.raises(ValueError,match='Fixed'):
        runner._check_metadata(SimpleNamespace(metadata=changed))


def test_qualification_refines_all_numerical_transport_settings():
    medium=campaign().numerical_settings(320)
    fine=campaign().numerical_settings(640)
    for key in ('path_step','wall_step','arc_step'):
        assert medium[key]==.001
        assert fine[key]==.0005
    assert medium['quadrature']==128
    assert fine['quadrature']==256
