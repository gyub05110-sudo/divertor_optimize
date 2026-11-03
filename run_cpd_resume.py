"""Resume an immutable CPD checkpoint using only its unused computation budget."""
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import traceback
from time import perf_counter
import numpy as np

from run_cpd_campaign import (CampaignRunner, PatternProposals, SOURCES, SEARCH_COUNT,
    numerical_settings, baseline_reasons, candidate_reasons, qualification_verdict,
    atomic_json)


def remaining_budget(report):
    total=float(report['requested_seconds']);spent=float(report['campaign_elapsed_seconds'])
    if not math.isfinite(total) or not math.isfinite(spent) or not 0 <= spent < total:
        raise ValueError('Checkpoint must have finite positive unused budget')
    return total-spent


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replay_proposals(evaluator, baseline, seed, source, rows, evaluations, expected_best,
                     probabilities=None):
    """Reconstruct the exact RNG and adaptive scale, without re-running physics."""
    proposal=PatternProposals(evaluator,baseline,seed)
    if probabilities is not None:
        if not np.allclose(proposal.face_probabilities,probabilities,rtol=0,atol=1e-12):
            raise ValueError('Baseline proposal probabilities changed')
        proposal.face_probabilities=np.asarray(probabilities,float)
    selected=[r for r in rows if r['source']==source and r.get('stage')!='qualification']
    if len(selected)!=evaluations:
        raise ValueError('Candidate log count does not match checkpoint')
    best=np.zeros(proposal.dimension)
    for i,row in enumerate(selected,1):
        if row['evaluation']!=i or row['count']!=SEARCH_COUNT:
            raise ValueError('Candidate log sequence changed')
        candidate=proposal.next(best)
        if not np.array_equal(candidate,np.asarray(row['coefficients'])):
            raise ValueError(f'Logged proposal differs during replay: {source} {i}')
        if row['accepted']:best=candidate
        proposal.feedback(row['accepted'])
    if not np.array_equal(best,np.asarray(expected_best)):
        raise ValueError('Logged best differs from checkpoint coefficients')
    return proposal


def load_history(study_path, visited=None):
    study_path=Path(study_path).resolve();visited=set() if visited is None else visited
    if study_path in visited:raise ValueError('Checkpoint ancestry cycle')
    visited.add(study_path)
    report=json.loads(study_path.read_text(encoding='utf-8'))
    rows=[]
    if 'resume' in report:
        parent=Path(report['resume']['parent_study'])
        if sha256(parent)!=report['resume']['parent_study_sha256']:
            raise ValueError('Parent checkpoint changed')
        if sha256(parent.parent/'candidates.jsonl')!=report['resume']['parent_candidate_log_sha256']:
            raise ValueError('Parent candidate log changed')
        rows.extend(load_history(parent,visited)[1])
    logfile=study_path.parent/'candidates.jsonl'
    if logfile.exists():
        rows.extend(json.loads(line) for line in logfile.read_text(encoding='utf-8').splitlines() if line)
    return report,rows


class ResumeRunner(CampaignRunner):
    def __init__(self,factory,output,parent_study,**kwargs):
        self._resume_begin=kwargs.get('clock',perf_counter)()
        self._resume_begin_utc=datetime.now(timezone.utc).isoformat()
        self.parent_path=Path(parent_study).resolve()
        self.parent,self.history=load_history(self.parent_path)
        if sum(row.get('stage')!='qualification' for row in self.history)!=self.parent['search_evaluations']:
            raise ValueError('Candidate log count differs from total checkpoint evaluations')
        seconds=remaining_budget(self.parent)
        if seconds<=900:raise ValueError('Insufficient remaining budget for search plus 900-second qualification reserve')
        super().__init__(factory,output,seconds=seconds,reserve=900,
                         seed=self.parent['seed'],**kwargs)
        self.prior_consumed=float(self.parent.get('prior_consumed_seconds',0))+float(self.parent['campaign_elapsed_seconds'])
        self.report.update(resume=dict(parent_study=str(self.parent_path),
            parent_study_sha256=sha256(self.parent_path),
            parent_candidate_log_sha256=sha256(self.parent_path.parent/'candidates.jsonl'),
            reconstruction='Replay original seeded proposals and acceptance decisions; no physics rerun during log replay'),
            prior_consumed_seconds=self.prior_consumed,
            original_budget_seconds=self.parent.get('original_budget_seconds',self.parent['requested_seconds']))

    def _checkpoint(self):
        self.report['total_consumed_seconds']=self.prior_consumed+self._elapsed()
        for source,proposal in self.proposals.items():
            self.report['sources'][source]['proposal_state']=dict(
                rng_state=deepcopy(proposal.rng.bit_generator.state),iteration=proposal.iteration,
                stagnation=proposal.stagnation,scale=proposal.scale)
        super()._checkpoint()

    def _prepare(self):
        # Recovery, cache reconstruction and best replay count against remaining time.
        self.started=self._resume_begin
        self.report.update(campaign_started=True,campaign_started_utc=self._resume_begin_utc,
                           status='recovering_checkpoint',fixed_inputs=deepcopy(self.parent['fixed_inputs']),
                           equilibrium_sha256=self.parent['equilibrium_sha256'],
                           geometry_sha256=self.parent['geometry_sha256'])
        self._checkpoint()
        for index,source in enumerate(SOURCES):
            old=self.parent['sources'][source]
            if old.get('qualification_rows'):
                raise ValueError('Search continuation requires a checkpoint before final qualification')
            with np.load(old['best_search_artifact'],allow_pickle=False) as data:
                np.testing.assert_array_equal(data['coefficients'],old['best_coefficients'])
                if str(data['equilibrium_sha256'])!=self.parent['equilibrium_sha256'] or str(data['geometry_sha256'])!=self.parent['geometry_sha256']:
                    raise ValueError('Saved artifact input hash mismatch')
            evaluator=self._timed(source,160,'prepare',lambda:self.factory(source,160))
            self._check_metadata(evaluator)
            baseline=self._timed(source,160,'evaluate',lambda:evaluator.evaluate(None))
            desc=evaluator.describe(baseline)
            if baseline_reasons(desc):raise ValueError('Restored baseline failed acceptance')
            for key in ('peak_MW_m2','area_90_m2'):
                if not np.isclose(desc[key],old['baseline'][key],rtol=1e-10,atol=1e-12):
                    raise ValueError('Restored baseline drift: '+key)
            proposal=replay_proposals(evaluator,baseline,self.seed+index,source,self.history,
                old['evaluations'],old['best_coefficients'],old['proposal_face_probabilities'])
            best=self._timed(source,160,'evaluate',lambda:evaluator.evaluate(np.asarray(old['best_coefficients'])))
            best_desc=evaluator.describe(best)
            if candidate_reasons(best_desc,desc):raise ValueError('Restored best failed acceptance')
            for key in ('peak_MW_m2','area_90_m2'):
                if not np.isclose(best_desc[key],old['best_search'][key],rtol=1e-10,atol=1e-12):
                    raise ValueError('Restored best drift: '+key)
            data=deepcopy(old)
            data.update(metadata=deepcopy(evaluator.metadata),qualification=qualification_verdict({}),qualification_rows={},
                baseline_artifact=self._export(evaluator,f'baseline-{source}-160.npz',baseline),
                best_search_artifact=self._export(evaluator,f'restored-{source}-best.npz',best),
                restored_best_replay=best_desc)
            self.report['sources'][source]=data
            self.evaluators[source]=evaluator;self.proposals[source]=proposal
            self._checkpoint()
        self.report['search_evaluations']=self.parent['search_evaluations']
        self.report['recovery_seconds']=self._elapsed()
        return True

    def run(self):
        # Explicit run shell avoids CampaignRunner.run resetting the recovery clock.
        self.output.mkdir(parents=True,exist_ok=True)
        if (self.output/'study.json').exists():raise FileExistsError('Choose a new resumed study directory')
        self._checkpoint()
        try:
            self._prepare()
            self._search();self._qualify()
            for data in self.report['sources'].values():
                data['qualification']=qualification_verdict(data['qualification_rows'])
            complete=all(v['qualification']['completed'] for v in self.report['sources'].values())
            self.report['status']='completed' if complete else 'qualification_incomplete'
            self.report['deadline_overrun_seconds']=max(0.,self._elapsed()-self.seconds)
            if self.report['deadline_overrun_seconds']>0:self.report['status']='deadline_overrun'
        except KeyboardInterrupt:
            self.report['status']='interrupted'
        except Exception as exc:
            self.report['status']='resume_failed'
            self.report['failure']=f'{type(exc).__name__}: {exc}'
            self.report['traceback']=traceback.format_exc()
            traceback.print_exc()
        finally:self._checkpoint()
        return self.report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent-study',required=True);p.add_argument('--output',required=True)
    p.add_argument('--nodes',required=True);p.add_argument('--equilibrium',required=True)
    p.add_argument('--cache-dir',required=True);p.add_argument('--log-dir')
    a=p.parse_args()
    if a.log_dir:
        logs=Path(a.log_dir);logs.mkdir(parents=True,exist_ok=True)
        sys.stdout=(logs/'stdout.log').open('a',encoding='utf-8',buffering=1)
        sys.stderr=(logs/'stderr.log').open('a',encoding='utf-8',buffering=1)
        import faulthandler
        faulthandler.enable(file=sys.stderr,all_threads=True)
    from cpd_optimization.evaluator import CPDEvaluator
    def factory(source,count):
        return CPDEvaluator(a.nodes,a.equilibrium,source,count=count,**numerical_settings(count),
                            max_path=30.,cache_dir=a.cache_dir)
    runner=ResumeRunner(factory,a.output,a.parent_study,progress=True)
    report=runner.run()
    print(json.dumps(dict(status=report['status'],remaining_budget_seconds=runner.seconds,
                         total_consumed_seconds=report['total_consumed_seconds'])),flush=True)
    return 0 if report['status']=='completed' else 2


if __name__=='__main__':raise SystemExit(main())
