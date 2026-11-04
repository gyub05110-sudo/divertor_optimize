"""Qualify archived OMP milestones without optimizing or changing historical data.

Scores every eligible milestone at 1280, then fully audits in increasing fine-peak
order until the first passing improvement per radius. Restart requires the exact
same reference, runtime, source and checksummed evidence.
"""
import argparse,hashlib,json,math,msvcrt,os,shutil,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from qualification_selection import ranked_improvements,passes_improvement

HERE=Path(__file__).resolve().parent
REPO=Path('D:/DivertorOptimization/worktrees/d3d-report-correction-20261104')
CAMPAIGN=Path('D:/DivertorOptimization/studies/imp-omp-global-20261029')
sys.path.insert(0,str(REPO))
sys.path.append(str(CAMPAIGN/'snapshot/.local_dependencies/d3d-fixed-equilibrium'))
from d3d_optimization.boundary_evaluator import prepare_boundary,BoundaryEvaluator,export_boundary
from d3d_optimization.campaign_geometry import Candidate,evaluator,evaluate_candidate,GEOMETRY_POLICY
from d3d_optimization.campaign_validation import geometric_audit
from d3d_optimization.campaign_incumbents import QualifiedStore,context,QUALIFICATION_POLICY
from d3d_optimization.reference_transport import ReferencePoint
from verification.validate_boundary_transport import LEVELS,audit,verdict
from verification.run_candidate_requalification import consumed_code_paths,check_preservation
from run_global_campaign import runtime_identity


def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def fingerprint(x):return hashlib.sha256(json.dumps(x,sort_keys=True,allow_nan=False).encode()).hexdigest()
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_suffix(p.suffix+'.tmp')
    temp.write_text(json.dumps(x,ensure_ascii=False,allow_nan=False,indent=2),encoding='utf-8')
    temp.replace(p)
def utc():return datetime.now(timezone.utc).isoformat()
def checked_save(p,payload):save(p,dict(schema='checksummed-record.v1',sha256=fingerprint(payload),payload=payload))
def checked_read(p):
    record=read(p)
    if record.get('schema')!='checksummed-record.v1' or record.get('sha256')!=fingerprint(record.get('payload')):
        raise ValueError('cached JSON checksum mismatch: '+str(p))
    return record['payload']


class Run:
    def __init__(self,out):
        self.out=out;self.originals={};self.evaluators={};self.base_audits={};self.wall_original=None
        data=read(HERE/'dist/data.json');campaign=read(CAMPAIGN/'campaign.json')
        if data['reference']!=campaign['reference']:raise ValueError('viewer/campaign reference mismatch')
        self.root=Path(campaign['frozen_root']);self.reference=ReferencePoint(**campaign['reference'])
        self.items={}
        for series in data['series']:
            if series['run']!='campaign' or series['source']!='OMP' or series['rmin'] not in (20,50):continue
            for i,row in enumerate(series['stages'],1):
                c=Candidate.from_record(row['candidate'])
                if c.identity!=row['id'] or c.source!='OMP' or c.radius_mm!=series['rmin']:raise ValueError('candidate identity/case mismatch')
                if c.identity in self.items:raise ValueError('duplicate candidate identity')
                self.items[c.identity]=dict(candidate=c.record(),identity=c.identity,stage=i,series=series['id'],
                    coarse_peak=row['peak'],wall=row['wall'],origin=row['origin'],line=row['line'])
        counts={mm:sum(r['candidate']['radius_mm']==mm for r in self.items.values()) for mm in (20,50)}
        if counts!={20:19,50:34}:raise ValueError('unexpected finite candidate catalog '+str(counts))
        protected={str(HERE/'dist'/name):sha(HERE/'dist'/name) for name in
                   ('data.json','legs.json','stage-profiles.json','full-saved-profiles.json')}
        for p in (HERE/'dist/profiles').glob('*.json'):protected[str(p)]=sha(p)
        for rel,h in read(CAMPAIGN/'frozen_files.json').items():
            path=self.root/rel
            if sha(path)!=h:raise ValueError('frozen input mismatch '+str(path))
            protected[str(path)]=h
        for path in (CAMPAIGN/'campaign.json',CAMPAIGN/'frozen_files.json',CAMPAIGN/'final_verification.json'):
            protected[str(path)]=sha(path)
        # Cover all previously published campaign manifests/arrays and the current
        # geometry requalification, without copying or changing their directories.
        prior=read(CAMPAIGN/'final_verification.json')
        for key in ('campaign_artifact_sha256','preserved_prior_manifest_sha256'):
            for p,h in prior[key].items():
                if sha(p)!=h:raise ValueError('protected campaign artifact mismatch '+p)
                protected[p]=h
        for folder in [CAMPAIGN/'OMP',Path('D:/DivertorOptimization/studies/candidate-requalification-20261103')]:
            for p in folder.rglob('*'):
                if p.is_file() and (p.name in ('manifest.json','validation.json') or p.name.startswith('selected_refined')):
                    protected[str(p)]=sha(p)
        for p in {row['origin'] for row in self.items.values()}:protected[p]=sha(p)
        runtime=runtime_identity('gpu')
        code={str(p):sha(p) for p in consumed_code_paths(REPO)}
        for p in [Path(__file__),HERE/'qualification_selection.py']:code[str(p)]=sha(p)
        config=dict(reference=campaign['reference'],runtime=runtime,backend='gpu',
            geometry_policy=GEOMETRY_POLICY,qualification_policy=QUALIFICATION_POLICY,
            data_sha256=sha(HERE/'dist/data.json'),candidate_catalog_sha256=fingerprint(self.items),
            code_sha256=code,protected_sha256=protected,selection_policy='all-1280-ascending-until-strict-pass.v1')
        config['fingerprint']=fingerprint(config);self.config=config;self.fp=config['fingerprint']
        self.report=dict(schema='archived-omp-qualification.v1',created_utc=utc(),config=config,
            candidates=self.items,counts=counts,ranking={},attempts={},results={},timings=[],complete=False,
            new_optimizations=0,global_optimum_proven=False,
            scope='Best passing fixed archived milestone by cubic1280 peak per radius; no continuous/global proof.')
        existing=out/'summary.json'
        if existing.exists():
            old=read(existing)
            if old['config']!=config:raise ValueError('resume reference/runtime/code/input mismatch')
            self.report=old;self.report.pop('failure',None)
        else:
            for p,h in code.items():
                src=Path(p);relative=src.relative_to(REPO) if src.is_relative_to(REPO) else Path('driver')/src.name
                target=out/'source'/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,target)
                if sha(target)!=h:raise ValueError('source changed while copying')
            save(out/'candidate_catalog.json',self.items)
        self.data=data
        if not self.report.get('complete'):self.persist('initialized')

    def persist(self,phase,**extra):
        self.report['updated_utc']=utc();self.report['phase']=phase;self.report.update(extra)
        save(self.out/'summary.json',self.report)

    def original(self,index):
        if index not in self.originals:
            method,count,step,nodes,wall_step=LEVELS[index]
            self.persist('preparing_paths',resolution=f'{method}_{count}');print('PREPARE',method,count,flush=True)
            started=time.perf_counter()
            e=prepare_boundary(self.root,'OMP',self.reference,count=count,step=step,nodes=nodes,
                               wall_step=wall_step,method=method)
            self.originals[index]=e
            peak=e.describe(e.baseline)['peak_MW_m2']
            if index in (0,2):
                expected=self.data['baseline']['coarse' if index==0 else 'refined']['OMP']
                if not math.isclose(peak,expected,rel_tol=1e-7):raise ValueError('fixed baseline mismatch')
            folder=self.out/'baseline'/f'{method}_{count}';folder.mkdir(parents=True,exist_ok=True)
            export_boundary(folder/'baseline.npz',e,e.baseline)
            self.report['timings'].append(dict(phase='preparation',level=f'{method}_{count}',seconds=time.perf_counter()-started))
            print('READY',method,count,'seconds',time.perf_counter()-started,flush=True)
        return self.originals[index]

    def engine(self,index,mm,backend='gpu'):
        key=(index,mm,backend)
        if key not in self.evaluators:self.evaluators[key]=evaluator(self.original(index),mm,backend)
        return self.evaluators[key]

    def rank_all(self):
        rows=[]
        for identity,item in self.items.items():
            folder=self.out/'ranking'/identity;meta=folder/'metadata.json';npz=folder/'selected.npz'
            if meta.exists():
                row=checked_read(meta)
                if row['fingerprint']!=self.fp or row['candidate']!=item['candidate'] or row['identity']!=identity:raise ValueError('rank context mismatch')
                if row.get('npz_sha256') and sha(npz)!=row['npz_sha256']:raise ValueError('rank artifact changed')
                if row.get('npz_sha256'):
                    with np.load(npz) as z:
                        if not math.isclose(float(z['q_MW_m2'].max()),row['description']['peak_MW_m2'],rel_tol=2e-11,abs_tol=1e-12):raise ValueError('cached ranked peak/array mismatch')
            else:
                c=Candidate.from_record(item['candidate']);e=self.engine(2,c.radius_mm)
                self.persist('ranking',current_identity=identity,ranking_completed=len(rows))
                started=time.perf_counter();r=evaluate_candidate(e,c);desc=e.describe(r)
                folder.mkdir(parents=True,exist_ok=True)
                if 'deposition' in r:export_boundary(npz,e,r)
                row=dict(identity=identity,candidate=c.record(),fingerprint=self.fp,description=desc,
                    geometry_screen=r.get('geometry_screen'),seconds=time.perf_counter()-started,
                    npz_sha256=sha(npz) if npz.exists() else None)
                # Store the canonical screen also inside the ranked description.
                row['description']['geometry_screen']=r.get('geometry_screen',{})
                checked_save(meta,row)
            rows.append(row)
            print('RANK',len(rows),len(self.items),item['candidate']['radius_mm'],identity[:12],row['description'].get('peak_MW_m2'),flush=True)
        self.report['ranking']={str(mm):[r for r in rows if r['candidate']['radius_mm']==mm] for mm in (20,50)}
        self.persist('ranked_all',ranking_completed=len(rows));return rows

    def baseline_audit(self,index):
        if index not in self.base_audits:
            o=self.original(index);self.base_audits[index]=audit(o,o.baseline)
        return self.base_audits[index]

    def validate(self,row):
        c=Candidate.from_record(row['candidate']);folder=self.out/f'R{c.radius_mm:03d}'/'validation'/c.identity
        folder.mkdir(parents=True,exist_ok=True);target=folder/'validation.json'
        if target.exists():
            old=checked_read(target)
            if old['fingerprint']!=self.fp or old['candidate']!=c.record():raise ValueError('validation context mismatch')
            if old.get('complete'):
                for rel,h in old.get('artifacts',{}).items():
                    if sha(folder/rel)!=h:raise ValueError('validation artifact changed')
                if old['verdict']!=verdict(old['levels'],old.get('wall_only_refinement')):raise ValueError('cached verdict mismatch')
                if len(old['levels'])>=3 and not math.isclose(old['levels'][2]['rows']['selected']['peak_MW_m2'],row['description']['peak_MW_m2'],rel_tol=2e-11,abs_tol=1e-12):raise ValueError('cached ranking/validation mismatch')
                return old
        record=dict(candidate=c.record(),identity=c.identity,fingerprint=self.fp,qualification_policy=QUALIFICATION_POLICY,
                    geometry_policy=GEOMETRY_POLICY,reference=self.config['reference'],levels=[],complete=False,artifacts={})
        for index,settings in enumerate(LEVELS):
            self.persist('qualifying',current_identity=c.identity,current_level=f'{settings[0]}_{settings[1]}')
            started=time.perf_counter();original=self.original(index);e=self.engine(index,c.radius_mm)
            r=evaluate_candidate(e,c);selected=audit(e,r)
            geom=geometric_audit(e,c,r) if r['geometry_valid'] else dict(passed_independent=False)
            selected['audit_passed'] &= geom['passed_independent'];base=self.baseline_audit(index)
            record['levels'].append(dict(method=settings[0],tubes=settings[1],trace_step_m=settings[2],nodes=settings[3],
                wall_step_m=settings[4],rows=dict(baseline=base,selected=selected),geometry_audit=geom))
            if 'deposition' in r:
                name=f'{settings[0]}_{settings[1]}.npz';export_boundary(folder/name,e,r);record['artifacts'][name]=sha(folder/name)
            if index==2 and not math.isclose(selected.get('peak_MW_m2',float('nan')),row['description']['peak_MW_m2'],rel_tol=2e-11,abs_tol=1e-12):raise ValueError('1280 ranking/validation changed')
            if index==0 and selected['audit_passed'] and base['audit_passed']:
                if self.wall_original is None:
                    self.wall_original=BoundaryEvaluator(original.direct,original.metadata,self.reference,original.wall_step/2,thermal_paths=original.thermal_paths)
                    self.wall_base=audit(self.wall_original,self.wall_original.baseline)
                    if not self.wall_base['audit_passed']:raise ValueError('fixed baseline wall-refinement audit failed')
                    export_boundary(self.out/'baseline/wall_refinement.npz',self.wall_original,self.wall_original.baseline)
                ef=evaluator(self.wall_original,c.radius_mm,'gpu');rf=evaluate_candidate(ef,c);s=audit(ef,rf)
                gf=geometric_audit(ef,c,rf) if rf['geometry_valid'] else dict(passed_independent=False)
                s['audit_passed'] &= gf['passed_independent']
                record['wall_only_refinement']=dict(rows=dict(baseline=self.wall_base,selected=s),geometry_audit=gf,
                    max_peak_relative_change=max(abs(self.wall_base['peak_q']/base['peak_q']-1),abs(s.get('peak_q',1e9)/selected['peak_q']-1)))
                if 'deposition' in rf:
                    name='wall_refinement.npz';export_boundary(folder/name,ef,rf);record['artifacts'][name]=sha(folder/name)
            checked_save(target,record)
            print('AUDIT',c.radius_mm,c.identity[:12],settings[0],settings[1],selected['audit_passed'],selected.get('peak_MW_m2'),'seconds',time.perf_counter()-started,flush=True)
            if not base['audit_passed']:raise ValueError('fixed baseline independent audit failed')
            if not selected['audit_passed']:break
        record['verdict']=verdict(record['levels'],record.get('wall_only_refinement'));record['complete']=True
        checked_save(target,record);print('VERDICT',c.radius_mm,c.identity[:12],record['verdict'],flush=True)
        return record

    def promote(self,row,v):
        if not passes_improvement(v):raise ValueError('cannot promote unconfirmed improvement')
        c=Candidate.from_record(row['candidate']);e=self.engine(2,c.radius_mm);r=evaluate_candidate(e,c)
        cpu=self.engine(2,c.radius_mm,'cpu');rc=evaluate_candidate(cpu,c)
        if not r['feasible'] or not rc['feasible']:raise ValueError('parity geometry failed')
        parity=dict(passed=True)
        for key,atol in [('q',1e-12),('cell_power',1e-14)]:
            np.testing.assert_allclose(r['deposition'][key],rc['deposition'][key],rtol=2e-11,atol=atol)
            parity['max_absolute_'+key]=float(np.max(abs(r['deposition'][key]-rc['deposition'][key])))
        peak_difference=abs(r['peak_q']/rc['peak_q']-1)
        if not math.isfinite(peak_difference) or peak_difference>1e-10:raise ValueError('peak parity failed')
        protection=check_preservation(self.config['protected_sha256']);code=check_preservation(self.config['code_sha256'])
        if not protection['unchanged'] or not code['unchanged']:raise ValueError('inputs/code changed before qualified publication')
        store=QualifiedStore(self.out/f'R{c.radius_mm:03d}',context(self.config,'OMP',c.radius_mm))
        manifest=dict(**context(self.config,'OMP',c.radius_mm),status='qualified',candidate=c.record(),
            selected=e.describe(r),baseline=self.original(2).describe(self.original(2).baseline),validation=v,
            cpu_gpu_peak_relative_difference=peak_difference,cpu_gpu_distribution=parity,
            qualification_utc=utc(),origin=self.items[c.identity],global_optimum_proven=False,
            prepublication_preservation=dict(checked_utc=utc(),inputs_unchanged=True,code_unchanged=True),
            selection_scope=self.report['scope'])
        def writer(folder):
            export_boundary(folder/'selected_refined.npz',e,r)
            export_boundary(folder/'baseline_refined.npz',self.original(2),self.original(2).baseline)
            save(folder/'validation.json',v);save(folder/'source_candidate.json',self.items[c.identity])
        published=store.publish(manifest,writer)
        if store.best()['candidate']!=c.record():raise ValueError('published best candidate mismatch')
        save(self.out/f'R{c.radius_mm:03d}'/'manifest.json',published)
        return published

    def execute(self):
        started=time.perf_counter()
        try:
            if self.report.get('complete'):
                if (self.out/'COMPLETED.json').exists():
                    marker=read(self.out/'COMPLETED.json')
                    if (marker['fingerprint']!=self.fp or marker['summary_sha256']!=sha(self.out/'summary.json')
                            or marker['artifact_index_sha256']!=sha(self.out/'artifact_sha256.json')):raise ValueError('completion marker mismatch')
                    for rel,h in read(self.out/'artifact_sha256.json').items():
                        if sha(self.out/rel)!=h:raise ValueError('completed artifact changed: '+rel)
                else:self.commit_index()
                print('ALREADY_COMPLETE',self.out,flush=True);return
            rows=self.rank_all()
            for mm in (20,50):
                store=QualifiedStore(self.out/f'R{mm:03d}',context(self.config,'OMP',mm));prior=store.best()
                if prior:
                    v=prior['validation'];recomputed=verdict(v['levels'],v.get('wall_only_refinement'))
                    if (not passes_improvement(v) or recomputed!=v['verdict']
                            or prior.get('cpu_gpu_distribution',{}).get('passed') is not True):
                        raise ValueError('recovered bundle fails stricter improvement/distribution contract')
                ranked=ranked_improvements([r for r in rows if r['candidate']['radius_mm']==mm]);attempts=[];selected=None
                for rank,row in enumerate(ranked,1):
                    # Re-run no completed heat audits; validate() verifies saved context/array checksums.
                    v=self.validate(row);attempts.append(dict(identity=row['identity'],fine_rank=rank,verdict=v['verdict'],passed=passes_improvement(v)))
                    self.report['attempts'][str(mm)]=attempts;self.persist('selection')
                    if passes_improvement(v):
                        if prior and Candidate.from_record(prior['candidate']).identity==row['identity']:
                            selected=prior
                        else:selected=self.promote(row,v)
                        break
                result=dict(status='qualified_improvement' if selected else 'no_passing_archived_improvement',
                    ranked_improving=len(ranked),fully_attempted=len(attempts),identity=Candidate.from_record(selected['candidate']).identity if selected else None,
                    selected_peak=selected['selected']['peak_MW_m2'] if selected else None,
                    baseline_peak=selected['baseline']['peak_MW_m2'] if selected else None,
                    bundle=selected['artifact_bundle'] if selected else None,
                    remaining_reason='equal-or-higher 1280 peak than first fully passing candidate' if selected else 'finite eligible list exhausted')
                self.report['results'][str(mm)]=result;self.persist('case_complete')
                print('RESULT',mm,result,flush=True)
            protection=check_preservation(self.config['protected_sha256']);code=check_preservation(self.config['code_sha256'])
            self.report.update(preservation=protection,code_preservation=code,elapsed_seconds=time.perf_counter()-started)
            if not protection['unchanged'] or not code['unchanged']:raise ValueError('protected inputs or code changed')
            self.persist('complete',complete=True)
            self.commit_index()
        except BaseException:
            self.persist('failed',failure=traceback.format_exc());raise

    def commit_index(self):
        # Recovery may finish an interrupted commit, never bless altered cached evidence.
        if (self.out/'artifact_sha256.json').exists():
            for rel,h in read(self.out/'artifact_sha256.json').items():
                if sha(self.out/rel)!=h:raise ValueError('indexed artifact changed: '+rel)
        for mm in (20,50):
            cached=[]
            for identity,item in self.items.items():
                if item['candidate']['radius_mm']!=mm:continue
                folder=self.out/'ranking'/identity;row=checked_read(folder/'metadata.json')
                if row['fingerprint']!=self.fp or row['candidate']!=item['candidate']:raise ValueError('cached rank context mismatch')
                if row.get('npz_sha256') and sha(folder/'selected.npz')!=row['npz_sha256']:raise ValueError('rank arrays changed')
                cached.append(row)
            if cached!=self.report['ranking'][str(mm)]:raise ValueError('summary ranking mismatch')
            for attempt in self.report['attempts'][str(mm)]:
                folder=self.out/f'R{mm:03d}'/'validation'/attempt['identity'];v=checked_read(folder/'validation.json')
                if v['fingerprint']!=self.fp or v['verdict']!=attempt['verdict'] or passes_improvement(v)!=attempt['passed']:raise ValueError('summary attempt mismatch')
                for rel,h in v['artifacts'].items():
                    if sha(folder/rel)!=h:raise ValueError('validation arrays changed')
            best=QualifiedStore(self.out/f'R{mm:03d}',context(self.config,'OMP',mm)).best()
            if best and (not passes_improvement(best['validation']) or best.get('cpu_gpu_distribution',{}).get('passed') is not True):raise ValueError('invalid recovered strict bundle')
            if self.report['results'][str(mm)]['identity']!=(Candidate.from_record(best['candidate']).identity if best else None):raise ValueError('result/bundle mismatch')
        save(self.out/'artifact_sha256.json',{p.relative_to(self.out).as_posix():sha(p) for p in self.out.rglob('*')
            if p.is_file() and p.name not in ('artifact_sha256.json','COMPLETED.json','writer.lock','stdout.log','stderr.log')})
        save(self.out/'COMPLETED.json',dict(fingerprint=self.fp,summary_sha256=sha(self.out/'summary.json'),
            artifact_index_sha256=sha(self.out/'artifact_sha256.json')))


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output.resolve()
    if out.drive.lower()!='d:' or any(out.is_relative_to(x) for x in (CAMPAIGN,REPO,HERE)):raise ValueError('output must be a separate D: study')
    out.mkdir(parents=True,exist_ok=True)
    with (out/'writer.lock').open('a+b') as lock:
        lock.seek(0);lock.write(b'0');lock.flush();lock.seek(0)
        msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        try:Run(out).execute()
        finally:lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)


if __name__=='__main__':main()
