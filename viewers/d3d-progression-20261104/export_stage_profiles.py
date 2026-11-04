"""Display-only evaluation of archived milestones. Never optimizes or promotes a candidate."""
from pathlib import Path
import argparse,hashlib,json,sys,time,traceback
import numpy as np
from stage_profile_data import make_profile

HERE=Path(__file__).resolve().parent
REPO=Path('D:/DivertorOptimization/worktrees/d3d-report-correction-20261104')
CAMPAIGN=Path('D:/DivertorOptimization/studies/imp-omp-global-20261029')
sys.path.insert(0,str(REPO))
# Read parquet with the already frozen, hash-checked campaign dependency.
sys.path.append(str(CAMPAIGN/'snapshot/.local_dependencies/d3d-fixed-equilibrium'))
from d3d_optimization.boundary_evaluator import prepare_boundary,export_boundary
from d3d_optimization.campaign_geometry import Candidate,evaluator,evaluate_candidate
from d3d_optimization.reference_transport import ReferencePoint
from verification.run_candidate_requalification import consumed_code_paths
from run_global_campaign import runtime_identity

def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    temporary=p.with_suffix(p.suffix+'.tmp')
    temporary.write_text(json.dumps(x,ensure_ascii=False,allow_nan=False,separators=(',',':')),encoding='utf-8')
    temporary.replace(p)

def run(limit=None):
    started=time.perf_counter();dist=HERE/'dist';out=HERE/'stage-evaluations';out.mkdir(exist_ok=True)
    data=read(dist/'data.json');config=read(CAMPAIGN/'campaign.json');root=Path(config['frozen_root'])
    if data['reference']!=config['reference']:raise ValueError('viewer reference mismatch')
    protected={str(dist/'data.json'):sha(dist/'data.json'),str(dist/'legs.json'):sha(dist/'legs.json')}
    for p in (dist/'profiles').glob('*.json'):protected[str(p)]=sha(p)
    for relative,h in read(CAMPAIGN/'frozen_files.json').items():
        p=root/relative
        if sha(p)!=h:raise ValueError('frozen input mismatch '+str(p))
        protected[str(p)]=h
    protected[str(CAMPAIGN/'campaign.json')]=sha(CAMPAIGN/'campaign.json')
    code={str(p):sha(p) for p in consumed_code_paths(REPO)}
    code.update({str(p):sha(p) for p in [Path(__file__),HERE/'stage_profile_data.py']})
    context=dict(reference=config['reference'],data_sha256=sha(dist/'data.json'),
        frozen_files_sha256=sha(CAMPAIGN/'frozen_files.json'),code_sha256=code,
        resolution=dict(method='cubic',donors=320,trace_step_m=.0005,nodes=161,wall_step_m=.001),backend='gpu',runtime=runtime_identity('gpu'))
    fingerprint=hashlib.sha256(json.dumps(context,sort_keys=True).encode()).hexdigest()
    items={};origins={}
    for s in data['series']:
        for n,r in enumerate(s['stages'],1):
            c=Candidate.from_record(r['candidate'])
            if c.identity!=r['id']:raise ValueError('candidate identity mismatch')
            if r['id'] in items and abs(items[r['id']]['peak']/r['peak']-1)>1e-7:raise ValueError('duplicate peak mismatch')
            items[r['id']]=r;origins.setdefault(r['id'],[]).append(dict(series=s['id'],stage=n,origin=r['origin'],line=r['line'],peak=r['peak']))
    index=dict(schema='d3d-stage-profiles.v1',data_sha256=context['data_sha256'],fingerprint=fingerprint,cases={},complete=False)
    audit=dict(context=context,fingerprint=fingerprint,origins=origins,protected_sha256=protected,
        display_only=True,new_optimizations=0,qualification_changed=False,parity={},timings=[],complete=False)
    old=out/'audit.json'
    if old.exists():
        previous=read(old)
        if previous['fingerprint']!=fingerprint:raise ValueError('resume context mismatch')
        audit['parity']=previous['parity'];audit['timings']=previous['timings']
    save(old,audit)
    try:
        for source in ['IMP','OMP']:
            original=None;by_radius={};count=0
            for identity,row in items.items():
                if row['candidate']['source']!=source:continue
                c=Candidate.from_record(row['candidate']);target=out/identity;meta_path=target/'metadata.json'
                relative='stage-profiles/'+identity+'.json';profile_path=dist/relative
                if meta_path.exists():
                    meta=read(meta_path)
                    if meta['fingerprint']!=fingerprint or meta['identity']!=identity:raise ValueError('cached identity/context mismatch')
                    if sha(profile_path)!=meta['profile_sha256'] or sha(target/'selected.npz')!=meta['npz_sha256']:raise ValueError('cached artifact hash mismatch')
                else:
                    if original is None:
                        print('PREPARE',source,flush=True);t=time.perf_counter()
                        original=prepare_boundary(root,source,ReferencePoint(**config['reference']),count=320,step=.0005,nodes=161,wall_step=.001)
                        baseline_peak=original.describe(original.baseline)['peak_MW_m2']
                        if abs(baseline_peak/data['baseline']['coarse'][source]-1)>1e-7:raise ValueError('baseline model mismatch')
                        audit['timings'].append(dict(source=source,phase='prepare',seconds=time.perf_counter()-t))
                        print('READY',source,'seconds',time.perf_counter()-t,flush=True)
                    if c.radius_mm not in by_radius:by_radius[c.radius_mm]=evaluator(original,c.radius_mm,'gpu')
                    e=by_radius[c.radius_mm];t=time.perf_counter();r=evaluate_candidate(e,c)
                    if not r.get('feasible') or not r.get('geometry_screen',{}).get('passed'):raise ValueError('stage reevaluation failed '+identity+' '+str(r.get('errors')))
                    description=e.describe(r)
                    if source not in audit['parity']:
                        rc=evaluate_candidate(evaluator(original,c.radius_mm,'cpu'),c)
                        np.testing.assert_allclose(r['deposition']['q'],rc['deposition']['q'],atol=1e-12,rtol=2e-11)
                        np.testing.assert_allclose(r['deposition']['cell_power'],rc['deposition']['cell_power'],atol=1e-14,rtol=2e-11)
                        audit['parity'][source]=dict(identity=identity,passed=True,max_absolute_q=float(np.max(abs(r['deposition']['q']-rc['deposition']['q']))))
                    target.mkdir(exist_ok=True);export_boundary(target/'selected.npz',e,r)
                    with np.load(target/'selected.npz') as z:profile=make_profile(z,row['wall'],row['peak'],config['reference']['power_W'])
                    profile.update(identity=identity,source=source,donors=320,provenance='archived-stage-reevaluation',fingerprint=fingerprint)
                    save(profile_path,profile)
                    meta=dict(identity=identity,candidate=c.record(),fingerprint=fingerprint,description=description,
                        profile=relative,profile_sha256=sha(profile_path),npz_sha256=sha(target/'selected.npz'),
                        recorded_peak_relative_error=profile['recorded_peak_relative_error'],seconds=time.perf_counter()-t)
                    save(meta_path,meta)
                    audit['timings'].append(dict(identity=identity,phase='evaluate_export',seconds=meta['seconds']))
                index['cases'][identity]=dict(source=source,profile=relative,donors=320,sha256=meta['profile_sha256'],
                    peak=meta['description']['peak_MW_m2'],recorded_peak_relative_error=meta['recorded_peak_relative_error'])
                count+=1
                audit['completed']=len(index['cases']);save(old,audit);save(out/'index.partial.json',index)
                print('PROFILE',source,c.radius_mm,count,'total',len(index['cases']),identity[:12],flush=True)
                if limit is not None and len(index['cases'])>=limit:return
        if len(index['cases'])!=len(items):raise ValueError('incomplete index')
        for p,h in {**protected,**code}.items():
            if sha(p)!=h:raise ValueError('protected file changed '+p)
        index['complete']=True;save(dist/'stage-profiles.json',index)
        audit.update(complete=True,protected_unchanged=True,elapsed_seconds=time.perf_counter()-started)
        print('COMPLETE',len(items),'profiles',flush=True)
    except BaseException:
        audit['failure']=traceback.format_exc();raise
    finally:
        save(old,audit)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int);args=p.parse_args();run(args.limit)
