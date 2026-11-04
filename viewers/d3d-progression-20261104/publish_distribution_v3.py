"""Final reviewed v3 reassessment/publication; all historical data read-only."""
import copy,time,shutil,math
from pathlib import Path
import numpy as np
from export_qualified_finals import verified_study,STUDY,HERE
from qualify_omp_stages import read,save,sha,checked_read,fingerprint,utc,check_preservation,Run,Candidate,evaluate_candidate,export_boundary
from distribution_policy_v3_final import from_npz,distribution_compare,qualify,POLICY,LIMITS
from stage_profile_data import make_profile

OUT=Path('D:/DivertorOptimization/studies/omp-distribution-final-20261104')
EXT=Path('D:/DivertorOptimization/studies/omp-distribution-extension-20261104')
NAMES=['cubic_320','cubic_640','cubic_1280','quintic_640','wall_refinement']
PAIRS=[(1,2),(1,3),(2,3),(0,4)]

def assess(paths):
    p=[from_npz(path) for path in paths]
    checks=[dict(pair=[NAMES[i],NAMES[j]],**distribution_compare(p[i],p[j],p[2])) for i,j in PAIRS]
    return dict(passed=all(c['passed'] for c in checks),checks=checks,
                coarse320_diagnostic=distribution_compare(p[0],p[2],p[2]))

def run():
    started=time.perf_counter();old=verified_study();ext=read(EXT/'distribution_summary.json')
    if not ext['complete'] or not ext['preservation']['unchanged'] or not ext['code_preservation']['unchanged']:raise ValueError('remaining candidate audits incomplete or preservation failure')
    ext_config=read(EXT/'summary.json')['config']
    if ext_config!=old['config']:raise ValueError('extension numerical context differs from original')
    if OUT.exists():raise FileExistsError('use a new run folder; immutable final publication')
    data=read(HERE/'dist/data.json');oldindex=read(HERE/'dist/qualified-finals.json')
    protected=dict(old['config']['protected_sha256'])
    protected.update({str(STUDY/p):h for p,h in read(STUDY/'artifact_sha256.json').items()})
    for n in ['summary.json','artifact_sha256.json','COMPLETED.json']:protected[str(STUDY/n)]=sha(STUDY/n)
    for n in ['data.json','legs.json','qualified-finals.json','stage-profiles.json','full-saved-profiles.json']:
        protected[str(HERE/'dist'/n)]=sha(HERE/'dist'/n)
    for folder in ['qualified-profiles','qualified-evidence']:
        for p in (HERE/'dist'/folder).glob('*.json'):protected[str(p)]=sha(p)
    for p in EXT.rglob('*.npz'):protected[str(p)]=sha(p)
    for p in EXT.rglob('validation.json'):protected[str(p)]=sha(p)
    for n in ['distribution_summary.json','summary.json']:protected[str(EXT/n)]=sha(EXT/n)
    code=dict(old['config']['code_sha256'])
    for n in ['distribution_policy_v3_final.py','publish_distribution_v3.py','test_distribution_v3_adversarial.py','test_distribution_v3_review.py','stage_profile_data.py','export_qualified_finals.py']:
        code[str(HERE/n)]=sha(HERE/n)
    config=dict(policy=POLICY,limits=LIMITS,reference=data['reference'],protected_sha256=protected,code_sha256=code,
                previous_summary_sha256=sha(STUDY/'summary.json'),sample_power_error_limit=.001)
    fp=fingerprint(config);OUT.mkdir(parents=True);save(OUT/'config.json',dict(**config,fingerprint=fp))
    for path in code:
        p=Path(path);target=OUT/'source'/sha(p)/(p.name);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    baseline=assess([STUDY/'baseline'/n/'baseline.npz' if n!='wall_refinement' else STUDY/'baseline/wall_refinement.npz' for n in NAMES])
    rows=[];validations={};paths={}
    for mm in (20,50):
        for path in sorted((STUDY/f'R{mm:03d}/validation').glob('*/validation.json'))+sorted((EXT/f'R{mm:03d}/validation').glob('*/validation.json')):
            v0=checked_read(path);identity=v0['identity'];item=old['candidates'][identity]
            if identity in validations:raise ValueError('duplicate qualification')
            if v0['candidate']!=item['candidate'] or v0['reference']!=data['reference'] or v0['fingerprint']!=old['config']['fingerprint']:raise ValueError('candidate/reference/numerical context mismatch')
            for rel,h in v0['artifacts'].items():
                if sha(path.parent/rel)!=h:raise ValueError('audit/NPZ mismatch')
            selected=assess([path.parent/(n+'.npz') for n in NAMES])
            v=qualify(v0,dict(policy=POLICY,passed=baseline['passed'] and selected['passed'],baseline=baseline,selected=selected))
            v.update(fingerprint=fp,prior_validation_path=str(path),prior_validation_sha256=sha(path))
            validations[identity]=v;paths[identity]=path.parent
            checks=selected['checks'];maximum=max(c['linf_relative'] for c in checks)
            row=dict(identity=identity,series=item['series'],stage=item['stage'],radius_mm=mm,
                peak=v['levels'][2]['rows']['selected']['peak_MW_m2'],previous_pass=v0['verdict']['improvement_confirmed'],
                passed=v['verdict']['improvement_confirmed'],max_linf_relative=maximum,
                max_l1_relative=max(c['l1_relative'] for c in checks),max_a2_relative=max(c['a2_relative'] for c in checks),
                max_location_change_m=v0['verdict']['max_peak_location_change_m'],
                failures=[dict(pair=c['pair'],failures=c['failures']) for c in checks if not c['passed']])
            rows.append(row);save(OUT/'validations'/identity/'validation.json',v)
            print('FINAL_SCREEN',mm,item['stage'],row['passed'],'Linf',maximum,flush=True)
    if sum(r['radius_mm']==20 for r in rows)!=19:raise ValueError('all19 OMP20 candidates required')
    entries=[];results={};parity_worker=None
    for mm in (20,50):
        eligible=sorted([r for r in rows if r['radius_mm']==mm and r['passed']],key=lambda r:(r['peak'],r['identity']))
        best=eligible[0] if eligible else None;results[str(mm)]=dict(status='qualified_improvement' if best else 'no_archived_candidate_passed_v3',best=best)
        if not best:continue
        identity=best['identity'];v=validations[identity];item=old['candidates'][identity]
        npz=paths[identity]/'cubic_1280.npz';base_npz=STUDY/'baseline/cubic_1280/baseline.npz'
        old_entry=next((e for e in oldindex['entries'] if e['final']['id']==identity),None)
        if old_entry:
            m0=read(HERE/'dist'/old_entry['evidence']['path'])
            if sha(HERE/'dist'/old_entry['evidence']['path'])!=old_entry['evidence']['sha256'] or sha(npz)!=m0['artifact_sha256']['selected_refined.npz']:raise ValueError('reused parity does not bind exact arrays')
            parity=m0['cpu_gpu_distribution'];peak_parity=m0['cpu_gpu_peak_relative_difference']
            parity={**parity,'reused_exact_npz_sha256':sha(npz),'source_manifest_sha256':old_entry['evidence']['sha256']}
        else:
            if parity_worker is None:parity_worker=Run(OUT/'parity-runtime')
            c=Candidate.from_record(item['candidate']);eg=parity_worker.engine(2,mm);ec=parity_worker.engine(2,mm,'cpu')
            rg=evaluate_candidate(eg,c);rc=evaluate_candidate(ec,c)
            if not rg['feasible'] or not rc['feasible']:raise ValueError('parity feasibility failed')
            parity=dict(passed=True)
            with np.load(npz) as z:
                np.testing.assert_allclose(rg['deposition']['q'],z['q_MW_m2'],rtol=2e-11,atol=1e-12)
                np.testing.assert_allclose(rg['deposition']['cell_power'],z['cell_power'],rtol=2e-11,atol=1e-14)
            for key,atol in [('q',1e-12),('cell_power',1e-14)]:
                np.testing.assert_allclose(rg['deposition'][key],rc['deposition'][key],rtol=2e-11,atol=atol)
                parity['max_absolute_'+key]=float(np.max(abs(rg['deposition'][key]-rc['deposition'][key])))
            peak_parity=abs(rg['peak_q']/rc['peak_q']-1)
        if not parity['passed'] or not math.isfinite(peak_parity) or peak_parity>1e-10:raise ValueError('parity failure')
        bundle=OUT/'qualified'/identity;bundle.mkdir(parents=True)
        shutil.copyfile(npz,bundle/'selected_refined.npz');shutil.copyfile(base_npz,bundle/'baseline_refined.npz')
        save(bundle/'validation.json',v)
        manifest=dict(status='qualified',fingerprint=fp,reference=data['reference'],source='OMP',radius_mm=mm,
            candidate=item['candidate'],origin=item,qualification_policy=POLICY,qualification_utc=utc(),validation=v,
            selected=v['levels'][2]['rows']['selected'],baseline=v['levels'][2]['rows']['baseline'],
            cpu_gpu_distribution=parity,cpu_gpu_peak_relative_difference=peak_parity,
            prepublication_preservation=dict(inputs_unchanged=check_preservation(protected)['unchanged'],code_unchanged=check_preservation(code)['unchanged']),
            artifact_sha256={'selected_refined.npz':sha(npz),'baseline_refined.npz':sha(base_npz)},global_optimum_proven=False)
        if not all(manifest['prepublication_preservation'].values()):raise ValueError('preservation failed')
        save(bundle/'manifest.json',manifest);save(bundle/'checksums.json',{p.name:sha(p) for p in bundle.iterdir() if p.is_file()})
        series=next(s for s in data['series'] if s['id']==item['series']);stage=next(s for s in series['stages'] if s['id']==identity)
        peak=best['peak'];base=v['levels'][2]['rows']['baseline']['peak_MW_m2']
        with np.load(npz) as z:profile=make_profile(z,stage['wall'],peak,data['reference']['power_W'])
        profile.update(identity=identity,source='OMP',donors=1280,fingerprint=fp,npz_sha256=sha(npz),qualified_utc=manifest['qualification_utc'],provenance=POLICY)
        rel=f'distribution-profiles/{identity}.json';ev=f'distribution-evidence/{identity}.json'
        save(OUT/'viewer'/rel,profile);save(OUT/'viewer'/ev,manifest)
        f=copy.deepcopy(stage);f.update(peak=peak,baseline=base,ratio=peak/base,profile=rel,validation='current',additional_validation=True,
            distribution_policy=POLICY,qualified_utc=manifest['qualification_utc'],resolution='cubic · 1280 donor · 분포 수렴 v3',
            archived_stage=item['stage'],archived_origin=stage['origin'],origin=str(bundle/'manifest.json'),constraints=v['levels'][2]['geometry_audit'],qualification_fingerprint=fp)
        entries.append(dict(series=series['id'],final=f,evidence=dict(path=ev,sha256=sha(OUT/'viewer'/ev)),
            profile=dict(identity=identity,source='OMP',donors=1280,profile=rel,sha256=sha(OUT/'viewer'/rel),fingerprint=fp,npz_sha256=sha(npz),peak=peak,qualified_utc=manifest['qualification_utc'])))
    report=dict(complete=True,policy=POLICY,fingerprint=fp,created_utc=utc(),reference=data['reference'],baseline=baseline,
        candidates=rows,results=results,elapsed_seconds=time.perf_counter()-started,
        preservation=check_preservation(protected),code_preservation=check_preservation(code),
        sensitivity=dict(half_linf_pass=sum(r['max_linf_relative']<=.005 for r in rows),nominal_linf_pass=sum(r['max_linf_relative']<=.01 for r in rows),double_linf_pass=sum(r['max_linf_relative']<=.02 for r in rows)),
        global_optimum_proven=False,scope='all19 OMP20 archived milestones; OMP50 lowest1280 candidate passed, other33 have no lower stored peak')
    if not report['preservation']['unchanged'] or not report['code_preservation']['unchanged']:raise ValueError('final preservation failed')
    save(OUT/'summary.json',report)
    index=dict(schema='d3d-qualified-finals.v1',complete=True,data_sha256=sha(HERE/'dist/data.json'),reference=data['reference'],fingerprint=fp,
        qualification_study=str(OUT),summary_sha256=sha(OUT/'summary.json'),entries=entries,reassessment=rows,policy=POLICY)
    save(OUT/'viewer/distribution-finals.json',index)
    files={p.relative_to(OUT).as_posix():sha(p) for p in OUT.rglob('*') if p.is_file()}
    save(OUT/'artifact_sha256.json',files);save(OUT/'COMPLETED.json',dict(summary_sha256=sha(OUT/'summary.json'),artifact_index_sha256=sha(OUT/'artifact_sha256.json'),fingerprint=fp))
    # Publish only after the complete study and immutable bundles have committed.
    for p in (OUT/'viewer').rglob('*.json'):
        if p.name=='distribution-finals.json':continue
        target=HERE/'dist'/p.relative_to(OUT/'viewer');target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists() and sha(target)!=sha(p):raise ValueError('refuse overwriting versioned asset')
        shutil.copyfile(p,target)
    save(HERE/'dist/distribution-finals.json',index)
    print('PUBLISHED',results,flush=True)

if __name__=='__main__':run()
