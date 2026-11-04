"""Reassess immutable archived multiresolution evidence without changing physics."""
from pathlib import Path
import copy,itertools,time,shutil
from export_qualified_finals import verified_study,STUDY
from qualify_omp_stages import read,save,sha,checked_read,fingerprint,utc,check_preservation
from distribution_policy_v3 import from_npz,distribution_compare,POLICY,LIMITS,qualify

HERE=Path(__file__).resolve().parent
OUT=Path('D:/DivertorOptimization/studies/omp-distribution-v3-20261104')
NAMES=['cubic_320','cubic_640','cubic_1280','quintic_640','wall_refinement']
PAIRS=[(1,2),(1,3),(2,3),(0,4)]

def assess(paths):
    profiles=[from_npz(p) for p in paths]
    checks=[]
    for i,j in PAIRS:
        checks.append(dict(pair=[NAMES[i],NAMES[j]],**distribution_compare(profiles[i],profiles[j],profiles[2])))
    coarse=distribution_compare(profiles[0],profiles[2],profiles[2])
    return dict(passed=all(c['passed'] for c in checks),checks=checks,coarse320_diagnostic=coarse,
        A2_cell_m2=[float(sum(p['q']*p['area'])**2/sum(p['q']**2*p['area'])) for p in profiles])

def run():
    start=time.perf_counter();old=verified_study()
    protected={str(STUDY/p):h for p,h in read(STUDY/'artifact_sha256.json').items()}
    protected.update({str(STUDY/n):sha(STUDY/n) for n in ['summary.json','COMPLETED.json','artifact_sha256.json']})
    for p in [HERE/'dist/qualified-finals.json',HERE/'dist/data.json',HERE/'dist/legs.json']:
        protected[str(p)]=sha(p)
    code={str(HERE/n):sha(HERE/n) for n in ['distribution_policy_v3.py','requalify_distribution_v3.py','test_distribution_policy_v3.py']}
    config=dict(policy=POLICY,limits=LIMITS,continuous_reconstruction_power_limit=.001,
        prior_summary_sha256=sha(STUDY/'summary.json'),protected_sha256=protected,code_sha256=code,
        reference=old['config']['reference'],original_numerical_gates='unchanged anchor semantics except argmax diagnostic',
        comparison='same-support saved continuous piecewise-linear q; all fine pairs + wall refinement')
    fp=fingerprint(config)
    if (OUT/'COMPLETED.json').exists():raise FileExistsError('completed v3 study must not be overwritten')
    OUT.mkdir(parents=True,exist_ok=True)
    save(OUT/'config.json',dict(**config,fingerprint=fp))
    for path in code:shutil.copyfile(path,OUT/Path(path).name)
    baseline_paths=[STUDY/'baseline'/n/'baseline.npz' if n!='wall_refinement' else STUDY/'baseline/wall_refinement.npz' for n in NAMES]
    baseline=assess(baseline_paths);print('BASELINE',baseline['passed'],[(c['pair'],c['failures'],c['linf_relative']) for c in baseline['checks']],flush=True)
    report=dict(complete=False,created_utc=utc(),fingerprint=fp,policy=POLICY,baseline=baseline,candidates=[],results={},
        scope='Reevaluation of all17 previously audited archived candidates; no new geometry search or continuous optimum proof',
        original_53_ranking=old['ranking'],global_optimum_proven=False)
    for mm in (20,50):
        for path in sorted((STUDY/f'R{mm:03d}/validation').glob('*/validation.json')):
            original=checked_read(path);identity=original['identity']
            if original['fingerprint']!=old['config']['fingerprint']:raise ValueError('original context mismatch')
            for rel,h in original['artifacts'].items():
                if sha(path.parent/rel)!=h:raise ValueError('NPZ hash changed')
            selected=assess([path.parent/(n+'.npz') for n in NAMES])
            distribution=dict(policy=POLICY,passed=baseline['passed'] and selected['passed'],baseline=baseline,selected=selected)
            v=qualify(original,distribution);v.update(fingerprint=fp,prior_validation_path=str(path),prior_validation_sha256=sha(path))
            item=old['candidates'][identity]
            row=dict(identity=identity,radius_mm=mm,stage=item['stage'],peak=v['levels'][2]['rows']['selected']['peak_MW_m2'],
                previous_pass=original['verdict']['improvement_confirmed'],passed=v['verdict']['improvement_confirmed'],
                selected_distribution_passed=selected['passed'],baseline_distribution_passed=baseline['passed'],
                failures=[dict(pair=c['pair'],failures=c['failures']) for c in selected['checks'] if not c['passed']])
            save(OUT/'validations'/identity/'validation.json',v);report['candidates'].append(row)
            save(OUT/'summary.json',report);print('CANDIDATE',mm,item['stage'],row['passed'],row['failures'],flush=True)
        passing=sorted([r for r in report['candidates'] if r['radius_mm']==mm and r['passed']],key=lambda r:(r['peak'],r['identity']))
        report['results'][str(mm)]=dict(status='requires_publication_and_parity' if passing else 'no_candidate_passed_v3',best=passing[0] if passing else None)
    preservation=check_preservation(protected);codes=check_preservation(code)
    if not preservation['unchanged'] or not codes['unchanged']:raise ValueError('protected input/code changed')
    report.update(complete=True,elapsed_seconds=time.perf_counter()-start,preservation=preservation,code_preservation=codes)
    save(OUT/'summary.json',report)
    artifacts={p.relative_to(OUT).as_posix():sha(p) for p in OUT.rglob('*') if p.is_file() and p.name not in ['COMPLETED.json','artifact_sha256.json']}
    save(OUT/'artifact_sha256.json',artifacts);save(OUT/'COMPLETED.json',dict(fingerprint=fp,summary_sha256=sha(OUT/'summary.json'),artifact_index_sha256=sha(OUT/'artifact_sha256.json')))
    print('COMPLETE',report['results'],flush=True)

if __name__=='__main__':run()
