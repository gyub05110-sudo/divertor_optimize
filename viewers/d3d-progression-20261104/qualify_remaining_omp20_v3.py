"""Complete the three previously unaudited OMP20 milestones, no new search."""
import time
from pathlib import Path
from qualify_omp_stages import Run,read,save,utc,sha,check_preservation
from export_qualified_finals import STUDY,verified_study
from requalify_distribution_v3 import assess,NAMES
from distribution_policy_v3 import qualify,POLICY

OUT=Path('D:/DivertorOptimization/studies/omp-distribution-extension-20261104')
def main():
    started=time.perf_counter();old=verified_study();worker=Run(OUT)
    baseline=read(Path('D:/DivertorOptimization/studies/omp-distribution-v3-20261104/summary.json'))['baseline']
    results=[]
    for row in sorted(old['ranking']['20'],key=lambda r:r['description']['peak_MW_m2']):
        identity=row['identity']
        if (STUDY/f'R020/validation/{identity}/validation.json').exists():continue
        original=worker.validate(row)
        folder=OUT/f'R020/validation/{identity}'
        selected=assess([folder/(n+'.npz') for n in NAMES])
        v=qualify(original,dict(policy=POLICY,passed=baseline['passed'] and selected['passed'],baseline=baseline,selected=selected))
        save(folder/'distribution_v3.json',v)
        item=old['candidates'][identity];results.append(dict(identity=identity,stage=item['stage'],passed=v['verdict']['improvement_confirmed'],peak=row['description']['peak_MW_m2'],distribution=selected))
        print('V3_REMAINING',item['stage'],results[-1]['passed'],flush=True)
        save(OUT/'distribution_summary.json',dict(complete=False,results=results))
    preservation=check_preservation(worker.config['protected_sha256']);code=check_preservation(worker.config['code_sha256'])
    if not preservation['unchanged'] or not code['unchanged']:raise ValueError('input/code changed')
    save(OUT/'distribution_summary.json',dict(complete=True,results=results,created_utc=utc(),elapsed_seconds=time.perf_counter()-started,preservation=preservation,code_preservation=code,driver_sha256=sha(Path(__file__))))
    print('DONE',[(r['stage'],r['passed']) for r in results],flush=True)

if __name__=='__main__':main()
