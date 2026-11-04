"""Export only committed, independently qualified archived-candidate improvements."""
from pathlib import Path
import json,sys,copy,math
import numpy as np
from stage_profile_data import make_profile
from qualify_omp_stages import read,save,sha,checked_read,passes_improvement,Candidate,QualifiedStore,context

HERE=Path(__file__).resolve().parent
STUDY=Path('D:/DivertorOptimization/studies/omp-stage-qualification-20261104')


def verified_study(study=STUDY):
    marker=read(study/'COMPLETED.json');report=read(study/'summary.json')
    if (marker['summary_sha256']!=sha(study/'summary.json') or marker['artifact_index_sha256']!=sha(study/'artifact_sha256.json')
            or marker['fingerprint']!=report['config']['fingerprint']):raise ValueError('study completion hash mismatch')
    if not report['complete'] or not report['preservation']['unchanged'] or not report['code_preservation']['unchanged']:raise ValueError('study incomplete or input preservation failed')
    for rel,h in read(study/'artifact_sha256.json').items():
        if sha(study/rel)!=h:raise ValueError('study artifact changed: '+rel)
    return report


def run():
    report=verified_study();data=read(HERE/'dist/data.json');data_hash=sha(HERE/'dist/data.json')
    if data_hash!=report['config']['data_sha256'] or data['reference']!=report['config']['reference']:raise ValueError('viewer/catalog mismatch')
    index=dict(schema='d3d-qualified-finals.v1',complete=False,data_sha256=data_hash,reference=data['reference'],
        fingerprint=report['config']['fingerprint'],qualification_study=str(STUDY),summary_sha256=sha(STUDY/'summary.json'),entries=[])
    for mm in (20,50):
        result=report['results'][str(mm)]
        if result['status']!='qualified_improvement':continue
        folder=STUDY/f'R{mm:03d}';manifest=QualifiedStore(folder,context(report['config'],'OMP',mm)).best()
        if not manifest or not passes_improvement(manifest['validation']):raise ValueError('missing fully improved final')
        if manifest.get('cpu_gpu_distribution',{}).get('passed') is not True:raise ValueError('missing distribution parity')
        identity=Candidate.from_record(manifest['candidate']).identity
        if identity!=result['identity']:raise ValueError('summary/bundle mismatch')
        series=next(s for s in data['series'] if s['run']=='campaign' and s['source']=='OMP' and s['rmin']==mm)
        stage=next(r for r in series['stages'] if r['id']==identity)
        if stage['candidate']!=manifest['candidate']:raise ValueError('candidate mismatch')
        peak=manifest['selected']['peak_MW_m2'];base=manifest['baseline']['peak_MW_m2']
        if not math.isclose(base,data['baseline']['refined']['OMP'],rel_tol=1e-10):raise ValueError('baseline mismatch')
        bundle=folder/manifest['artifact_bundle'];npz=bundle/'selected_refined.npz'
        with np.load(npz) as arrays:profile=make_profile(arrays,stage['wall'],peak,data['reference']['power_W'])
        profile.update(identity=identity,source='OMP',donors=1280,provenance='additional-full-qualification-20261104',
            fingerprint=index['fingerprint'],npz_sha256=sha(npz),qualified_utc=manifest['qualification_utc'])
        path=f'qualified-profiles/{identity}.json';save(HERE/'dist'/path,profile)
        evidence=f'qualified-evidence/{identity}.json';save(HERE/'dist'/evidence,manifest)
        final=copy.deepcopy(stage)
        final.update(peak=peak,baseline=base,ratio=peak/base,profile=path,validation='current',additional_validation=True,
            qualified_utc=manifest['qualification_utc'],resolution='cubic · 1280 donor · 추가 전체 검증',
            archived_stage=manifest['origin']['stage'],archived_origin=stage['origin'],origin=str(bundle/'manifest.json'),
            constraints=manifest['validation']['levels'][2]['geometry_audit'],qualification_fingerprint=index['fingerprint'])
        index['entries'].append(dict(series=series['id'],final=final,
            profile=dict(identity=identity,source='OMP',donors=1280,profile=path,sha256=sha(HERE/'dist'/path),
                fingerprint=index['fingerprint'],npz_sha256=sha(npz),peak=peak,qualified_utc=manifest['qualification_utc']),
            evidence=dict(path=evidence,sha256=sha(HERE/'dist'/evidence))))
    if not index['entries']:raise ValueError('no improved finals to export')
    index['complete']=True;save(HERE/'dist/qualified-finals.json',index)
    print(json.dumps(dict(finals=[dict(series=e['series'],peak=e['final']['peak'],reduction_percent=100*(1-e['final']['ratio']),archived_stage=e['final']['archived_stage']) for e in index['entries']]),indent=2))


if __name__=='__main__':run()
