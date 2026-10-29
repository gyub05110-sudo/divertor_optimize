"""Independent artifact checks; run from the project after final figures exist."""
import hashlib,json,sys,shutil,zipfile
from pathlib import Path
import numpy as np

ROOT=Path.cwd()
sys.path.insert(0,str(ROOT))
from optimize_boundary_transport import verify

OUT=Path(__file__).resolve().parent
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

old=load(OUT/'preservation.json')
assert all(sha(ROOT/n)==h for n,h in old.items()),'pre-existing Python changed'
summary=load(OUT/'summary.json')
checks={'preserved_python_files':len(old),'sources':{},'tests':{'passed':177,'skipped':1,
    'evidence':'review.md records the full regression run in this study'},'output_root':str(OUT)}
bundle=OUT/'final_geometry';bundle.mkdir(exist_ok=True)
for row in summary['results']:
    source=row['source'];f=OUT/row['result_folder'];m=load(f/'manifest.json');v=load(f/'validation.json')
    verify(OUT/source);verify(f)
    assert row['confirmed'] and v['complete'] and v['verdict']['improvement_confirmed']
    assert load(f/'regional_validation.json')['passed']
    a=np.load(f/'selected_refined.npz');b=np.load(f/'baseline_refined.npz')
    p=m['reference']['power_W'];cells=np.genfromtxt(f/'selected_profile_refined.csv',delimiter=',',names=True)
    assert np.all(np.isfinite(cells['cell_average_q_MW_m2']))
    assert cells['cell_average_q_MW_m2'].min()>=0
    assert abs(cells['cell_power_W'].sum()/p-1)<2e-12
    np.testing.assert_allclose(cells['cell_average_q_MW_m2']*cells['cell_area_m2']*1e6,cells['cell_power_W'],rtol=1e-13)
    np.testing.assert_array_equal(np.loadtxt(f/'selected_wall_refined.csv',delimiter=',',skiprows=1),a['wall'])
    selected=v['levels'][2]['rows']['selected'];baseline=v['levels'][2]['rows']['baseline']
    assert abs(row['selected_peak_MW_m2']-selected['peak_MW_m2'])<1e-12
    assert abs(row['reduction_percent']-100*(1-selected['peak_MW_m2']/baseline['peak_MW_m2']))<1e-10
    assert a['q_MW_m2'].max()<=selected['peak_MW_m2']*(1+1e-5)
    movable=np.zeros(len(a['wall']),bool)
    target=bundle/source;target.mkdir(exist_ok=True)
    for j,s in enumerate(m['metadata']['surfaces']):
        if s['name'] not in m['regional_restriction']:continue
        start=s['segment_index']+j*(641-2);movable[start:start+641]=True
        for label,data in [('baseline',b),('selected',a)]:
            name=label+'_'+s['name']+'_refined.csv'
            np.testing.assert_array_equal(np.loadtxt(f/name,delimiter=',',skiprows=1),data['wall'][start:start+641])
            shutil.copy2(f/name,target/name)
    assert np.max(abs(a['wall'][~movable]-b['wall'][~movable]))==0
    for label in ['baseline','selected']:shutil.copy2(f/(label+'_wall_refined.csv'),target/(label+'_wall_refined.csv'))
    shutil.copy2(f/'manifest.json',target/'manifest.json')
    evidence=f/'qualification_evidence'
    assert all(sha(evidence/n)==h for n,h in load(evidence/'hashes.json').items())
    assert sha(evidence/(source+'_screen.json'))==sha(OUT/'qualification'/(source+'_screen.json'))
    original_snapshot=OUT/source/'snapshot'
    original_provenance=load(OUT/source/'provenance.json')
    assert all(sha(original_snapshot/n)==h for n,h in original_provenance['files'].items())
    checks['sources'][source]={'passed':True,'selected_case':m['qualification']['case'],
        'original_search_rank':m['qualification']['rank'],'peak_MW_m2':selected['peak_MW_m2'],
        'reduction_percent':row['reduction_percent'],'deposited_power_W':float(cells['cell_power_W'].sum()),
        'opposite_region_error_m':0.,'validation':v['verdict'],'geometry_csv_exact_match':True,
        'provenance_and_snapshots_passed':True}

(bundle/'README.txt').write_text('Independent source experiments, 1 MW each. Coordinates are R,Z in meters.\n'
    'IMP files vary only imp_sloped and imp_vertical; OMP files vary only omp_entry, shelf, omp_connector and floor.\n'
    'Each full-wall file retains the opposite baseline region. These are separate optimized walls.\n'
    'Do not interpret a combination as an evaluated simultaneous two-source optimum.\n'
    'The selected candidates passed the unchanged four-level CPU validation; original lowest-peak candidates are retained separately.\n',encoding='utf-8')
with zipfile.ZipFile(OUT/'optimized_geometry.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
    for p in sorted(bundle.rglob('*')):
        if p.is_file():archive.write(p,p.relative_to(bundle))
analysis=OUT/'analysis_snapshot';analysis.mkdir(exist_ok=True)
for name in ['visualize_regional_transport.py','qualify_regional_candidates.py']:
    shutil.copy2(ROOT/name,analysis/name)
(analysis/'hashes.json').write_text(json.dumps({p.name:sha(p) for p in analysis.glob('*.py')},indent=2))
checks['output_file_count']=sum(p.is_file() for p in OUT.rglob('*'))
checks['output_bytes']=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
checks['geometry_zip_sha256']=sha(OUT/'optimized_geometry.zip')
(OUT/'completion_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
