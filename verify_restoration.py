"""Verify historical-copy fidelity, restored math, and reviewed topology."""
from pathlib import Path
import hashlib,json
import numpy as np
from shapely.geometry import Point,LineString,Polygon
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'RESTORE_MANIFEST.json').read_text(encoding='utf-8'))
overrides={'catalog.json','quality_catalog.json','start_viewer.ps1'}
count=0
for row in r['copied_files']:
    p=ROOT/row['path']
    if p.parent.name=='viewer3d' and p.name in overrides:continue
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
    count+=1
for old,new in [('frozen_field.py','frozen_field.py'),('frozen_mapping.py','frozen_mapping.py'),('shape.py','shape.py'),
                ('d3d_frozen_case.py','case.py'),('d3d_multisurface.py','multisurface.py')]:
    original=(ROOT/'archive_math'/old).read_text(encoding='utf-8').replace('.d3d_frozen_case','.case').replace('divertor_optimization.','d3d_optimization.')
    assert original==(ROOT/'d3d_optimization'/new).read_text(encoding='utf-8')
assert not (ROOT/'heat_unified').exists() and not (ROOT/'native_bo').exists()
assert not list(ROOT.glob('crossfield*.py'))
v=ROOT/'input/D3D/optimized/d3d-divertor-floor-shelf/viewer3d'
cat=json.loads((v/'catalog.json').read_text(encoding='utf-8'));assert len(cat['cases'])==1
case=json.loads((v/cat['cases'][0]['bundle_url']).read_text(encoding='utf-8'))
assert case['retained_case']==6758 and len(case['stages'])==29
assert abs(case['stages'][-1]['score']-.6442257832483174)<1e-12
q=json.loads((v/'quality/cae740322035-1400ms-a395bb40.json').read_text(encoding='utf-8'))
bp=q['models']['baseline']['summary']['3']['all_wall_peak_q_over_P']
cp=q['models']['candidate']['summary']['3']['all_wall_peak_q_over_P']
assert f'{bp:.3f}'=='2.909' and f'{cp:.3f}'=='1.756'
assert f'{100*(1-cp/bp):.1f}'=='39.6'
d=json.loads((ROOT/'baseline_review/baseline_review.json').read_text(encoding='utf-8'))
wall=LineString(d['wall']);xp=np.asarray(d['xpoint_RZ_m'])
assert set(d['legs'])=={'inner','outer'}
for name,leg in d['legs'].items():
    assert np.linalg.norm(np.array(leg['points'][0])-xp)<1e-12
    assert wall.distance(Point(leg['strike_RZ_m']))<1e-9
    assert leg['max_psiN_residual']<1e-6 and leg['strike_refinement_delta_m']<1e-5
assert d['legs']['inner']['segment']==69 and d['legs']['outer']['segment']==68
assert set(d['sources']['IMP']['hit_counts'])=={'69','70'}
assert set(d['sources']['OMP']['hit_counts'])=={'53','54','55','58','59','60','68'}
assert all(s['unresolved_paths']==0 for s in d['sources'].values())
assert Polygon(d['pfr_polygon']).is_valid
for row in d['footprints']:
    assert all(wall.distance(Point(p))<1e-8 for p in row['points'])
assert {r['segment'] for r in d['optimization_footprints'] if r['source']=='OMP'}=={58,59,60,68}
assert all(max(p[1] for p in r['points'])<=-1.15 for r in d['optimization_footprints'] if r['source']=='IMP')
assert any(abs(p[1]+1.15)<1e-12 for r in d['optimization_footprints'] if r['segment']==70 for p in r['points'])
assert d['optimization_policy']['excluded_wall_heat_remains_in_accounting']
summary=dict(historical_files_verified=count,intentional_viewer_overrides=sorted(overrides),restored_math_verified=True,
             screenshot_metrics=dict(baseline_q_over_P=bp,candidate_q_over_P=cp,reduction_percent=100*(1-cp/bp)),
             history_stages=29,legacy_retained_case=6758,leg_contour_residuals={n:l['max_psiN_residual'] for n,l in d['legs'].items()},
             sampled_topology_verified=True,continuous_footprint_boundaries_verified=False,independent_optimization_checked_by_this_script=False)
(ROOT/'baseline_review/restoration_checks.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
