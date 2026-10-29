"""Export only verified, recorded geometry and heat; no optimization execution."""
import csv,json,hashlib,urllib.request
from pathlib import Path
import numpy as np
from shapely.geometry import LineString,Point
p=Path(__file__).resolve().parent;base=p.parent
a=np.load(base/'case_arrays.npz');m=json.loads((base/'manifest.json').read_text())
assert m['numerically_verified_feasible']
samples=list(csv.DictReader((base/'comparison/heat_samples.csv').open(encoding='utf-8-sig')))
summary=list(csv.DictReader((base/'width_scan.csv').open()))
result={'widths':[1,2,3,5,10],'time_ms':m['time_ms'],'retained_case':m['retained_case'],
        'units':'q/P_outer (m^-2)','lcfs':a['lcfs'].tolist(),'models':{},
        'source_geometry_sha256':m['config']['geometry_sha256']}
for key,label in [('baseline','기준 형상'),('candidate','최적화 형상')]:
 wall=a[key+'_wall'];line=LineString(wall)
 model={'wall':wall.tolist(),'floor':a[key+'_floor'].tolist(),'shelf':a[key+'_shelf'].tolist(),'heat':{},'summary':{}}
 for mm in result['widths']:
  rows=[r for r in samples if r['geometry']==label and int(r['lambda_q_mm'])==mm]
  projected=[]
  for row in rows:
   if row['component']=='미포착':continue
   pt=Point(float(row['R_m']),float(row['Z_m']))
   assert line.distance(pt)<1e-7
   projected.append([line.project(pt),float(row['q_over_P_outer_per_m2']),row['component']])
  projected.sort();model['heat'][str(mm)]=projected
  row=next(r for r in summary if r['geometry']==key and int(r['lambda_q_mm'])==mm)
  peak=max(row[1] for row in projected);assert abs(peak-float(row['all_wall_peak_q_over_P']))<1e-10
  model['summary'][str(mm)]={k:float(row[k]) for k in ['all_wall_peak_q_over_P','all_wall_peak_ratio','peak_R_m','peak_Z_m','floor_power','shelf_power','fixed_surface_power']}
 result['models'][key]=model
(p/'data.json').write_text(json.dumps(result,ensure_ascii=False,separators=(',',':')),encoding='utf-8')
vendor=p/'vendor';vendor.mkdir(exist_ok=True)
files={'three.module.js':'build/three.module.js','three.core.js':'build/three.core.js','OrbitControls.js':'examples/jsm/controls/OrbitControls.js','LICENSE':'LICENSE'}
hashes={}
for name,remote in files.items():
 target=vendor/name
 if not target.exists():
  url='https://unpkg.com/three@0.180.0/'+remote
  with urllib.request.urlopen(url,timeout=30) as response:target.write_bytes(response.read())
 hashes[name]=hashlib.sha256(target.read_bytes()).hexdigest()
(vendor/'provenance.json').write_text(json.dumps({'package':'three','version':'0.180.0','source':'https://unpkg.com/three@0.180.0/','sha256':hashes},indent=2))
print('Exported verified models and offline Three.js assets',flush=True)
