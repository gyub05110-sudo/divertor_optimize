"""User-selected design portions, separate from physical first-hit coverage."""
import copy
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
POLICY=ROOT/'baseline_review/optimization_region_policy.json'

def apply_design_policy(result):
    policy=json.loads(POLICY.read_text(encoding='utf-8'))
    selected=[];excluded=[]
    cutoff=policy['imp_maximum_Z_m']
    for row in result['footprints']:
        a,b=row['points'];item=copy.deepcopy(row)
        if row['source']=='OMP' and row['segment'] in policy['excluded_omp_segments']:
            item['exclusion_reason']='User excluded the outer sloped green wall region'
            excluded.append(item);continue
        if row['source']=='IMP':
            if a[1]>cutoff and b[1]>cutoff:
                item['exclusion_reason']='IMP wall above Z=-1.15 m';excluded.append(item);continue
            if (a[1]-cutoff)*(b[1]-cutoff)<0:
                t=(cutoff-a[1])/(b[1]-a[1]);cross=[a[0]+t*(b[0]-a[0]),cutoff]
                low=a if a[1]<=cutoff else b;high=b if a[1]<=cutoff else a
                item['points']=[low,cross]
                upper=copy.deepcopy(row);upper['points']=[cross,high]
                upper['exclusion_reason']='IMP wall above Z=-1.15 m';excluded.append(upper)
        selected.append(item)
    result['optimization_policy']=policy
    result['optimization_footprints']=selected
    result['excluded_optimization_footprints']=excluded
    result['requested_design']='floor+shelf union selected IMP portions at Z<=-1.15 and OMP first hits excluding segments 53/54/55'
    return result

if __name__=='__main__':
    path=ROOT/'baseline_review/baseline_review.json'
    data=apply_design_policy(json.loads(path.read_text(encoding='utf-8')))
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    assert all(max(p[1] for p in r['points'])<=-1.15 for r in data['optimization_footprints'] if r['source']=='IMP')
    assert not any(r['source']=='OMP' and r['segment'] in [53,54,55] for r in data['optimization_footprints'])
    assert any(abs(p[1]+1.15)<1e-12 for r in data['optimization_footprints'] if r['segment']==70 for p in r['points'])
    print(json.dumps({'selected_segments':{s:[r['segment'] for r in data['optimization_footprints'] if r['source']==s] for s in ['IMP','OMP']},'IMP_cutoff_Z_m':-1.15,'excluded_OMP_segments':[53,54,55],'physical_first_hits_preserved':True}))
