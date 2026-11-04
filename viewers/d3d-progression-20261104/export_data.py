"""Export recorded D3D improvements; no optimizer or thermal solver is run."""
from pathlib import Path
import sys,json,hashlib,math,time
import numpy as np
from shapely.geometry import LineString

HERE=Path(__file__).resolve().parent
REPO=Path('D:/DivertorOptimization/worktrees/d3d-report-correction-20261104')
STUDIES=Path('D:/DivertorOptimization/studies')
sys.path.insert(0,str(REPO))
from verification.test_campaign_consistency import geometry_example
from d3d_optimization.campaign_geometry import Candidate,geometry_screen,assembled_wall

class RunningBest:
    """Original running records, before later geometry screening."""
    def __init__(self):self.best=1.
    def consider(self,ratio):
        if not isinstance(ratio,(int,float)) or not math.isfinite(ratio) or not 0<ratio<self.best-1e-9:return False
        self.best=ratio
        return True

def heat_cells(z,power_W):
    """Display bins preserve deposited power/area; never join unordered q samples."""
    wall=z['receiving_wall'];area=np.diff(z['area_edges']);power=z['cell_power']
    parent=z['cell_parent'];segments=[];i=0
    if len(wall)!=len(power)+1 or np.any(area<=0):raise ValueError('invalid receiving mesh')
    while i<len(power):
        j=i+1
        while j<min(i+4,len(power)) and parent[j]==parent[i]:j+=1
        a,b=wall[i],wall[j]
        if min(a[1],b[1])<-1.05 and min(a[0],b[0])<1.85:
            q=float(power[i:j].sum()*power_W/1e6/area[i:j].sum())
            segments.append([float(a[0]),float(a[1]),float(b[0]),float(b[1]),q])
        i=j
    return segments

def same_shape(a,b,tol=2e-5):
    return LineString(a).hausdorff_distance(LineString(b))<=tol

def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8',newline='\n')

def run():
    start=time.perf_counter();dist=HERE/'dist';dist.mkdir(exist_ok=True)
    evidence={};audit=[];screen_cache={};evaluators={}
    def tracked(p):
        p=Path(p);evidence[str(p)]=sha(p);return read(p)
    fixture=tracked(REPO/'verification/fixtures/campaign_geometry_20261029.json')
    campaign=STUDIES/'imp-omp-global-20261029'
    reference=tracked(campaign/'campaign.json')['reference']
    base_coarse={};base_fine={}
    for source in ['IMP','OMP']:
        m=tracked(campaign/source/'R050/manifest.json')
        base_coarse[source]=m['validation']['levels'][0]['rows']['baseline']['peak_MW_m2']
        base_fine[source]=m['validation']['levels'][2]['rows']['baseline']['peak_MW_m2']
    def evaluator(source,rad):
        key=(source,rad)
        if key not in evaluators:
            e,_=geometry_example('IMP/R020' if source=='IMP' else 'OMP/R020',nodes=161)
            e.rmin=rad/1000;evaluators[key]=e
        return evaluators[key]
    def inspect(c):
        if c.identity not in screen_cache:
            e=evaluator(c.source,c.radius_mm);coeff=c.bind(e)
            constraints=e.constraints(coeff);screen=geometry_screen(e,coeff)
            lengths,_=e.precise_lengths(coeff)
            good=constraints['passed'] and screen['passed'] and max(lengths.values())<=1.2+1e-10
            constraints={**constraints,'length_ratios':lengths}
            screen_cache[c.identity]=(good,constraints,screen)
        return screen_cache[c.identity]
    def geometry(c):
        e=evaluator(c.source,c.radius_mm);coeff=c.bind(e)
        return np.round(assembled_wall(e,coeff),10).tolist()
    def active_curves(c):
        e=evaluator(c.source,c.radius_mm);coeff=c.bind(e)
        curves=e.curves(coeff)
        names=['imp_sloped','imp_vertical'] if c.source=='IMP' else ['omp_entry','shelf','omp_connector','floor']
        return {name:np.round(curves[name].point(np.linspace(0,1,161)),10).tolist() for name in names}
    def candidate(row,source,rad):
        return Candidate.from_record(row['candidate']) if 'candidate' in row else Candidate(source,rad,tuple(row['coefficients_m'][:16]))
    def profile(path,name,wall=None):
        path=Path(path)
        if not path.exists():return None
        evidence[str(path)]=sha(path)
        with np.load(path) as z:
            if wall is not None and not same_shape(wall,z['wall']):raise ValueError('profile shape mismatch '+str(path))
            p=z['points'];q=z['q_MW_m2'];imax=int(np.argmax(q))
            segments=heat_cells(z,reference['power_W'])
            payload={'segments':segments,'peak':float(q[imax]),'peak_point':p[imax].tolist(),
                     'power_fraction':float(z['cell_power'].sum()),
                     'note':'area-weighted average of up to four consecutive saved cells within each parent segment; separate refined peak'}
            if not np.all(np.isfinite(segments)):raise ValueError('nonfinite profile')
        rel='profiles/'+name+'.json';save(dist/rel,payload);return rel
    qualification=STUDIES/'candidate-requalification-20261103'
    base_profiles={s:profile(qualification/f'baseline_{s}/cubic_1280/selected.npz','baseline-'+s,fixture['wall']) for s in ['IMP','OMP']}
    qualified=tracked(qualification/'summary.json')['candidates']
    qmap={};qpaths={}
    for name,entry in qualified.items():
        p=qualification/name/'validation.json'
        if p.exists():
            validation=tracked(p)
            if validation['identity']!=entry['identity']:raise ValueError('validation identity mismatch')
            qmap[entry['identity']]=validation;qpaths[entry['identity']]=p
    series=[]
    specs=[]
    for s in ['IMP','OMP']:
        for rad in [20,50,100]:specs.append(('campaign',s,rad,campaign/s/f'R{rad:03d}','72시간 다중 시작 탐색 · 10/29–11/01'))
    for rad in [20,50,100]:
        specs.append(('tangent','OMP',rad,REPO/f'results/omp-tangent-curvature-20261029/R{rad:03d}','접선·곡률 비교 연구 · 10/29'))
    for rad in [20,50,100]:
        specs.append(('pilot','IMP',rad,STUDIES/f'imp-global-pilot-20261029-v3/IMP/R{rad:03d}','IMP 접선 연결 초기 탐색 · 10/29'))
    for runid,source,rad,folder,label in specs:
        log=folder/'evaluations.jsonl'
        if not log.exists():continue
        print('Reading',runid,source,rad,flush=True)
        m=tracked(folder/'manifest.json') if (folder/'manifest.json').exists() else {}
        stats={'lines':0,'feasible_logged':0,'raw_records':0,'excluded_geometry':0,'excluded_duplicate':0}
        milestones=[];records=RunningBest();h=hashlib.sha256();seen=set()
        with log.open('rb') as f:
            for line in f:
                h.update(line);stats['lines']+=1
                if b'"feasible": true' not in line and b'"feasible":true' not in line:continue
                r=json.loads(line);stats['feasible_logged']+=1
                ratio=r.get('objective',r.get('peak_ratio',math.inf))
                if not records.consider(ratio):continue
                # These are actual logged record milestones; screen after selection.
                stats['raw_records']+=1
                c=candidate(r,source,rad)
                if c.identity in seen:stats['excluded_duplicate']+=1;continue
                seen.add(c.identity)
                good,constraints,screen=inspect(c)
                if not good:stats['excluded_geometry']+=1;continue
                peak=float(r.get('peak_MW_m2',ratio*base_coarse[source]))
                actual_ratio=peak/base_coarse[source]
                if abs(actual_ratio-ratio)>1e-7:raise ValueError('inconsistent baseline ratio')
                milestones.append({'id':c.identity,'candidate':c.record(),'ratio':actual_ratio,'peak':peak,'baseline':base_coarse[source],
                    'sequence':r.get('sequence',r.get('case',stats['lines'])),'epoch':r.get('epoch'),
                    'method':r.get('stage','search'),'wall':geometry(c),'active_curves':active_curves(c),'constraints':constraints,'geometry_screen':screen,
                    'validation':'geometry','profile':None,'resolution':'cubic · 320 donor',
                    'origin':str(log),'line':stats['lines']})
        evidence[str(log)]=h.hexdigest()
        final=None
        cr=m.get('candidate')
        if cr is None and 'coefficients_m' in m:cr=Candidate(source,rad,tuple(m['coefficients_m'][:16])).record()
        if cr:
            c=Candidate.from_record(cr);good,constraints,screen=inspect(c)
            validation=qmap.get(c.identity) or m.get('validation')
            validation_origin='current' if c.identity in qmap else 'historical'
            if validation is None and (folder/'validation.json').exists():validation=tracked(folder/'validation.json')
            if good and validation and validation.get('verdict',{}).get('resolution_screen_passed'):
                level=next((v for v in validation.get('levels',[]) if v.get('tubes')==1280),None)
                if level:
                    peak=level['rows']['selected']['peak_MW_m2'];base=level['rows']['baseline']['peak_MW_m2']
                    if peak<base:
                        wall=geometry(c)
                        final={'id':c.identity,'candidate':c.record(),'ratio':peak/base,'peak':peak,'baseline':base,'sequence':None,
                               'method':'final','wall':wall,'active_curves':active_curves(c),'constraints':constraints,'geometry_screen':screen,
                               'validation':'current' if validation_origin=='current' else 'historical',
                               'resolution':'cubic · 1280 donor','origin':str(qpaths.get(c.identity,folder/'manifest.json')),
                               'profile':profile(qpaths[c.identity].parent/'cubic_1280/selected.npz' if c.identity in qpaths else folder/'selected_refined.npz',f'{runid}-{source}-{rad}',wall)}
        series.append({'id':f'{runid}-{source}-{rad}','source':source,'rmin':rad,'label':label,
                       'run':runid,'stages':milestones,'final':final,'statistics':stats})
        audit.append({'id':series[-1]['id'],**stats,'displayed':len(milestones),'final':final is not None})
        print(audit[-1],flush=True)
    payload={'schema':'d3d-progression.v1','created':'2026-11-04','reference':reference,
             'baseline':{'wall':fixture['wall'],'lcfs':fixture['lcfs'],'coarse':base_coarse,'refined':base_fine,'profiles':base_profiles},
             'surfaces':fixture['surfaces'],'series':series,
             'selection':'original running-best records, then current geometry screen; ratio decrease > 1e-9',
             'heat':'saved profile only; intermediate recorded peaks do not imply a stored profile'}
    save(dist/'data.json',payload)
    # Hash again after export; never mutate historical evidence.
    for p,h in evidence.items():
        if sha(p)!=h:raise ValueError('input changed: '+p)
    code_files=[Path(__file__),REPO/'verification/test_campaign_consistency.py',*(REPO/'d3d_optimization').glob('*.py')]
    save(HERE/'data-audit.json',{'studies':audit,'input_sha256':evidence,'seconds':time.perf_counter()-start,
         'code_sha256':{str(p):sha(p) for p in code_files},'inputs_unchanged_verified':True,
         'old_viewer_research_data_used':False,'new_thermal_solves':0,'new_optimizations':0})
    print('DONE',round(time.perf_counter()-start,2),'s',len(series),'series',sum(len(s['stages']) for s in series),'stages',flush=True)

if __name__=='__main__':run()
