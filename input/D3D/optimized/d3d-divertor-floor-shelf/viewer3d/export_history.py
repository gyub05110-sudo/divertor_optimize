"""Register finished multisurface runs for equilibrium-specific history inspection."""
import argparse,csv,hashlib,json,importlib,sys,types
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Polygon,LineString,Point
from history_data import core_surfaces,normalize_psi,truncate_path,select_history_records

p=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--case-dir',type=Path,default=p.parent)
parser.add_argument('--code-dir',type=Path,help='archived math directory used by this legacy run')
args=parser.parse_args();base=args.case_dir.resolve()
# The original code moved with the results. Import its preserved math under an
# isolated namespace; do not depend on, or modify, the active optimization package.
archives=sorted(p.parent.parent.glob('migration-*/d3d_multisurface.py'))
code_dir=args.code_dir.resolve() if args.code_dir else (archives[-1].parent if archives else None)
if code_dir:
    package=types.ModuleType('_d3d_viewer_math');package.__path__=[str(code_dir)]
    sys.modules['_d3d_viewer_math']=package
    case_code=importlib.import_module('_d3d_viewer_math.d3d_frozen_case')
    MultiSurfaceEvaluator=importlib.import_module('_d3d_viewer_math.d3d_multisurface').MultiSurfaceEvaluator
    load_case,separatrix,pretrace,sha256,WIDTHS_MM=[getattr(case_code,n) for n in ['load_case','separatrix','pretrace','sha256','WIDTHS_MM']]
else:
    from divertor_optimization.d3d_frozen_case import load_case,separatrix,pretrace,sha256,WIDTHS_MM
    from divertor_optimization.d3d_multisurface import MultiSurfaceEvaluator
m=json.loads((base/'manifest.json').read_text(encoding='utf-8'))
assert m['numerically_verified_feasible']
assert sha256(m['shot_path'])==m['shot_sha256']
assert sha256(base/'DIIID_geom.json')==m['config']['geometry_sha256']
assert [s['name'] for s in m['config']['surfaces']]==['floor','shelf']
field,lcfs,axis=load_case(m['shot_path'],m['time_ms'],'cubic')
r_sep,xp,psi_sep=separatrix(field,lcfs,axis)
g=json.loads((base/'DIIID_geom.json').read_text());wall=np.array(g['limiter']);wall=np.vstack([wall,wall[0]])
edges,paths,inside,_,drift=pretrace(field,r_sep,axis[1],wall,320,.002)
e=MultiSurfaceEvaluator(field,lcfs,wall,paths,inside,edges,r_sep,axis[1],m['config']['surfaces'],
                        vessel_solids=shapely.union_all([Polygon(v['pts']) for v in g['vv']]))
core,normalization=core_surfaces(field,axis,psi_sep)
records=list(map(json.loads,(base/'evaluations.jsonl').read_text(encoding='utf-8').splitlines()))
safe={int(row['case']):row['vessel_safe_geometry_valid']=='True' for row in csv.DictReader((base/'vessel_geometry_audit.csv').open())}
verification=json.loads((base/'candidate_verification.json').read_text())
checked={row['case']:row['passed'] for row in verification}
improving=select_history_records(records,safe)
run_hash=hashlib.sha256(base.name.encode()).hexdigest()[:8]
case_id=m['shot_sha256'][:12]+f"-{m['time_ms']:g}ms-"+run_hash
directory=p/'cases'/case_id;stages_dir=directory/'stages';stages_dir.mkdir(parents=True,exist_ok=True)
sol_indices=[]
centers=(edges[:-1]+edges[1:])/2
for offset in [1,2,3,5,10,20,40,60,80]:
    idx=int(np.argmin(abs(centers-offset/1000)))
    if inside[idx] and idx not in sol_indices:sol_indices.append(idx)

def pack(result):
    line=LineString(result['wall']);heat={};summary={}
    for mm,h,bpeak in zip(WIDTHS_MM,result['maps'],e.baseline_peaks):
        samples=[]
        for i,seg in enumerate(h.segment):
            if seg<0:continue
            component='fixed_surface'
            for name,ids in result['surface_segments'].items():
                if seg in ids:component=name
            samples.append([line.project(Point(h.points[i])),float(h.q[i]),component])
        samples.sort();heat[str(mm)]=samples;idx=int(np.argmax(h.q))
        summary[str(mm)]={'all_wall_peak_q_over_P':float(max(h.q)),
                          'all_wall_peak_ratio':float(max(h.q)/bpeak),'peak_R_m':float(h.points[idx,0]),
                          'peak_Z_m':float(h.points[idx,1]),'floor_power':result['surface_power']['floor'][list(WIDTHS_MM).index(mm)],
                          'shelf_power':result['surface_power']['shelf'][list(WIDTHS_MM).index(mm)]}
    sol=[]
    for idx in sol_indices:
        h=result['maps'][0]
        if h.segment[idx]<0:continue
        curve=truncate_path(paths[idx],h.points[idx])
        pn=normalize_psi(field.values(curve),normalization['psi_axis'],psi_sep)
        if max(pn)-min(pn)>.001:raise ValueError('SOL psi drift exceeds display tolerance')
        sol.append({'source_offset_mm':float(centers[idx]*1000),'psi_normalized':float(pn[0]),'points':curve,
                    'hit_original_segment':int(result['original_segment'][h.segment[idx]])})
    return {'wall':result['wall'].tolist(),'floor':result['targets']['floor'].tolist(),
            'shelf':result['targets']['shelf'].tolist(),'heat':heat,'summary':summary,'sol':sol}

descriptors=[];baseline=None;final=None
final_row=next(row for row in records if row['case']==m['retained_case'])
all_rows=[(row,'baseline' if j==0 else 'improvement') for j,row in enumerate(improving)]+[(final_row,'verified_final')]
for j,(row,kind) in enumerate(all_rows):
    r=e.evaluate(row['coefficients_m']);assert r['geometry_valid']
    if kind!='baseline':assert r['feasible']
    computed=max(float(max(h.q)/b) for h,b in zip(r['maps'],e.baseline_peaks))
    assert abs(computed-row['worst_all_wall_peak_ratio'])<1e-10
    key='verified_final' if kind=='verified_final' else f"case_{row['case']}"
    descriptor={'id':key,'order':j,'kind':kind,'case':row['case'],'search_stage':row['stage'],
                'score':computed,'verified':kind=='verified_final','constraint_feasible':bool(r['feasible']),
                'refinement_rejected':checked.get(row['case']) is False,
                'label':('기준 형상' if kind=='baseline' else '최종 검증 형상' if kind=='verified_final'
                         else f"{'제약 충족 개선' if not improving[0]['feasible'] else '개선'} {j:02d} · case {row['case']} · {row['stage']}"),
                'detail_url':f'cases/{case_id}/stages/{key}.json'}
    detail=pack(r);detail['stage']=descriptor
    (stages_dir/(key+'.json')).write_text(json.dumps(detail,separators=(',',':'),ensure_ascii=False,allow_nan=False),encoding='utf-8')
    descriptors.append(descriptor)
    if kind=='baseline':baseline=detail
    if kind=='verified_final':final=detail
    print('export stage',j+1,'/',len(all_rows),key,flush=True)
bundle={'schema_version':2,'id':case_id,'time_ms':m['time_ms'],'retained_case':m['retained_case'],
        'shot_name':Path(m['shot_path']).name,'shot_sha256':m['shot_sha256'],'run_name':base.name,
        'widths':list(map(int,WIDTHS_MM)),'lcfs':lcfs.tolist(),'flux_surfaces':core,'normalization':normalization,
        'sol_scope':'lower outboard poloidal SOL branch: upstream midplane to first material hit',
        'models':{'baseline':baseline,'candidate':final},'stages':descriptors,
        'source_geometry_sha256':m['config']['geometry_sha256'],'max_path_psi_drift':drift}
(directory/'case.json').write_text(json.dumps(bundle,separators=(',',':'),ensure_ascii=False,allow_nan=False),encoding='utf-8')
catalog_path=p/'catalog.json'
catalog=json.loads(catalog_path.read_text(encoding='utf-8')) if catalog_path.exists() else {'schema_version':2,'cases':[]}
catalog['cases']=[c for c in catalog['cases'] if c['id']!=case_id]+[{'id':case_id,'shot_name':bundle['shot_name'],
    'shot_sha256':m['shot_sha256'],'time_ms':m['time_ms'],'run_name':base.name,'stage_count':len(descriptors),
    'bundle_url':f'cases/{case_id}/case.json'}]
catalog_path.write_text(json.dumps(catalog,indent=2,ensure_ascii=False),encoding='utf-8')
(p/'history_export_checks.json').write_text(json.dumps({'case_id':case_id,'stage_count':len(descriptors),
     'core_levels':[c['psi_normalized'] for c in core],'normalization':normalization,
     'max_contour_residual':max(c['max_psiN_residual'] for c in core),'first_hit_sol_paths_per_stage':len(final['sol']),
     'recorded_scores_verified':True,'sol_truncated_to_material':True},indent=2),encoding='utf-8')
print('Registered',case_id,'stages',len(descriptors),flush=True)
