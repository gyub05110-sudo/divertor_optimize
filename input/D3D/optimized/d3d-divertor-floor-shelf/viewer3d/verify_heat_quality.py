"""Additional tracing and independent surface quadrature for registered results.

Writes separate viewer quality artifacts, never modifies original run records.
"""
import argparse,hashlib,importlib,json,sys,types,time
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import LineString,Point,Polygon
from scipy.optimize import brentq
from heat_quality import refinement_cells,connection_status,compare_levels
from quality_result import complete_report
from source_domain import physical_source_limit,require_physical_extension

P=Path(__file__).resolve().parent
WIDTHS=[1,2,3,5,10]
def fingerprint(model):
    return hashlib.sha256(json.dumps(model['wall'],separators=(',',':')).encode()).hexdigest()
def save(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case-id');parser.add_argument('--no-register',action='store_true');parser.add_argument('--source-max-mm',type=float,help='Explicit source range; must end inside public wall and EFIT grid');args=parser.parse_args()
    catalog=json.loads((P/'catalog.json').read_text(encoding='utf-8'))
    entry=next(c for c in catalog['cases'] if not args.case_id or c['id']==args.case_id)
    bundle=json.loads((P/entry['bundle_url']).read_text(encoding='utf-8'))
    case_base=Path(bundle['case_directory']) if bundle.get('case_directory') else P.parent.parent/entry['run_name']
    manifest=json.loads((case_base/'manifest.json').read_text(encoding='utf-8'))
    diffusion_ratio=bundle.get('sigma_over_lambda_q',0)
    if diffusion_ratio:
        sys.path.insert(0,str(P.parents[4]))
        from d3d_optimization import case, frozen_mapping as mapping, frozen_field as field_module
        from d3d_optimization.transverse_diffusion import ReflectingProfile
    else:
        archives=sorted(P.parent.parent.glob('migration-*/d3d_multisurface.py'))
        if not archives:raise RuntimeError('Preserved case math required')
        code_dir=case_base/'case_math' if (case_base/'case_math').is_dir() else archives[-1].parent
        package=types.ModuleType('_heat_quality_case');package.__path__=[str(code_dir)]
        sys.modules[package.__name__]=package
        case=importlib.import_module(package.__name__+'.d3d_frozen_case')
        mapping=importlib.import_module(package.__name__+'.frozen_mapping')
        field_module=importlib.import_module(package.__name__+'.frozen_field')
    if case.sha256(manifest['shot_path'])!=bundle['shot_sha256']:raise ValueError('Shot hash mismatch')
    field,lcfs,axis=case.load_case(manifest['shot_path'],bundle['time_ms'],'cubic')
    r_sep,_,_=case.separatrix(field,lcfs,axis)
    models=bundle['models'];walls={k:np.asarray(v['wall']) for k,v in models.items()}
    source_domain={k:physical_source_limit(w,r_sep,float(axis[1]),float(field.r[-1])) for k,w in walls.items()}
    physical_max=min(d['physical_max_offset_m'] for d in source_domain.values())
    source_max=.1 if args.source_max_mm is None else args.source_max_mm/1000
    if diffusion_ratio and source_max!=.1:raise ValueError('diffusion audit must retain the searched 100 mm source band')
    profiles={mm:ReflectingProfile(field,r_sep,axis[1],mm/1000,diffusion_ratio*mm/1000) for mm in WIDTHS} if diffusion_ratio else None
    if args.source_max_mm is not None:require_physical_extension(source_max,physical_max)
    # The original 100 mm source is retained only for legacy comparison. Its
    # outside-wall starts remain explicitly unmapped; new ranges must be physical.
    lines={k:LineString(v) for k,v in walls.items()};polygons={k:Polygon(v) for k,v in walls.items()}
    tags={}
    for k,model in models.items():
        middle=(walls[k][:-1]+walls[k][1:])/2
        tags[k]=np.array(['fixed_surface']*len(middle),dtype=object)
        for name in ['floor','shelf']:
            target=LineString(model[name])
            tags[k][[target.distance(Point(p))<1e-8 for p in middle]]=name
    def trace(x,step=.001):
        starts=np.column_stack([r_sep+np.asarray(x),np.full(len(x),axis[1])])
        paths=[LineString(p) for p in mapping.trace_paths(field,starts,step=step)]
        return starts,paths
    def evaluate(edges):
        centers=(edges[:-1]+edges[1:])/2;starts,paths=trace(centers)
        sources=[profiles[mm].heat_source(edges) if profiles else field_module.exponential_source(field,r_sep,axis[1],edges,mm/1000) for mm in WIDTHS]
        results={};brief={}
        for k,wall in walls.items():
            h=mapping.first_hits(paths,wall)
            inside=shapely.contains(polygons[k],shapely.points(starts))
            h.segment[~inside]=-1;h.points[~inside]=np.nan;h.distance[~inside]=np.inf
            maps=[mapping.map_heat(field,h,s.weights,s.h,wall) for s in sources]
            arcs=np.array([lines[k].project(Point(p)) if i>=0 else np.nan for p,i in zip(h.points,h.segment)])
            components=[tags[k][i] if i>=0 else 'unmapped' for i in h.segment]
            results[k]={'hits':h,'maps':maps,'arcs':arcs,'components':components}
            imax=[int(np.argmax(m.q)) for m in maps]
            brief[k]={'peaks':[float(m.q[i]) for m,i in zip(maps,imax)],
                      'locations':[h.points[i].tolist() for i in imax],
                      'unmapped_power':[m.unmapped_power for m in maps],
                      'source_accounting_error':[m.balance_error for m in maps],
                      'component_power':{name:[float(np.sum(m.segment_power[tags[k]==name])) for m in maps]
                                         for name in ['floor','shelf','fixed_surface']},
                      'shelf_min_angle_deg':float(np.min(maps[0].angles[np.array(components)=='shelf']))}
        ratios=(np.array(brief['candidate']['peaks'])/np.array(brief['baseline']['peaks'])).tolist()
        for k in brief:brief[k]['ratios']=ratios
        return results,brief,sources
    levels=[];previous=None
    for count in [320,640,1280]:
        edges=source_max*np.linspace(0,1,count+1)**2
        results,brief,sources=evaluate(edges)
        checks={k:compare_levels(previous[k],brief[k]) for k in brief} if previous else None
        levels.append({'kind':'uniform','source_cells':count,'models':brief,'comparison':checks})
        previous=brief;print('uniform',count,'ratios',brief['candidate']['ratios'],flush=True)
    for iteration in range(3):
        ids=set()
        for r in results.values():
            ids.update(refinement_cells(edges,r['arcs'],r['components'],np.array([m.q for m in r['maps']])))
        if not ids:break
        edges=np.sort(np.r_[edges,[(edges[i]+edges[i+1])/2 for i in ids]])
        results,brief,sources=evaluate(edges)
        checks={k:compare_levels(previous[k],brief[k]) for k in brief}
        levels.append({'kind':'adaptive','source_cells':len(edges)-1,'models':brief,'comparison':checks})
        print('adaptive',iteration+1,len(edges)-1,'converged',all(c['passed'] for c in checks.values()),flush=True)
        previous=brief
        if iteration>=1 and all(c['passed'] for c in checks.values()):break
    # Fine surface quadrature is independent of source power assignment. Each
    # midpoint is psi-inverted, then forward-traced to certify branch and first hit.
    up_x=np.linspace(0,source_max,2001);up_psi=field.values(np.column_stack([r_sep+up_x,np.full(len(up_x),axis[1])]))
    if not np.all(np.diff(up_psi)>0) and not np.all(np.diff(up_psi)<0):raise ValueError('Nonmonotonic upstream psi')
    psi_min,psi_max=float(min(up_psi)),float(max(up_psi))
    surface_results={};out_models={};surface_peak_checks={}
    for k,model in models.items():
        patches=[]
        for name in ['floor','shelf']:
            target=np.asarray(model[name])
            for a,b in zip(target[:-1],target[1:]):
                n=max(1,int(np.ceil(np.linalg.norm(b-a)/.0005)))
                for j in range(n):patches.append({'a':(a+(b-a)*j/n).tolist(),'b':(a+(b-a)*(j+1)/n).tolist(),'component':name})
        def inspect(items):
            points=np.array([(np.array(p['a'])+p['b'])/2 for p in items]);psi=field.values(points)
            in_scope=(psi>=psi_min)&(psi<=psi_max);xs=np.full(len(points),np.nan)
            for i in np.flatnonzero(in_scope):
                xs[i]=brentq(lambda x:float(field.values([[r_sep+x,axis[1]]])[0]-psi[i]),0,source_max,xtol=1e-12)
            valid=np.flatnonzero(in_scope);starts,paths=trace(xs[valid],.0005)
            hits=mapping.first_hits(paths,walls[k]);inside=shapely.contains(polygons[k],shapely.points(starts))
            cert={}
            for j,i in enumerate(valid):
                reached=bool(inside[j] and paths[j].distance(Point(points[i]))<=.0001)
                first=bool(hits.segment[j]>=0 and np.linalg.norm(hits.points[j]-points[i])<=.0001)
                cert[i]=(reached,first)
            for i,(p,point) in enumerate(zip(items,points)):
                reached,first=cert.get(i,(False,False));status=connection_status(bool(in_scope[i]),reached,first)
                p['status']=status;p['q']={};p['source_offset_mm']=float(xs[i]*1000) if in_scope[i] else None
                if status=='outside_source_scope':
                    sep_side=bool(psi[i]<psi_min) if up_psi[-1]>up_psi[0] else bool(psi[i]>psi_max)
                    p['scope_reason']='separatrix_side' if sep_side else 'beyond_source_end'
                d=np.asarray(p['b'])-p['a'];normal=np.array([-d[1],d[0]])/np.linalg.norm(d)
                p['area_m2']=float(2*np.pi*point[0]*np.linalg.norm(d))
                p['angle_deg']=float(field.incidence([point],[normal])[0]) if status=='direct_hit' else None
                for mm,src in zip(WIDTHS,sources):
                    if status!='direct_hit':p['q'][str(mm)]=None;continue
                    bup=np.linalg.norm(field.magnetic([[r_sep+xs[i],axis[1]]])[0])
                    h=float(profiles[mm].h_at([xs[i]])[0]) if profiles else np.exp(-xs[i]/(mm/1000))/bup/src.normalization
                    p['q'][str(mm)]=float(h*abs(np.dot(field.magnetic([point])[0,:2],normal)))
            return items
        patches=inspect(patches)
        spatial_checks=[]
        for refinement in range(3):
            integrals={name:[sum(p['area_m2']*(p['q'][str(mm)] or 0) for p in patches if p['component']==name) for mm in WIDTHS] for name in ['floor','shelf']}
            residual={name:((np.asarray(integrals[name])-brief[k]['component_power'][name])/np.maximum(brief[k]['component_power'][name],1e-5)).tolist() for name in integrals}
            spatial_peaks=[max((p['q'][str(mm)] or 0) for p in patches) for mm in WIDTHS]
            spatial_checks.append({'patch_count':len(patches),'integrals':integrals,'relative_residuals':residual,'sampled_peaks':spatial_peaks})
            print('surface',k,refinement,len(patches),'max integral residual',max(abs(x) for a in residual.values() for x in a),flush=True)
            if refinement==2:break
            poor={name for name,r in residual.items() if max(abs(x) for x in r)>.01}
            peak={str(mm):max((p['q'][str(mm)] or 0) for p in patches) for mm in WIDTHS}
            add=[];keep=[]
            for i,p in enumerate(patches):
                boundary=any(p['component']==patches[j]['component'] and p['status']!=patches[j]['status'] for j in [i-1,i+1] if 0<=j<len(patches))
                nearpeak=any(peak[str(mm)]>0 and (p['q'][str(mm)] or 0)>.8*peak[str(mm)] for mm in WIDTHS)
                if (p['component'] in poor or boundary or nearpeak) and np.linalg.norm(np.array(p['b'])-p['a'])>.0001:
                    mid=((np.array(p['a'])+p['b'])/2).tolist()
                    add.extend([{'a':p['a'],'b':mid,'component':p['component']},{'a':mid,'b':p['b'],'component':p['component']}])
                else:keep.append(p)
            patches=keep+inspect(add)
            patches.sort(key=lambda p:lines[k].project(Point(p['a'])))
        audit_peaks=[];audit_positions=[]
        for mm,index in zip(WIDTHS,range(len(WIDTHS))):
            best=max(patches,key=lambda p:p['q'][str(mm)] or 0)
            if (best['q'][str(mm)] or 0)>brief[k]['peaks'][index]:
                audit_peaks.append(best['q'][str(mm)]);audit_positions.append(((np.array(best['a'])+best['b'])/2).tolist())
            else:audit_peaks.append(brief[k]['peaks'][index]);audit_positions.append(brief[k]['locations'][index])
        status_lengths={s:sum(np.linalg.norm(np.array(p['b'])-p['a']) for p in patches if p['status']==s) for s in ['direct_hit','blocked_first_hit','outside_source_scope','trace_unresolved']}
        surface_peak_checks[k]={'relative_increase':(np.asarray(audit_peaks)/brief[k]['peaks']-1).tolist(),
                               'passed':bool(np.max(np.asarray(audit_peaks)/brief[k]['peaks']-1)<=.01)}
        surface_results[k]={'quadrature_levels':spatial_checks,'status_length_m':status_lengths}
        # Retain fine forward samples for other walls. Audited target patches
        # supply their own colors and do not interpolate across missing coverage.
        heat={};summary={}
        for index,mm in enumerate(WIDTHS):
            r=results[k];h=r['maps'][index]
            heat[str(mm)]=sorted([[float(r['arcs'][i]),float(h.q[i]),r['components'][i]] for i in range(len(h.q)) if h.segment[i]>=0])
            summary[str(mm)]={**model['summary'][str(mm)],'all_wall_peak_q_over_P':audit_peaks[index],
                              'peak_R_m':audit_positions[index][0],'peak_Z_m':audit_positions[index][1]}
        out_models[k]={'shape_sha256':fingerprint(model),'reference_wall':model['wall'],'patches':patches,'heat':heat,'summary':summary}
    last_checks=levels[-1]['comparison']
    integral_ok=all(max(abs(x) for r in s['quadrature_levels'][-1]['relative_residuals'].values() for x in r)<=.01 for s in surface_results.values())
    source_ok=all(max(b['unmapped_power'])<=.001 and max(b['source_accounting_error'])<=1e-4 for b in brief.values())
    unresolved=sum(s['status_length_m']['trace_unresolved'] for s in surface_results.values())
    ratios=[out_models['candidate']['summary'][str(mm)]['all_wall_peak_q_over_P']/out_models['baseline']['summary'][str(mm)]['all_wall_peak_q_over_P'] for mm in WIDTHS]
    for k in out_models:
        for mm,ratio in zip(WIDTHS,ratios):out_models[k]['summary'][str(mm)]['all_wall_peak_ratio']=ratio if k=='candidate' else 1.
    report={'schema_version':1,'case_id':bundle['id'],'shot_sha256':bundle['shot_sha256'],
            'configured_source_max_mm':source_max*1000,'source_domain_bounds':source_domain,
            'geometry_sha256':bundle['source_geometry_sha256'],'scope':'baseline and retained final; fixed stored polygon geometry, lower outer SOL only',
            'criteria':{'peak_relative':.01,'location_m':.002,'ratio_relative':.01,'integral_relative':.01,'unmapped_power':.001},
            'levels':levels,'surfaces':surface_results,'final_ratios':ratios,
            'forward_converged':all(c['passed'] for c in last_checks.values()),'surface_integral_passed':integral_ok,
            'source_accounting_passed':source_ok,'unresolved_target_length_m':unresolved,
            'surface_peak_crosscheck':surface_peak_checks,
            'verified':bool(all(c['passed'] for c in last_checks.values()) and integral_ok and source_ok and unresolved<1e-6 and all(v['passed'] for v in surface_peak_checks.values())),
            'limitations':['fixed polygon geometry resolution not refined',f'finite outer SOL source interval 0–{source_max*1000:g} mm',
                           'reduced upstream-equivalent diffusion; no volume transport, radiation or measured power calibration' if profiles else 'no radiation, cross-field transport or measured power calibration'],
            'models':out_models}
    if profiles:report.update(sigma_over_lambda_q=diffusion_ratio,source_model=manifest['diffusion'],source_profile_power_errors={str(mm):p.mass_error for mm,p in profiles.items()})
    directory=P/'quality';directory.mkdir(exist_ok=True)
    complete_report(report,bundle)
    save(directory/(bundle['id']+'.json'),report)
    if not args.no_register:
        quality_catalog=json.loads((P/'quality_catalog.json').read_text(encoding='utf-8')) if (P/'quality_catalog.json').exists() else {'cases':{}}
        quality_catalog['cases'][bundle['id']]='quality/'+bundle['id']+'.json'
        save(P/'quality_catalog.json',quality_catalog)
    print('quality saved; verified=',report['verified'],'ratios=',ratios,flush=True)

if __name__=='__main__':main()
