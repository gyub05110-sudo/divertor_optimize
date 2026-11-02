"""Bounded continuous-path validation study; leaves fast optimization unchanged."""
import argparse
import copy
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import BSpline,bisplev
from d3d_optimization.continuous_path_bounds import trace
from d3d_optimization.continuous_transport_certificate import ValidatedPaths,certify_box,cover,replay,replay_cover,first_hits
from d3d_optimization.enclosure_math import I,D,DOWN,interval
from d3d_optimization.geometry_transport_bounds import Surface
from verification.run_optimality_bounds import CAMPAIGN,REPO,digest,save,snapshot_code

PREVIOUS=Path('D:/DivertorOptimization/studies/geometry-transport-enclosure-20261102-v3')


def numeric_audit(recorded):
    """Independent DOP853 chart propagation; finite tolerance, not a proof."""
    m=recorded['inputs']['field'];kr,kz=m['degree_r'],m['degree_z']
    tr,tz=np.asarray(m['r_knots'],float),np.asarray(m['z_knots'],float)
    tck=(tz,tr,np.asarray(m['coefficients'],float).ravel(),kz,kr)
    if min(kr,kz)==1:
        # FITPACK bisplev excludes derivative order == degree. Explicit
        # SciPy basis derivatives support the valid linear-field case.
        coeff=np.asarray(m['coefficients'],float)
        br=BSpline(tr,np.eye(len(tr)-kr-1),kr);bz=BSpline(tz,np.eye(len(tz)-kz-1),kz)
        dr=br.derivative();dz=bz.derivative()
    p=np.asarray(recorded['inputs']['guide'][0],float);K=float(m['toroidal_K'])
    def field(p):
        if not tr[kr]<=p[0]<=tr[-kr-1] or not tz[kz]<=p[1]<=tz[-kz-1]:
            raise AssertionError('numerical audit left spline domain')
        if min(kr,kz)==1:g=np.array([bz(p[1])@coeff@dr(p[0]),dz(p[1])@coeff@br(p[0])])
        else:g=np.array([bisplev(p[1],p[0],tck,dy=1),bisplev(p[1],p[0],tck,dx=1)])
        return g,float(np.hypot(np.linalg.norm(g),K)/p[0])
    _,b0=field(p);grad_sep=float(recorded['inputs']['grad_sep'][0]);L=J=0.;checks=0;max_excess=0.
    def check(x,bound):
        nonlocal checks,max_excess
        lo,hi=map(float,bound);excess=max(lo-x,x-hi,0.)
        max_excess=max(max_excess,excess)
        if excess>2e-8*max(1.,abs(x)):raise AssertionError((x,lo,hi,excess))
        checks+=1
    for cell in recorded['cells']:
        axis=cell['axis'];other=1-axis;sign=cell['sign'];a=p[axis]
        target=float(cell['end'][axis][0]);h=sign*(target-a)
        if h<=0:raise AssertionError('independent chart no longer forward')
        def rhs(t,state):
            point=np.empty(2);point[axis]=a+sign*t;point[other]=state[0]
            g,B=field(point);metric=np.hypot(np.linalg.norm(g),K)/abs(g[other])
            return [-sign*g[axis]/g[other],metric,B/b0*metric]
        sol=solve_ivp(rhs,(0,h),[p[other],L,J],method='DOP853',rtol=3e-13,atol=1e-14,dense_output=True)
        if not sol.success:raise AssertionError(sol.message)
        for t in [0,h/2,h]:
            state=sol.sol(t);q=np.empty(2);q[axis]=a+sign*t;q[other]=state[0]
            for k in range(2):check(q[k],cell['tube'][k])
            check(state[1],[cell['length_start'][0],cell['length_end'][1]])
            check(state[2],[cell['integral_start'][0],cell['integral_end'][1]])
            g,B=field(q);check(np.linalg.norm(g)/grad_sep*b0/B,cell['mapping'])
        p=q;L,J=state[1:]
        for k in range(2):check(p[k],cell['end'][k])
        check(L,cell['length_end']);check(J,cell['integral_end'])
    return {'passed':True,'checks':checks,'max_raw_excess':max_excess,'is_proof':False,
            'tolerance':'2e-8 * max(1,abs(value))','rtol':3e-13,'atol':1e-14,
            'method':'independent floating DOP853; no invariant projection'}


def refined_guide(path,selected,factor=8):
    """Propose more coordinate planes near possible hits; never a proof step."""
    if not isinstance(factor,int) or not 2<=factor<=32:raise ValueError('refinement factor 2..32 required')
    def midpoint(bounds):
        x=interval(bounds);return DOWN.divide(DOWN.add(x.lo,x.hi),D(2))
    guide=[copy.deepcopy(path['inputs']['guide'][0])]
    for index,cell in enumerate(path['cells']):
        a=[midpoint(x) for x in cell['start']];b=[midpoint(x) for x in cell['end']]
        n=factor if index in selected else 1
        for j in range(1,n+1):
            guide.append([str(DOWN.add(x,DOWN.divide(DOWN.multiply(DOWN.subtract(y,x),D(j)),D(n)))) for x,y in zip(a,b)])
    return guide


def run(output,stride=20,splits=2,reuse=None):
    if not isinstance(stride,int) or not 1<=stride<=100:raise ValueError('guide stride 1..100 required')
    if not isinstance(splits,int) or not 0<=splits<=8:raise ValueError('study split budget 0..8 required')
    out=Path(output).resolve()
    if out.drive.lower()!='d:' or out.is_relative_to(CAMPAIGN.resolve()) or out.is_relative_to(PREVIOUS.resolve()):
        raise ValueError('new D: study outside preserved studies required')
    out.mkdir(parents=True,exist_ok=False);started=perf_counter()
    settings=json.loads((PREVIOUS/'settings.json').read_text(encoding='utf-8'))
    protected=dict(settings['protected_sha256'])
    protected.update({str(PREVIOUS/p):h for p,h in json.loads((PREVIOUS/'artifact_sha256.json').read_text()).items()})
    if reuse:
        reuse=Path(reuse).resolve()
        for source in ['IMP','OMP']:
            for p in (reuse/source).glob('path_???.json'):protected[str(p)]=digest(p)
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError('protected input mismatch: '+p)
    save(out/'settings.json',{'stride':stride,'splits':splits,'reuse':str(reuse) if reuse else None,'protected_sha256':protected,
                             'code_sha256':snapshot_code(REPO,out/'source')})
    summary=[]
    for source in ['IMP','OMP']:
        source_out=out/source;source_out.mkdir()
        model=json.loads((PREVIOUS/f'{source}_wide.json').read_text())['model']
        save(source_out/'model.json',model);records=[];timings=[];audits=[]
        surfaces=[Surface(s,list(map(interval,model['domains'])),8) for s in model['surfaces']]
        for i,donor in enumerate(model['donors']):
            guide=donor['path'][::stride]
            if guide[-1]!=donor['path'][-1]:guide.append(donor['path'][-1])
            t=perf_counter();cached=reuse/source/f'path_{i:03d}.json' if reuse else None
            reused=bool(cached and cached.exists());rounds=[]
            if reused:
                path=json.loads(cached.read_text())
                if path['inputs']['field']!=model['field'] or path['inputs']['guide']!=guide:
                    raise ValueError('reused initial path input mismatch')
            else:path=trace(model['field'],guide,model['grad_sep'],max_cells=2048,max_depth=8)
            for attempt in range(2):
                hit=first_hits(path,surfaces)
                if hit['guaranteed_hit']:break
                indices={h['cell'] for h in hit['possible_hits']}
                selected={j for c in indices for j in (c-1,c,c+1) if 0<=j<len(path['cells'])}
                if not selected:break
                save(source_out/f'path_{i:03d}_before_refinement_{attempt}.json',path)
                guide=refined_guide(path,selected)
                rounds.append({'selected_cells':sorted(selected),'factor':8})
                path=trace(model['field'],guide,model['grad_sep'],max_cells=2048,max_depth=8)
            elapsed=perf_counter()-t;records.append(path);save(source_out/f'path_{i:03d}.json',path)
            audit=numeric_audit(path);audits.append({'donor':i,**audit})
            timings.append({'donor':i,'seconds':elapsed,'cells':len(path['cells']),
                            'initial_cache_reused':reused,'refinements':rounds,
                            'complete':path['complete'],'unresolved':path['unresolved']})
            save(source_out/'path_timings.json',timings);save(source_out/'numeric_audit.json',audits)
            print(source,i,'cells',len(path['cells']),'complete',path['complete'],f'{elapsed:.2f}s',flush=True)
        t=perf_counter();cache=ValidatedPaths(model,records);replay_seconds=perf_counter()-t
        for label in ['local','wide']:
            box=json.loads((PREVIOUS/f'{source}_{label}.json').read_text())['box']
            t=perf_counter();certificate=certify_box(model,box,cache,area_panels=8,thermal_panels=64)
            elapsed=perf_counter()-t;save(source_out/f'{label}.json',certificate)
            replay(certificate,cache)
            row={'case':source+'_'+label,'continuous_transport_proven':certificate['transport_parameter_enclosure_proven'],
                 'resolved_donors':sum(d['status']=='enclosed' for d in certificate['donors']),
                 'donors':len(model['donors']),'reused_box_seconds':elapsed,
                 'path_preparation_seconds':sum(t['seconds'] for t in timings),'path_replay_seconds':replay_seconds,
                 'peak_bracket_MW_m2':certificate.get('peak_bracket_MW_m2'),'unresolved':certificate['unresolved']}
            summary.append(row);save(out/'summary.json',summary);print(row,flush=True)
        full=copy.deepcopy(model)
        full['domains']=[['0','1'] if source=='IMP' and j==8 else ['-0.1','0.1'] for j in range(len(model['domains']))]
        t=perf_counter();tree=cover(full,cache,splits=splits,area_panels=4,thermal_panels=16)
        save(source_out/'full_cover.json',tree)
        summary.append({'case':source+'_full_domain','seconds':perf_counter()-t,'leaves':len(tree['leaves']),
                        'resolved_leaves':sum(e['certificate']['transport_parameter_enclosure_proven'] for e in tree['leaves']),
                        'split_history':tree['history'],
                        'unresolved':[{ 'id':e['id'],'reasons':e['certificate']['unresolved']} for e in tree['leaves']]})
        save(out/'summary.json',summary);replay_cover(tree,cache)
        save(source_out/'replay.json',{'paths':len(records),'geometry':['local','wide','full_cover'],'passed':True})
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError('protected input changed: '+p)
    save(out/'completion.json',{'seconds':perf_counter()-started,'protected_files':len(protected),
                               'geometry_global_optimum_proven':False,'certified_geometry_optimality_gap':None})
    save(out/'artifact_sha256.json',{p.relative_to(out).as_posix():digest(p) for p in out.rglob('*') if p.is_file()})
    print('Completed study',perf_counter()-started,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    parser.add_argument('--stride',type=int,default=20);parser.add_argument('--splits',type=int,default=2)
    parser.add_argument('--reuse')
    args=parser.parse_args();run(args.output,args.stride,args.splits,args.reuse)
