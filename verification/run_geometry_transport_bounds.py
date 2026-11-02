"""Bounded IMP/OMP inclusion study; never resumes the expired campaign."""
import argparse
from dataclasses import asdict
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
from time import perf_counter
import numpy as np
from d3d_optimization.campaign_geometry import Candidate,evaluator,evaluate_candidate
from d3d_optimization.imp_tangent import junction_bounds
from d3d_optimization.enclosure_math import I,interval,record
from d3d_optimization.geometry_transport_certificate import certify_box,replay,cover
from d3d_optimization.transport_peak_bounds import certify as certify_kernel
from verification.run_optimality_bounds import CAMPAIGN,REPO,digest,save,snapshot_code


def freeze_model(e,candidate,halfwidth,junction_halfwidth):
    coefficients=candidate.bind(e);curves=e.curves(coefficients)
    active=['imp_sloped','imp_vertical'] if candidate.source=='IMP' else ['omp_entry','shelf','omp_connector','floor']
    by_segment={s['segment_index']:s['name'] for s in e.surfaces}
    surfaces=[]
    for j,(a,b) in enumerate(zip(e.wall,e.wall[1:])):
        if j in by_segment:
            name=by_segment[j];q=curves[name];controls=q.spline.c.tolist()
            if name in active:
                for k in range(4):controls[k+2]={'variable':4*active.index(name)+k}
            if candidate.source=='IMP' and name in active:
                lo,hi=junction_bounds(curves['imp_sloped'].v,curves['imp_vertical'].v)
                at=6 if name=='imp_sloped' else 1
                controls[at]={'junction':{'variable':8,'angle_low':float(lo),'angle_high':float(hi),
                                         'sign':-1 if name=='imp_sloped' else 1}}
            surfaces.append({'name':name,'a':a.tolist(),'b':b.tolist(),'normal':q.n.tolist(),
                             'controls':controls,'degree':3,'knots':q.spline.t.tolist()})
        else:
            v=b-a;n=np.array([-v[1],v[0]])/np.linalg.norm(v)
            surfaces.append({'name':f'fixed_{j}','a':a.tolist(),'b':b.tolist(),'normal':n.tolist(),
                             'controls':[0,0],'degree':1,'knots':[0,0,1,1],'exact_straight_area':True})
    domains=[]
    for j,value in enumerate(candidate.vector):
        h=I(junction_halfwidth if candidate.source=='IMP' and j==8 else halfwidth)
        x=I(value)+I(h.hi.copy_negate(),h.hi)
        lo,hi=(0,1) if candidate.source=='IMP' and j==8 else (-.1,.1)
        domains.append(record(I(max(x.lo,I(lo).lo),min(x.hi,I(hi).hi))))
    field=e.field;tz,tr=field.spline.get_knots();kz,kr=field.spline.degrees
    field_model={'r_knots':tr.tolist(),'z_knots':tz.tolist(),'degree_r':int(kr),'degree_z':int(kz),
        'coefficients':field.spline.get_coeffs().reshape(len(tz)-kz-1,len(tr)-kr-1).tolist(),
        'toroidal_K':float(2e-7*144*field.tf_current_ka*1000)}
    donors=[];thermal=e.thermal_paths
    def finite(values):return [float(x) if np.isfinite(x) else None for x in values]
    for i,path in enumerate(thermal.lengths.paths):
        b,spatial,_=thermal.data[i];metric=thermal.lengths.metric[i]
        donors.append({'path':path.tolist(),'metric':finite(metric),
            'integrand':finite(b/b[0]*metric),'mapping':finite(spatial/metric*b[0]/b),
            'b0':float(b[0]),'qu':float(thermal.qu[i]),
            'weight_MW':float(e.source.weights[i]*e.reference.power_W/1e6)})
    return {'schema':'geometry-transport-input.v1','source':candidate.source,'candidate':candidate.record(),
        'domains':domains,'surfaces':surfaces,'field':field_model,'reference':asdict(e.reference),
        'grad_sep':float(thermal.grad_sep),'donors':donors,
        'definitions':{'path':'binary64 vertices interpreted as exact constants',
            'metric':'exact frozen nodal metric, piecewise-linear between exact path vertices',
            'mapping':'exact frozen positive nodal mapping, piecewise-linear',
            'spreading':'continuous transformed-coordinate integral, enclosed by positive range',
            'field':'exact exported spline coefficients; SciPy construction error not included',
            'source':'finite frozen donor weights/qu; continuous source quadrature not certified',
            'normal':'frozen binary64 normal; geometry uses exact v=b-a and v dot v',
            'wall':'continuous normal-graph spline, not production sampled polyline'}}


def run(output,count=16):
    from d3d_optimization.boundary_evaluator import prepare_boundary
    from d3d_optimization.reference_transport import ReferencePoint
    if not isinstance(count,int) or not 2<=count<=320:raise ValueError('count 2..320 required')
    out=Path(output).resolve()
    if out.drive.lower()!='d:' or out.is_relative_to(CAMPAIGN.resolve()):
        raise ValueError('new D: study outside expired campaign required')
    out.mkdir(parents=True,exist_ok=False);start=perf_counter()
    campaign=json.loads((CAMPAIGN/'campaign.json').read_text());root=Path(campaign['frozen_root'])
    final=json.loads((CAMPAIGN/'final_verification.json').read_text())
    protected={**final['campaign_artifact_sha256'],**final['preserved_prior_manifest_sha256']}
    protected.update({str(root/p):h for p,h in json.loads((CAMPAIGN/'frozen_files.json').read_text()).items()})
    for name in ['campaign.json','frozen_files.json','final_verification.json']:protected[str(CAMPAIGN/name)]=digest(CAMPAIGN/name)
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError(f'protected input mismatch: {p}')
    settings={'created_utc':datetime.now(timezone.utc).isoformat(),'count':count,'protected_sha256':protected,
              'code_sha256':snapshot_code(REPO,out/'source'),'python':sys.version,'executable':sys.executable}
    save(out/'settings.json',settings);rows=[]
    for source,mm in [('IMP',20),('OMP',50)]:
        candidate=Candidate.from_record(json.loads((CAMPAIGN/f'{source}/R{mm:03d}/manifest.json').read_text())['candidate'])
        print(f'Preparing {source}, {count} frozen paths',flush=True);t=perf_counter()
        original=prepare_boundary(root,source,ReferencePoint(**campaign['reference']),count=count,step=.0005,nodes=161,wall_step=.001)
        e=evaluator(original,mm,'cpu');preparation=perf_counter()-t
        for label,width,jwidth in [('local','0.0000005','0.000005'),('wide','0.0005','0.0005')]:
            model=freeze_model(e,candidate,width,jwidth);name=f'{source}_{label}'
            print(f'Enclosing {name}, all {len(model["domains"])} variables',flush=True);t=perf_counter()
            certificate=certify_box(model,model['domains'],root_depth=5,area_panels=8,thermal_panels=64)
            elapsed=perf_counter()-t;save(out/f'{name}.json',certificate)
            row={'case':name,'variables':len(model['domains']),'donors':count,'seconds':elapsed,
                 'preparation_seconds':preparation,'resolved':certificate['transport_parameter_enclosure_proven'],
                 'resolved_donors':sum(d.get('status')=='enclosed' for d in certificate['donors']),
                 'unresolved':certificate['unresolved'],'total_area_m2':certificate.get('total_area_m2')}
            if certificate['kernel'] is not None:
                # Small spatial budget, separate from exact geometry inclusion.
                universal=certificate['peak_bracket_MW_m2']
                row['frozen_family_peak_bracket']=universal
                try:
                    bounds=certify_kernel(certificate['kernel'],tolerance='.01',max_splits=16)
                    save(out/f'{name}_peak.json',bounds)
                    other=bounds['peak_bracket_MW_m2']
                    row['frozen_family_peak_bracket']=[str(max(interval(universal).lo,interval(other).lo)),
                        str(min(interval(universal).hi,interval(other).hi))]
                except ValueError as exc:
                    if 'image-series budget exceeded' not in str(exc):raise
                    row['spatial_refinement_status']=str(exc)
            rows.append(row);save(out/'summary.json',rows)
            print(f'{name}: {row["resolved_donors"]}/{count} enclosed in {elapsed:.1f}s',flush=True)
        # Full conventional variable domain: retain unresolved outcome honestly.
        model=freeze_model(e,candidate,'0.0005','0.0005')
        model['domains']=[['0','1'] if source=='IMP' and j==8 else [str(I(-.1).lo),str(I(.1).hi)] for j in range(len(candidate.vector))]
        whole=cover(model,splits=0,root_depth=2,area_panels=4)
        save(out/f'{source}_full_domain.json',whole)
        rows.append({'case':f'{source}_full_domain','resolved':whole['all_boxes_resolved'],
            'unresolved':whole['leaves'][0]['certificate']['unresolved']});save(out/'summary.json',rows)
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError(f'protected input changed: {p}')
    save(out/'completion.json',{'seconds':perf_counter()-start,'protected_files':len(protected),
        'continuous_fieldline_enclosure_proven':False,'geometry_global_optimum_proven':False})
    save(out/'artifact_sha256.json',{p.relative_to(out).as_posix():digest(p) for p in out.rglob('*') if p.is_file()})
    print(f'Completed in {perf_counter()-start:.1f}s',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output');parser.add_argument('--replay');parser.add_argument('--count',type=int,default=16)
    args=parser.parse_args()
    if args.replay:
        from d3d_optimization.geometry_transport_certificate import replay_cover
        r=json.loads(Path(args.replay).read_text())
        (replay_cover if r['schema']=='geometry-transport-cover.v1' else replay)(r)
        print('Geometry-to-transport certificate replay passed.')
    elif args.output:run(args.output,args.count)
    else:parser.error('--output or --replay required')
