"""Bounded one-variable IMP pilot, separate from the expired global campaign.

Run as python -m verification.run_optimality_bounds --output D:/.../new-folder.
Replay only: python -m verification.run_optimality_bounds --replay bounds.json.
"""
import argparse
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from d3d_optimization import optimality_bounds as ob
from d3d_optimization.campaign_geometry import Candidate,evaluator,evaluate_candidate

REPO=Path(__file__).resolve().parents[1]
CAMPAIGN=Path('D:/DivertorOptimization/studies/imp-omp-global-20261029')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def snapshot_code(root,destination):
    """Preserve exact source bytes, refusing an existing destination."""
    root=Path(root);destination=Path(destination)
    destination.mkdir(parents=True,exist_ok=False)
    hashes={}
    for directory in ['d3d_optimization','verification']:
        (destination/directory).mkdir()
        for source in sorted((root/directory).glob('*.py')):
            relative=source.relative_to(root);expected=digest(source)
            shutil.copyfile(source,destination/relative)
            if digest(destination/relative)!=expected:raise ValueError('source changed during snapshot')
            hashes[relative.as_posix()]=expected
    return hashes


def freeze_family(e,candidate,index,power_MW):
    """Freeze binary64 constants; exact wall endpoints define the continuous model."""
    count=8 if candidate.source=='IMP' else 16
    if not isinstance(index,int) or not 0<=index<count:
        raise ValueError('choose an offset coefficient; junction variation is unsupported')
    coefficients=candidate.bind(e)
    curves=e.curves(coefficients)
    model={'wall':e.wall.tolist(),'curves':[],'power_MW':power_MW}
    active_names=(['imp_sloped','imp_vertical'] if candidate.source=='IMP' else
                  ['omp_entry','shelf','omp_connector','floor'])
    active=active_names[index//4];variable=None
    for surface in e.surfaces:
        name=surface['name'];q=curves[name];segment=int(surface['segment_index'])
        if not np.array_equal(e.wall[segment:segment+2],np.array([q.a,q.b])):
            raise ValueError('spline/wall anchor mismatch')
        model['curves'].append({'segment':segment,'normal':q.n.tolist(),
                                'knots':q.spline.t.tolist(),'controls':q.spline.c.tolist()})
        if name==active:variable=[segment,2+index%4]
    ob.validate_model(model)
    if variable is None:raise ValueError('active surface missing')
    return model,variable


def run(output):
    from d3d_optimization.boundary_evaluator import prepare_boundary
    from d3d_optimization.reference_transport import ReferencePoint
    output=Path(output).resolve()
    if output.drive.lower()!='d:':raise ValueError('pilot data must be saved on D:')
    if output.is_relative_to(CAMPAIGN.resolve()):raise ValueError('expired campaign is read-only')
    output.mkdir(parents=True,exist_ok=False)
    start=perf_counter()
    config_path=CAMPAIGN/'campaign.json';manifest_path=CAMPAIGN/'IMP/R020/manifest.json'
    config=json.loads(config_path.read_text());manifest=json.loads(manifest_path.read_text())
    candidate=Candidate.from_record(manifest['candidate'])
    reference=ReferencePoint(**config['reference'])
    root=Path(config['frozen_root'])
    hashes=json.loads((CAMPAIGN/'frozen_files.json').read_text())
    # Read-only source inputs and old selected/manifests are checked again at completion.
    protected={str(config_path):digest(config_path),str(manifest_path):digest(manifest_path)}
    protected.update({str(root/p):h for p,h in hashes.items()})
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError(f'preserved input hash mismatch: {p}')
    index=2;halfwidth=.0005;center=candidate.vector[index]
    domain=[max(-.1,center-halfwidth),min(.1,center+halfwidth)]
    code_hashes=snapshot_code(REPO,output/'source')
    settings={'created_utc':datetime.now(timezone.utc).isoformat(),'seed_candidate':candidate.record(),
              'coefficient_index':index,'coefficient_domain_m':domain,'reference':asdict(reference),
              'tubes':320,'step_m':.0005,'wall_step_m':.001,'target_nodes':161,
              'bound_splits':16,'bound_subcells':4,'numerical_probes':9,'backend':'cpu',
              'scope':'one-dimensional restricted family; numerical probes are not validated upper bounds',
              'code_sha256':code_hashes,
              'protected_sha256':protected}
    save(output/'settings.json',settings)
    print('Preparing 320-tube IMP evaluator (new pilot).',flush=True)
    original=prepare_boundary(root,'IMP',reference,count=320,step=.0005,nodes=161,wall_step=.001)
    e=evaluator(original,20,'cpu')
    model,variable=freeze_family(e,candidate,index,reference.power_W/1e6)
    print('Computing exact rational area lower bounds.',flush=True)
    bounds=ob.build_bounds(model,variable,domain,splits=16,subcells=4)
    save(output/'bounds.json',bounds)
    observations=[]
    # Include original center exactly, independently of linspace rounding.
    values=np.linspace(*domain,9);values[4]=center
    for value in values:
        vector=list(candidate.vector);vector[index]=float(value)
        probe=Candidate('IMP',20,tuple(vector));result=evaluate_candidate(e,probe)
        row={'candidate':probe.record(),'identity':probe.identity,'coefficient_m':float(value),
             'feasible_numerical':bool(result['feasible']),'geometry_valid':bool(result['geometry_valid']),
             'errors':result.get('errors',[]),'geometry_screen':result.get('geometry_screen'),
             'shape_constraints':result.get('shape_constraints')}
        if 'deposition' in result:
            description=e.describe(result)
            row['peak_MW_m2']=description['peak_MW_m2'];row['ledger']=description['ledger']
        observations.append(row)
        save(output/'observations.json',observations)
        print(f'probe {len(observations)}/9: feasible={row["feasible_numerical"]}, peak={row.get("peak_MW_m2")}',flush=True)
    feasible=[r for r in observations if r['feasible_numerical'] and 'peak_MW_m2' in r]
    selected=min(feasible,key=lambda r:r['peak_MW_m2']) if feasible else None
    summary={'scope':settings['scope'],'lower_bound_scope':bounds['scope'],'observations':len(observations),
             'feasible_numerical_count':len(feasible),'selected_numerical':selected,
             'root_lower_bound_MW_m2':float(ob.F(bounds['history'][0]['lower_bound_MW_m2'])),
             'final_lower_bound_MW_m2':float(ob.F(bounds['history'][-1]['lower_bound_MW_m2'])),
             'global_optimum_proven':False,'elapsed_seconds':perf_counter()-start}
    if selected:summary['comparison']=ob.compare_numerical(bounds,selected['peak_MW_m2'])
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError(f'preserved input changed: {p}')
    for p,h in code_hashes.items():
        if digest(REPO/p)!=h:raise ValueError(f'running source changed: {p}')
    summary['protected_files_unchanged']=len(protected)
    summary['bound_replay_passed']=ob.replay(bounds)
    save(output/'summary.json',summary)
    save(output/'artifact_sha256.json',{p.name:digest(p) for p in output.glob('*.json')})
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output',type=Path);group.add_argument('--replay',type=Path)
    args=parser.parse_args()
    if args.replay:
        ob.replay(json.loads(args.replay.read_text(encoding='utf-8')))
        print('Exact frozen-model lower-bound replay passed; no production optimality certificate.')
    else:run(args.output)
