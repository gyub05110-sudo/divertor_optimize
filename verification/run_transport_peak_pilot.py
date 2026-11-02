"""Coupled 3-coefficient IMP sampling and frozen-transport peak verification."""
import argparse
import copy
from dataclasses import asdict
from datetime import datetime,timezone
from itertools import product
import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path
from time import perf_counter
import numpy as np
from d3d_optimization import transport_peak_bounds as pb
from d3d_optimization.campaign_geometry import Candidate,evaluator,evaluate_candidate
from verification.run_optimality_bounds import CAMPAIGN,REPO,digest,save,snapshot_code


def candidates(seed):
    if seed.source!='IMP':raise ValueError('this bounded study is IMP only')
    for offsets in product(np.linspace(-.0005,.0005,5),repeat=3):
        vector=list(seed.vector)
        for i,value in enumerate(offsets):vector[i]+=float(value)
        yield Candidate(seed.source,seed.radius_mm,tuple(vector))


def freeze_kernel(e,r):
    dep=r['deposition'];period=float(dep['total_area_m2']);rows=[]
    for center,sigma,weight in zip(dep['donor_area_m2'],dep['sigma_area_m2'],e.source.weights):
        if weight<=0 or not np.isfinite(center):continue
        # These binary64 donor powers become exact constants of this verifier.
        power=float(weight*e.reference.power_W/1e6)
        rows.append({'weight_MW':[power,power],'center_m2':[float(center),float(center)],
                     'sigma_m2':[float(sigma),float(sigma)]})
    return {'period_m2':[period,period],'donors':rows}


def uncertainty_example(model):
    """An explicitly assumed kernel family, NEVER a proved geometry-box image."""
    result=copy.deepcopy(model);period=pb.I(*model['period_m2'])
    for row in result['donors']:
        width=pb.I(*row['sigma_m2']);center=pb.I(*row.pop('center_m2'))/period
        radius=pb.I('.01')*width/period
        center_range=center+pb.I(-1,1)*radius
        row['center_fraction']=[str(max(pb.D(0),center_range.lo)),str(min(pb.D(1),center_range.hi))]
        spread=width*pb.I('.99','1.01');row['sigma_m2']=[str(spread.lo),str(spread.hi)]
    return result


def run(output):
    from d3d_optimization.boundary_evaluator import prepare_boundary
    from d3d_optimization.reference_transport import ReferencePoint
    output=Path(output).resolve()
    if output.drive.lower()!='d:' or output.is_relative_to(CAMPAIGN.resolve()):
        raise ValueError('use a new D: folder outside the expired campaign')
    output.mkdir(parents=True,exist_ok=False);start=perf_counter()
    campaign=json.loads((CAMPAIGN/'campaign.json').read_text())
    manifest_path=CAMPAIGN/'IMP/R020/manifest.json'
    seed=Candidate.from_record(json.loads(manifest_path.read_text())['candidate'])
    root=Path(campaign['frozen_root']);reference=ReferencePoint(**campaign['reference'])
    verification=json.loads((CAMPAIGN/'final_verification.json').read_text())
    protected={**verification['campaign_artifact_sha256'],**verification['preserved_prior_manifest_sha256']}
    protected.update({str(root/p):h for p,h in json.loads((CAMPAIGN/'frozen_files.json').read_text()).items()})
    protected.update({str(CAMPAIGN/name):digest(CAMPAIGN/name)
                      for name in ['campaign.json','frozen_files.json','final_verification.json']})
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError(f'protected input changed: {p}')
    hashes=snapshot_code(REPO,output/'source')
    settings={'created_utc':datetime.now(timezone.utc).isoformat(),'seed':seed.record(),
              'reference':asdict(reference),'indices':[0,1,2],'halfwidth_m':.0005,'levels_per_index':5,
              'count':125,'geometry_global_optimum_proven':False,'code_sha256':hashes,'protected_sha256':protected,
              'runtime':{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
                         'packages':{name:version(name) for name in ['numpy','scipy','shapely']}},
              'backend':'cpu','reason':'small pilot; directed-decimal verification is CPU only'}
    save(output/'settings.json',settings)
    print('Prepare shared 320-tube IMP evaluator',flush=True)
    prepared=perf_counter();original=prepare_boundary(root,'IMP',reference,count=320,step=.0005,nodes=161,wall_step=.001)
    e=evaluator(original,20,'cpu');preparation_seconds=perf_counter()-prepared
    rows=[];best=None;best_peak=np.inf
    for index,candidate in enumerate(candidates(seed),1):
        r=evaluate_candidate(e,candidate)
        row={'candidate':candidate.record(),'identity':candidate.identity,'feasible':bool(r['feasible']),
             'errors':r.get('errors',[]),'shape_constraints':r.get('shape_constraints'),
             'geometry_screen':r.get('geometry_screen')}
        if 'deposition' in r:row['description']=e.describe(r)
        rows.append(row)
        if r['feasible'] and r['peak_q']<best_peak:
            best=candidate;best_peak=r['peak_q']
            save(output/'best_numerical.json',row)
        if index%25==0:
            save(output/'observations.json',rows)
            print(f'{index}/125 sampled, feasible={sum(row["feasible"] for row in rows)}',flush=True)
    save(output/'observations.json',rows)
    if best is None:raise RuntimeError('no numerically feasible candidate; observations preserved')
    resolution=[]
    for count in (320,640,1280):
        if count!=320:
            print(f'Prepare {count}-tube evaluation',flush=True)
            original=prepare_boundary(root,'IMP',reference,count=count,step=.0005,nodes=161,wall_step=.001)
            e=evaluator(original,20,'cpu')
        row={'tubes':count,'baseline_peak_MW_m2':float(e.baseline_peak*reference.power_W/1e6)}
        for name,candidate in [('seed',seed),('selected',best)]:
            r=evaluate_candidate(e,candidate);row[name]=e.describe(r)
            if not r['feasible']:raise RuntimeError(f'{name} failed at {count}; partial outputs preserved')
            if name=='selected':
                model=freeze_kernel(e,r);save(output/f'kernel_{count}.json',model)
                print(f'Certifying frozen {count}-donor peak',flush=True)
                cert_start=perf_counter();certificate=pb.certify(model,tolerance='.001',max_splits=512)
                save(output/f'certificate_{count}.json',certificate)
                row['frozen_kernel_peak_bracket']=certificate['peak_bracket_MW_m2']
                row['peak_certificate_tolerance_met']=certificate['tolerance_met']
                row['certificate_seconds']=perf_counter()-cert_start
                print(f'Peak bracket {certificate["peak_bracket_MW_m2"]}, tolerance_met={certificate["tolerance_met"]}',flush=True)
        resolution.append(row);save(output/'resolution.json',resolution)
    # Conditional transport-aware lower bound for a stated kernel-input family.
    family=uncertainty_example(json.loads((output/'kernel_1280.json').read_text()))
    print('Bound explicitly assumed arrival/width uncertainty family',flush=True)
    robust=pb.certify(family,tolerance='.001',max_splits=64)
    save(output/'conditional_family_certificate.json',robust)
    summary={'scope':'restricted three-coefficient numerical search plus independent frozen-kernel peak enclosure',
             'sample_count':len(rows),'numerically_feasible_count':sum(row['feasible'] for row in rows),
             'selected_candidate':best.record(),'seed_candidate':seed.record(),'resolution':resolution,
             'conditional_family_peak_bracket':robust['peak_bracket_MW_m2'],
             'family_assumptions':'fixed donor powers and period; each sigma +/-1%; center +/-1% of its nominal sigma (clipped to [0,period])',
             'family_covers_geometry_box_proven':False,'geometry_global_optimum_proven':False,
             'certified_geometry_optimality_gap':None,'promoted_to_qualified_store':False,
             'preparation_seconds_320':preparation_seconds,'elapsed_seconds':perf_counter()-start}
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError(f'protected input changed: {p}')
    for p,h in hashes.items():
        if digest(REPO/p)!=h:raise ValueError(f'running code changed: {p}')
    summary['preserved_files_unchanged']=len(protected)
    save(output/'summary.json',summary)
    save(output/'artifact_sha256.json',{p.relative_to(output).as_posix():digest(p) for p in output.rglob('*') if p.is_file()})
    print(f'Complete: {summary["elapsed_seconds"]:.1f}s; no geometry optimum claimed.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output',type=Path);group.add_argument('--replay',type=Path)
    args=parser.parse_args()
    if args.replay:
        pb.replay(json.loads(args.replay.read_text()));print('Frozen/conditional kernel certificate replay passed.')
    else:run(args.output)
