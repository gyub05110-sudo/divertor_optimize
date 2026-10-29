"""Fixed-input fast divertor search with a conservative boundary footprint."""
import argparse,json,time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
from scipy.stats import qmc
from d3d_optimization.boundary_evaluator import prepare_boundary,export_boundary
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.broadening_study import save,load_json,digest
from d3d_optimization.run_provenance import REQUIRED
from compare_broadening import HISTORICAL,select_warm_start

ROOT=Path(__file__).resolve().parent
PRESERVATION='docs/superpowers/plans/2026-10-28-boundary-preservation.json'
MODEL_FILES=('d3d_optimization/boundary_transport.py','d3d_optimization/boundary_evaluator.py',
             'optimize_boundary_transport.py','verification/validate_boundary_transport.py',PRESERVATION)


def preserved(root=ROOT):
    original=load_json(root/PRESERVATION)
    changed=[n for n,v in original.items() if not (root/n).is_file() or digest(root/n)!=v]
    if changed:raise ValueError(f'original files changed: {changed}')
    return dict(files=len(original),unchanged=True)


def capture(folder,extra=()):
    preserved();names=sorted(set(REQUIRED)|set(MODEL_FILES)|set(load_json(ROOT/PRESERVATION))|set(extra))
    save(folder/'provenance.json',dict(schema='boundary-study-provenance.v1',files={n:digest(ROOT/n) for n in names}),exclusive=True)


def verify(folder):
    preserved();p=load_json(Path(folder)/'provenance.json')
    expected=set(REQUIRED)|set(MODEL_FILES)|set(load_json(ROOT/PRESERVATION))
    if p.get('schema')!='boundary-study-provenance.v1' or not expected.issubset(p.get('files',{})):
        raise ValueError('incomplete boundary provenance')
    if any(not (ROOT/n).is_file() or digest(ROOT/n)!=v for n,v in p['files'].items()):
        raise ValueError('study code/input changed since capture')


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=['preflight','optimize','verify'],default='optimize')
    p.add_argument('--source',choices=['IMP','OMP','both'],default='OMP')
    for name,default in [('power-mw',1.),('density',1e19),('chi',.3),('lambda-mm',3.),('gamma',7.),('kappa0',2000.),
                         ('step-mm',.5),('wall-step-mm',1.),('bound-mm',100.)]:
        p.add_argument('--'+name,type=float,default=default)
    for name,default in [('tubes',320),('nodes',161),('quadrature',64),('sobol-power',6),('population',24),('generations',6),('seed',20261028)]:
        p.add_argument('--'+name,type=int,default=default)
    p.add_argument('--smoke',action='store_true')
    a=p.parse_args(argv)
    if (not all(np.isfinite(getattr(a,k)) and getattr(a,k)>0 for k in
        ['power_mw','density','chi','lambda_mm','gamma','kappa0','step_mm','wall_step_mm','bound_mm'])
        or a.tubes<8 or a.nodes<9 or a.quadrature<8 or a.population<5 or a.generations<0
        or not 0<=a.sobol_power<=20 or not 0<=a.seed<2**32):p.error('positive finite inputs and supported counts required')
    if a.smoke:a.tubes,a.nodes,a.sobol_power,a.population,a.generations=48,81,3,8,1
    return a


def search(e,folder,args,warm):
    """Exact proposal memoization; baseline retained; fixed-budget exploration."""
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    if not e.baseline['feasible']:return dict(status='baseline_evaluation_failed',best_coefficients_m=None,evaluated=0)
    start=time.perf_counter();seen={};pool=[];best=None;stage='initial'
    with (folder/'evaluations.jsonl').open('x',encoding='utf-8') as log:
        def objective(c):
            nonlocal best
            c=np.asarray(c,float);key=tuple(c)
            if key in seen:return seen[key]
            r=e.evaluate(c);row=dict(case=len(seen)+1,stage=stage,coefficients_m=c.tolist(),**e.describe(r))
            valid=r['feasible'] and np.isscalar(r.get('peak_q')) and np.isfinite(r['peak_q'])
            # Scaling by the fixed baseline preserves the absolute-peak ordering.
            value=float(r['peak_q']/e.baseline_peak) if valid else 1e6+float(r.get('penalty',1e3))
            if valid:
                pool.append((value,c.copy()))
                if best is None or value<best[0]:
                    best=value,c.copy();save(folder/'best_checkpoint.json',row)
            seen[key]=value;row['objective']=value;log.write(json.dumps(row,allow_nan=False)+'\n')
            if len(seen)%16==0:
                log.flush();print(e.metadata['source'],stage,len(seen),'best ratio',best[0] if best else None,flush=True)
            return value
        objective(np.zeros(e.dimension))
        for c in warm:objective(c)
        stage='sobol';samples=2*qmc.Sobol(e.dimension,scramble=True,seed=args.seed).random_base2(args.sobol_power)-1
        local=min(.005,.1*min(e.lengths.values()),e.bound)
        for i,sample in enumerate(samples):
            center=warm[i%len(warm)] if warm and i%2 else np.zeros(e.dimension)
            amplitude=local if i%4 else min(.02,e.bound)
            objective(np.clip(center+amplitude*sample,-e.bound,e.bound))
        init=[c for _,c in sorted(pool,key=lambda v:v[0])[:args.population]]
        rng=np.random.default_rng(args.seed);center=best[1]
        while len(init)<args.population:init.append(np.clip(center+rng.uniform(-local,local,e.dimension),-e.bound,e.bound))
        stage='differential_evolution'
        opt=differential_evolution(objective,[(-e.bound,e.bound)]*e.dimension,init=np.asarray(init),seed=args.seed,
            maxiter=args.generations,tol=0.,atol=0.,polish=False)
    result=dict(status='candidate_found',best_coefficients_m=best[1].tolist(),best_peak_ratio=best[0],
        evaluated=len(seen),feasible_evaluations=len(pool),elapsed_seconds=time.perf_counter()-start,
        optimizer_message=str(opt.message),global_optimum_proven=False)
    save(folder/'search.json',result);return result


def run_source(args,source):
    folder=args.output.resolve()/source
    if folder.exists():raise FileExistsError(f'refusing to overwrite study: {folder}')
    capture(folder,[HISTORICAL[source]])
    ref=ReferencePoint(args.power_mw*1e6,args.density,args.chi,args.lambda_mm/1000,args.gamma,args.kappa0)
    start=time.perf_counter()
    e=prepare_boundary(ROOT,source,ref,count=args.tubes,step=args.step_mm/1000,nodes=args.nodes,
        wall_step=args.wall_step_mm/1000,quadrature=args.quadrature,bound=args.bound_mm/1000)
    historical=load_json(ROOT/HISTORICAL[source])
    if any(historical[k]!=e.metadata[k] for k in ['shot_sha256','geometry_sha256','surfaces']):raise ValueError('warm geometry mismatch')
    warm,status=select_warm_start(historical,e.dimension,e.bound);warm=[] if warm is None else [warm]
    baseline=e.describe(e.baseline)
    manifest=dict(schema='boundary-transport-study.v1',source=source,reference=asdict(ref),metadata=e.metadata,
        settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        setup_seconds=time.perf_counter()-start,baseline=baseline,candidates={'baseline':[0.]*e.dimension},
        historical_warm_start_status=status,validation_status='not_run',measurement_agreement='outside_user_goal',
        model='one-way two-point thermal width; area-measure periodic boundary footprint; no grazing projection',
        global_optimum_proven=False)
    export_boundary(folder/'baseline.npz',e,e.baseline)
    save(folder/'manifest.json',manifest)
    print(source,'baseline',baseline.get('peak_MW_m2'),'MW/m2; feasible',baseline['feasible'],flush=True)
    if args.mode=='optimize':
        result=search(e,folder/'search',args,warm);manifest['search']=result
        if result.get('best_coefficients_m') is not None:
            c=result['best_coefficients_m'];r=e.evaluate(c);manifest['candidates']['selected']=c
            manifest['selected']=e.describe(r);export_boundary(folder/'selected.npz',e,r)
    manifest['elapsed_seconds']=time.perf_counter()-start;save(folder/'manifest.json',manifest)
    print(source,'saved',folder,'seconds',round(manifest['elapsed_seconds'],2),flush=True)
    return manifest


def main(argv=None):
    args=parse_args(argv)
    if args.mode=='verify':
        from verification.validate_boundary_transport import validate_study
        return validate_study(args.output,args.source)
    for source in ['IMP','OMP'] if args.source=='both' else [args.source]:run_source(args,source)


if __name__=='__main__':main()
