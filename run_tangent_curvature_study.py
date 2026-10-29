"""Reproducible OMP-only G1/radius/opening constrained heat-flux search."""
import argparse,json,time,shutil
from pathlib import Path
from dataclasses import asdict
import numpy as np
from scipy.optimize import minimize,differential_evolution
from scipy.interpolate import BSpline
from d3d_optimization.tangent_geometry import KNOTS
from d3d_optimization.boundary_evaluator import prepare_boundary
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.tangent_evaluator import constrained_evaluator,ACTIVE
from d3d_optimization.broadening_study import save,load_json,digest
from optimize_regional_transport import write_arrays
from optimize_boundary_transport import capture,verify

ROOT=Path(__file__).resolve().parent
FILES=['run_tangent_curvature_study.py','d3d_optimization/tangent_geometry.py',
       'd3d_optimization/tangent_direct.py','d3d_optimization/tangent_evaluator.py',
       'd3d_optimization/accelerated_transport.py','d3d_optimization/gpu_footprint.py',
       'optimize_regional_transport.py','results/regional-20261028/qualified/OMP/manifest.json']


def expand(x):
    return np.r_[np.asarray(x,float),np.zeros(8)]


def project_floor(e,target):
    def curve(x): return e.direct.base.curve('floor',x)
    # Smooth signed sampled curvature constraints make seed generation fast.
    # Sampling NEVER decides acceptance: exact span extrema remain mandatory.
    base=curve(np.zeros(4));u=np.linspace(0,1,1001); basis=BSpline(KNOTS,np.eye(8),3)
    b0=basis(u)[:,2:6];b1=basis.derivative()(u)[:,2:6];b2=basis.derivative(2)(u)[:,2:6]
    d1=base.d1(u);d2=base.d2(u);points=base.point(u)
    def margins(y):
        x=y*.03; first=d1+b1@x; second=d2+b2@x
        curvature=base.length*second/(base.length**2+first**2)**1.5
        xy=points+(b0@x)[:,None]*base.n
        gap2=np.sum((xy-e.lip)**2,axis=1)
        length=np.trapezoid(np.sqrt(base.length**2+first**2),u)/base.length
        return np.r_[.995-e.rmin*curvature,.995+e.rmin*curvature,
                     (gap2-(e.gap_min+1e-6)**2)/.03**2,1.1999-length]
    target=np.asarray(target,float)
    r=minimize(lambda y:float(np.sum((y-target/.03)**2)),np.zeros(4),method='SLSQP',
        bounds=[(-e.bound/.03,e.bound/.03)]*4,constraints=[dict(type='ineq',fun=margins)],
        options=dict(maxiter=100,ftol=1e-9))
    x=r.x*.03;q=curve(x)
    passed=q.max_curvature()<=1/e.rmin+1e-7 and q.minimum_distance(e.lip)[0]>=e.gap_min-1e-10 and q.length_ratio()<=1.2
    return x if passed else None


def search(e,folder,args,warm,cpu):
    rng=np.random.default_rng(args.seed+round(e.rmin*1000)); folder.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter(); scales=np.concatenate([np.full(4,min(.012,e.lengths[n]**2/(100*e.rmin))) for n in ACTIVE])
    scales[12:16]=.025
    seeds=[]; old=np.array(load_json(ROOT/'results/regional-20261028/qualified/OMP/manifest.json')['candidates']['selected'])
    targets=[np.zeros(4),old[12:16],-old[12:16]]
    targets += [rng.uniform(-.05,.05,4) for _ in range(args.population)]
    for i,target in enumerate(targets):
        c=project_floor(e,target)
        if c is None: continue
        x=np.zeros(16); x[12:16]=c
        if e.constraints(expand(x))['passed']:
            seeds.append(x)
        if i%8==0: print('radius',e.rmin,'geometry seeds',len(seeds),flush=True)
    for x in warm:
        if e.constraints(expand(x))['passed']:seeds.append(np.asarray(x))
    if not seeds:return dict(status='no_analytic_seed_found',global_optimum_proven=False)
    seen={}; pool=[]; best=None; stage='seed'
    with (folder/'evaluations.jsonl').open('x',encoding='utf-8') as log:
        def objective(x):
            nonlocal best
            x=np.asarray(x,float); key=tuple(x)
            if key in seen:return seen[key]
            c=expand(x); check=e.constraints(c)
            if not check['passed']:
                value=10+max(0,e.rmin*max(check['curvature_max_1_m'].values())-1)+max(0,(e.gap_min-check['gap_m'])/.03)
                row=dict(feasible=False,shape_constraints=check)
            else:
                result=e.evaluate(c); row=e.describe(result)
                value=float(result['peak_q']/e.baseline_peak) if result['feasible'] else 20+min(result.get('penalty',1000),1000)/1000
                if result['feasible']:
                    if best is None:
                        cr=cpu.evaluate(c); parity=abs(result['peak_q']/cr['peak_q']-1)
                        if parity>1e-10 or not cr['feasible']:raise ValueError('initial candidate CPU/GPU parity failed')
                        save(folder/'backend_preflight.json',dict(peak_relative_difference=parity,passed=True,coefficients_m=c.tolist()))
                    pool.append((value,x.copy()))
                    if best is None or value<best[0]:
                        best=(value,x.copy());save(folder/'best_checkpoint.json',dict(coefficients_m=c.tolist(),**row))
            row.update(case=len(seen)+1,stage=stage,coefficients_m=c.tolist(),objective=value)
            seen[key]=value;log.write(json.dumps(row,allow_nan=False)+'\n');log.flush()
            if len(seen)%24==0:print('radius',e.rmin,stage,len(seen),'feasible',len(pool),'best',best[0] if best else None,flush=True)
            return value
        for x in seeds: objective(x)
        # Include feasible perturbations on the other OMP components, both signs.
        stage='regional_seeds'
        for i in range(args.population):
            x=seeds[i%len(seeds)].copy(); x[:12]+=rng.uniform(-1,1,12)*scales[:12]
            objective(x)
        if best is None:return dict(status='no_wall_physics_feasible_seed_found',global_optimum_proven=False)
        ranked=sorted(pool,key=lambda a:a[0]); init=[x for _,x in ranked[:args.population]]
        while len(init)<args.population:
            init.append(np.clip(best[1]+rng.uniform(-1,1,16)*scales,-e.bound,e.bound))
        # Retain global bounds: upward/downward are both allowed. The initial
        # population is geometry-informed; no baseline feasibility exemption.
        stage='differential_evolution'
        differential_evolution(objective,[(-e.bound,e.bound)]*16,init=np.asarray(init),
            seed=args.seed,maxiter=args.generations,tol=0,atol=0,polish=False,
            mutation=(.3,.9),recombination=.75)
        # Deterministic multiscale coordinate polishing of the current best.
        stage='coordinate_polish'
        for factor in [.35,.12,.04]:
            for j in range(16):
                for sign in [-1,1]:
                    x=best[1].copy();x[j]=np.clip(x[j]+sign*factor*scales[j],-e.bound,e.bound);objective(x)
    pool=sorted(pool,key=lambda a:a[0])
    result=dict(evaluated=len(seen),feasible_evaluations=len(pool),elapsed_seconds=time.perf_counter()-start,
        best_coefficients_m=expand(best[1]).tolist(),best_peak_ratio=best[0],global_optimum_proven=False,
        feasible_pool=[dict(ratio=v,coefficients_m=expand(x).tolist()) for v,x in pool])
    save(folder/'search.json',result);return result


def run(args):
    if args.output.exists():raise FileExistsError('refusing to overwrite study')
    capture(args.output,FILES)
    p=load_json(args.output/'provenance.json')
    for name in p['files']:
        dest=args.output/'snapshot'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    start=time.perf_counter();reference=ReferencePoint()
    original=prepare_boundary(ROOT,'OMP',reference)
    save(args.output/'study.json',dict(reference=asdict(reference),baseline=original.describe(original.baseline),
        settings=vars(args)|{'output':str(args.output)},metadata=original.metadata,
        comparison_baseline_satisfies_new_tangent=False,setup_seconds=time.perf_counter()-start))
    write_arrays(args.output,original,original.baseline,'baseline')
    warm=[]
    for mm in [100,50,20]:
        e=constrained_evaluator(original,mm/1000,args.backend); folder=args.output/f'R{mm:03d}'
        cpu=constrained_evaluator(original,mm/1000,'cpu')
        outcome=search(e,folder,args,warm,cpu)
        if 'best_coefficients_m' not in outcome:
            save(folder/'manifest.json',dict(radius_mm=mm,search=outcome,selected=None));continue
        c=np.array(outcome['best_coefficients_m']);r=e.evaluate(c);rcpu=cpu.evaluate(c)
        parity=abs(r['peak_q']/rcpu['peak_q']-1)
        if parity>1e-10 or r['feasible']!=rcpu['feasible']:raise ValueError('backend parity failed')
        m=dict(radius_mm=mm,reference=asdict(reference),metadata=e.metadata,selected=e.describe(r),
            baseline=original.describe(original.baseline),coefficients_m=c.tolist(),search=outcome,
            cpu_gpu_peak_relative_difference=parity,global_optimum_proven=False)
        save(folder/'manifest.json',m);write_arrays(folder,e,r,'selected');warm=[c[:16]]
        print('FINISHED',mm,'peak',r['peak_q'],'ratio',r['peak_q']/e.baseline_peak,flush=True)
    verify(args.output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--population',type=int,default=24);p.add_argument('--generations',type=int,default=12)
    p.add_argument('--seed',type=int,default=20261029);p.add_argument('--backend',choices=['cpu','gpu'],default='gpu')
    a=p.parse_args()
    if a.population<5 or a.generations<0:p.error('population >=5 and nonnegative generations required')
    run(a)
