"""Independent actual IMP/OMP searches using the restored Sobol + DE recipe."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.stats import qmc
from scipy.optimize import differential_evolution
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.case import sha256

ROOT=Path(__file__).resolve().parent

def save(path, value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

def summary(e,r):
    row=dict(geometry_valid=r['geometry_valid'],feasible=r['feasible'])
    if r['geometry_valid']:
        peaks=r.get('peak_q',np.array([max(m.q) for m in r['maps']]))
        ratios=[float(q/b) for q,b in zip(peaks,e.baseline_peaks)]
        row.update(peak_ratios=ratios,worst_all_wall_peak_ratio=max(ratios),angles=r['angles'],
                   surface_power=r['surface_power'],unmapped_power=r['missing'],power_balance_error=r['power_balance_error'])
        if 'boundary_witness' in r:row['boundary_witness']=r['boundary_witness']
    return row

def export(e,r,p,manifest=None):
    np.savez_compressed(p/'case_arrays.npz',baseline_wall=e.wall,candidate_wall=r['wall'],lcfs=e.lcfs,
                        **{f'{label}_q_{mm}mm':m.q for label,result in [('baseline',e.baseline),('candidate',r)]
                           for mm,m in zip([1,2,3,5,10],result['maps'])},
                        **{f'{label}_hit_points':result['maps'][0].points for label,result in [('baseline',e.baseline),('candidate',r)]})
    heat={}
    for mm,m in zip([1,2,3,5,10],r['maps']):
        keep=(m.segment>=0)&np.isfinite(m.q)&np.isfinite(m.points).all(axis=1)
        heat[str(mm)]=dict(points=m.points[keep].tolist(),q_over_P=m.q[keep].tolist())
    witness=r.get('boundary_witness')
    if witness and witness['resolved']:
        for mm,q in zip([1,2,3,5,10],witness['q']):
            heat[str(mm)]['points'].append(witness['point'])
            heat[str(mm)]['q_over_P'].append(q)
    save(p/'viewer_result.json',dict(schema='independent-midplane-viewer.v1',wall=r['wall'].tolist(),
         heat=heat,source=manifest['source'] if manifest else None,smoke_test=manifest.get('smoke_test',False) if manifest else False,
         numerically_verified_feasible=manifest.get('numerically_verified_feasible',False) if manifest else False,
         qualified_search=manifest.get('qualified_search',False) if manifest else False,
         boundary_witness=witness,design_surfaces=manifest.get('surfaces') if manifest else None,
         boundary_mode=manifest.get('boundary_mode','legacy') if manifest else 'legacy',
         imp_buffer_m=manifest.get('imp_buffer_m',0.) if manifest else 0.,
         validity_note=manifest.get('validity_note','아직 다중 해상도 검증 전입니다.') if manifest else None,
         validity_audit_summary=manifest.get('validity_audit_summary') if manifest else None))

def verify(p, max_attempts):
    manifest=json.loads((p/'manifest.json').read_text(encoding='utf-8'))
    if manifest.get('boundary_mode')=='continuous':
        from verification.validate_continuous_run import validate_run
        return validate_run(p)['verification_status']=='passed'
    if manifest['mode']!='optimize': raise ValueError('verification requires a saved optimization')
    evaluators=[]
    for method,count,step,nodes in [('cubic',320,.002,161),('cubic',640,.001,321),('quintic',640,.001,321)]:
        e,meta=prepare_evaluator(ROOT,manifest['source'],count,step,nodes,method,manifest['bound_m'])
        if any(meta[k]!=manifest[k] for k in ['shot_sha256','geometry_sha256','policy','surfaces']):
            raise ValueError('saved inputs/design policy differ from current inputs')
        evaluators.append(e)
    rows=[json.loads(x) for x in (p/'evaluations.jsonl').read_text(encoding='utf-8').splitlines()]
    ranked=sorted((r for r in rows if r['feasible']),key=lambda r:r['worst_all_wall_peak_ratio'])
    attempts=[]; chosen=None
    for row in ranked[:max_attempts]:
        checks=[summary(e,e.evaluate(row['coefficients_m'])) for e in evaluators]
        passed=all(c['feasible'] and c.get('worst_all_wall_peak_ratio',2)<=1+1e-10 for c in checks)
        attempts.append(dict(case=row['case'],passed=passed,checks=checks))
        if passed: chosen=row;break
    manifest['numerically_verified_feasible']=chosen is not None
    manifest['verification_status']='passed' if chosen else 'no_candidate_passed_within_attempt_limit'
    manifest['verification_attempt_limit']=max_attempts
    if chosen:
        manifest['verified_coefficients_m']=chosen['coefficients_m']
        manifest['verified_case']=chosen['case']
        export(evaluators[0],evaluators[0].evaluate(chosen['coefficients_m']),p,manifest)
    save(p/'candidate_verification.json',attempts);save(p/'manifest.json',manifest)
    print(manifest['source'],manifest['verification_status'],flush=True)
    return chosen is not None

def run(args,source):
    p=args.output.resolve()/source
    if args.mode=='verify': return verify(p,args.verification_attempts)
    if (p/'manifest.json').exists() or (p/'evaluations.jsonl').exists():
        raise ValueError(f'refusing to overwrite result: {p}')
    p.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    e,meta=prepare_evaluator(ROOT,source,count=48 if args.smoke else 320,
                           nodes=41 if args.smoke else 161,bound=args.bound_mm/1000,
                           boundary_mode=args.boundary_mode,imp_buffer_m=args.imp_buffer_mm/1000)
    manifest=dict(meta,mode=args.mode,bound_m=e.bound,smoke_test=args.smoke,
                  baseline_peak=e.baseline_peaks.tolist(),baseline=summary(e,e.baseline),
                  objective='worst all-wall peak ratio across five widths, independent unit-power source',
                  numerically_verified_feasible=False,certified_heat_optimum=False,
                  historical_recipe='Sobol 1024/512/512 then two DE runs with 64 initial members and 40 generations')
    manifest['implementation_sha256']={name:sha256(ROOT/name) for name in [
        'd3d_optimization/independent_sources.py','d3d_optimization/continuous_boundary.py',
        'optimize_independent_sources.py']}
    if args.mode=='preflight':
        manifest['elapsed_seconds']=time.monotonic()-start;save(p/'manifest.json',manifest)
        export(e,e.baseline,p,manifest);print(source,'preflight',manifest['baseline'],flush=True);return True
    counts=dict(evaluated=0,geometry_valid=0,feasible=0);seen={};pool=[];best=None;stage='baseline'
    with (p/'evaluations.jsonl').open('w',encoding='utf-8') as log:
        def evaluate(c):
            nonlocal best
            c=np.asarray(c);key=tuple(np.round(c,12))
            if key in seen:return seen[key]
            r=e.evaluate(c);counts['evaluated']+=1
            row=dict(case=counts['evaluated'],stage=stage,coefficients_m=c.tolist(),**summary(e,r))
            if r['geometry_valid']:
                counts['geometry_valid']+=1;score=row['worst_all_wall_peak_ratio'];value=score+r['penalty']
                if r['feasible']:
                    counts['feasible']+=1;pool.append((score,c.tolist()))
                    if best is None or score<best[0]:
                        best=(score,c.copy(),r);save(p/'best_checkpoint.json',row)
            else:value=1000+r['penalty']
            row['search_objective']=float(value);seen[key]=float(value)
            log.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n')
            if counts['evaluated']%128==0:
                log.flush();print(source,stage,counts,'best',best[0] if best else None,flush=True)
            return float(value)
        evaluate(np.zeros(e.dimension))
        for i,(limit,power) in enumerate([(min(.02,e.bound),10),(min(.05,e.bound),9),(e.bound,9)]):
            stage=f'Sobol_{limit*1000:g}mm'
            for c in limit*(2*qmc.Sobol(e.dimension,scramble=True,seed=20261003+i).random_base2(3 if args.smoke else power)-1):evaluate(c)
        minimum=16 if args.smoke else args.min_cases
        local_limit=min(.02,e.bound,.15*min(e.lengths.values()))
        manifest['additional_valid_sampling_bound_m']=local_limit
        while counts['geometry_valid']<minimum:
            if counts['evaluated']>=args.max_cases:raise RuntimeError('geometry-valid minimum not reached; evaluation log retained')
            stage=f'extra_valid_{local_limit*1000:g}mm'
            for c in local_limit*(2*qmc.Sobol(e.dimension,scramble=True,seed=counts['evaluated']).random_base2(5 if args.smoke else 9)-1):evaluate(c)
        runs=[]
        for seed in [42,93]:
            stage=f'DE_{seed}';rng=np.random.default_rng(seed);ranked=sorted(pool)
            init=np.vstack([np.asarray([x[1] for x in ranked[:32]]).reshape(-1,e.dimension),rng.uniform(-e.bound,e.bound,(32,e.dimension))])
            if best:init[0]=best[1]
            opt=differential_evolution(evaluate,[(-e.bound,e.bound)]*e.dimension,init=init,seed=seed,
                                       maxiter=1 if args.smoke else args.generations,tol=0,atol=0,polish=False)
            runs.append(dict(seed=seed,generations=int(opt.nit),nfev=int(opt.nfev),message=str(opt.message),population=len(init)))
    manifest.update(counts=counts,optimizer=runs,elapsed_seconds=time.monotonic()-start,
                    qualified_search=not args.smoke and counts['geometry_valid']>=args.min_cases)
    if best:
        manifest.update(coefficients_m=best[1].tolist(),worst_peak_ratio=best[0],candidate=summary(e,best[2]))
        export(e,best[2],p,manifest)
    else:manifest['status']='no_feasible_candidate'
    save(p/'manifest.json',manifest)
    print(source,'DONE',counts,'elapsed',manifest['elapsed_seconds'],'best',best[0] if best else None,flush=True)
    return best is not None

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source',choices=['IMP','OMP','both'],default='both')
    parser.add_argument('--mode',choices=['preflight','optimize','verify'],default='optimize')
    parser.add_argument('--smoke',action='store_true',help='small integration check; not a qualified optimization')
    parser.add_argument('--boundary-mode',choices=['legacy','continuous'],default='continuous')
    parser.add_argument('--imp-buffer-mm',type=float,default=0.,help='explicit IMP design extension toward PFR; default strict scope')
    parser.add_argument('--min-cases',type=int,default=1000)
    parser.add_argument('--generations',type=int,default=40)
    parser.add_argument('--bound-mm',type=float,default=100.)
    parser.add_argument('--max-cases',type=int,default=100000,help='cap on extra geometry-valid sampling')
    parser.add_argument('--verification-attempts',type=int,default=64)
    args=parser.parse_args()
    if args.min_cases<1000 or args.generations<0 or not 0<args.bound_mm<=1000 or args.verification_attempts<1:
        parser.error('invalid minimum, generations, bound or verification attempts')
    passed=[]
    for source in ['IMP','OMP'] if args.source=='both' else [args.source]:passed.append(run(args,source))
    if not all(passed):raise SystemExit(1)

if __name__=='__main__':main()
