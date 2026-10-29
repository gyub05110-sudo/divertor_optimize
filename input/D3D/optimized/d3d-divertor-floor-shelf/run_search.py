"""Corrected eight-variable lower divertor search; no shelf-only shortcut."""
import csv
import json
from pathlib import Path
import shutil
import time
import numpy as np
import shapely
from shapely.geometry import Polygon
from scipy.optimize import differential_evolution
from scipy.stats import qmc
from divertor_optimization.d3d_multisurface import MultiSurfaceEvaluator
from divertor_optimization.d3d_frozen_case import load_case, separatrix, pretrace, save_json, select_target, sha256

p=Path(__file__).resolve().parent
workspace=p.parent.parent
base=p.parent/"d3d-frozen-heat-target-tagged"
old=json.loads((base/"manifest.json").read_text(encoding="utf-8"))
config=json.loads((workspace/"configs/d3d_lower_divertor_floor_shelf.json").read_text())
assert sha256(old["shot_path"]) == old["shot_sha256"]
for name in ["DIIID_geom.json","geometry_provenance.json"]: shutil.copyfile(base/name,p/name)
assert sha256(p/"DIIID_geom.json") == config["geometry_sha256"]
geometry=json.loads((p/"DIIID_geom.json").read_text())
wall=np.asarray(geometry["limiter"],dtype=float)
wall=np.vstack([wall,wall[0]])
for s in config["surfaces"]: select_target(wall,config["geometry_sha256"],dict(s,geometry_sha256=config["geometry_sha256"]))
field,lcfs,axis=load_case(old["shot_path"],1400.)
r_sep,xp,_=separatrix(field,lcfs,axis)
edges,paths,inside,_,drift=pretrace(field,r_sep,axis[1],wall,320,.002)
vessel=shapely.union_all([Polygon(v['pts']) for v in geometry['vv']])
e=MultiSurfaceEvaluator(field,lcfs,wall,paths,inside,edges,r_sep,axis[1],config["surfaces"],vessel_solids=vessel)
peak_ids=[int(np.argmax(m.q)) for m in e.baseline["maps"]]
design,angle=e.hit_masks(e.baseline)
assert all(design[i] and not angle[i] for i in peak_ids), "main peak must be on designed floor without shelf floor"
print("Baseline primary peaks",e.baseline_peaks.tolist(),"surface angles",e.baseline["angles"],flush=True)
start=time.monotonic(); stage="baseline"; records=[]; seen={}; best=None; pool=[]
counts={"evaluated":0,"geometry_valid":0,"feasible":0}
log=(p/"evaluations.jsonl").open("w",encoding="utf-8")


def evaluate(coefficients):
    global best
    coefficients=np.asarray(coefficients)
    key=tuple(np.round(coefficients,12))
    if key in seen: return seen[key]
    r=e.evaluate(coefficients)
    counts["evaluated"]+=1
    row={"case":counts["evaluated"],"stage":stage,"coefficients_m":coefficients.tolist(),
         "geometry_valid":r["geometry_valid"],"feasible":r["feasible"]}
    if r["geometry_valid"]:
        counts["geometry_valid"]+=1
        ratios=[float(max(h.q)/b) for h,b in zip(r["maps"],e.baseline_peaks)]
        score=max(ratios)
        row.update(worst_all_wall_peak_ratio=score,peak_ratios=ratios,angles=r["angles"],
                   surface_power=r["surface_power"],unmapped_power=r["missing"],
                   angle_penalty=r["angle_penalty"],
                   max_displacement_m={n:float(np.max(np.linalg.norm(t-e.targets[n],axis=1))) for n,t in r["targets"].items()})
        if r["feasible"]:
            counts["feasible"]+=1; pool.append((score,coefficients.tolist()))
            if best is None or score<best[0]:
                best=(score,coefficients.copy(),r)
                save_json(p/"best_checkpoint.json",row)
        value=score+r["penalty"]
    else: value=1000+r["penalty"]
    row["search_objective"]=float(value); seen[key]=float(value); records.append(row)
    log.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+"\n")
    if counts["evaluated"]%128==0:
        log.flush(); print(stage,counts,"best ratio",best[0] if best else None,"elapsed",round(time.monotonic()-start,1),flush=True)
    return float(value)


evaluate(np.zeros(8))
for i,(limit,power) in enumerate([(.02,10),(.05,9),(.10,9)]):
    stage=f"Sobol_{limit*1000:g}mm"
    for c in limit*(2*qmc.Sobol(8,scramble=True,seed=20261003+i).random_base2(power)-1): evaluate(c)
while counts["geometry_valid"]<1000:
    stage="extra_valid_20mm"
    for c in .02*(2*qmc.Sobol(8,scramble=True,seed=counts["evaluated"]).random_base2(9)-1): evaluate(c)
runs=[]
for seed in [42,93]:
    stage=f"DE_{seed}"
    rng=np.random.default_rng(seed)
    ranked=sorted(pool)
    init=np.vstack([np.array([x[1] for x in ranked[:32]]),rng.uniform(-.1,.1,(32,8))])
    init[0]=best[1]
    generations=[0]
    def callback(x,convergence):
        generations[0]+=1
        if generations[0]%5==0: print(stage,"generation",generations[0],counts,"best ratio",best[0],flush=True)
    opt=differential_evolution(evaluate,[(-.1,.1)]*8,init=init,seed=seed,maxiter=40,tol=0,atol=0,
                               polish=False,callback=callback)
    runs.append({"seed":seed,"generations":int(opt.nit),"nfev":int(opt.nfev),"message":str(opt.message)})
log.close()
assert counts["geometry_valid"]>=1000 and best is not None
manifest={"shot_path":old["shot_path"],"shot_sha256":old["shot_sha256"],"geometry":old["geometry"],
          "config":config,"time_ms":1400.,"tubes":320,"step_m":.002,"max_path_psi_drift":drift,
          "widths_mm":[1,2,3,5,10],"outer_branch_power":1.,"angle_sensitivity":False,
          "objective":"worst all-wall peak ratio across five widths",
          "coefficients_m":best[1].tolist(),"worst_peak_ratio":best[0],"counts":counts,"optimizer":runs,
          "baseline_angles":e.baseline["angles"],"baseline_surface_power":e.baseline["surface_power"],
          "candidate_angles":best[2]["angles"],"numerically_verified_feasible":False,
          "certified_heat_optimum":False,"elapsed_seconds":time.monotonic()-start}
save_json(p/"manifest.json",manifest)
save_json(p/"numerical_checks.json",[])
print("DONE",counts,"best ratio",best[0],"coefficients",best[1].tolist(),flush=True)
