"""Independent IMP/OMP regional search; all study data saved to chosen output."""
import argparse,os,time,shutil
from pathlib import Path
from dataclasses import asdict
import numpy as np
from d3d_optimization.boundary_evaluator import prepare_boundary,export_boundary
from d3d_optimization.accelerated_transport import AcceleratedEvaluator
from d3d_optimization.regional_transport import RegionalEvaluator,REGIONS
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.broadening_study import load_json,save,digest
from optimize_boundary_transport import ROOT,capture,search,verify
from compare_broadening import HISTORICAL,select_warm_start

FILES=['optimize_regional_transport.py','d3d_optimization/gpu_footprint.py',
       'd3d_optimization/accelerated_transport.py','d3d_optimization/regional_transport.py']


def write_arrays(folder,e,r,name):
    export_boundary(folder/(name+'.npz'),e,r)
    if 'deposition' not in r:return
    np.savetxt(folder/(name+'_wall.csv'),r['wall'],delimiter=',',header='R_m,Z_m',comments='')
    for surface in e.surfaces:
        n=surface['name'];np.savetxt(folder/(name+'_'+n+'.csv'),r['targets'][n],delimiter=',',header='R_m,Z_m',comments='')
    d=r['deposition'];areas=d['areas'];wall=d['wall'];mid=(wall[:-1]+wall[1:])/2
    np.savetxt(folder/(name+'_cells.csv'),np.c_[mid,d['cell_parent'],areas,d['cell_power']*e.reference.power_W,
        d['cell_power']/areas*e.reference.power_W/1e6],delimiter=',',
        header='R_m,Z_m,parent_segment,area_m2,power_W,cell_average_q_MW_m2',comments='')
    np.savetxt(folder/(name+'_heat_samples.csv'),np.c_[d['points'],d['q']*e.reference.power_W/1e6],delimiter=',',
        header='R_m,Z_m,q_MW_m2',comments='')
    np.savetxt(folder/(name+'_donors.csv'),np.c_[e.source.x,e.source.weights*e.reference.power_W,r['hits'].points,
        r['thermal']['length_m'],r['thermal']['Tu_eV'],r['thermal']['Tt_eV'],d['sigma_physical_m']],delimiter=',',
        header='source_offset_m,power_W,hit_R_m,hit_Z_m,L_parallel_m,Tu_eV,Tt_eV,sigma_physical_m',comments='')


def run(args,source):
    folder=args.output/source
    if folder.exists():raise FileExistsError(folder)
    extra=FILES+[HISTORICAL[source]]
    previous='input/D3D/optimized/boundary-optimization-20261028/OMP/manifest.json'
    if source=='OMP':extra.append(previous)
    capture(folder,extra)
    provenance=load_json(folder/'provenance.json')
    for n in provenance['files']:
        dest=folder/'snapshot'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/n,dest)
    start=time.perf_counter();ref=ReferencePoint(args.power_mw*1e6,args.density,args.chi,args.lambda_mm/1000)
    cpu=prepare_boundary(ROOT,source,ref,count=args.tubes,nodes=args.nodes,step=args.step_mm/1000,
        wall_step=args.wall_step_mm/1000,quadrature=args.quadrature,bound=args.bound_mm/1000)
    full=AcceleratedEvaluator(cpu.direct,cpu.metadata,ref,cpu.wall_step,thermal_paths=cpu.thermal_paths) if args.backend=='gpu' else cpu
    parity=abs(full.baseline['peak_q']/cpu.baseline['peak_q']-1)
    if parity>1e-10 or full.baseline['feasible']!=cpu.baseline['feasible']:raise ValueError('baseline backend parity failed')
    e=RegionalEvaluator(full,source)
    historical=load_json(ROOT/HISTORICAL[source])
    if any(historical[k]!=e.metadata[k] for k in ['shot_sha256','geometry_sha256','surfaces']):raise ValueError('warm input mismatch')
    w,_=select_warm_start(historical,full.dimension,full.bound);warm=[] if w is None else [e.project(w)]
    if source=='OMP':warm.append(e.project(load_json(ROOT/previous)['candidates']['selected']))
    m=dict(schema='boundary-transport-study.v1',source=source,reference=asdict(ref),metadata=e.metadata,
        settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        baseline=e.describe(e.baseline),candidates={'baseline':[0.]*full.dimension},
        regional_restriction=list(e.active_names),backend=args.backend,baseline_backend_peak_difference=parity,
        independent_source=True,simultaneous_two_source_optimization=False,global_optimum_proven=False,
        validation_status='not_run',setup_seconds=time.perf_counter()-start)
    write_arrays(folder,full,e.baseline,'baseline');save(folder/'manifest.json',m)
    print(source,'baseline',m['baseline']['peak_MW_m2'],'active',e.active_names,'backend',args.backend,flush=True)
    outcome=search(e,folder/'search',args,warm);m['search']=outcome
    if outcome.get('best_coefficients_m') is not None:
        regional=outcome['best_coefficients_m'];c=e.expand(regional);r=full.evaluate(c)
        m['candidates']['selected']=c.tolist();m['selected']=e.describe(r)
        m['selected_regional_coefficients_m']=regional
        assert m['selected']['opposite_region_displacement_m']==0.
        write_arrays(folder,full,r,'selected')
    m['elapsed_seconds']=time.perf_counter()-start;save(folder/'manifest.json',m)
    print(source,'SEARCH COMPLETE',m.get('selected',{}).get('peak_ratio'),'seconds',m['elapsed_seconds'],flush=True)
    return m


def check_region(folder):
    folder=Path(folder);m=load_json(folder/'manifest.json');v=load_json(folder/'validation.json')
    a=np.load(folder/'baseline_refined.npz');b=np.load(folder/'selected_refined.npz');nodes=641
    movable=np.zeros(len(a['wall']),bool)
    for j,s in enumerate(m['metadata']['surfaces']):
        if s['name'] in m['regional_restriction']:
            start=s['segment_index']+j*(nodes-2);movable[start:start+nodes]=True
    error=float(np.max(abs(a['wall'][~movable]-b['wall'][~movable])))
    frozen=[j for j,s in enumerate(m['metadata']['surfaces']) if s['name'] not in m['regional_restriction']]
    coefficients=np.asarray(m['candidates']['selected']).reshape(-1,4)
    ok=error==0. and np.all(coefficients[frozen]==0.)
    save(folder/'regional_validation.json',dict(opposite_and_fixed_wall_max_error_m=error,
        opposite_coefficients_zero=bool(np.all(coefficients[frozen]==0.)),passed=bool(ok)))
    if not ok:raise ValueError('regional invariant failed')
    for name,data in [('baseline',a),('selected',b)]:
        np.savetxt(folder/(name+'_wall_refined.csv'),data['wall'],delimiter=',',header='R_m,Z_m',comments='')
        for j,s in enumerate(m['metadata']['surfaces']):
            if s['name'] in m['regional_restriction']:
                start=s['segment_index']+j*(nodes-2)
                np.savetxt(folder/(name+'_'+s['name']+'_refined.csv'),data['wall'][start:start+nodes],
                    delimiter=',',header='R_m,Z_m',comments='')
    return v['verdict']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source',choices=['IMP','OMP','both'],default='both');p.add_argument('--mode',choices=['optimize','verify'],default='optimize')
    p.add_argument('--backend',choices=['cpu','gpu'],default='gpu')
    for n,d in [('tubes',320),('nodes',161),('quadrature',64),('sobol-power',6),('population',24),('generations',10),('seed',20261028)]:p.add_argument('--'+n,type=int,default=d)
    for n,d in [('power-mw',1.),('density',1e19),('chi',.3),('lambda-mm',3.),('step-mm',.5),('wall-step-mm',1.),('bound-mm',100.)]:p.add_argument('--'+n,type=float,default=d)
    a=p.parse_args();a.output=a.output.resolve();a.smoke=False
    if (a.tubes<8 or a.nodes<9 or a.quadrature<8 or a.population<5 or a.generations<0 or not 0<=a.sobol_power<=20
        or not 0<=a.seed<2**32 or not all(np.isfinite(getattr(a,k)) and getattr(a,k)>0 for k in
        ['power_mw','density','chi','lambda_mm','step_mm','wall_step_mm','bound_mm'])):p.error('invalid numerical settings')
    for s in ['IMP','OMP'] if a.source=='both' else [a.source]:
        if a.mode=='optimize':run(a,s)
        else:
            from verification.validate_boundary_transport import validate_source
            validate_source(a.output/s);print(s,check_region(a.output/s),flush=True)


if __name__=='__main__':main()
