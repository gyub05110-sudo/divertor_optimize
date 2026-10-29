"""Fixed-baseline legacy comparison and hypothetical-input robustness screen."""
import argparse,time,sys
from pathlib import Path
from dataclasses import replace,asdict
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from d3d_optimization.boundary_evaluator import prepare_boundary,BoundaryEvaluator,export_boundary
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.transport_evaluator import TransportEvaluator,summarize
from d3d_optimization.broadening_study import load_json,save,digest
from compare_broadening import HISTORICAL,select_warm_start
from optimize_boundary_transport import capture


def run(output):
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    capture(output,[str(Path(__file__).relative_to(ROOT)).replace('\\','/'),*HISTORICAL.values()])
    record=dict(scope='fixed real baseline and hypothetical inputs; no measurement comparison',models={},input_sweep={})
    reference=ReferencePoint()
    for source in ['OMP','IMP']:
        start=time.perf_counter();e=prepare_boundary(ROOT,source,reference)
        record['models'][source]=dict(setup_seconds=time.perf_counter()-start,new=e.describe(e.baseline))
        export_boundary(output/(source+'_baseline.npz'),e,e.baseline)
        start=time.perf_counter();old=TransportEvaluator(e.direct,e.metadata,reference,mode='thermal',wall_step_m=.001,
                                                        thermal_paths=e.thermal_paths,quadrature_order=64)
        record['models'][source]['legacy']=summarize(old,old.baseline)
        record['models'][source]['legacy_baseline_seconds']=time.perf_counter()-start
        if 'deposition' in old.baseline:
            d=old.baseline['deposition'];np.savez_compressed(output/(source+'_legacy.npz'),points=d['points'],q_MW_m2=d['q']*reference.power_W/1e6)
        print(source,'new baseline',record['models'][source]['new']['feasible'],flush=True)
        if source=='OMP':
            warm,_=select_warm_start(load_json(ROOT/HISTORICAL[source]),e.dimension,e.bound)
            variants={}
            for key,attr in [('power','power_W'),('density','density_m3'),('chi','chi_perp_m2_s')]:
                for factor in [.5,2.]:variants[f'{key}_{factor:g}x']=replace(reference,**{attr:getattr(reference,attr)*factor})
            for factor in [.01,10.]:variants[f'chi_{factor:g}x']=replace(reference,chi_perp_m2_s=reference.chi_perp_m2_s*factor)
            variants['lambda_1mm']=replace(reference,lambda_q_m=.001)
            variants['lambda_10mm']=replace(reference,lambda_q_m=.010)
            for label,ref in variants.items():
                start=time.perf_counter()
                other=(BoundaryEvaluator(e.direct,e.metadata,ref,.001,64) if ref.lambda_q_m==reference.lambda_q_m
                       else prepare_boundary(ROOT,source,ref,count=320,nodes=161))
                base=other.describe(other.baseline);selected=other.describe(other.evaluate(warm))
                record['input_sweep'][label]=dict(reference=asdict(ref),baseline=base,historical_shape=selected,
                    elapsed_seconds=time.perf_counter()-start,scope='fixed candidate screen, not a new optimum or convergence claim')
                print(label,'baseline',base['feasible'],'candidate',selected['feasible'],'ratio',selected.get('peak_ratio'),flush=True)
                save(output/'benchmark.json',record)
    record['complete']=True;record['benchmark_script_sha256']=digest(Path(__file__))
    save(output/'benchmark.json',record)
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.output)
