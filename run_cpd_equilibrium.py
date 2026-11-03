"""Prepare CPD-5 equilibrium; this does not start the optimization clock."""
import argparse
from pathlib import Path
import json
from cpd_optimization.equilibrium import reference_inputs,target_shape,solve_pressure,save_solution,measure,compare_refinement,verify_reference_source
from cpd_optimization.geometry import load_geometry

def main():
    p=argparse.ArgumentParser();p.add_argument('--nodes',default='input/CPD_5/CPD_5_nodes.yaml');p.add_argument('--output',required=True);p.add_argument('--n',type=int,default=65);p.add_argument('--paxis',type=float,default=20000.);p.add_argument('--fit-iterations',type=int,default=5);p.add_argument('--compare-to',help='Coarse output directory; validate and mark final output only after comparison')
    a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    inputs=reference_inputs();target=target_shape();geometry=load_geometry(a.nodes)
    source=Path(__file__).resolve().parent/'input/D3D/d3d/d3d_shot_1a3155370c.parquet'
    provenance=verify_reference_source(source)
    (out/'reference_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    (out/'requested_inputs.json').write_text(json.dumps(dict(inputs=inputs,target=target,n=a.n,initial_pressure_Pa=a.paxis),indent=2))
    history=[];eq=None;paxis=a.paxis
    try:
        for iteration in range(a.fit_iterations):
            eq,history=solve_pressure(a.nodes,n=a.n,paxis=paxis,inputs=inputs,target=target,eq=eq,history=history)
            report=save_solution(eq,a.nodes,out,inputs,target,history,geometry)
            report['pressure_fit_iteration']=iteration
            print(json.dumps(report['metrics']),flush=True)
            if report['metrics']['betaN_relative_error']<.002: break
            paxis*=inputs['betaN']/report['metrics']['betaN']
        report['pressure_fit_converged']=report['metrics']['betaN_relative_error']<.002
        (out/'equilibrium.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        if a.compare_to:
            validation=compare_refinement(a.compare_to,out)
            print(json.dumps(validation),flush=True)
            if validation['failures']: raise RuntimeError('Equilibrium refinement validation failed')
    except Exception as e:
        import traceback
        (out/'failure.json').write_text(json.dumps(dict(error=str(e),traceback=traceback.format_exc(),history=history),indent=2))
        raise
if __name__=='__main__': main()
