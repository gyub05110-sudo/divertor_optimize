"""Evaluate the CPD reference geometry; no optimization is performed."""
import argparse
import json
from pathlib import Path
from cpd_optimization.evaluator import CPDEvaluator


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes',default='input/CPD_5/CPD_5_nodes.yaml')
    parser.add_argument('--equilibrium',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--count',type=int,default=160)
    parser.add_argument('--path-step',type=float,default=.002)
    parser.add_argument('--wall-step',type=float,default=.003)
    parser.add_argument('--arc-step',type=float,default=.003)
    parser.add_argument('--max-path',type=float,default=30.)
    parser.add_argument('--quadrature',type=int,default=64)
    parser.add_argument('--allow-diagnostic',action='store_true')
    args=parser.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    summary={}
    for source in ('IMP','OMP'):
        evaluator=CPDEvaluator(args.nodes,args.equilibrium,source,count=args.count,
            path_step=args.path_step,wall_step=args.wall_step,arc_step=args.arc_step,
            max_path=args.max_path,quadrature=args.quadrature,
            cache_dir=out.parent/'path_cache',allow_diagnostic=args.allow_diagnostic)
        result=evaluator.evaluate()
        summary[source]=dict(metadata=evaluator.metadata,result=evaluator.describe(result))
        evaluator.export(out/f'{source}-baseline.npz',result)
        (out/'baseline.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
        print(json.dumps({source:summary[source]}),flush=True)
    return summary


if __name__=='__main__':main()
