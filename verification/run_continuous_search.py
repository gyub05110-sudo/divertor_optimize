"""Capture complete provenance before launching the unchanged search entry point."""
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from d3d_optimization.run_provenance import capture
from optimize_independent_sources import main

if __name__=='__main__':
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source',choices=['IMP','OMP','both'],default='both')
    parser.add_argument('--mode',default='optimize')
    args,_=parser.parse_known_args()
    if args.mode!='optimize':raise ValueError('this launcher is only for new searches')
    for source in ['IMP','OMP'] if args.source=='both' else [args.source]:
        capture(ROOT,args.output.resolve()/source)
    main()
