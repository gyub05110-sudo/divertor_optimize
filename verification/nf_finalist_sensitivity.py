"""Same-donor interpolation stress test of the finite benchmark winners and their baselines."""
import argparse
import gc
import sys
from pathlib import Path
from time import perf_counter
from d3d_optimization.broadening_study import load_json,save
from d3d_optimization.boundary_evaluator import export_boundary
from d3d_optimization.campaign_geometry import Candidate,evaluator,evaluate_candidate
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.operating_sweep import reuse_boundary
from d3d_optimization.physical_profile import read_physical_profile
from verification.nf_search_benchmark import ROOT,prep,measurement
sys.path.insert(0,str(ROOT/'viewers/d3d-progression-20261104'))
from distribution_policy_v3_final import from_npz,distribution_compare


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--grid',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    design=load_json(a.grid/'design.json');report={}
    for src in ('OMP','IMP'):
        profiles={};stats={}
        for method,level in [('cubic',1),('quintic',3)]:
            initial=None
            for spec in design['references']:
                ref=ReferencePoint(**spec);label=f'{src}_chi{ref.chi_perp_m2_s:g}'
                grid=load_json(a.grid/label/'result.json')
                if not grid['selected']:
                    report[label]={'status':'no_eligible_finalist'};continue
                i=grid['selected'][0]
                if initial is None:
                    initial,_=prep(ROOT,src,ref,level);current=initial
                else:current=reuse_boundary(initial,ref)
                e=evaluator(current,design['radius_mm'],design['backend'])
                c=Candidate(src,design['radius_mm'],tuple(design['families'][src][i]))
                for name,result in [('baseline',current.baseline),('candidate',evaluate_candidate(e,c))]:
                    stats[label,method,name]=measurement(result,ref)
                    if 'deposition' not in result:continue
                    path=a.output/f'{label}_{method}_{name}.npz'
                    export_boundary(path,current,result)
                    physical=read_physical_profile(path,ref.power_W)
                    prof=from_npz(path,power_MW=physical['assigned_power_W']/1e6)
                    profiles[label,method,name]=prof
                del current,e
            del initial
            gc.collect()
        for spec in design['references']:
            label=f"{src}_chi{spec['chi_perp_m2_s']:g}"
            if label in report:continue
            row={}
            for name in ('baseline','candidate'):
                saved=[stats[label,m,name] for m in ('cubic','quintic')]
                row[name]={'cubic':saved[0],'quintic':saved[1]}
                if all((label,m,name) in profiles for m in ('cubic','quintic')):
                    diff=distribution_compare(profiles[label,'cubic',name],profiles[label,'quintic',name])
                    row[name]['comparison']=diff
            report[label]=row
            save(a.output/(label+'.json'),row,exclusive=True)
        print(src,'interpolation completed',flush=True)
    save(a.output/'summary.json',report,exclusive=True)

if __name__=='__main__':main()

