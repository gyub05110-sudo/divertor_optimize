"""Re-optimize the locked finite domain with each donor's baseline thermal width.
Candidate hit locations, gradients, area, operating inputs and engineering tests still change.
"""
import argparse
from pathlib import Path
import numpy as np
from d3d_optimization.broadening_study import load_json,save
from d3d_optimization.campaign_geometry import Candidate,evaluator,evaluate_candidate
from d3d_optimization.operating_sweep import reuse_boundary
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_search_benchmark import ROOT,prep,shortlist,measurement


class BaselineWidths:
    def __init__(self, thermal, widths):
        self.thermal=thermal
        self.widths=np.array(widths,dtype=float,copy=True)
        if self.widths.ndim!=1 or not np.isfinite(self.widths).all() or np.any(self.widths<=0):
            raise ValueError('the comparison requires an arrived baseline for every donor')
        self.widths.flags.writeable=False
        self.grad_sep=thermal.grad_sep

    def evaluate(self,hits,**kwargs):
        resolved=self.thermal.evaluate(hits,**kwargs)
        if resolved['sigma_m'].shape!=self.widths.shape:
            raise ValueError('baseline and candidate donor samples differ')
        return {**resolved,'sigma_m':self.widths}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    grid=load_json(a.grid/'design.json')
    for src in ('OMP','IMP'):
        initial=None
        for spec in grid['references']:
            ref=ReferencePoint(**spec);label=f'{src}_chi{ref.chi_perp_m2_s:g}'
            results=load_json(a.grid/label/'result.json')
            if initial is None:
                initial,_=prep(ROOT,src,ref,0);original=initial
            else:original=reuse_boundary(initial,ref)
            e=evaluator(original,grid['radius_mm'],grid['backend'])
            e.thermal_paths=BaselineWidths(original.thermal_paths,original.baseline['thermal']['sigma_m'])
            rows=[]
            for i,v in enumerate(grid['families'][src]):
                c=Candidate(src,grid['radius_mm'],tuple(v))
                r=evaluate_candidate(e,c)
                rows.append(dict(i=i,id=c.identity,**measurement(r,ref)))
                assert results['coarse'][i]['id']==c.identity
            ids=shortlist(rows,3)
            refined=shortlist([results['fine'][i] for i in ids],1)
            overlapping=[(x,y) for x,y in zip(rows,results['coarse']) if x['eligible'] and y['eligible']]
            max_delta=max((abs(x['peak']/y['peak']-1) for x,y in overlapping),default=None)
            save(a.output/(label+'.json'),dict(frozen_coarse=rows,frozen_shortlist=ids,
                 dynamic_shortlist=results['shortlisted'],frozen_top1=ids[:1],dynamic_coarse_top1=results['shortlisted'][:1],
                 refined_choice_from_frozen=refined,refined_choice_from_dynamic=results['selected'],
                 jointly_eligible=len(overlapping),max_coarse_relative_peak_change=max_delta),exclusive=True)
            print(label,ids,refined,flush=True)

if __name__=='__main__':main()

