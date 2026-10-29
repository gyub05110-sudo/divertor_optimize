"""Ranked fixed-candidate screen; never changes the acceptance thresholds."""
import argparse,json,time,shutil
from pathlib import Path
import numpy as np
from d3d_optimization.boundary_evaluator import prepare_boundary
from d3d_optimization.regional_transport import RegionalEvaluator
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.broadening_study import load_json,save,digest
from optimize_boundary_transport import ROOT,verify
from optimize_regional_transport import write_arrays


def screen(root,source):
    target=Path(root)/'qualified'/source
    destination=Path(root)/'qualification'
    for existing in [target,destination/(source+'_screen.json')]:
        if existing.exists():raise FileExistsError(existing)
    folder=Path(root)/source;verify(folder);m=load_json(folder/'manifest.json')
    history=[json.loads(l) for l in (folder/'search/evaluations.jsonl').read_text().splitlines()]
    ranked=sorted([r for r in history if r['feasible'] and r.get('peak_ratio',1)<1],key=lambda r:r['peak_ratio'])
    indices=sorted(set(min(i,len(ranked)-1) for i in [0,1,2,4,7,10,15,22,32,45,65,90]))
    selected=[ranked[i] for i in indices];record=dict(source=source,screening_only=True,
        candidates=[dict(rank=i+1,case=r['case'],regional_coefficients_m=r['coefficients_m'],rows={}) for i,r in zip(indices,selected)])
    destination=Path(root)/'qualification';destination.mkdir(exist_ok=True)
    for method in ['cubic','quintic']:
        e=prepare_boundary(ROOT,source,ReferencePoint(**m['reference']),count=640,step=.00025,nodes=321,
            wall_step=.0005,method=method,quadrature=m['settings']['quadrature'],bound=m['settings']['bound_mm']/1000)
        re=RegionalEvaluator(e,source);baseline=e.describe(e.baseline)
        record[method+'_baseline']=baseline
        for item in record['candidates']:
            r=re.evaluate(item['regional_coefficients_m']);item['rows'][method]=re.describe(r)
            print(source,method,'rank',item['rank'],'ratio',item['rows'][method].get('peak_ratio'),flush=True)
            save(destination/(source+'_screen.json'),record)
    passing=[]
    for item in record['candidates']:
        a,b=item['rows']['cubic'],item['rows']['quintic']
        if not (a['feasible'] and b['feasible']):continue
        ratio_error=abs(a['peak_ratio']-b['peak_ratio']);peak_error=abs(a['peak_q']/b['peak_q']-1)
        item['ratio_difference']=ratio_error;item['peak_relative_difference']=peak_error
        item['screen_passed']=bool(ratio_error<=.002 and peak_error<=.005)
        if item['screen_passed']:passing.append(item)
    if not passing:raise ValueError('no screened improved candidate passes; retain original evidence')
    chosen=min(passing,key=lambda item:item['rows']['cubic']['peak_q']);record['chosen_case']=chosen['case']
    save(destination/(source+'_screen.json'),record)
    e=prepare_boundary(ROOT,source,ReferencePoint(**m['reference']),count=m['settings']['tubes'],nodes=m['settings']['nodes'],
        step=m['settings']['step_mm']/1000,wall_step=m['settings']['wall_step_mm']/1000,
        quadrature=m['settings']['quadrature'],bound=m['settings']['bound_mm']/1000)
    re=RegionalEvaluator(e,source);c=re.expand(chosen['regional_coefficients_m']);r=e.evaluate(c)
    target=Path(root)/'qualified'/source
    if target.exists():raise FileExistsError(target)
    target.mkdir(parents=True)
    p=load_json(folder/'provenance.json');p['files'][Path(__file__).name]=digest(Path(__file__));save(target/'provenance.json',p)
    m['original_search_best']=dict(coefficients_m=m['candidates']['selected'],result=m['selected'],validation_path=str(folder/'validation.json'))
    m['qualification']=dict(case=chosen['case'],rank=chosen['rank'],screen_path=str(destination/(source+'_screen.json')),
        reason='original lowest-peak candidate failed interpolation ratio tolerance; thresholds unchanged')
    m['candidates']['selected']=c.tolist();m['selected']=re.describe(r);m['selected_regional_coefficients_m']=chosen['regional_coefficients_m']
    m['validation_status']='not_run';save(target/'manifest.json',m)
    shutil.copytree(folder/'search',target/'search')
    write_arrays(target,e,e.baseline,'baseline');write_arrays(target,e,r,'selected')
    print('QUALIFIED CANDIDATE, full validation still required:',target,chosen['case'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--source',choices=['IMP','OMP'],required=True)
    a=p.parse_args();screen(a.output,a.source)
