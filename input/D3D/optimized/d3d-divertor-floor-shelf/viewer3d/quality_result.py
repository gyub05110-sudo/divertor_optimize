"""Validate and complete stored quality artifacts without re-running traces."""
import hashlib,json
import numpy as np

def complete_report(report,bundle):
    if report['case_id']!=bundle['id'] or report['shot_sha256']!=bundle['shot_sha256']:
        raise ValueError('Audit and registered case differ')
    last=report['levels'][-1];crosschecks={};angles={};spatial={}
    for key in ['baseline','candidate']:
        model=report['models'][key];reference=bundle['models'][key]['wall']
        digest=hashlib.sha256(json.dumps(reference,separators=(',',':')).encode()).hexdigest()
        if digest!=model['shape_sha256']:raise ValueError('Shape changed since audit')
        model['reference_wall']=reference
        actual=[model['summary'][str(mm)]['all_wall_peak_q_over_P'] for mm in [1,2,3,5,10]]
        coords=[[model['summary'][str(mm)]['peak_R_m'],model['summary'][str(mm)]['peak_Z_m']] for mm in [1,2,3,5,10]]
        change=np.asarray(actual)/last['models'][key]['peaks']-1
        location=np.linalg.norm(np.asarray(coords)-last['models'][key]['locations'],axis=1)
        crosschecks[key]={'relative_increase':change.tolist(),'location_change_m':location.tolist(),
                          'passed':bool(np.max(np.abs(change))<=.01 and np.max(location)<=.002)}
        for p in model['patches']:
            values=list(p['q'].values())
            if p['status']=='direct_hit':
                if any(q is None or not np.isfinite(q) or q<0 for q in values):raise ValueError('Invalid direct-hit heat')
            elif any(q is not None for q in values):raise ValueError('Missing connectivity must retain null heat')
        minimum=min(p['angle_deg'] for p in model['patches'] if p['component']=='shelf' and p['status']=='direct_hit')
        minimum=min(minimum,last['models'][key]['shelf_min_angle_deg'])
        angles[key]={'floor_lower_bound_deg':None,'shelf_lower_bound_deg':1.5,'shelf_min_angle_deg':minimum,'passed':bool(minimum>=1.5)}
        levels=report['surfaces'][key]['quadrature_levels'];previous=levels[-2]['integrals'];latest=levels[-1]['integrals']
        changes={name:((np.asarray(latest[name])-previous[name])/np.maximum(last['models'][key]['component_power'][name],1e-5)).tolist() for name in latest}
        spatial[key]={'integral_relative_change':changes,'passed':all(max(abs(x) for x in a)<=.01 for a in changes.values())}
    report['surface_peak_crosscheck']=crosschecks;report['angle_scope_checks']=angles;report['surface_quadrature_convergence']=spatial
    report['source_accounting_passed']=all(max(m['unmapped_power'])<=.001 and max(m['source_accounting_error'])<=1e-4 for m in last['models'].values())
    report['baseline_constraint_feasible']=angles['baseline']['passed']
    report['verified']=bool(report['forward_converged'] and report['surface_integral_passed'] and report['source_accounting_passed'] and
                            report['unresolved_target_length_m']<1e-6 and all(x['passed'] for x in crosschecks.values()) and
                            all(x['passed'] for x in spatial.values()) and angles['candidate']['passed'])
    return report

if __name__=='__main__':
    from pathlib import Path
    p=Path(__file__).resolve().parent
    catalog=json.loads((p/'catalog.json').read_text(encoding='utf-8'))
    for entry in catalog['cases']:
        path=p/'quality'/(entry['id']+'.json')
        if not path.exists():continue
        report=json.loads(path.read_text(encoding='utf-8'));bundle=json.loads((p/entry['bundle_url']).read_text(encoding='utf-8'))
        complete_report(report,bundle)
        path.write_text(json.dumps(report,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')
        print(entry['id'],'verified',report['verified'],'crosscheck',report['surface_peak_crosscheck'],'angle',report['angle_scope_checks'],flush=True)
