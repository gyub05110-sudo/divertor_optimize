"""Measured finalist-only refinement; all original studies remain read-only."""
import argparse
import copy
import json
import math
from pathlib import Path
from time import perf_counter
from d3d_optimization.enclosure_math import I,interval,record,intersect
from d3d_optimization.continuous_transport_certificate import ValidatedPaths,certify_box,replay
from d3d_optimization.selective_transport import trace_tight,refine_donor,select_cells,next_donor
from d3d_optimization.selective_peak import certify as peak_certify,replay as peak_replay
from verification.run_continuous_path_study import numeric_audit
from verification.run_optimality_bounds import digest,save,snapshot_code,REPO


PREVIOUS=Path('D:/DivertorOptimization/studies/continuous-path-validation-20261102-v2')


def verified_previous(previous):
    previous=Path(previous).resolve();manifest=previous/'artifact_sha256.json'
    recorded=json.loads(manifest.read_text());protected={str(manifest):digest(manifest)}
    for relative,expected in recorded.items():
        path=(previous/relative).resolve()
        if not path.is_relative_to(previous):raise ValueError('recorded path outside previous study')
        if digest(path)!=expected:raise ValueError('previous recorded hash mismatch: '+str(path))
        protected[str(path)]=expected
    settings=json.loads((previous/'settings.json').read_text())
    for path,expected in settings['protected_sha256'].items():
        if digest(path)!=expected:raise ValueError('inherited recorded hash mismatch: '+path)
        protected[path]=expected
    return protected


def diagnostics(certificate):
    return [{k:d[k] for k in ('index','weight_MW','status','length_m','integral_Bratio_m',
             'sigma_flux_m','sigma_physical_m','sigma_area_m2','center_fraction','reason') if k in d}
            for d in certificate['donors']]


def run(output,previous=PREVIOUS,sources=('IMP','OMP'),seconds=180,peak_seconds=5,max_rounds=8,case='point'):
    if not math.isfinite(seconds) or seconds<0 or not math.isfinite(peak_seconds) or peak_seconds<0:
        raise ValueError('finite nonnegative budgets required')
    if not isinstance(max_rounds,int) or not 0<=max_rounds<=64:raise ValueError('round budget 0..64 required')
    if case not in ('point','local','wide') or not sources or any(s not in ('IMP','OMP') for s in sources):
        raise ValueError('supported source and case required')
    out=Path(output).resolve();previous=Path(previous).resolve()
    if out.drive.lower()!='d:' or out.is_relative_to(previous):raise ValueError('new D: directory outside preserved inputs required')
    out.mkdir(parents=True,exist_ok=False);started=perf_counter()
    protected=verified_previous(previous)
    code=snapshot_code(REPO,out/'source')
    save(out/'settings.json',dict(previous=str(previous),sources=list(sources),seconds=seconds,peak_seconds=peak_seconds,
        max_rounds=max_rounds,case=case,protected_sha256=protected,source_sha256=code,
        budget_scope='refinement only; initial cache replay, diagnostics and final independent replay timed separately'))
    summary=[]
    for source in sources:
        folder=out/source;folder.mkdir()
        saved=json.loads((previous/source/('wide.json' if case=='wide' else 'local.json')).read_text())
        model=saved['model'];box=saved['box'] if case!='point' else [record(I(x)) for x in model['candidate']['vector']]
        paths=[json.loads((previous/source/f'path_{j:03d}.json').read_text()) for j in range(len(model['donors']))]
        t=perf_counter();cache=ValidatedPaths(model,paths);cache_seconds=perf_counter()-t
        print(source,'initial cache replay',cache_seconds,flush=True)
        history=[];visits={};peak_interval=None;cert=None;current_peak=None

        def evaluate(label,mode):
            nonlocal peak_interval,cert,current_peak
            t=perf_counter();cert=certify_box(model,box,cache,field_enclosure=mode)
            cost=perf_counter()-t;save(folder/f'{label}_transport.json',cert)
            row=dict(label=label,transport_seconds=cost,donors=diagnostics(cert),unresolved=cert['unresolved'])
            if cert['kernel']:
                current_peak=peak_certify(cert['kernel'],seconds=peak_seconds,tolerance='.01',max_splits=256)
                save(folder/f'{label}_peak.json',current_peak)
                bracket=interval(current_peak['peak_bracket_MW_m2'])
                peak_interval=bracket if peak_interval is None else intersect(peak_interval,bracket)
                if peak_interval is None:raise ArithmeticError('same-geometry independent bounds disagree')
                row.update(analytic_peak=cert['peak_bracket_MW_m2'],spatial_peak=record(bracket),
                           retained_peak=record(peak_interval),peak_seconds=current_peak['elapsed_seconds'])
            history.append(row);save(folder/'history.json',history)
            print(source,label,row.get('retained_peak'),cost,flush=True)

        evaluate('legacy_point' if case=='point' else 'legacy_'+case,'legacy')
        evaluate('centered_arrival','centered-v1')
        refinement_start=perf_counter();attempts=[]
        for round_index in range(max_rounds):
            if perf_counter()-refinement_start>=seconds:break
            donor=next_donor(cert['donors'],visits)
            if donor is None:break
            row=cert['donors'][donor];path=cache._paths[donor];inputs=path['inputs'];round_started=perf_counter();t=round_started
            visits[donor]=visits.get(donor,0)+1
            attempt=dict(round=round_index,donor=donor,visit=visits[donor],accepted=False)
            try:
                if visits[donor]==1:
                    refined=trace_tight(inputs['field'],inputs['guide'],inputs['grad_sep'])
                    attempt['action']='centered-field-and-picard-contraction'
                else:
                    hits=[h['cell'] for h in row.get('first_hit',{}).get('possible_hits',[])]
                    cells=select_cells(path,hits,limit=24)
                    refined=refine_donor(path,cells,factor=4)
                    attempt.update(action='selected-hit-and-integral-cells',cells=sorted(cells),factor=4)
                attempt['trace_seconds']=perf_counter()-t
                save(folder/f'path_{donor:03d}_round_{round_index:02d}.json',refined)
                if not refined['complete']:raise ValueError('refined original horizon incomplete')
                t=perf_counter();new_cache=cache.replacing(model,donor,refined)
                attempt['path_replay_seconds']=perf_counter()-t
                t=perf_counter();attempt['numeric_audit']=numeric_audit(refined);attempt['audit_seconds']=perf_counter()-t
                cache=new_cache;attempt['accepted']=True
                evaluate(f'round_{round_index:02d}','centered-v1')
            except ValueError as exc:attempt['reason']=str(exc)
            attempt['atomic_round_seconds']=perf_counter()-round_started
            attempts.append(attempt);save(folder/'attempts.json',attempts)
        refinement_elapsed=perf_counter()-refinement_start
        for j,path in enumerate(cache._paths):save(folder/f'final_path_{j:03d}.json',path)
        t=perf_counter();replay(cert,cache)
        if current_peak:peak_replay(current_peak)
        final_replay_seconds=perf_counter()-t
        row=dict(source=source,case=case,donors=len(paths),cache_replay_seconds=cache_seconds,
                 initial_diagnostic_seconds=sum(x['transport_seconds']+x.get('peak_seconds',0) for x in history[:2]),
                 refinement_seconds=refinement_elapsed,refinement_budget_seconds=seconds,
                 budget_overrun_seconds=max(0,refinement_elapsed-seconds),final_replay_seconds=final_replay_seconds,
                 accepted_rounds=sum(a['accepted'] for a in attempts),retained_peak=record(peak_interval) if peak_interval else None,
                 initial_peak=history[0].get('retained_peak'),global_optimum_proven=False,
                 continuous_source_enclosure_proven=False,geometry_feasibility_proven=False)
        summary.append(row);save(out/'summary.json',summary)
    for p,h in protected.items():
        if digest(p)!=h:raise ValueError('preserved input changed: '+p)
    save(out/'completion.json',dict(seconds=perf_counter()-started,protected_files_unchanged=len(protected),
        geometry_global_optimum_proven=False))
    save(out/'artifact_sha256.json',{p.relative_to(out).as_posix():digest(p) for p in out.rglob('*') if p.is_file()})
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True)
    p.add_argument('--previous',default=str(PREVIOUS));p.add_argument('--source',choices=['IMP','OMP','both'],default='both')
    p.add_argument('--seconds',type=float,default=180);p.add_argument('--peak-seconds',type=float,default=5)
    p.add_argument('--rounds',type=int,default=8);p.add_argument('--case',choices=['point','local','wide'],default='point')
    a=p.parse_args();run(a.output,a.previous,('IMP','OMP') if a.source=='both' else (a.source,),a.seconds,a.peak_seconds,a.rounds,a.case)
