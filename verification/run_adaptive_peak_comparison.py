"""Shared, bounded finalist ranking of previously qualified frozen kernels."""
import argparse
import json
from pathlib import Path
from time import perf_counter
from verification.run_optimality_bounds import digest,save,snapshot_code,REPO
from verification.run_selective_peak_comparison import (
    problem_identity,require_accepted_report,verify_code_scope)


def load_inputs(study,level,code_scope_audit=None):
    study=Path(study).resolve();manifest=study/'artifact_sha256.json'
    recorded=json.loads(manifest.read_text(encoding='utf-8'))
    protected={str(manifest):digest(manifest)}
    def read(relative):
        path=(study/relative).resolve()
        if (not path.is_relative_to(study) or relative not in recorded
            or digest(path)!=recorded[relative]):
            raise ValueError('recorded input hash mismatch: '+relative)
        protected[str(path)]=recorded[relative]
        return json.loads(path.read_text(encoding='utf-8'))
    report=read('summary.json')
    if code_scope_audit:protected.update(verify_code_scope(study,report,code_scope_audit))
    else:require_accepted_report(report)
    models={};ids={};pairs=[]
    for source in ('IMP','OMP'):
        baseline='baseline_'+source
        labels=[name for name,c in report['candidates'].items()
                if c['candidate']['source']==source and c['verdict']['resolution_screen_passed']]
        if not labels:continue
        pid=problem_identity(report,source,level)
        for label in [baseline]+labels:
            models[label]=read(f'{label}/{level}/kernel.json');ids[label]=pid
        pairs.extend((label,baseline) for label in labels)
    if all(k in models for k in ('historical_OMP50','campaign_OMP50')):
        pairs.append(('historical_OMP50','campaign_OMP50'))
    if not pairs:raise ValueError('qualified candidates required')
    return dict(models=models,problem_ids=ids,pairs=pairs),protected


def run(study,output,level='cubic_320',seconds=136,max_operations=4096,code_scope_audit=None):
    from d3d_optimization.adaptive_peak import rank
    start=perf_counter();study=Path(study).resolve();out=Path(output).resolve()
    if out.drive.lower()!='d:' or out.is_relative_to(study):
        raise ValueError('new D: output outside input study required')
    inputs,protected=load_inputs(study,level,code_scope_audit)
    out.mkdir(parents=True,exist_ok=False)
    code=snapshot_code(REPO,out/'source')
    save(out/'settings.json',dict(study=str(study),level=level,seconds=seconds,max_operations=max_operations,
        pairs=inputs['pairs'],problem_ids=inputs['problem_ids'],protected_sha256=protected,source_sha256=code,
        code_scope_audit=str(code_scope_audit) if code_scope_audit else None,
        scope='frozen numerically qualified donor kernels; no continuous-source or geometry proof'))
    preparation=perf_counter()-start
    print('START',level,'kernels',len(inputs['models']),'pairs',len(inputs['pairs']),'preparation',preparation,flush=True)
    result=rank(**inputs,seconds=seconds,max_operations=max_operations)
    save(out/'ranking.json',result)
    for label,certificate in result['certificates'].items():save(out/f'{label}.json',certificate)
    for path,h in protected.items():
        if digest(path)!=h:raise ValueError('source input changed: '+path)
    for relative,h in code.items():
        if digest(REPO/relative)!=h or digest(out/'source'/relative)!=h:
            raise ValueError('running source changed: '+relative)
    summary={k:v for k,v in result.items() if k not in ('certificates','schedule')}
    summary.update(level=level,preparation_seconds=preparation,elapsed_seconds=perf_counter()-start,
        rows=[dict(label=name,peak=cert['peak_bracket_MW_m2'],operations=len(cert['operations']),
                   splits=cert['completed_splits']) for name,cert in result['certificates'].items()],
        source_and_inputs_unchanged=True,continuous_transport_enclosure=False,
        timing_note='elapsed includes input verification, source snapshot, generation, mandatory final replay and result serialization; final manifest write excluded')
    save(out/'summary.json',summary)
    save(out/'artifact_sha256.json',{p.relative_to(out).as_posix():digest(p) for p in out.rglob('*') if p.is_file()})
    print('COMPLETE',level,summary['elapsed_seconds'],[(c['a'],c['b'],c['status']) for c in result['comparisons']],flush=True)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',required=True);p.add_argument('--output',required=True)
    p.add_argument('--level',choices=['cubic_320','cubic_640','cubic_1280','quintic_640'],default='cubic_320')
    p.add_argument('--seconds',type=float,default=136);p.add_argument('--max-operations',type=int,default=4096)
    p.add_argument('--code-scope-audit')
    a=p.parse_args();run(a.study,a.output,a.level,a.seconds,a.max_operations,a.code_scope_audit)
