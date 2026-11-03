"""Budgeted ranking of requalified finite frozen distributions, not ODE proofs."""
import argparse
import json
import math
from pathlib import Path
from time import perf_counter
from d3d_optimization.geometry_transport_certificate import digest as identity
from d3d_optimization.selective_peak import certify,compare
from verification.run_optimality_bounds import digest,save


def problem_identity(report,source,level):
    rows=[r for r in report['baselines'].get(source,[]) if r['level']==level]
    if len(rows)!=1:raise ValueError('one shared baseline per source and level required')
    return identity({'reference':report['reference'],'source':source,'level':level,
                     'metadata':rows[0]['metadata'],'definition':'frozen production finite periodic Gaussian mixture'})


def require_accepted_report(report,*,verified_code_scope=False):
    if (report.get('complete') is not True or report.get('failure')
        or report.get('preservation',{}).get('unchanged') is not True
        or (report.get('code_preservation',{}).get('unchanged') is not True and not verified_code_scope)):
        raise ValueError('completed, preserved requalification required')


def verify_code_scope(study,report,audit_path):
    """Verify an explicit immutable-source exception; never waive changed inputs."""
    study=Path(study).resolve();audit_path=Path(audit_path).resolve()
    audit=json.loads(audit_path.read_text());checked={str(audit_path):digest(audit_path)}
    def check(path,expected):
        path=Path(path)
        if digest(path)!=expected:raise ValueError('code scope hash mismatch: '+str(path))
        checked[str(path)]=expected
    if (audit.get('schema')!='candidate-requalification-code-scope-audit.v1'
        or audit.get('accepted') is not True or not audit.get('historical_preservation_unchanged')
        or audit.get('static_changed_module_references')!=[]
        or audit.get('import_probe',{}).get('exit_code')!=0):
        raise ValueError('accepted unused-code audit required')
    check(study/'summary.json',audit['summary_sha256'])
    check(study/'artifact_sha256.json',audit['artifact_index_sha256'])
    check(study/'supplemental_snapshot_dependencies.json',audit['supplemental_dependencies_sha256'])
    allowed={'d3d_optimization/continuous_path_bounds.py','d3d_optimization/continuous_transport_certificate.py'}
    changed=report['code_preservation']['changed']
    if (audit['changed_code_paths']!=changed or set(audit['allowed_unused_changed_names'])!=allowed
        or len(changed)!=2 or {Path(p).name for p in changed}!={Path(p).name for p in allowed}):
        raise ValueError('unexpected changed code')
    consumed=audit['consumed_snapshot_sha256']
    if not consumed or allowed.intersection(consumed):raise ValueError('changed code was consumed')
    snapshot=(study/'source').resolve();probe=audit['import_probe']
    imported=[Path(p).resolve() for p in probe.get('imported_local_paths',[])]
    if (Path(probe.get('cwd','')).resolve()!=snapshot
        or any(not p.is_relative_to(snapshot) for p in imported)
        or {p.relative_to(snapshot).as_posix() for p in imported}!=set(consumed)):
        raise ValueError('import inventory mismatch')
    for path,h in report['code_sha256'].items():
        matches=[r for r in audit['snapshot_code_sha256'] if Path(path).as_posix().endswith('/'+r)]
        if len(matches)!=1 or audit['snapshot_code_sha256'][matches[0]]!=h:
            raise ValueError('original guarded snapshot omitted or changed')
    if audit['consumed_code_sha256']!=audit['consumed_code_after_sha256']:
        raise ValueError('consumed code changed')
    for relative,h in audit['snapshot_code_sha256'].items():
        path=(study/'source'/relative).resolve()
        if not path.is_relative_to(study/'source'):raise ValueError('snapshot outside source')
        check(path,h)
    evidence=audit['source_evidence']
    if {row['path'] for row in evidence}!=set(consumed):raise ValueError('incomplete source evidence')
    for row in evidence:
        if row['kind']=='original_code_guard':
            before=report['code_sha256'][row['before_path']]
            after=report['code_preservation']['after_sha256'][row['after_path']]
        elif row['kind']=='guarded_campaign_snapshot_supplement':
            before=report['protected_sha256'][row['before_path']]
            after=report['preservation']['after_sha256'][row['after_path']]
            check(row['before_path'],before)
        else:raise ValueError('unknown source evidence')
        if not before==after==consumed[row['path']]==audit['snapshot_code_sha256'][row['path']]:
            raise ValueError('consumed snapshot differs from guarded source')
    for path,h in report['protected_sha256'].items():check(path,h)
    require_accepted_report(report,verified_code_scope=True)
    return checked


def run(study,output,level='cubic_320',seconds=30,baseline_seconds=8,code_scope_audit=None):
    if any(not math.isfinite(x) or x<0 for x in (seconds,baseline_seconds)):raise ValueError('nonnegative finite budgets required')
    study=Path(study).resolve();out=Path(output).resolve()
    if out.drive.lower()!='d:' or out.is_relative_to(study):raise ValueError('new D: output outside source study required')
    report=json.loads((study/'summary.json').read_text())
    recorded=json.loads((study/'artifact_sha256.json').read_text());protected={}
    def read(relative):
        path=(study/relative).resolve()
        if not path.is_relative_to(study) or relative not in recorded or digest(path)!=recorded[relative]:
            raise ValueError('recorded input hash mismatch: '+relative)
        protected[str(path)]=recorded[relative]
        return json.loads(path.read_text())
    report=read('summary.json')
    if code_scope_audit:
        protected.update(verify_code_scope(study,report,code_scope_audit))
    else:require_accepted_report(report)
    out.mkdir(parents=True,exist_ok=False);start=perf_counter();rows=[];certs={};comparisons=[]
    for source in ('IMP','OMP'):
        labels=[f'baseline_{source}']+[name for name,c in report['candidates'].items()
            if c['candidate']['source']==source and c['verdict']['resolution_screen_passed']]
        pid=problem_identity(report,source,level)
        for label in labels:
            kernel=read(f'{label}/{level}/kernel.json')
            result=certify(kernel,seconds=baseline_seconds if label.startswith('baseline_') else seconds,
                           max_splits=1024,tolerance='.01',problem_id=pid)
            save(out/f'{label}.json',result);certs[label]=result
            rows.append(dict(label=label,source=source,level=level,peak=result['peak_bracket_MW_m2'],
                seconds=result['elapsed_seconds'],splits=result['completed_splits'],stop=result['stop_reason']))
            save(out/'progress.json',rows);print(label,rows[-1],flush=True)
        for label in labels[1:]:
            comparison=compare(certs[label],certs[labels[0]],problem_id_a=pid,problem_id_b=pid)
            comparisons.append(dict(a=label,b=labels[0],**comparison))
        if source=='OMP' and all(k in certs for k in ('historical_OMP50','campaign_OMP50')):
            comparisons.append(dict(a='historical_OMP50',b='campaign_OMP50',
                **compare(certs['historical_OMP50'],certs['campaign_OMP50'],problem_id_a=pid,problem_id_b=pid)))
    for path,h in protected.items():
        if digest(path)!=h:raise ValueError('source input changed')
    summary=dict(level=level,rows=rows,comparisons=comparisons,elapsed_seconds=perf_counter()-start,
        protected_sha256=protected,code_scope_audit=str(code_scope_audit) if code_scope_audit else None,
        scope='frozen donor locations, widths and powers from numerically qualified candidates',
        continuous_transport_enclosure=False,geometry_global_optimum_proven=False,
        note='time includes new spatial evidence and independent replay; source numerical preparation/audits are in source study')
    save(out/'summary.json',summary)
    save(out/'artifact_sha256.json',{p.name:digest(p) for p in out.glob('*.json')})
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',required=True);p.add_argument('--output',required=True)
    p.add_argument('--level',choices=['cubic_320','cubic_640','cubic_1280','quintic_640'],default='cubic_320')
    p.add_argument('--seconds',type=float,default=30);p.add_argument('--baseline-seconds',type=float,default=8)
    p.add_argument('--code-scope-audit',help='explicit immutable source audit for an unused-code guard failure')
    a=p.parse_args();run(a.study,a.output,a.level,a.seconds,a.baseline_seconds,a.code_scope_audit)
