"""Read-only historical candidate requalification under the current wall policy.

Run with python -m verification.run_candidate_requalification --output D:/.../new.
All thermal originals and baseline audits are shared by source and resolution.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
from time import perf_counter
import traceback
import numpy as np
from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.boundary_evaluator import prepare_boundary, BoundaryEvaluator, export_boundary
from d3d_optimization.campaign_geometry import Candidate, evaluator, evaluate_candidate, GEOMETRY_POLICY
from d3d_optimization.campaign_incumbents import QUALIFICATION_POLICY
from d3d_optimization.campaign_validation import Qualifier, geometric_audit
from d3d_optimization.reference_transport import ReferencePoint
from verification.validate_boundary_transport import LEVELS, audit, verdict
from verification.run_transport_peak_pilot import freeze_kernel

REPO=Path(__file__).resolve().parents[1]
CAMPAIGN=Path('D:/DivertorOptimization/studies/imp-omp-global-20261029')
HISTORICAL=Path('D:/DivertorOptimization/studies/omp-tangent-curvature-20261029-v2')


def load_candidate(path, source, radius_mm, reference):
    path=Path(path).resolve();m=load_json(path)
    if m.get('reference')!=reference:raise ValueError('candidate operating reference mismatch')
    if 'candidate' in m:c=Candidate.from_record(m['candidate']);conversion=None
    else:
        coefficients=m.get('coefficients_m',[])
        if source!='OMP' or len(coefficients)!=24 or any(v!=0 for v in coefficients[16:]):
            raise ValueError('historical candidate must preserve the opposite geometry')
        c=Candidate('OMP',m['radius_mm'],tuple(coefficients[:16]))
        conversion='exact first sixteen OMP coefficients; eight fixed IMP coefficients verified zero'
    if c.source!=source or c.radius_mm!=radius_mm:raise ValueError('candidate source/radius mismatch')
    return dict(candidate=c.record(),identity=c.identity,origin=str(path),origin_sha256=digest(path),
                historical_peak_MW_m2=m.get('selected',{}).get('peak_MW_m2'),
                geometry_conversion=conversion,requires_current_validation=True)


def new_output(path, protected_roots):
    path=Path(path).resolve()
    if any(path.is_relative_to(Path(p).resolve()) for p in protected_roots):
        raise ValueError('output overlaps a protected historical study')
    path.mkdir(parents=True,exist_ok=False)
    return path


def check_preservation(hashes):
    after={p:digest(p) if Path(p).is_file() else None for p in hashes}
    changed=[p for p in hashes if hashes[p]!=after[p]]
    return dict(unchanged=not changed,checked_files=len(hashes),changed=changed,after_sha256=after)


def consumed_code_paths(root, modules=None):
    """Snapshot imported local code, including indirect script dependencies."""
    root=Path(root).resolve();modules=sys.modules if modules is None else modules
    paths=set()
    for module in list(modules.values()):
        filename=getattr(module,'__file__',None)
        if filename:
            path=Path(filename).resolve()
            if path.suffix=='.py' and path.is_relative_to(root):paths.add(path)
    return sorted(paths)


def export_evidence(folder,e,r):
    folder.mkdir(parents=True,exist_ok=True)
    if 'deposition' not in r:return
    export_boundary(folder/'selected.npz',e,r)
    save(folder/'kernel.json',freeze_kernel(e,r))


def run_source(output,source,rows,root,config,report):
    timings=report['timings'];start=perf_counter();settings=LEVELS[0]
    original=prepare_boundary(root,source,ReferencePoint(**config['reference']),count=settings[1],
        step=settings[2],nodes=settings[3],wall_step=settings[4],method=settings[0])
    timings.append(dict(source=source,level='cubic_320',phase='preparation',seconds=perf_counter()-start,
                        baseline_evaluation_seconds=original.baseline['elapsed_seconds']))
    qualifier=Qualifier(root,source,config,original)
    records={}
    for label,row in rows.items():
        records[label]=dict(**row,reference=config['reference'],qualification_policy=QUALIFICATION_POLICY,
            geometry_policy=GEOMETRY_POLICY,levels=[],complete=False,eligible=None)
    report['candidates'].update(records)
    for index,settings in enumerate(LEVELS):
        method,count,step,nodes,wall_step=settings;level=f'{method}_{count}'
        if index:
            start=perf_counter();original=qualifier.original(index)
            timings.append(dict(source=source,level=level,phase='preparation',seconds=perf_counter()-start,
                                baseline_evaluation_seconds=original.baseline['elapsed_seconds']))
        start=perf_counter();base=audit(original,original.baseline)
        timings.append(dict(source=source,level=level,phase='baseline_audit',seconds=perf_counter()-start))
        report['baselines'].setdefault(source,[]).append(dict(level=level,description=base,metadata=original.metadata))
        export_evidence(output/f'baseline_{source}'/level,original,original.baseline)
        refined=None;refined_base=None
        if index==0:
            start=perf_counter();refined=BoundaryEvaluator(original.direct,original.metadata,qualifier.reference,
                original.wall_step/2,thermal_paths=original.thermal_paths)
            timings.append(dict(source=source,level='wall_refinement',phase='preparation',seconds=perf_counter()-start,
                                baseline_evaluation_seconds=refined.baseline['elapsed_seconds']))
            start=perf_counter();refined_base=audit(refined,refined.baseline)
            timings.append(dict(source=source,level='wall_refinement',phase='baseline_audit',seconds=perf_counter()-start))
            report['baseline_wall_refinement'][source]=refined_base
            export_evidence(output/f'baseline_{source}'/'wall_refinement',refined,refined.baseline)
        for label,record in records.items():
            if record['eligible'] is False:continue
            candidate=Candidate.from_record(record['candidate'])
            start=perf_counter();e=evaluator(original,candidate.radius_mm,config['backend'])
            timings.append(dict(candidate=label,level=level,phase='evaluator_setup',seconds=perf_counter()-start))
            start=perf_counter();r=evaluate_candidate(e,candidate)
            timings.append(dict(candidate=label,level=level,phase='candidate_evaluation',seconds=perf_counter()-start))
            if index==0:record['eligible']=bool(r['geometry_valid'] and r.get('geometry_screen',{}).get('passed',False))
            start=perf_counter();selected=audit(e,r)
            geom=geometric_audit(e,candidate,r) if r['geometry_valid'] else dict(passed_independent=False)
            selected['audit_passed'] &= geom['passed_independent']
            timings.append(dict(candidate=label,level=level,phase='candidate_audit',seconds=perf_counter()-start))
            record['levels'].append(dict(method=method,tubes=count,trace_step_m=step,nodes=nodes,wall_step_m=wall_step,
                rows=dict(baseline=base,selected=selected),geometry_audit=geom,geometry_screen=r.get('geometry_screen')))
            export_evidence(output/label/level,e,r)
            if index==0 and record['eligible']:
                start=perf_counter();ef=evaluator(refined,candidate.radius_mm,config['backend'])
                timings.append(dict(candidate=label,level='wall_refinement',phase='evaluator_setup',seconds=perf_counter()-start))
                start=perf_counter();rf=evaluate_candidate(ef,candidate)
                timings.append(dict(candidate=label,level='wall_refinement',phase='candidate_evaluation',seconds=perf_counter()-start))
                start=perf_counter();s=audit(ef,rf)
                gf=geometric_audit(ef,candidate,rf) if rf['geometry_valid'] else dict(passed_independent=False)
                s['audit_passed'] &= gf['passed_independent']
                timings.append(dict(candidate=label,level='wall_refinement',phase='candidate_audit',seconds=perf_counter()-start))
                record['wall_only_refinement']=dict(rows=dict(baseline=refined_base,selected=s),geometry_audit=gf,
                    max_peak_relative_change=max(abs(refined_base['peak_q']/base['peak_q']-1),
                        abs(s.get('peak_q',1e9)/selected.get('peak_q',1e-30)-1)))
                export_evidence(output/label/'wall_refinement',ef,rf)
            save(output/label/'validation.json',record);save(output/'summary.json',report)
            print('CANDIDATE',label,level,'eligible',record['eligible'],'audit',selected['audit_passed'],
                  'peak',selected.get('peak_MW_m2'),'geometry',r.get('geometry_screen'),flush=True)
        save(output/'summary.json',report)
    for label,record in records.items():
        record['verdict']=verdict(record['levels'],record.get('wall_only_refinement'))
        record['complete']=True
        save(output/label/'validation.json',record)
        print('VERDICT',label,record['verdict'],flush=True)
    save(output/'summary.json',report)


def run(output,backend='cpu'):
    output=Path(output).resolve()
    if output.drive.lower()!='d:':raise ValueError('save new studies on D:')
    output=new_output(output,[CAMPAIGN,HISTORICAL]);started=perf_counter();protected={}
    report=dict(schema='saved-candidate-requalification.v1',complete=False,candidates={},baselines={},
        baseline_wall_refinement={},timings=[],promoted_to_qualified_store=False,global_optimum_proven=False)
    try:
        campaign=load_json(CAMPAIGN/'campaign.json');root=Path(campaign['frozen_root'])
        rows={}
        for source,mm in [('IMP',20),('IMP',50),('OMP',50)]:
            label=f'campaign_{source}{mm}'
            rows[label]=load_candidate(CAMPAIGN/source/f'R{mm:03d}'/'manifest.json',source,mm,campaign['reference'])
        rows['historical_OMP50']=load_candidate(HISTORICAL/'R050/manifest.json','OMP',50,campaign['reference'])
        verification=load_json(CAMPAIGN/'final_verification.json')
        expected={**verification['campaign_artifact_sha256'],**verification['preserved_prior_manifest_sha256']}
        expected.update({str(root/p):h for p,h in load_json(CAMPAIGN/'frozen_files.json').items()})
        for p,h in expected.items():
            if digest(p)!=h:raise ValueError(f'historical protected hash mismatch: {p}')
        protected.update(expected)
        # Include every file consumed from historical selected cases and all manifests.
        for folder in [CAMPAIGN/'IMP/R020',CAMPAIGN/'IMP/R050',CAMPAIGN/'OMP/R050',HISTORICAL/'R050']:
            for p in folder.rglob('*'):
                if p.is_file() and (p.suffix in ['.json','.npz']):protected[str(p)]=digest(p)
        for p in [CAMPAIGN/'campaign.json',CAMPAIGN/'frozen_files.json',CAMPAIGN/'final_verification.json']:
            protected[str(p)]=digest(p)
        files=consumed_code_paths(REPO)
        code_hashes={}
        for p in files:
            relative=p.relative_to(REPO);target=output/'source'/relative;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(p,target);code_hashes[str(p)]=digest(p)
            if digest(target)!=code_hashes[str(p)]:raise ValueError('code changed while copying')
        config=dict(reference=campaign['reference'],backend=backend,deadline_epoch=float('inf'))
        report.update(reference=config['reference'],geometry_policy=GEOMETRY_POLICY,qualification_policy=QUALIFICATION_POLICY,
            backend=backend,created_utc=datetime.now(timezone.utc).isoformat(),python=sys.version,
            reference_origin=str(CAMPAIGN/'campaign.json'),input_root=str(root),candidate_inputs=rows,
            code_sha256=code_hashes,protected_sha256=protected,
            preparation_timing_note='Preparation includes original baseline evaluation, separately reported as a subcomponent.')
        save(output/'summary.json',report)
        for source in ['IMP','OMP']:
            run_source(output,source,{k:v for k,v in rows.items() if v['candidate']['source']==source},root,config,report)
        report['complete']=True
    except BaseException as exc:
        report['failure']=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
        raise
    finally:
        report['elapsed_seconds']=perf_counter()-started
        preservation=check_preservation(protected);report['preservation']=preservation
        if 'code_sha256' in report:report['code_preservation']=check_preservation(report['code_sha256'])
        report['timing_totals_seconds']={phase:sum(t['seconds'] for t in report['timings'] if t['phase']==phase)
            for phase in sorted({t['phase'] for t in report['timings']})}
        save(output/'summary.json',report)
        save(output/'preservation.json',preservation)
        save(output/'artifact_sha256.json',{p.relative_to(output).as_posix():digest(p) for p in output.rglob('*')
            if p.is_file() and p.name!='artifact_sha256.json'})
    if not preservation['unchanged'] or not report['code_preservation']['unchanged']:
        raise ValueError('protected inputs or running code changed; inspect preservation evidence')
    print('COMPLETE',output,'seconds',report['elapsed_seconds'],flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--backend',choices=['cpu','gpu'],default='cpu')
    args=parser.parse_args();run(args.output,args.backend)
