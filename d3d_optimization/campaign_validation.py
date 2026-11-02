"""Fine-grid selection and independent qualification of immutable candidates."""
import time
from pathlib import Path
import numpy as np
from .boundary_evaluator import prepare_boundary,BoundaryEvaluator
from .campaign_geometry import Candidate,evaluator,evaluate_candidate,geometry_screen
from .campaign_state import should_stop
from .reference_transport import ReferencePoint
from .broadening_study import load_json
from .campaign_state import public_save as save
from .campaign_search import scales
from .campaign_incumbents import QualifiedStore,context,protected_candidates,QUALIFICATION_POLICY
from optimize_regional_transport import write_arrays
from verification.validate_boundary_transport import LEVELS,audit,verdict


def diverse_shortlist(rows,limit,scale):
    if not rows:return []
    result=[rows[0]]
    for row in rows[1:]:
        x=np.array(row['candidate']['vector'])
        if all(np.linalg.norm((x-np.array(r['candidate']['vector']))/scale)>.2 for r in result):result.append(row)
        if len(result)>=limit:break
    if len(result)<limit:
        for row in rows:
            if row not in result:result.append(row)
            if len(result)>=limit:break
    return result


def geometric_audit(e,candidate,r):
    c=candidate.bind(e);check=e.constraints(c);dense={};u=np.linspace(0,1,20001)
    for name in check['curvature_max_1_m']:
        q=e.curves(c)[name];points=q.point(u)
        first=np.gradient(points,u,axis=0,edge_order=2);second=np.gradient(first,u,axis=0,edge_order=2)
        k=abs(first[:,0]*second[:,1]-first[:,1]*second[:,0])/np.linalg.norm(first,axis=1)**3
        dense[name]=float(max(k[2:-2]))
    opposite=('omp_entry','shelf','omp_connector','floor') if candidate.source=='IMP' else ('imp_sloped','imp_vertical')
    error=max(float(np.max(abs(r['targets'][n]-e.baseline['targets'][n]))) for n in opposite)
    excess=max(dense[n]-check['curvature_max_1_m'][n] for n in dense)
    screen=geometry_screen(e,c);continuous_screen=screen['passed']
    return dict(**check,dense_curvature_1_m=dense,opposite_coordinate_error_m=error,
                chord_tolerance_m=screen['chord_tolerance_m'],tessellated_full_wall_valid=continuous_screen,
                geometry_screen=screen,
                passed_independent=bool(check['passed'] and excess<.02 and error<1e-12 and continuous_screen))


class Qualifier:
    def __init__(self,root,source,config,coarse):
        self.root,self.source,self.config=root,source,config
        self.reference=ReferencePoint(**config['reference']);self.originals={0:coarse};self.base_audits={}

    def original(self,index):
        if index not in self.originals:
            started=time.perf_counter()
            method,count,step,nodes,wall_step=LEVELS[index]
            self.originals[index]=prepare_boundary(self.root,self.source,self.reference,count=count,step=step,
                                                 nodes=nodes,wall_step=wall_step,method=method)
            print('VALIDATION_SETUP',self.source,method,count,'seconds',time.perf_counter()-started,flush=True)
        return self.originals[index]

    def shortlist(self,case,inherited=(),limit=8,excluded=()):
        eligible=[row for row in case.state['archive'] if Candidate.from_record(row['candidate']).identity not in excluded]
        candidates=[Candidate.from_record(row['candidate']) for row in diverse_shortlist(eligible,limit,scales(case.e))]
        for row in inherited:
            candidates.append(Candidate(self.source,case.mm,tuple(row['vector'])))
        candidates=list({c.identity:c for c in candidates if c.identity not in excluded}.values());rows=[]
        if not candidates:return rows
        fine=self.original(1);e=evaluator(fine,case.mm,self.config['backend'])
        for candidate in candidates:
            if should_stop(self.config,search=False):break
            r=evaluate_candidate(e,candidate);row=dict(candidate=candidate.record(),identity=candidate.identity,**e.describe(r))
            rows.append(row)
        rows.sort(key=lambda r:r.get('peak_ratio',1e9) if r['feasible'] else 1e9)
        save(case.folder/'shortlist_640.json',dict(epoch=time.time(),rows=rows,complete=True,global_optimum_proven=False))
        return [r for r in rows if r['feasible']]

    def validate(self,case,candidate):
        folder=case.folder/'validation'/candidate.identity;folder.mkdir(parents=True,exist_ok=True)
        target=folder/'validation.json'
        if target.exists():
            existing=load_json(target)
            if (existing.get('complete') and existing.get('fingerprint')==self.config['fingerprint']
                    and existing.get('qualification_policy')==QUALIFICATION_POLICY
                    and existing.get('candidate')==candidate.record()):return existing
        record=dict(candidate=candidate.record(),identity=candidate.identity,levels=[],complete=False,
                    fingerprint=self.config['fingerprint'],qualification_policy=QUALIFICATION_POLICY)
        for index,settings in enumerate(LEVELS):
            if should_stop(self.config,search=False):return record
            original=self.original(index);e=evaluator(original,case.mm,self.config['backend'])
            r=evaluate_candidate(e,candidate);selected=audit(e,r)
            if r['geometry_valid']:
                geom=geometric_audit(e,candidate,r);selected['audit_passed'] &= geom['passed_independent']
            else:geom={'passed_independent':False}
            if index not in self.base_audits:self.base_audits[index]=audit(original,original.baseline)
            base=self.base_audits[index]
            record['levels'].append(dict(method=settings[0],tubes=settings[1],trace_step_m=settings[2],nodes=settings[3],
                wall_step_m=settings[4],rows=dict(baseline=base,selected=selected),geometry_audit=geom))
            if index==0 and selected['audit_passed']:
                refined=BoundaryEvaluator(original.direct,original.metadata,self.reference,original.wall_step/2,
                                          thermal_paths=original.thermal_paths)
                ef=evaluator(refined,case.mm,self.config['backend']);rf=evaluate_candidate(ef,candidate)
                b=audit(refined,refined.baseline);s=audit(ef,rf)
                record['wall_only_refinement']=dict(rows=dict(baseline=b,selected=s),
                    max_peak_relative_change=max(abs(b['peak_q']/base['peak_q']-1),abs(s.get('peak_q',1e9)/selected['peak_q']-1)))
            save(target,record)
            if not selected['audit_passed'] or not base['audit_passed']:break
        record['verdict']=verdict(record['levels'],record.get('wall_only_refinement'))
        record['complete']=True;save(target,record)
        return record

    def finalize(self,cases):
        shortlists={};manifests={};stores={}
        for case in cases:
            stores[case.mm]=QualifiedStore(case.folder,context(self.config,self.source,case.mm))
            recovered=stores[case.mm].best()
            manifests[case.mm]=dict(**context(self.config,self.source,case.mm),global_optimum_proven=False,
                qualification_attempts=[],status='no_qualified_candidate_found',search_best=case.state['best'])
            path=case.folder/'manifest.json'
            if path.exists():
                previous=load_json(path)
                for key in ['fingerprint','source','radius_mm','reference']:
                    if previous.get(key)!=manifests[case.mm][key]:raise ValueError('finalization manifest belongs to a different case')
                if previous.get('status')=='qualified' and not previous.get('validation',{}).get('verdict',{}).get('resolution_screen_passed'):
                    raise ValueError('qualified incumbent has no passing validation')
                if (previous.get('status')=='qualified' and not should_stop(self.config,search=False)
                        and previous.get('qualification_policy')!=QUALIFICATION_POLICY):
                    raise ValueError('legacy qualified manifest: start a new campaign and import for revalidation')
                if (previous.get('status')=='qualified' and previous.get('qualification_policy')==QUALIFICATION_POLICY
                        and recovered is None):
                    raise ValueError('qualified manifest has no durable bundle')
                manifests[case.mm]=previous
            if recovered is not None:
                attempts=manifests[case.mm].get('qualification_attempts',[])
                manifests[case.mm]=dict(recovered,qualification_attempts=attempts,
                    search_best=case.state['best'])
            save(case.folder/'manifest.json',manifests[case.mm])

        # Protected historical geometries are independent of coarse/fine top-K.
        # Revalidate them first, round-robin across radii, before using the
        # remaining budget on newly discovered coarse finalists.
        for case in cases:
            protected=protected_candidates(self.root,self.source,case.mm) if hasattr(self,'root') else []
            for other in cases:
                if other.mm>=case.mm:
                    m=manifests[other.mm]
                    if m.get('status')=='qualified' and 'candidate' in m:
                        c=Candidate(self.source,case.mm,tuple(m['candidate']['vector']))
                        protected.append(dict(candidate=c.record(),origin='current-qualified-bundle'))
            shortlists[case.mm]=list({Candidate.from_record(r['candidate']).identity:r for r in protected}.values())

        def attempt(case,row):
            candidate=Candidate.from_record(row['candidate']);manifest=manifests[case.mm]
            # Same case/policy already has fully qualified evidence and arrays.
            if any(m['candidate']==candidate.record() for m in stores[case.mm].records()):return
            record=self.validate(case,candidate)
            manifest.setdefault('qualification_attempts',[]).append(dict(identity=candidate.identity,
                origin=row.get('origin','search-shortlist'),verdict=record.get('verdict'),complete=record.get('complete',False)))
            qualified=record.get('verdict',{}).get('resolution_screen_passed')
            if qualified and record['levels'][2]['rows']['selected']['peak_MW_m2']<manifest.get('selected',{}).get('peak_MW_m2',float('inf')):
                original=self.original(2);e=evaluator(original,case.mm,self.config['backend']);r=evaluate_candidate(e,candidate)
                cpu=evaluator(original,case.mm,'cpu');rc=evaluate_candidate(cpu,candidate)
                if not r['feasible'] or not rc['feasible']:raise ValueError('final CPU/GPU parity geometry failed')
                np.testing.assert_allclose(r['deposition']['q'],rc['deposition']['q'],atol=1e-12,rtol=2e-11)
                promoted=dict(manifest,qualification_policy=QUALIFICATION_POLICY,status='qualified',candidate=candidate.record(),
                    selected=e.describe(r),baseline=original.describe(original.baseline),validation=record,
                    cpu_gpu_peak_relative_difference=abs(r['peak_q']/rc['peak_q']-1))
                manifest=stores[case.mm].publish(promoted,lambda folder:write_arrays(folder,e,r,'selected_refined'))
                manifests[case.mm]=manifest
                base_folder=case.folder.parent
                if not (base_folder/'baseline_refined.npz').exists():write_arrays(base_folder,original,original.baseline,'baseline_refined')
            save(case.folder/'manifest.json',manifest)
            print('QUALIFIED',self.source,case.mm,manifest['status'],manifest.get('selected',{}).get('peak_MW_m2'),flush=True)

        seen={case.mm:set() for case in cases}
        for rank in range(max((len(v) for v in shortlists.values()),default=0)):
            for case in cases:
                if should_stop(self.config,search=False):return
                if rank>=len(shortlists[case.mm]):continue
                row=shortlists[case.mm][rank];seen[case.mm].add(Candidate.from_record(row['candidate']).identity)
                attempt(case,row)

        def propagate_qualified():
            # Promotions may differ from coarse state.best. Re-read the actual
            # qualified winners until no untested same/looser-radius pair remains.
            while not should_stop(self.config,search=False):
                pending=[]
                for case in cases:
                    for other in cases:
                        m=manifests[other.mm]
                        if other.mm<case.mm or m.get('status')!='qualified' or 'candidate' not in m:continue
                        c=Candidate(self.source,case.mm,tuple(m['candidate']['vector']))
                        if c.identity not in seen[case.mm]:
                            seen[case.mm].add(c.identity)
                            pending.append((case,dict(candidate=c.record(),origin='new-qualified-stricter-case')))
                if not pending:return
                for case,row in pending:
                    if should_stop(self.config,search=False):return
                    attempt(case,row)

        propagate_qualified()
        for case in cases:
            inherited=[other.state['best']['candidate'] for other in cases if other.mm>=case.mm and other.state['best']]
            rows=self.shortlist(case,inherited,limit=12,excluded=seen[case.mm]) if not should_stop(self.config,search=False) else []
            eligible=[r for r in rows if Candidate.from_record(r['candidate']).identity not in seen[case.mm]]
            shortlists[case.mm]=diverse_shortlist(eligible,4,scales(case.e))
        for rank in range(4):
            for case in cases:
                if should_stop(self.config,search=False):return
                if rank<len(shortlists[case.mm]):
                    row=shortlists[case.mm][rank];seen[case.mm].add(Candidate.from_record(row['candidate']).identity)
                    attempt(case,row)
        propagate_qualified()
