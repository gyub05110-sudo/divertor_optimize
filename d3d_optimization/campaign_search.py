"""Regional multi-start search with exact geometric screening and checkpoints."""
import json,time
from pathlib import Path
from collections import OrderedDict
import numpy as np
from .campaign_state import atomic_save,read_state,new_population,propose,accept,should_stop
from .campaign_geometry import Candidate,evaluator,evaluate_candidate
from .imp_tangent import project_imp
from .broadening_study import load_json
from .campaign_state import public_save as save
from .campaign_incumbents import protected_candidates
from run_tangent_curvature_study import project_floor


def bounds(source):return np.array([[-.1,.1]]*(8 if source=='IMP' else 16)+([[0,1]] if source=='IMP' else []))


def scales(e):
    if e.metadata['source']=='IMP':return np.r_[np.full(8,.01),.15]
    result=np.concatenate([np.full(4,min(.012,e.lengths[n]**2/(100*e.rmin))) for n in ('omp_entry','shelf','omp_connector','floor')])
    result[12:]=.025;return result


def warm_vectors(root,source,mm):
    for row in protected_candidates(root,source,mm):
        yield np.array(row['candidate']['vector'])
    pilot=root/'campaign_warm_candidates.json'
    if pilot.exists():
        for candidate in load_json(pilot)['candidates']:
            if candidate['source']==source and candidate['radius_mm']>=mm:
                yield np.array(candidate['vector'])
    if source=='OMP':
        for other in [100,50,20]:
            if other<mm:continue
            path=root/'results/omp-tangent-curvature-20261029'/f'R{other:03d}'/'manifest.json'
            if path.exists():yield np.array(load_json(path)['coefficients_m'][:16])
    else:
        path=root/'results/regional-20261028/qualified/IMP/manifest.json'
        if path.exists():yield np.r_[load_json(path)['candidates']['selected'][16:],.5]


def seed_population(e,root,config,seed,archive=(),inherited=()):
    rng=np.random.default_rng(seed); source=e.metadata['source']; mm=round(e.rmin*1000)
    vectors=[np.array(row['candidate']['vector']) for row in archive[:12]]
    vectors+=list(inherited)+list(warm_vectors(root,source,mm))
    for j in range(config['projected_seeds']):
        if should_stop(config):break
        if source=='IMP':
            target=np.zeros(8) if j==0 else rng.uniform(-.025,.025,8)
            x=project_imp(e,target,.5 if j==0 else rng.uniform(.3,.7),return_best=True)
        else:
            target=np.zeros(4) if j==0 else rng.uniform(-.08,.08,4)
            floor=project_floor(e,target); x=None if floor is None else np.r_[np.zeros(12),floor]
        if x is not None:vectors.append(x)
    if not vectors:vectors=[np.r_[np.zeros(8),.5] if source=='IMP' else np.zeros(16)]
    # Independent global restarts remain within the explicit bounds. Projected
    # seeds and signed perturbations explore useful thin feasible components.
    initial=list(vectors[:config['population']]); scale=scales(e); b=bounds(source)
    while len(initial)<config['population']:
        i=len(initial); center=vectors[i%len(vectors)]
        amplitude=[.05,.2,.7,2.][i%4]
        x=center+rng.uniform(-1,1,len(scale))*scale*amplitude
        if i%13==0:x=rng.uniform(b[:,0],b[:,1])
        initial.append(np.clip(x,b[:,0],b[:,1]))
    return new_population(initial,seed)


def retain(archive,row,scale):
    by_id={r['identity']:r for r in archive};by_id[row['identity']]=row
    ordered=sorted(by_id.values(),key=lambda r:r['objective']);result=ordered[:16]
    for item in ordered[16:]:
        x=np.array(item['candidate']['vector'])
        if all(np.linalg.norm((x-np.array(r['candidate']['vector']))/scale)>.2 for r in result):result.append(item)
        if len(result)>=48:break
    return result


class SearchCase:
    def __init__(self,original,mm,folder,root,config):
        self.original,self.mm,self.folder,self.root,self.config=original,mm,Path(folder),root,config
        self.source=original.metadata['source'];self.e=evaluator(original,mm,config['backend'])
        self.cpu=evaluator(original,mm,'cpu');self.memo=OrderedDict();self.path=self.folder/'state.json'
        self._checkpoint_time=0.
        self.folder.mkdir(parents=True,exist_ok=True)
        if self.path.exists() or self.path.with_suffix('.json.previous').exists():
            self.state=read_state(self.path)
            if self.state['fingerprint']!=config['fingerprint'] or self.state['source']!=self.source or self.state['radius_mm']!=mm:
                raise ValueError('checkpoint belongs to different inputs or case')
        else:
            self.state=dict(fingerprint=config['fingerprint'],source=self.source,radius_mm=mm,
                restart=0,total_proposals=0,physics_evaluations=0,feasible_evaluations=0,archive=[],best=None,
                best_infeasible=None,active_seconds=0.,last_improvement=0,parity=None,optimizer=None)

    def initialize(self,inherited=()):
        started=time.monotonic()
        seed=self.config['seed']+self.mm+(0 if self.source=='IMP' else 10000)+100003*self.state['restart']
        self.state['optimizer']=seed_population(self.e,self.root,self.config,seed,self.state['archive'],inherited)
        self.state['active_seconds']+=time.monotonic()-started
        self.state['last_improvement']=self.state['total_proposals'];self.checkpoint()

    def checkpoint(self,force=True):
        if not force and time.monotonic()-self._checkpoint_time<2:return
        self.state['updated_epoch']=time.time();atomic_save(self.path,self.state)
        self._checkpoint_time=time.monotonic()
        # Human-readable lightweight progress, separate from the restart state.
        save(self.folder/'progress.json',{k:self.state[k] for k in ['source','radius_mm','restart','total_proposals',
             'physics_evaluations','feasible_evaluations','active_seconds','updated_epoch','best','best_infeasible','parity']}|
             dict(generation=self.state['optimizer']['generation'] if self.state['optimizer'] else 0,global_optimum_proven=False))

    def score(self,x,stage):
        candidate=Candidate(self.source,self.mm,tuple(map(float,x)));key=candidate.identity
        if key in self.memo:
            row=dict(self.memo[key],cached=True,stage=stage)
            self.memo.move_to_end(key);return row
        start=time.perf_counter();c=candidate.bind(self.e);check=self.e.constraints(c);thermal=False
        if not check['passed']:
            violation=max(0.,self.e.rmin*max(check['curvature_max_1_m'].values())-1)
            if self.source=='OMP':violation+=max(0.,(self.e.gap_min-check['gap_m'])/.03)
            lengths,_=self.e.precise_lengths(c);violation+=max(0.,max(lengths.values())-1.2)
            row=dict(feasible=False,objective=10+violation,shape_constraints=check,errors=['continuous_shape_constraint_failed'])
        else:
            r=evaluate_candidate(self.e,candidate);thermal='deposition' in r
            row=self.e.describe(r)
            if r['feasible']:row['objective']=float(r['peak_q']/self.e.baseline_peak)
            elif not r['geometry_valid']:
                # Length is a continuous geometric residual, not a jump to a
                # worse tier than an almost curvature-feasible proposal.
                length_excess=max(0.,max(r['precise_length_ratios'].values())-1.2)
                collision=0. if 'material_length_exceeds_limit' in r['errors'] else 1.
                row['objective']=10+length_excess+collision
            else:row['objective']=30+min(r['penalty'],1000)/1000
            if r['feasible'] and self.state['parity'] is None:
                # Compare the actual GPU density array, not only the CPU-refined peak.
                t=time.perf_counter();cr=evaluate_candidate(self.cpu,candidate);cpu_seconds=time.perf_counter()-t
                if not cr['feasible']:raise ValueError('CPU candidate failed GPU parity gate')
                np.testing.assert_allclose(r['deposition']['points'],cr['deposition']['points'],atol=1e-13,rtol=0)
                np.testing.assert_allclose(r['deposition']['q'],cr['deposition']['q'],atol=1e-12,rtol=2e-11)
                peak=abs(r['peak_q']/cr['peak_q']-1)
                if peak>1e-10:raise ValueError('CPU/GPU peak mismatch')
                density=float(np.max(abs(r['deposition']['q']-cr['deposition']['q']))/max(cr['deposition']['q']))
                quadrature=abs(r['deposition']['quadrature_integrated_power']-cr['deposition']['quadrature_integrated_power'])
                self.state['parity']=dict(candidate=candidate.record(),peak_relative_difference=peak,
                    density_scaled_max_difference=density,quadrature_absolute_difference=quadrature,
                    cpu_seconds=cpu_seconds,gpu_seconds=r['elapsed_seconds'],passed=True)
                save(self.folder/'backend_preflight.json',self.state['parity'])
        row.update(candidate=candidate.record(),identity=key,cached=False,stage=stage,
                   physics_evaluated=thermal,evaluation_seconds=time.perf_counter()-start)
        self.memo[key]=row
        while len(self.memo)>2000:self.memo.popitem(last=False)
        return row

    def step(self):
        s=self.state;opt=s['optimizer'];b=bounds(self.source)
        # Interleave local signed polishing and independent full-bound proposals
        # without disturbing DE's saved target index and RNG sequence.
        stage='differential_evolution';is_de=True
        if s['best'] and s['total_proposals']%19==0:
            stage='coordinate_polish';is_de=False;x=np.array(s['best']['candidate']['vector'])
            j=(s['total_proposals']//19)%(2*len(x));factor=[.25,.08,.02,.005][(s['total_proposals']//(38*len(x)))%4]
            x[j//2]+=(-1 if j%2 else 1)*factor*scales(self.e)[j//2];x=np.clip(x,b[:,0],b[:,1])
        else:x=propose(opt,b)
        row=self.score(x,stage)
        if is_de:accept(opt,x,row['objective'])
        s['total_proposals']+=1
        if not row['cached']:
            s['physics_evaluations']+=int(row['physics_evaluated'])
            s['feasible_evaluations']+=int(row['feasible'])
        if row['feasible']:
            s['archive']=retain(s['archive'],row,scales(self.e))
            if s['best'] is None or row['objective']<s['best']['objective']:
                s['best']=row;s['last_improvement']=s['total_proposals']
                save(self.folder/'best_checkpoint.json',row)
        elif s['best_infeasible'] is None or row['objective']<s['best_infeasible']['objective']:
            s['best_infeasible']=row
        compact={k:row[k] for k in ['identity','candidate','objective','feasible','cached','stage','evaluation_seconds']}
        compact.update(sequence=s['total_proposals'],restart=s['restart'],epoch=time.time(),physics_evaluated=row['physics_evaluated'])
        if not row['feasible']:
            compact['errors']=row.get('errors',[])
            if 'geometry_screen' in row:
                compact['geometry_errors']=row['geometry_screen']['errors']
        with (self.folder/'evaluations.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(compact,allow_nan=False)+'\n')
        self.checkpoint(force=False)

    def run_slice(self,seconds,inherited=()):
        start=time.monotonic()
        if self.state['optimizer'] is None:self.initialize(inherited)
        while time.monotonic()-start<seconds and not should_stop(self.config) and not (self.folder.parents[1]/'STOP').exists():
            step_start=time.monotonic();self.step();self.state['active_seconds']+=time.monotonic()-step_start
            opt=self.state['optimizer'];stalled=self.state['total_proposals']-self.state['last_improvement']
            if opt['generation']>=100 or (opt['generation']>=15 and stalled>25*self.config['population']):
                self.state['restart']+=1
                if self.state['best_infeasible'] and not self.state['archive']:
                    inherited=list(inherited)+[self.state['best_infeasible']['candidate']['vector']]
                self.initialize(inherited)
        self.checkpoint()
