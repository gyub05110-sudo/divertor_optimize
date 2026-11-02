"""Launch/resume a frozen, two-source GPU search with an absolute 72-hour cap."""
import argparse,hashlib,json,os,shutil,sys,time,subprocess,threading,traceback,platform
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import numpy as np
import scipy,shapely
from d3d_optimization.broadening_study import load_json,digest
from d3d_optimization.campaign_state import public_save as save
from d3d_optimization.boundary_evaluator import prepare_boundary
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.campaign_state import process_lock,should_stop,keep_awake
from d3d_optimization.campaign_search import SearchCase
from d3d_optimization.campaign_validation import Qualifier
from d3d_optimization.campaign_incumbents import collect_protected,QUALIFICATION_POLICY
from optimize_boundary_transport import capture
from optimize_regional_transport import write_arrays

ROOT=Path(__file__).resolve().parent


def runtime_identity(backend):
    identity=dict(executable=str(Path(sys.executable).resolve()),python=platform.python_version(),
                  numpy=np.__version__,scipy=scipy.__version__,shapely=shapely.__version__)
    if backend=='gpu':
        from d3d_optimization.gpu_footprint import cuda
        cp=cuda();properties=cp.cuda.runtime.getDeviceProperties(0)
        name=properties['name'];name=name.decode() if isinstance(name,bytes) else name
        identity.update(cupy=cp.__version__,cuda_runtime=cp.cuda.runtime.runtimeGetVersion(),
                        cuda_driver=cp.cuda.runtime.driverGetVersion(),gpu=name)
    return identity


def snapshot(output,settings):
    if output.exists():raise FileExistsError('new campaign requires a new output folder')
    extra=[str(p.relative_to(ROOT)).replace('\\','/') for pattern in ['*.py','d3d_optimization/*.py','verification/*.py'] for p in ROOT.glob(pattern)]
    extra+=['AGENTS.md','docs/superpowers/plans/2026-10-29-imp-global-search.md']
    extra+=['results/regional-20261028/qualified/IMP/manifest.json']
    extra += [f'results/omp-tangent-curvature-20261029/R{mm:03d}/manifest.json' for mm in [20,50,100]]
    capture(output,extra);files=load_json(output/'provenance.json')['files'];frozen=output/'snapshot'
    for name in files:
        destination=frozen/name;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,destination)
    deps=ROOT/'.local_dependencies/d3d-fixed-equilibrium'
    if deps.exists():
        # Preservation hashes already copied the pure-Python dependency files.
        # Add its binary/runtime files into that same new snapshot directory.
        shutil.copytree(deps,frozen/'.local_dependencies/d3d-fixed-equilibrium',dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for p in (frozen/'.local_dependencies').rglob('*'):
            if p.is_file():files[str(p.relative_to(frozen)).replace('\\','/')]=digest(p)
    if settings['warm_study']:
        warm=Path(settings['warm_study']);rows=[];origins={}
        for source in ['IMP','OMP']:
            for mm in [20,50,100]:
                path=warm/source/f'R{mm:03d}'/'progress.json'
                if path.exists():
                    record=load_json(path)
                    if record.get('best'):
                        rows.append(record['best']['candidate']);origins[str(path)]=digest(path)
        destination=frozen/'campaign_warm_candidates.json'
        save(destination,dict(candidates=rows,origins=origins,source_campaign_sha256=digest(warm/'campaign.json'),
                              note='warm geometry only; all constraints and physics are evaluated again'))
        files[destination.name]=digest(destination)
    protected=frozen/'campaign_protected_candidates.json'
    save(protected,dict(candidates=collect_protected(ROOT,settings['warm_study']),
         note='protected geometry candidates; current inputs and policy always require revalidation'))
    files[protected.name]=digest(protected)
    save(output/'frozen_files.json',files)
    now=time.time();duration=settings['hours']*3600
    reserve=min(settings['validation_reserve_seconds'],duration*.25)
    config=dict(schema='regional-global-campaign.v2',qualification_policy=QUALIFICATION_POLICY,created_epoch=now,deadline_epoch=now+duration,
        search_deadline_epoch=now+duration-reserve,reference=asdict(ReferencePoint()),**settings,
        frozen_root=str(frozen),original_root=str(ROOT),frozen_files_sha256=digest(output/'frozen_files.json'),
        executable=sys.executable,python=platform.python_version(),numpy=np.__version__,
        objective='minimum global wall peak at fixed independent source power',global_optimum_proven=False)
    config['runtime']=runtime_identity(settings['backend'])
    if settings['search_seconds'] is not None:
        config['search_deadline_epoch']=min(config['search_deadline_epoch'],now+settings['search_seconds'])
    config['fingerprint']=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    save(output/'campaign.json',config)
    return config


def verify_frozen(output):
    config=load_json(output/'campaign.json');identity=config.pop('fingerprint')
    if hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()!=identity:raise ValueError('campaign configuration changed')
    config['fingerprint']=identity
    if runtime_identity(config['backend'])!=config['runtime']:raise ValueError('runtime changed; cached energies cannot be resumed')
    if ROOT.resolve()!=Path(config['frozen_root']).resolve():raise ValueError('worker must execute from frozen snapshot')
    if digest(output/'frozen_files.json')!=config['frozen_files_sha256']:raise ValueError('frozen manifest changed')
    for name,value in load_json(output/'frozen_files.json').items():
        if not (ROOT/name).is_file() or digest(ROOT/name)!=value:raise ValueError('frozen file changed: '+name)
    return config


def watchdog(output,source,deadline):
    while time.time()<deadline:time.sleep(min(30,max(.01,deadline-time.time())))
    try:save(output/source/'hard_deadline.json',dict(epoch=time.time(),status='hard_deadline_stop',pid=os.getpid()))
    finally:os._exit(124)


def worker(output,source,finalize_only=False):
    config=verify_frozen(output);folder=output/source;folder.mkdir(parents=True,exist_ok=True)
    with process_lock(folder/'worker.lock'),keep_awake():
        if should_stop(config,search=False):raise RuntimeError('campaign deadline already expired')
        threading.Thread(target=watchdog,args=(output,source,config['deadline_epoch']),daemon=True).start()
        save(folder/'worker_status.json',dict(status='starting',pid=os.getpid(),epoch=time.time(),deadline_epoch=config['deadline_epoch']))
        original=prepare_boundary(ROOT,source,ReferencePoint(**config['reference']))
        write_arrays(folder,original,original.baseline,'baseline')
        save(folder/'baseline.json',original.describe(original.baseline))
        cases=[SearchCase(original,mm,folder/f'R{mm:03d}',ROOT,config) for mm in [100,50,20]]
        qualifier=Qualifier(ROOT,source,config,original);last_qualification=time.time()
        for case in cases:
            if should_stop(config) or finalize_only:break
            if case.state['optimizer'] is None:case.initialize()
        try:
            while not finalize_only and not should_stop(config) and not (output/'STOP').exists():
                for case in sorted(cases,key=lambda c:c.state['active_seconds']):
                    if should_stop(config) or (output/'STOP').exists():break
                    inherited=[other.state['best']['candidate']['vector'] for other in cases
                               if other.mm>=case.mm and other.state['best']]
                    save(folder/'worker_status.json',dict(status='searching',pid=os.getpid(),radius_mm=case.mm,
                         epoch=time.time(),deadline_epoch=config['deadline_epoch']))
                    case.run_slice(config['slice_seconds'],inherited)
                    print('SLICE',source,case.mm,'proposals',case.state['total_proposals'],
                          'best',case.state['best']['objective'] if case.state['best'] else None,flush=True)
                if time.time()-last_qualification>config['qualify_interval_seconds'] and not should_stop(config):
                    for case in cases:
                        if should_stop(config):break
                        qualifier.shortlist(case,limit=6)
                    last_qualification=time.time()
            save(folder/'worker_status.json',dict(status='validating',pid=os.getpid(),epoch=time.time()))
            qualifier.finalize(cases)
            save(folder/'worker_status.json',dict(status='completed',pid=os.getpid(),epoch=time.time(),
                 reason='stop_file' if (output/'STOP').exists() else 'search_budget_finished'))
        except BaseException as exc:
            save(folder/'worker_status.json',dict(status='failed',pid=os.getpid(),epoch=time.time(),error=repr(exc),traceback=traceback.format_exc()))
            raise


def launch(output,sources=('IMP','OMP')):
    config=load_json(output/'campaign.json');children=[]
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
        TEMP='D:/DivertorOptimization/tmp',TMP='D:/DivertorOptimization/tmp',PYTHONUNBUFFERED='1',
        CUPY_CACHE_DIR='D:/DivertorOptimization/runtime/cupy-cache')
    for source in sources:
        folder=output/source;folder.mkdir(parents=True,exist_ok=True)
        with (folder/'stdout.log').open('ab') as out,(folder/'stderr.log').open('ab') as err:
            p=subprocess.Popen([config['executable'],str(Path(config['frozen_root'])/'run_global_campaign.py'),
                '--mode','worker','--output',str(output),'--source',source],cwd=config['frozen_root'],env=env,
                stdout=out,stderr=err,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        children.append(dict(source=source,pid=p.pid,launched_epoch=time.time()))
    save(output/f'launch_{int(time.time())}.json',dict(children=children,deadline_epoch=config['deadline_epoch']))
    print(json.dumps(dict(output=str(output),children=children,
        deadline_utc=datetime.fromtimestamp(config['deadline_epoch'],timezone.utc).isoformat())),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=['create','launch','worker','finalize'],default='create')
    p.add_argument('--source',choices=['IMP','OMP','both'],default='both');p.add_argument('--hours',type=float,default=72.)
    p.add_argument('--population',type=int,default=64);p.add_argument('--projected-seeds',type=int,default=12)
    p.add_argument('--slice-seconds',type=float,default=300.);p.add_argument('--seed',type=int,default=20261029)
    p.add_argument('--backend',choices=['cpu','gpu'],default='gpu');p.add_argument('--qualify-interval-seconds',type=float,default=7200.)
    p.add_argument('--search-seconds',type=float,default=None,help='optional shorter search budget for end-to-end pilots')
    p.add_argument('--validation-reserve-seconds',type=float,default=3600.)
    p.add_argument('--warm-study',default=None,help='previous study: coarse seeds plus protected qualified geometries for revalidation')
    a=p.parse_args();a.output=a.output.resolve()
    if not 0<a.hours<=72 or not 5<=a.population<=512 or not 1<=a.projected_seeds<=64 or not 1<=a.slice_seconds<=3600:
        p.error('positive duration <=72h, population 5..512, projected seeds 1..64, slice 1..3600s required')
    if (not np.isfinite(a.validation_reserve_seconds) or a.validation_reserve_seconds<=0 or
        not np.isfinite(a.qualify_interval_seconds) or a.qualify_interval_seconds<=0 or
        (a.search_seconds is not None and (not np.isfinite(a.search_seconds) or a.search_seconds<=0)) or
        not 0<=a.seed<2**32):p.error('positive finite timing settings and a 32-bit nonnegative seed required')
    if a.mode=='create':
        snapshot(a.output,{k:v for k,v in vars(a).items() if k not in ('output','source','mode')});launch(a.output)
    elif a.mode=='launch':launch(a.output,('IMP','OMP') if a.source=='both' else (a.source,))
    else:
        if a.source=='both':p.error('one worker source required')
        worker(a.output,a.source,finalize_only=a.mode=='finalize')
