"""Crash-safe state and deterministic, resumable immediate-update DE."""
from contextlib import contextmanager
import hashlib,json,os,time
from pathlib import Path
import numpy as np


def replace_retry(source,target):
    for attempt in range(9):
        try:os.replace(source,target);return
        except PermissionError:
            if attempt==8:raise
            time.sleep(.01*2**attempt)


def public_save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    payload=json.dumps(value,allow_nan=False,indent=2)
    temp=path.with_suffix(path.suffix+'.tmp')
    for attempt in range(9):
        try:
            with temp.open('w',encoding='utf-8') as f:f.write(payload);f.flush();os.fsync(f.fileno())
            break
        except PermissionError:
            if attempt==8:raise
            time.sleep(.01*2**attempt)
    replace_retry(temp,path)


@contextmanager
def keep_awake():
    if os.name!='nt':
        yield;return
    import ctypes
    request=ctypes.windll.kernel32.SetThreadExecutionState
    request.argtypes=[ctypes.c_uint];request.restype=ctypes.c_uint
    if not request(0x80000001):raise OSError('could not inhibit automatic system sleep')
    try:yield
    finally:request(0x80000000)  # Display sleep remains allowed; no power-plan edit.


def atomic_save(path,value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    payload=json.dumps(value,allow_nan=False,separators=(',',':'))
    envelope=json.dumps(dict(sha256=hashlib.sha256(payload.encode()).hexdigest(),payload=payload))
    temp=path.with_suffix(path.suffix+'.tmp')
    with temp.open('w',encoding='utf-8') as f: f.write(envelope); f.flush(); os.fsync(f.fileno())
    if path.exists():
        try: _read(path)
        except (ValueError,OSError,KeyError): pass
        else: replace_retry(path,path.with_suffix(path.suffix+'.previous'))
    replace_retry(temp,path)


def _read(path):
    envelope=json.loads(Path(path).read_text(encoding='utf-8'))
    if hashlib.sha256(envelope['payload'].encode()).hexdigest()!=envelope['sha256']:
        raise ValueError('checkpoint checksum mismatch')
    return json.loads(envelope['payload'])


def read_state(path):
    path=Path(path)
    try: return _read(path)
    except (ValueError,OSError,KeyError): return _read(path.with_suffix(path.suffix+'.previous'))


@contextmanager
def process_lock(path):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    f=path.open('a+b')
    if path.stat().st_size==0: f.write(b'0'); f.flush()
    f.seek(0)
    try:
        if os.name=='nt':
            import msvcrt
            try: msvcrt.locking(f.fileno(),msvcrt.LK_NBLCK,1)
            except OSError as exc: raise RuntimeError('another process owns this campaign source') from exc
        else:
            import fcntl
            try: fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except OSError as exc: raise RuntimeError('another process owns this campaign source') from exc
        yield
    finally:
        f.close()  # OS releases the advisory lock even after an unclean exit.


def should_stop(config,now=None,search=True):
    return (time.time() if now is None else now)>=config['search_deadline_epoch' if search else 'deadline_epoch']


def new_population(population,seed):
    p=np.asarray(population,float)
    if p.ndim!=2 or len(p)<5 or not np.all(np.isfinite(p)): raise ValueError('five finite population members required')
    rng=np.random.default_rng(seed)
    return dict(population=p.tolist(),energies=[None]*len(p),index=0,generation=0,
                calls=0,rng=rng.bit_generator.state)


def propose(state,bounds):
    p=np.array(state['population']); i=state['index']; scores=state['energies']
    if scores[i] is None:return p[i].copy()
    rng=np.random.default_rng(); rng.bit_generator.state=state['rng']
    choices=np.delete(np.arange(len(p)),i); a,b,c=rng.choice(choices,3,replace=False)
    scale=rng.uniform(.35,.95); best=p[int(np.argmin(scores))]
    mutant=p[a]+scale*(p[b]-p[c]) if rng.random()<.6 else p[i]+scale*(best-p[i])+scale*(p[a]-p[b])
    cross=rng.random(p.shape[1])<rng.uniform(.5,.95); cross[rng.integers(p.shape[1])]=True
    trial=np.where(cross,mutant,p[i]); bounds=np.asarray(bounds)
    outside=(trial<bounds[:,0])|(trial>bounds[:,1])
    trial[outside]=rng.uniform(bounds[outside,0],bounds[outside,1])
    state['rng']=rng.bit_generator.state
    return trial


def accept(state,x,value):
    if not np.isfinite(value):raise ValueError('finite objective required')
    i=state['index']
    if state['energies'][i] is None or value<=state['energies'][i]:
        state['population'][i]=np.asarray(x,float).tolist(); state['energies'][i]=float(value)
    state['calls']+=1; state['index']=(i+1)%len(state['population'])
    if state['index']==0:state['generation']+=1
