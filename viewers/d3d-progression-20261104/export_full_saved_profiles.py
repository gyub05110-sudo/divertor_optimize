"""Extend existing baseline/refined profiles to full-wall coverage, without a solve."""
from pathlib import Path
import json,hashlib
import numpy as np
from stage_profile_data import make_profile

HERE=Path(__file__).resolve().parent;DIST=HERE/'dist'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    data=read(DIST/'data.json');audit=read(HERE/'data-audit.json')
    inputs=[Path(p) for p in audit['input_sha256'] if p.endswith('.npz')]
    records=[(data['baseline']['profiles'][s],s,None,data['baseline']['wall'],data['baseline']['refined'][s]) for s in ['IMP','OMP']]
    records += [(s['final']['profile'],s['source'],s['final']['id'],s['final']['wall'],s['final']['peak']) for s in data['series'] if s['final']]
    if len(records)!=len(inputs):raise ValueError('source inventory mismatch')
    entries={}
    for (old,source,identity,wall,peak),path in zip(records,inputs):
        if sha(path)!=audit['input_sha256'][str(path)]:raise ValueError('saved source changed')
        # Shape + peak + source name also check the mapping from original exporter order.
        if source.lower() not in str(path).lower():raise ValueError('source mapping mismatch')
        with np.load(path) as z:p=make_profile(z,wall,peak,data['reference']['power_W'])
        p.update(source=source,identity=identity,donors=1280,provenance='saved-refined-evaluation')
        relative='full-saved-profiles/'+Path(old).name;target=DIST/relative;target.parent.mkdir(exist_ok=True)
        target.write_text(json.dumps(p,separators=(',',':'),allow_nan=False),encoding='utf-8')
        entries[old]=dict(profile=relative,source=source,identity=identity,sha256=sha(target),input=str(path),input_sha256=sha(path))
    index=dict(data_sha256=sha(DIST/'data.json'),profiles=entries)
    (DIST/'full-saved-profiles.json').write_text(json.dumps(index,separators=(',',':')),encoding='utf-8')
    print('Full saved profiles:',len(entries))
if __name__=='__main__':run()
