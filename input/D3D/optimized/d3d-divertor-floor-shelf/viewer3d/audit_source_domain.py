"""Audit proposed 100/150/200 mm sources without launching paths outside material."""
import argparse,importlib,json,sys,types
from pathlib import Path
import numpy as np
from source_domain import physical_source_limit,scope_reason

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case-id');args=parser.parse_args()
    p=Path(__file__).resolve().parent
    catalog=json.loads((p/'catalog.json').read_text(encoding='utf-8'))
    entry=next(c for c in catalog['cases'] if c['id']==args.case_id) if args.case_id else next(c for c in catalog['cases'] if c['run_name']==p.parent.name)
    case_base=p.parent.parent/entry['run_name']
    archive=sorted(p.parent.parent.glob('migration-*/d3d_multisurface.py'))[-1].parent
    code_dir=case_base/'case_math' if (case_base/'case_math').is_dir() else archive
    package=types.ModuleType('_source_domain_case');package.__path__=[str(code_dir)];sys.modules[package.__name__]=package
    case=importlib.import_module(package.__name__+'.d3d_frozen_case')
    math=importlib.import_module(package.__name__+'.frozen_field')
    m=json.loads((case_base/'manifest.json').read_text(encoding='utf-8'))
    bundle=json.loads((p/entry['bundle_url']).read_text(encoding='utf-8'))
    quality_path=p/'quality'/(entry['id']+'.json');quality=json.loads(quality_path.read_text(encoding='utf-8'))
    if case.sha256(m['shot_path'])!=bundle['shot_sha256']:raise ValueError('Shot mismatch')
    field,lcfs,axis=case.load_case(m['shot_path'],bundle['time_ms'],'cubic');r_sep,_,_=case.separatrix(field,lcfs,axis)
    domain={k:physical_source_limit(v['wall'],r_sep,float(axis[1]),float(field.r[-1])) for k,v in bundle['models'].items()}
    physical=min(d['physical_max_offset_m'] for d in domain.values())
    upper=min(.2,float(field.r[-1])-r_sep)
    # Use small intervals with 4-point Gaussian quadrature. A single wide source
    # cell badly under-resolves the 1 mm profile and must never normalize the tail.
    edges=np.unique(np.r_[upper*np.linspace(0,1,4097)**2,physical,.1,.15])
    edges=edges[(edges>=0)&(edges<=upper)]
    mass={}
    for mm in [1,2,3,5,10]:
        source=math.exponential_source(field,r_sep,axis[1],edges,mm/1000)
        x=source.x;weights=source.weights
        mass[str(mm)]={'inside_public_wall':float(sum(weights[x<physical])),
                       'wall_to_100mm':float(sum(weights[(x>=physical)&(x<.1)])),
                       '100_to_150mm':float(sum(weights[(x>=.1)&(x<.15)])),
                       '150_to_200mm':float(sum(weights[x>=.15])),
                       'normalization_physical_domain_to_200mm':source.normalization}
    psi_start,psi_end=field.values([[r_sep,axis[1]],[r_sep+.1,axis[1]]])
    reasons={}
    for key,model in quality['models'].items():
        length={'separatrix_side':0.,'beyond_source_end':0.,'within_source_range':0.}
        for patch in model['patches']:
            if patch['status']!='outside_source_scope':continue
            point=(np.array(patch['a'])+patch['b'])/2
            reason=scope_reason(float(field.values([point])[0]),psi_start,psi_end)
            patch['scope_reason']=reason
            length[reason]+=float(np.linalg.norm(np.array(patch['b'])-patch['a']))
        reasons[key]=length
    audit={'case_id':bundle['id'],'shot_sha256':bundle['shot_sha256'],'geometry_sha256':bundle['source_geometry_sha256'],
           'configured_source_max_mm':100.,'physical_source_max_mm':physical*1000,'domains':domain,
           'proposed_extensions_mm':{str(mm):{'physically_valid':bool(mm/1000<=physical),'reason':'outside_public_wall' if mm/1000>physical else 'within_wall'} for mm in [100,150,200]},
           'tail_quadrature_upper_mm':upper*1000,'hypothetical_profile_mass_fractions':mass,
           'outside_target_length_m_by_reason':reasons,
           'interpretation':'target scope reasons and the physical wall limit are recorded separately; outside-wall source starts are excluded',
           'source_cutoff_flagged_target_length_m':sum(x['beyond_source_end'] for x in reasons.values()),
           'invalid_extension_paths_launched':0}
    quality['source_domain_audit']=audit
    quality_path.write_text(json.dumps(quality,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')
    (p/'quality'/(bundle['id']+'-source-domain.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
