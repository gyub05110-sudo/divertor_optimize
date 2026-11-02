"""Append-only qualified bundles, separate from disposable search rankings.

One campaign source owns the writer lock. A directory becomes visible only
after all artifacts and checksums have been flushed. Legacy records provide
candidate geometry, never authority to skip current qualification.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import uuid
from .campaign_geometry import Candidate, GEOMETRY_POLICY
from .campaign_state import public_save

QUALIFICATION_POLICY = 'regional-multiresolution-shared-geometry.v2'
CONTEXT_KEYS = ('fingerprint','source','radius_mm','reference','qualification_policy')


def digest_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def digest_record(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def context(config,source,mm):
    return dict(fingerprint=config['fingerprint'],reference=config['reference'],source=source,
                radius_mm=mm,qualification_policy=QUALIFICATION_POLICY)


def validate_manifest(manifest,expected):
    if any(manifest.get(k)!=expected[k] for k in CONTEXT_KEYS):
        raise ValueError('qualified bundle context mismatch')
    c=Candidate.from_record(manifest['candidate']);v=manifest.get('validation',{})
    if c.source!=expected['source'] or c.radius_mm!=expected['radius_mm']:
        raise ValueError('qualified candidate case mismatch')
    levels=v.get('levels',[])
    passed=(manifest.get('status')=='qualified' and v.get('complete') is True
            and v.get('fingerprint')==expected['fingerprint']
            and v.get('qualification_policy')==QUALIFICATION_POLICY
            and v.get('candidate')==c.record() and v.get('identity')==c.identity
            and v.get('verdict',{}).get('resolution_screen_passed') is True
            and [(x.get('method'),x.get('tubes')) for x in levels]==
                [('cubic',320),('cubic',640),('cubic',1280),('quintic',640)])
    for level in levels:
        g=level.get('geometry_audit',{});s=g.get('geometry_screen',{})
        passed &= (g.get('passed_independent') is True and s.get('passed') is True
                   and s.get('policy')==GEOMETRY_POLICY
                   and all(level.get('rows',{}).get(k,{}).get('audit_passed') is True for k in ('selected','baseline')))
    wall=v.get('wall_only_refinement',{})
    passed &= all(wall.get('rows',{}).get(k,{}).get('audit_passed') is True for k in ('selected','baseline'))
    for key in ('selected','baseline'):
        peak=manifest.get(key,{}).get('peak_MW_m2',float('nan'))
        passed &= math.isfinite(peak) and peak>0
        if len(levels)==4:
            validated_peak=levels[2].get('rows',{}).get(key,{}).get('peak_MW_m2',float('nan'))
            passed &= math.isclose(peak,validated_peak,rel_tol=2e-11,abs_tol=1e-12)
    parity=manifest.get('cpu_gpu_peak_relative_difference',float('nan'))
    passed &= math.isfinite(parity) and 0<=parity<=1e-10
    if not passed:raise ValueError('incomplete or incompatible qualification evidence')


class QualifiedStore:
    def __init__(self,folder,expected):
        self.folder=Path(folder);self.root=self.folder/'qualified';self.context=expected

    def records(self):
        records=[]
        if not self.root.exists():return records
        for folder in sorted(self.root.iterdir()):
            if folder.name.startswith('.'):continue  # Unpublished transaction.
            if not folder.is_dir():continue
            try:
                checks=json.loads((folder/'checksums.json').read_text(encoding='utf-8'))
                if digest_record(checks)!=folder.name:raise ValueError('bundle index checksum mismatch')
                for name,value in checks.items():
                    if Path(name).name!=name or digest_file(folder/name)!=value:
                        raise ValueError('bundle artifact checksum mismatch')
                if 'manifest.json' not in checks or 'selected_refined.npz' not in checks:
                    raise ValueError('incomplete bundle checksum index')
                m=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
                validate_manifest(m,self.context)
            except (OSError,KeyError,json.JSONDecodeError) as exc:
                raise ValueError('qualified bundle missing/corrupt: '+str(folder)) from exc
            m['artifact_bundle']=folder.relative_to(self.folder).as_posix()
            m['artifact_sha256']=checks
            records.append(m)
        return records

    def best(self):
        return min(self.records(),key=lambda m:m['selected']['peak_MW_m2'],default=None)

    def publish(self,manifest,writer):
        validate_manifest(manifest,self.context)
        self.root.mkdir(parents=True,exist_ok=True)
        staging=self.root/('.pending-'+uuid.uuid4().hex);staging.mkdir()
        writer(staging)
        if not (staging/'selected_refined.npz').is_file():raise ValueError('missing selected arrays')
        record={k:v for k,v in manifest.items() if k not in ('artifact_bundle','artifact_sha256')}
        public_save(staging/'manifest.json',record)
        files=[p for p in staging.iterdir() if p.is_file()]
        for path in files:
            with path.open('r+b') as f:os.fsync(f.fileno())
        checks={p.name:digest_file(p) for p in files}
        public_save(staging/'checksums.json',checks)
        destination=self.root/digest_record(checks)
        if not destination.exists():os.replace(staging,destination)
        # If interrupted after the directory publish, restart discovers the
        # complete bundle even when the convenience manifest was not updated.
        return next(m for m in self.records() if m['artifact_bundle']==destination.relative_to(self.folder).as_posix())


def protected_candidates(root,source,mm):
    path=Path(root)/'campaign_protected_candidates.json'
    if not path.exists():return []
    result={}
    for row in json.loads(path.read_text(encoding='utf-8'))['candidates']:
        old=Candidate.from_record(row['candidate'])
        if old.source!=source or old.radius_mm<mm:continue
        candidate=Candidate(source,mm,old.vector)
        result.setdefault(candidate.identity,dict(candidate=candidate.record(),origin=row['origin'],
            origin_sha256=row['origin_sha256'],historical_peak=row.get('historical_peak'),
            geometry_conversion=row.get('geometry_conversion'),
            requires_current_validation=True))
    return sorted(result.values(),key=lambda r:r['historical_peak'] if r['historical_peak'] is not None else float('inf'))


def collect_protected(root,warm_study=None):
    """Freeze qualified historical geometries even if progress.best lost them."""
    paths=[];root=Path(root)
    if warm_study:
        paths += [Path(warm_study)/s/f'R{mm:03d}'/'manifest.json' for s in ('IMP','OMP') for mm in (20,50,100)]
        paths += sorted(Path(warm_study).glob('*/R*/qualified/*/manifest.json'))
    paths += sorted((root/'results').glob('*/IMP/R*/manifest.json'))
    paths += sorted((root/'results').glob('*/OMP/R*/manifest.json'))
    paths += sorted((root/'results').glob('*/*/R*/qualified/*/manifest.json'))
    paths += sorted((root/'results/omp-tangent-curvature-20261029').glob('R*/manifest.json'))
    paths += [root/'results/regional-20261028/qualified/IMP/manifest.json']
    rows=[]
    for path in dict.fromkeys(paths):
        if not path.exists():continue
        m=json.loads(path.read_text(encoding='utf-8'))
        legacy_imp=(path.parent.name=='IMP' and path.parent.parent.name=='qualified'
                    and 'selected' in m.get('candidates',{}))
        if m.get('status')!='qualified' and not ('coefficients_m' in m and 'qualification' in m) and not legacy_imp:continue
        candidate=m.get('candidate')
        conversion=None
        if legacy_imp:
            candidate=dict(source='IMP',radius_mm=100,vector=list(m['candidates']['selected'][16:])+[.5])
            # All radii receive this migration seed. Its curve changes when G1
            # endpoints are imposed; the old peak is NOT attached to this shape.
            conversion='legacy normal-displacement coefficients mapped to G1 with junction fraction 0.5; unvalidated at every radius'
        if candidate is None and 'coefficients_m' in m:
            candidate=dict(source='OMP',radius_mm=int(path.parent.name[1:]),vector=m['coefficients_m'][:16])
        if candidate is None:continue
        c=Candidate.from_record(candidate)
        rows.append(dict(candidate=c.record(),origin=str(path),origin_sha256=digest_file(path),
                         historical_peak=None if conversion else m.get('selected',{}).get('peak_MW_m2'),
                         geometry_conversion=conversion))
    return rows
