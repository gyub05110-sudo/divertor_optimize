"""Fail-closed fingerprints of the evaluator and its actual geometry inputs."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

REQUIRED=(
    'd3d_optimization/__init__.py', 'd3d_optimization/case.py',
    'd3d_optimization/frozen_field.py', 'd3d_optimization/frozen_mapping.py',
    'd3d_optimization/multisurface.py', 'd3d_optimization/shape.py',
    'd3d_optimization/independent_sources.py', 'd3d_optimization/continuous_boundary.py',
    'optimize_independent_sources.py', 'baseline_review/baseline_review.json',
    'input/D3D/d3d/d3d_shot_1a3155370c.parquet',
    'input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json',
    'input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json',
)

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def capture(root, run, phase='before_search'):
    root,run=Path(root),Path(run)
    if phase not in ('before_search','during_existing_search'):
        raise ValueError('explicit provenance capture phase required')
    destination=run/'input_provenance.json'
    if destination.exists():
        verify(root,run)
        return
    if phase=='before_search' and ((run/'manifest.json').exists() or (run/'evaluations.jsonl').exists()):
        raise ValueError('cannot label existing results as captured before search')
    record=dict(schema='complete-direct-mapping-inputs.v1',capture_phase=phase,
                captured_at_utc=datetime.now(timezone.utc).isoformat(),
                files={name:digest(root/name) for name in REQUIRED})
    run.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(record,indent=2),encoding='utf-8')

def verify(root, run):
    root,run=Path(root),Path(run)
    path=run/'input_provenance.json'
    if not path.exists():raise ValueError('missing complete provenance; launch searches via verification/run_continuous_search.py')
    record=json.loads(path.read_text(encoding='utf-8'))
    if record.get('schema')!='complete-direct-mapping-inputs.v1' or set(record.get('files',{}))!=set(REQUIRED):
        raise ValueError('incomplete evaluator/input provenance')
    for name,expected in record['files'].items():
        if not (root/name).is_file() or digest(root/name)!=expected:
            raise ValueError(f'evaluator or actual input changed: {name}')
    return record
