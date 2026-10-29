"""Verify publication copies, restored dependencies and saved study fingerprints."""
import hashlib,json,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
def load(path):return json.loads(path.read_text(encoding='utf-8'))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

copy=load(ROOT/'docs/publication/source_copy_hashes.json')
for name,expected in copy['files'].items():
    path=ROOT/('docs/README_LEGACY.md' if name=='README.md' else name)
    assert path.is_file() and digest(path)==expected,f'published source differs: {name}'

preserved=load(ROOT/'results/regional-20261028/preservation.json')
for name,expected in preserved.items():
    path=ROOT/name
    assert path.is_file() and digest(path)==expected,f'preserved dependency/source differs: {name}'

from optimize_boundary_transport import verify
for source in ['IMP','OMP']:
    verify(ROOT/'results/regional-20261028'/source)
    folder=ROOT/'results/regional-20261028/qualified'/source
    verify(folder)
    assert load(folder/'validation.json')['verdict']['improvement_confirmed']
    assert load(folder/'regional_validation.json')['passed']
print(json.dumps({'copied_project_files_verified':len(copy['files']),
    'preserved_python_files_verified':len(preserved),'study_provenance_verified':4,
    'saved_qualified_verdicts_passed':2,'scope':'fingerprints and saved verdicts; not a fresh numerical refinement'}))
