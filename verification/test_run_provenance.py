import importlib.util
import json
from pathlib import Path
import pytest

def test_provenance_detects_changed_wall_and_preserved_kernel(tmp_path):
    spec=importlib.util.find_spec('d3d_optimization.run_provenance')
    assert spec is not None, 'complete evaluator/input provenance missing'
    from d3d_optimization.run_provenance import capture, verify, REQUIRED
    root=tmp_path/'project';run=tmp_path/'run'
    for name in REQUIRED:
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('original')
    capture(root,run,'before_search')
    verify(root,run)
    for name in ['baseline_review/baseline_review.json','d3d_optimization/frozen_mapping.py']:
        p=root/name;p.write_text('changed')
        with pytest.raises(ValueError,match='changed'):verify(root,run)
        p.write_text('original')
    record=json.loads((run/'input_provenance.json').read_text())
    del record['files']['d3d_optimization/shape.py']
    (run/'input_provenance.json').write_text(json.dumps(record))
    with pytest.raises(ValueError,match='incomplete'):verify(root,run)

def test_validation_requires_preexisting_provenance(tmp_path):
    spec=importlib.util.find_spec('d3d_optimization.run_provenance')
    assert spec is not None, 'complete evaluator/input provenance missing'
    from d3d_optimization.run_provenance import verify
    with pytest.raises(ValueError,match='missing'):verify(tmp_path,tmp_path)
