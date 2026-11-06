import importlib.util
import json
from pathlib import Path
import sys
import zipfile

import pytest


def builder():
    path = Path(__file__).resolve().parents[1] / 'results/nf-anytime-fullspace-20261106/build_evidence_archive.py'
    spec = importlib.util.spec_from_file_location('anytime_evidence_builder', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding='utf-8')


def invoke(module, monkeypatch, study, out):
    if study.drive.lower() != 'd:':
        pytest.skip('The study and its test scratch directory must be on D')
    monkeypatch.setattr(sys, 'argv', ['bundle', '--study', str(study), '--out', str(out)])
    module.main()


def test_missing_rejection_phase_refuses_archive(tmp_path, monkeypatch):
    study, out = tmp_path/'study', tmp_path/'out'
    for phase in ('reference', 'uncertainty'):
        put(study/phase/'summary.json', {})
    with pytest.raises(ValueError, match='screen-rejections'):
        invoke(builder(), monkeypatch, study, out)
    assert not (out/'evidence.zip').exists()


def test_bundle_retains_raw_bytes_and_refuses_overwrite(tmp_path, monkeypatch):
    study, out = tmp_path/'study', tmp_path/'out'
    prior, source = tmp_path/'prior', tmp_path/'source'
    for phase in ('reference', 'uncertainty', 'screen-rejections'):
        put(study/phase/'summary.json', {})
    put(study/'frozen.json', {'source_study': str(source)})
    put(study/'amendment.json', {'prior_study': str(prior)})
    put(prior/'interruption.json', {'retained': True})
    for name in ('selection.json', 'downloads.json', 'preflight.json', 'protocol.ko.md'):
        put(source/name, {})
    raw = study/'live/raw.jsonl'
    raw.parent.mkdir(parents=True)
    raw.write_bytes(b'{"x":1}\r\n{"x":2}\r\n')
    put(study/'posttrial-runtime/cache.json', {'omit': True})
    module = builder()
    invoke(module, monkeypatch, study, out)
    with zipfile.ZipFile(out/'evidence.zip') as bundle:
        name = 'studies/study/live/raw.jsonl'
        assert bundle.read(name) == raw.read_bytes()
        manifest = json.loads(bundle.read('evidence-manifest.json'))
        assert manifest['entries'][name]['sha256'] == module.sha(raw.read_bytes())
        assert not any('posttrial-runtime/' in member for member in bundle.namelist())
    before = (out/'evidence.zip').read_bytes()
    with pytest.raises(FileExistsError):
        invoke(module, monkeypatch, study, out)
    assert (out/'evidence.zip').read_bytes() == before
