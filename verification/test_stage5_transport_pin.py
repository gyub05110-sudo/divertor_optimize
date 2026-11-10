import hashlib
import json

import pytest

from verification import stage5_native as native


def test_corrected_transport_reads_sealed_library_and_rejects_changed_bytes(tmp_path):
    rel = "native-physics/cap-sync-v1/libFIREFLY.so"
    p = tmp_path / rel
    p.parent.mkdir(parents=True)
    p.write_bytes(b"transport")
    expected = hashlib.sha256(p.read_bytes()).hexdigest()
    (tmp_path/"native-physics/active-library.txt").write_text(rel+"\n")
    env = {"active_transport_library": "/stage5/"+rel, "native_files_sha256": {"/stage5/"+rel:expected}}
    (tmp_path/"execution-environment.json").write_text(json.dumps(env))
    assert native.active_transport(tmp_path)==dict(active_transport_library="/stage5/"+rel,active_transport_sha256=expected)
    p.write_bytes(b"altered")
    with pytest.raises(ValueError,match="thermal library"):
        native.active_transport(tmp_path)


def test_library_must_be_separately_versioned(tmp_path):
    (tmp_path/"native-physics").mkdir()
    (tmp_path/"native-physics/active-library.txt").write_text("../old/libFIREFLY.so")
    with pytest.raises(ValueError,match="thermal library path"):
        native.active_transport(tmp_path)
