"""Stage immutable held-out DIII-D inputs using the existing material-wall checks."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil

import numpy as np
import pyarrow.parquet as pq

from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.case import load_case
from verification import nf_equilibrium_times as geometry
from verification.nf_search_benchmark import ROOT

SHOT = geometry.SHOT
REVIEW = geometry.REVIEW
GEOMETRY = geometry.GEOMETRY
MANIFEST = geometry.MANIFEST


def nearest_time(times, target=1400.):
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or not len(times) or not np.isfinite(times).all():
        raise ValueError("need finite EFIT timestamps")
    index = int(np.lexsort((times, abs(times-target)))[0])
    time = float(times[index])
    if abs(time-target) >= 100:
        raise ValueError("nearest stored EFIT time is at least 100 ms from 1400 ms")
    if np.count_nonzero(times == time) != 1:
        raise ValueError("selected EFIT time must appear once")
    return index, time


def original_fingerprints():
    record = pq.read_table(ROOT/SHOT, columns=["efit_psirz"]).to_pylist()[0]
    return [np.asarray(x, dtype=float) for x in record["efit_psirz"]]


def prepare(study):
    selections = load_json(study/"selection.json")
    downloads = load_json(study/"downloads.json")
    if selections["revision"] != downloads["revision"]:
        raise ValueError("download revision mismatch")
    if selections["protocol_sha256"] != digest(study/"protocol.ko.md"):
        raise ValueError("protocol changed after selection")
    if [row["path"] for row in downloads["files"]] != selections["selected"]:
        raise ValueError("source order changed")
    original = load_json(ROOT/REVIEW)
    old_flux = original_fingerprints()
    all_flux = list(old_flux)
    outputs = []
    selected_count = 0
    for item in downloads["files"]:
        case_id = Path(item["path"]).stem.removeprefix("d3d_shot_")
        staged = study/"cases"/case_id/"input-root"
        record = dict(case_id=case_id, status="excluded", source=item, root=str(staged.resolve()))
        try:
            if "error" in item:
                raise ValueError("source download failed: "+item["error"])
            raw = Path(item["local"])
            if digest(raw) != item["sha256"]:
                raise ValueError("download hash mismatch")
            data = pq.read_table(raw, columns=["source", "efit_times", "efit_psirz"]).to_pylist()
            if len(data) != 1 or data[0]["source"] != "DIII-D":
                raise ValueError("expected exactly one DIII-D discharge")
            row = data[0]
            i, time = nearest_time(row["efit_times"])
            flux = np.asarray(row["efit_psirz"][i], dtype=float)
            record.update(time_ms=time, frame_index=i)
            if any(np.array_equal(flux, old) for old in all_flux):
                raise ValueError("identical equilibrium frame to a previous input")
            all_flux.append(flux)
            staged.mkdir(parents=True, exist_ok=False)
            for source, dest in ((raw, SHOT), (ROOT/GEOMETRY, GEOMETRY)):
                (staged/dest).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, staged/dest)
            save(staged/REVIEW, original, exclusive=True)
            manifest = load_json(ROOT/MANIFEST)
            manifest.update(time_ms=time, shot_sha256=digest(staged/SHOT))
            save(staged/MANIFEST, manifest, exclusive=True)
            # This old helper already checks exact OMP endpoints and both leg traces.
            # Its ROOT is scoped to this serial call and restored even on failure.
            previous_root = geometry.ROOT
            try:
                geometry.ROOT = staged
                details, review, reason = geometry.geometry_at(time)
            finally:
                geometry.ROOT = previous_root
            record.update(details)
            if reason:
                raise ValueError(reason)
            save(staged/REVIEW, review)
            np.testing.assert_array_equal(review["wall"], original["wall"])
            field, lcfs, axis = load_case(staged/SHOT, time)
            axis_flux = float(field.values([axis])[0])
            xp_flux = float(field.values([details["X_RZ_m"]])[0])
            flux_scale = xp_flux-axis_flux
            values = (field.values(lcfs)-axis_flux)/flux_scale
            record.update(tf_current_kA=field.tf_current_ka, grid_shape=list(field.psi.shape),
                          lower_X_flux=xp_flux, lcfs_normalized_flux_quantiles=np.quantile(values,[0,.05,.5,.95,1]).tolist(),
                          geometry_passed=True, inputs={str(p):digest(staged/p) for p in (SHOT,REVIEW,GEOMETRY,MANIFEST)})
            if selected_count < 3:
                record["status"] = "passed"
                selected_count += 1
            else:
                record["status"] = "reserve_geometry_passed"
        except (ValueError, RuntimeError, OSError, AssertionError) as exc:
            record["reason"] = f"{type(exc).__name__}: {exc}"
        save(study/"cases"/case_id/"preflight.json", record, exclusive=True)
        outputs.append(record)
        print(case_id,record["status"],record.get("time_ms"),record.get("reason"),flush=True)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  protocol_sha256=digest(study/"protocol.ko.md"),
                  selection_sha256=digest(study/"selection.json"),
                  downloads_sha256=digest(study/"downloads.json"),
                  original_sha256={str(p):digest(ROOT/p) for p in (SHOT,REVIEW,GEOMETRY,MANIFEST)},
                  driver_sha256=digest(__file__), geometry_driver_sha256=digest(geometry.__file__),
                  cases=outputs)
    save(study/"preflight.json",report,exclusive=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,required=True)
    args=parser.parse_args()
    prepare(args.study.resolve())

