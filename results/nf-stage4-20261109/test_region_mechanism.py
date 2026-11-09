"""Post-study regional attribution checks on synthetic native exports."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


def load():
    p=Path(__file__).with_name("region_mechanism.py")
    assert p.is_file(), "standalone regional analysis helper is required"
    spec=importlib.util.spec_from_file_location("regional_mechanism",p)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def put(p, value):
    p.parent.mkdir(parents=True,exist_ok=True)
    b=value if isinstance(value,bytes) else json.dumps(value).encode()
    p.write_bytes(b)
    return hashlib.sha256(b).hexdigest()


def make_case(study, predicted="floor-only", missing=None):
    names=("hotspot","full","floor-only","floor-entry","floor-shelf","floor-connector")
    rows={}
    cid="0123456789"
    for i, arm in enumerate(names):
        rel=f"runs/{cid}/30713/{arm}"
        root=study/rel
        peak=(8.,7.,9.,7.,7.,8.)[i]
        def measure(kind, q, a2, power):
            if kind=="baseline":
                surfaces=dict(floor=dict(peak_MW_m2=10.,deposited_W=70.,first_hit_W=50.),
                              shelf=dict(peak_MW_m2=10.,deposited_W=30.,first_hit_W=25.))
            else:
                surfaces=dict(floor=dict(peak_MW_m2=q-1,deposited_W=power,first_hit_W=power-10),
                              shelf=dict(peak_MW_m2=q,deposited_W=100.-power,first_hit_W=80.-power))
            files={name:put(root/kind/name,(kind+arm).encode())
                   for name in ("geometry.csv","profile.csv.gz","donors.csv.gz","boundary.npz")}
            return dict(physics_measurable=True,candidate_eligible=kind=="selected",within_budget=True,
                        peak_MW_m2=q,peak_region="floor" if kind=="baseline" else "shelf",
                        A2_m2=a2,surfaces=surfaces,directory=kind,sha256=files)
        cand=None if arm==missing else dict(id=arm)
        data=dict(inputs=dict(equilibrium_sha256="shot",time_ms=1400),arguments=dict(seed=30713,family=arm),
                  candidate=cand,baseline=measure("baseline",10.,2.,70.),
                  selected=None if cand is None else measure("selected",peak,3.+i,40.-i))
        status="qualified" if cand else "no_design"
        sha=put(root/"result.json",dict(data,status=status))
        rows[arm]=dict(case_id=cid,seed=30713,arm=arm,status=status,receipt_status="completed",result_sha256=sha,
                       relative_output=rel,candidate=cand,checkpoints={"600":dict(id=arm if cand else None,
                       peak_MW_m2=peak if cand else None)})
    pred=dict(case_id=cid,seed=30713,predicted_family=predicted,abstained=predicted=="full",
              scorable=predicted!="full" and missing is None,all_controls=missing is None,
              best_families=["floor-entry","floor-shelf"] if missing is None else [],
              competitive_families=["floor-entry","floor-shelf"] if missing is None else [],
              selected_competitive=False if predicted!="full" and missing is None else None)
    case=dict(id=cid,rank=1,time_ms=1400.,sha256="shot",dataset_path="d3d_shot_0123456789.parquet",status="ready")
    return case,pred,rows


def test_tracks_preexisting_hotspots_even_when_winner_peak_moves_and_includes_all_strict_ties(tmp_path):
    mod=load()
    case,pred,rows=make_case(tmp_path)
    result=mod.case_mechanism(tmp_path,case,pred,rows)
    assert result["classification"]=="miss"
    assert result["baseline"]["tied_peak_regions"]==["floor","shelf"]
    assert result["winners"]["hotspot"]["peak_region"]=="shelf"
    assert result["winners"]["hotspot"]["peak_region_power"]["deposited_W"]==60
    assert result["winners"]["hotspot"]["baseline_hotspot_regions"]["floor"]["first_hit_W"]==30
    assert result["winners"]["hotspot"]["change_from_baseline"]["baseline_hotspot_regions"]["floor"]["first_hit_W"]==-20
    assert result["product_minus_full"]["global_peak_MW_m2"]==1
    assert set(result["predicted_minus_strict_best"])=={"floor-entry","floor-shelf"}
    assert result["predicted_minus_strict_best"]["floor-entry"]["global_peak_MW_m2"]==2
    assert result["winners"]["floor-entry"]["result_sha256"]==rows["floor-entry"]["result_sha256"]
    (tmp_path/rows["floor-entry"]["relative_output"]/"selected/donors.csv.gz").write_bytes(b"modified")
    with pytest.raises(ValueError,match="digest"):
        mod.case_mechanism(tmp_path,case,pred,rows)


@pytest.mark.parametrize("predicted,missing,label",[("full",None,"abstained"),("floor-only","floor-shelf","incomplete_controls")])
def test_retains_abstentions_and_incomplete_panels_without_an_inferred_reference(tmp_path,predicted,missing,label):
    case,pred,rows=make_case(tmp_path,predicted,missing)
    result=load().case_mechanism(tmp_path,case,pred,rows)
    assert result["classification"]==label and result["predicted_minus_strict_best"]=={}
    assert result["product_minus_full"]["global_peak_MW_m2"]==1


def test_helper_reads_no_run_before_completed_final_report(tmp_path):
    put(tmp_path/"stage4_report.json",dict(complete=False))
    with pytest.raises(ValueError,match="complete"):
        load().analyze(tmp_path)


def test_complete_report_keeps_all_files_and_hashes_the_native_inputs(tmp_path):
    mod=load()
    case,pred,first=make_case(tmp_path)
    cases=[case]
    preds=[pred]
    rows=list(first.values())
    for n in range(2,25):
        cid=f"{n:010d}"
        cases.append(dict(case,id=cid,rank=n,status="missing_time"))
        preds.append(dict(pred,case_id=cid,scorable=False,predicted_family=None,abstained=None,
                          selected_competitive=None,all_controls=False,best_families=[],competitive_families=[]))
        for arm in mod.ARMS:
            rows.append(dict(case_id=cid,seed=30713,arm=arm,status="registration_excluded",receipt_status="registration_excluded",
                             relative_output=f"runs/{cid}/30713/{arm}",result_sha256=None))
    for item in cases:
        rows.extend(dict(case_id=item["id"],seed=s,arm=a,status="no_design") for s in (30731,30757) for a in mod.PRODUCT)
    reg=put(tmp_path/"inputs/registration.json",dict(cases=cases))
    plan=put(tmp_path/"expected.json",{"synthetic":True})
    digest=put(tmp_path/"stage4_report.json",dict(schema="nf-stage4-forecast-report.v1",complete=True,
               registration_sha256=reg,plan_sha256=plan,predictions=preds,runs=rows))
    result=mod.analyze(tmp_path)
    assert len(result["cases"])==24 and result["cases"][0]["seed"]==30713
    assert result["cases"][-1]["baseline"] is None
    assert result["report_sha256"]==digest and result["classification_counts"]==dict(miss=1,incomplete_controls=23)
