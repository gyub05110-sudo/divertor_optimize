"""Synthetic publication checks; no production heat data required."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest


CASES = ("3264fdd71f","7ce879f1b3","70481becc0","5d7f95a730","1255ab8051","d880609264","f317e77cf7")


def fixture():
    def peak(q,region):
        return {"peak_MW_m2":q,"peak_surface":region,"peak_ci95_MW_m2":[.92*q,1.1*q]}
    def pair(name,ref,val):
        return dict(target=name,references=ref,stage4_pct=val+.1,model_pct=val+.2,firefly=dict(
            estimate_pct=val,ci95_pct=[max(0.,val-2.) if len(ref)>1 else val-2.,val+3.],best_references=ref[:1],
            uncertainty_pp=3.,ci_only_1pct="unresolved",precision_1pct=False,verdict_1pct="unresolved",direction="unresolved"))
    cases=[]
    for i,cid in enumerate(CASES):
        family="floor-shelf" if i in (0,3,5) else "floor-only"
        forced=["floor-only","floor-entry","floor-shelf","floor-connector"]
        rows=[]; roles=[dict(role="baseline",seed=30713,wall_id="base",
                             model=peak(4+i/10,"floor"),stage4=peak(4.1+i/10,"shelf"),firefly=peak(4.2+i/10,"shelf"))]
        for s in ([30713,30757] if cid==CASES[2] else [30713]):
            p=dict(optimizer_seed=s,hotspot_vs_full=pair("hotspot",["full"],i-.4),
                   hotspot_vs_baseline=pair("hotspot",["baseline"],-8.-i),
                   full_vs_baseline=pair("full",["baseline"],-10.-i))
            names=["hotspot","full"]
            if s==30713:
                p.update(predicted_forced_vs_best_forced=pair(family,forced,.1*i),
                         floor_only_vs_best_forced=pair("floor-only",forced,i+.1),
                         predicted_forced_vs_floor_only=pair(family,["floor-only"],-2.5 if family=="floor-shelf" else 0.))
                names+=forced
            for j,r in enumerate(names):
                roles.append(dict(role=r,seed=s,wall_id=r+str(s),stage4=peak(3+j/3,"floor"),
                                  model=peak(3+j/3,"floor"),firefly=peak(3.1+j/3,"shelf")))
            rows.append(p)
        cases.append(dict(case_id=cid,predicted_family=family,roles=roles,comparisons=rows,
                          baseline_hotspot=dict(stage4="shelf",model="floor",firefly="shelf",predicted_region="shelf")))
    return dict(schema="nf-stage5-whole-wall-metrics.v1",eligible_files=7,role_count=51,distinct_case_walls=46,
                ordered_bootstrap_draws=46656,unique_bootstrap_counts=462,original_registration_files=24,cases=cases)


def load_summary():
    spec=importlib.util.spec_from_file_location("summarize_stage5",Path(__file__).with_name("summarize.py"))
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_two_full_precision_tables_and_report_figure(tmp_path):
    m=load_summary()
    metrics=fixture()
    (tmp_path/"metrics.json").write_text(json.dumps(metrics))
    forecast,products=m.publish(tmp_path,tmp_path)
    assert len(forecast)==7 and len(products)==8
    assert forecast[0]["predicted_vs_floor_pct"] == -2.5
    assert products[2]["optimizer_seed"]==30713 and products[3]["optimizer_seed"]==30757
    assert forecast[0]["model_baseline_region"]=="floor" and forecast[0]["firefly_baseline_region"]=="shelf"
    text=m.markdown(forecast,products)
    assert "95%" in text and "70481becc0" in text and "30757" in text
    for name,n in (("firefly-forecast.csv",7),("firefly-products.csv",8)):
        with (tmp_path/name).open(encoding="utf-8-sig",newline="") as f:
            assert len(list(csv.DictReader(f)))==n
    assert (tmp_path/"firefly-comparisons.pdf").read_bytes().startswith(b"%PDF")
    assert (tmp_path/"firefly-comparisons.png").read_bytes()[:8]==b"\x89PNG\r\n\x1a\n"


def test_incomplete_metrics_cannot_generate_summary(tmp_path):
    m=load_summary()
    metrics=fixture(); metrics["cases"][2]["comparisons"].pop()
    (tmp_path/"metrics.json").write_text(json.dumps(metrics))
    with pytest.raises(ValueError,match="seed"):
        m.publish(tmp_path,tmp_path)
    assert not list(tmp_path.glob("*.csv"))
