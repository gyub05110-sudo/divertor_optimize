"""Cross-check Stage 5 paired decisions using the published per-seed ring watts."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path

import numpy as np


SEEDS=(51031,51047,51059,51071,51089,51109)
FORCED=("floor-only","floor-entry","floor-shelf","floor-connector")


def ordered_counts():
    # Enumerate the ordered Cartesian product independently of the production
    # multinomial coefficient implementation; coalesce identical draws.
    histogram=Counter(tuple(draw.count(i) for i in range(6)) for draw in product(range(6),repeat=6))
    return np.array(list(histogram),dtype=float),np.array(list(histogram.values()),dtype=int)


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def peaks(path,expected_sha,counts):
    if digest(path)!=expected_sha:
        raise ValueError("published profile checksum differs: "+str(path))
    with gzip.open(path,"rt",newline="",encoding="utf8") as f:
        cells=csv.DictReader(f)
        powers=np.array([[float(row["area_m2"]),*[float(row[f"seed_{s}_W"]) for s in SEEDS]] for row in cells])
    q=powers[:,1:]/powers[:,:1]/1e6
    if not np.all(np.isfinite(q)) or np.any(q<0) or np.any(powers[:,0]<=0):
        raise ValueError("invalid published cell area or watts")
    return float(q.mean(axis=1).max()),np.max(q@counts.T,axis=0)/6


def compare(candidate,references,frequency):
    point,draws=candidate
    best=min(r[0] for r in references)
    sampled=np.stack([r[1] for r in references]).min(axis=0)
    gap=100*(point/best-1)
    low,high=np.quantile(np.repeat(100*(draws/sampled-1),frequency),[.025,.975])
    width=max(gap-low,high-gap)
    interval="supported" if high<=1+1e-12 else "contradicted" if low>1+1e-12 else "unresolved"
    return dict(estimate_pct=gap,ci95_pct=[float(low),float(high)],uncertainty_pp=float(width),
                ci_only_1pct=interval,precision_1pct=bool(width<=1+1e-12),
                verdict_1pct=interval if width<=1+1e-12 else "unresolved",
                direction="lower" if high<0 else "higher" if low>0 else "unresolved")


def run(study):
    metrics_path=study/"metrics.json"
    metrics=json.loads(metrics_path.read_bytes())
    reg=json.loads((study/"wall-registration.json").read_bytes())
    counts,frequency=ordered_counts()
    cases=[]
    tally={}
    pairs=0
    for c in metrics["cases"]:
        cid=c["case_id"]
        saved={p["optimizer_seed"]:p for p in c["comparisons"]}
        data={w:peaks(study/g["profile"],g["profile_sha256"],counts) for w,g in c["geometries"].items()}
        roles={(r["optimizer_seed"],r["role"]):data[r["wall_id"]] for r in reg["roles"] if r["case_id"]==cid}
        for seed in ((30713,30757) if cid=="70481becc0" else (30713,)):
            pairs+=1
            hotspot=roles[seed,"hotspot"]
            full=roles[seed,"full"]
            baseline=roles[None,"baseline"]
            tests=dict(hotspot_vs_full=(hotspot,[full]),hotspot_vs_baseline=(hotspot,[baseline]),full_vs_baseline=(full,[baseline]))
            if seed==30713:
                selected=roles[seed,reg["predictions"][cid]]
                forced=[roles[seed,f] for f in FORCED]
                tests.update(predicted_forced_vs_best_forced=(selected,forced),
                             floor_only_vs_best_forced=(roles[seed,"floor-only"],forced),
                             predicted_forced_vs_floor_only=(selected,[roles[seed,"floor-only"]]))
            for name,(candidate,refs) in tests.items():
                row=compare(candidate,refs,frequency)
                recorded=saved[seed][name]["firefly"]
                for key in ("estimate_pct","ci95_pct","uncertainty_pp"):
                    if not np.allclose(recorded[key],row[key],rtol=1e-11,atol=1e-11):
                        raise ValueError(f"different {key}: {cid}/{seed}/{name}")
                for key in ("ci_only_1pct","precision_1pct","verdict_1pct","direction"):
                    if recorded[key]!=row[key]:
                        raise ValueError(f"different {key}: {cid}/{seed}/{name}")
                hit=tally.setdefault(name,dict(verdict=Counter(),direction=Counter(),ci_only=Counter()))
                hit["verdict"].update([row["verdict_1pct"]])
                hit["direction"].update([row["direction"]])
                hit["ci_only"].update([row["ci_only_1pct"]])
                cases.append(dict(case_id=cid,seed=seed,contrast=name,**row))
    if len(metrics["cases"])!=7 or pairs!=8 or sum(frequency)!=6**6 or len(cases)!=45:
        raise ValueError("incomplete 7-file/8-product independent audit")
    return dict(metrics_sha256=digest(metrics_path),ordered_draws=int(sum(frequency)),unique_counts=len(counts),
                forecast_files=7,product_pairs=8,by_contrast=tally,rows=cases)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--study",required=True,type=Path)
    cli.add_argument("--out",type=Path)
    args=cli.parse_args()
    result=run(args.study)
    if args.out:
        data=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
        if args.out.exists() and args.out.read_bytes()!=data:
            raise ValueError("existing cross-check differs")
        if not args.out.exists():
            args.out.write_bytes(data)
    print(json.dumps(result["by_contrast"],indent=2))
