"""Register fixed Stage 4 walls and validate the native-input panel before FIREFLY heating."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

from verification.stage3_campaign import public_save, read, require, sha
from verification.stage4_campaign import FOCUSED

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT/"results/nf-stage5-firefly-20261110"
STAGE4 = Path("D:/DivertorOptimization/studies/nf-stage4-20261109")
STUDY = Path("D:/DivertorOptimization/studies/nf-stage5-firefly-20261110")
CASES = ("3264fdd71f","7ce879f1b3","70481becc0","5d7f95a730","1255ab8051","d880609264","f317e77cf7")
EXCEPTION = "70481becc0"
SEEDS = (51031,51047,51059,51071,51089,51109)
ARMS = (*FOCUSED,"hotspot","full")
REPORT_SHA = "ea51d468e649f5e1cdef8015167c47c923564b0d958da37b1aca5df3e3520b29"


def jobs(roles):
    expected = {(case,role,30713) for case in CASES for role in ARMS}
    expected.update((case,"baseline",None) for case in CASES)
    expected.update((EXCEPTION,role,30757) for role in ("hotspot","full"))
    keys = [(r["case_id"],r["role"],r["optimizer_seed"]) for r in roles]
    require(len(keys)==len(expected) and set(keys)==expected, "missing or duplicate registered case/role/seed")
    tasks = []
    for case in CASES:
        walls = dict.fromkeys(r["wall_id"] for r in roles if r["case_id"]==case)
        for seed in SEEDS:
            for wall in walls:
                require(len(wall)==64 and all(ch in "0123456789abcdef" for ch in wall),"invalid wall SHA-256")
                tasks.append((case,wall,seed))
    return tasks


def verify_panel(roles,panel):
    expected = {(r["case_id"],r["role"],r["optimizer_seed"]):r["wall_id"] for r in roles}
    actual = {}
    for case in panel["cases"]:
        for value in case["roles"]:
            role=value["role"]
            key=(case["case_id"],role,None if role=="baseline" else value["seed"])
            require(key not in actual,"duplicate input-panel role")
            actual[key]=value["wall_id"]
    require(actual==expected,"native input-panel wall differs from independently registered geometry")


def archive_role(root,row,role):
    folder = root/row["relative_output"]
    for name,key in (("result.json","result_sha256"),("events.jsonl","events_sha256"),("receipt.json","receipt_sha256")):
        require(sha(folder/name)==row[key],"archived Stage 4 run differs: "+str(folder/name))
    result = read(folder/"result.json")
    require(row["status"]==result["status"]=="qualified" and result["arguments"]["seed"]==row["seed"]
            and result["arguments"]["family"]==row["arm"],"wrong original winner")
    source = "baseline" if role=="baseline" else "selected"
    measurement = result[source]
    for name,digest in measurement["sha256"].items():
        require(sha(folder/source/name)==digest,"archived wall or measurement changed")
    with np.load(folder/source/"boundary.npz",allow_pickle=False) as raw:
        wall=np.ascontiguousarray(raw["wall"],dtype="<f8")
        require(wall.ndim==2 and wall.shape[1]==2 and np.array_equal(wall[0],wall[-1]) and
                np.isfinite(wall).all() and np.all(wall[:,0]>0),"invalid physical wall vertices")
        require(float(raw["q_MW_m2"].max())==measurement["peak_MW_m2"],"native peak differs from samples")
    return dict(case_id=row["case_id"],role=role,optimizer_seed=None if role=="baseline" else row["seed"],
                wall_id=hashlib.sha256(wall.tobytes()).hexdigest(),vertex_count=len(wall),
                stage4_directory=row["relative_output"]+"/"+source,result_sha256=row["result_sha256"],
                events_sha256=row["events_sha256"],receipt_sha256=row["receipt_sha256"],
                equilibrium_sha256=result["inputs"]["equilibrium_sha256"],peak_MW_m2=measurement["peak_MW_m2"],
                measurement_sha256=measurement["sha256"],candidate_id=None if role=="baseline" else row["candidate"]["id"])


def make_registration(stage4=STAGE4):
    path=stage4/"stage4_report.json"
    require(sha(path)==REPORT_SHA,"wrong Stage 4 results")
    report=read(path)
    predictions=[p for p in report["predictions"] if p["scorable"]]
    require(tuple(p["case_id"] for p in predictions)==CASES and report["complete"],"changed eligibility")
    lookup={(r["case_id"],r["seed"],r["arm"]):r for r in report["runs"]}
    roles=[]
    for case in CASES:
        roles.append(archive_role(stage4,lookup[case,30713,"hotspot"],"baseline"))
        for arm in ARMS:
            roles.append(archive_role(stage4,lookup[case,30713,arm],arm))
        if case==EXCEPTION:
            for arm in ("hotspot","full"):
                roles.append(archive_role(stage4,lookup[case,30757,arm],arm))
    tasks=jobs(roles)
    return dict(schema="nf-stage5-fixed-wall-registration.v1",stage4_report_sha256=REPORT_SHA,
                stage4_source_commit="427de9e1180dbf45de4ba50904a8272c135eeb9b",
                protocol_sha256=sha(PAPER/"protocol.ko.md"),predictions={p["case_id"]:p["predicted_family"] for p in predictions},
                seeds=SEEDS,histories_per_seed=4096,mpi_ranks=4,primary_ring_m=.002,roles=roles,jobs=tasks)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("register","verify-panel"))
    parser.add_argument("--study",type=Path,default=STUDY)
    args=parser.parse_args()
    registration=args.study/"wall-registration.json"
    record=make_registration()
    if registration.exists():
        require(read(registration)==json.loads(json.dumps(record)),"registration already fixed with other inputs")
    else:
        require(args.action=="register" and not (args.study/"native/results").exists(),"missing pre-heat registration")
        public_save(registration,record)
    if args.action=="verify-panel":
        verify_panel(record["roles"],read(args.study/"inputs/panel.json"))
    print(len(record["roles"]),"roles,",len(record["jobs"])//len(SEEDS),"distinct case/walls,",
          len(record["jobs"]),"fixed case/wall/seed jobs. SHA-256",sha(registration))


if __name__=="__main__":
    main()
