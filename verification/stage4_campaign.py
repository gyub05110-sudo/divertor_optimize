"""Sealed serial experiment for the supplied-equilibrium hotspot policy and forced families."""
import argparse
from collections import Counter
from pathlib import Path
import re
import subprocess
import sys

from verification import stage3_campaign as runner


ROOT = runner.ROOT
PAPER = ROOT / "results/nf-stage4-20261109"
SEEDS = (30713, 30731, 30757)
FOCUSED = ("floor-only", "floor-entry", "floor-shelf", "floor-connector")
PRODUCT = ("hotspot", "full")
ALL = PRODUCT + FOCUSED
PROTOCOL_SHA = "506285e6e9d7c36aead7f5a38a096bea23aeda836d43444de76f0e975fdcf47d"


def build_tasks(cases):
    ids = [case["id"] for case in cases]
    runner.require(len(cases)==24 and len(set(ids))==24 and all(re.fullmatch(r"[\w-]+",i,re.ASCII) for i in ids),
                   "expected 24 unique safe case identifiers")
    runner.require([c["rank"] for c in cases]==list(range(1,25)) and all(c["time_ms"]==1400 for c in cases),
                   "case registration rank or time changed")
    jobs = []
    for i, case in enumerate(cases):
        runner.require(bool(case.get("path")) == (case["status"] != "acquisition_failed"),
                       "registered acquisition and raw path disagree")
        rotation = i % len(ALL)
        families = (ALL[rotation:] + ALL[:rotation], PRODUCT if i%2 == 0 else PRODUCT[::-1],
                    PRODUCT[::-1] if i%2 == 0 else PRODUCT)
        for seed, order in zip(SEEDS, families):
            for family in order:
                jobs.append(dict(case_id=case["id"],input_rank=case["rank"],registration_status=case["status"],
                                 registration_reason=case.get("reason"),equilibrium=case.get("path"),sha256=case["sha256"],
                                 time_ms=case["time_ms"],seed=seed,arm=family,budget_seconds=600,family=family,
                                 method="screened",relative_output=f"runs/{case['id']}/{seed}/{family}"))
    return jobs


def source_files():
    files = set(runner.study_files())
    files.update((ROOT/"verification").glob("stage4_*.py"))
    files.update((PAPER/"protocol.ko.md", PAPER/"execution-environment.json", PAPER/"inputs/registration.json"))
    return {p.relative_to(ROOT).as_posix():runner.sha(p) for p in sorted(files)}


def invariants(study):
    registry = study/"inputs/registration.json"
    reg = runner.read(registry)
    runner.require(runner.sha(registry)==runner.sha(PAPER/"inputs/registration.json") and
                   reg["protocol_sha256"]==PROTOCOL_SHA==runner.sha(PAPER/"protocol.ko.md"),
                   "registration or prospective protocol changed")
    return dict(schema="nf-stage4-hotspot.v1",files=source_files(),registration_sha256=runner.sha(registry),
                settings=runner.SETTINGS,seeds=list(SEEDS),arms=list(ALL),tasks=build_tasks(reg["cases"]),
                python=dict(executable=str(Path(sys.executable).resolve()),sha256=runner.sha(sys.executable),
                            version=sys.version),worker_environment=runner.worker_env_settings(study))


def make_plan(study):
    return dict(**invariants(study),registered_utc=runner.utc(),
                revision=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip())


def check_plan(study, plan):
    for key, value in invariants(study).items():
        runner.require(plan[key]==value, "stage 4 seal mismatch: "+key)


def load_seal(study):
    return runner.load_seal(study,plan_checker=check_plan)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("seal","run","worker","check"))
    parser.add_argument("--study",required=True,type=Path)
    parser.add_argument("--index",type=int)
    args=parser.parse_args()
    study=args.study.resolve()
    if args.action=="seal":
        plan=runner.seal(study,plan_factory=make_plan,plan_checker=check_plan)
        print("Sealed",len(plan["tasks"]),"slots",runner.sha(study/"expected.json"))
    elif args.action=="worker":
        if args.index is None or not 0 <= args.index < 240:
            parser.error("worker requires an index in the sealed 240-slot schedule")
        runner.worker(study,args.index,plan_loader=load_seal)
    elif args.action=="run":
        runner.run(study,plan_loader=load_seal,worker_module="verification.stage4_campaign")
    else:
        plan=load_seal(study)
        counts=Counter()
        for job in plan["tasks"]:
            receipt=runner.verified_receipt(study/job["relative_output"],job,plan)
            counts[receipt["status"]]+=1
        print("Verified",sum(counts.values()),"terminal receipts",dict(counts))


if __name__=="__main__":
    main()
