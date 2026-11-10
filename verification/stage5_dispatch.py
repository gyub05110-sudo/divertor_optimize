"""Run the registered fixed-wall FIREFLY jobs in isolated native processes and verify resume evidence."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter

from d3d_optimization.campaign_state import keep_awake, process_lock
from verification.stage3_campaign import public_save, read, require, sha
from verification import stage5_campaign as campaign


def utc():
    return datetime.now(timezone.utc).isoformat()


def result_paths(study,job):
    case,wall,seed=job
    raw=study/"native/results"/case/wall/f"{seed}.nc"
    return raw,raw.with_suffix(".json")


def fingerprints(study,job):
    raw,ledger=result_paths(study,job)
    return dict(raw_path=raw.relative_to(study).as_posix(),raw_sha256=sha(raw),
                native_receipt_path=ledger.relative_to(study).as_posix(),native_receipt_sha256=sha(ledger))


def validate_receipt(study,job,seal,receipt):
    require(receipt["job"]==list(job) and receipt["execution_sha256"]==seal,"wrong dispatch job/seal")
    require(all(receipt[k]==v for k,v in fingerprints(study,job).items()),"native artifacts changed")


def receipt_path(study,job):
    case,wall,seed=job
    return study/"receipts"/case/wall/f"{seed}.json"


def bindings(study):
    verify=campaign.ROOT/"verification"
    code=[verify/x for x in ("stage5_campaign.py","stage5_inputs.py","stage5_dispatch.py",
                              "stage5_bridge_report.py","stage5_metrics.py")]
    code.extend(p for p in sorted(verify.glob("stage5_native*")) if p.is_file())
    snapshots={p.relative_to(campaign.ROOT).as_posix():sha(p) for p in code}
    inputs=study/"inputs"
    fixed={p.relative_to(study).as_posix():sha(p) for p in sorted(inputs.rglob("*")) if p.is_file()}
    from verification.stage5_native import wsl_path
    env_path=study/"execution-environment.json"
    fixed["execution-environment.json"]=sha(env_path)
    environment=read(env_path)
    library_paths=list(environment["native_files_sha256"])
    command=["wsl","-d","Ubuntu","-u","root","--","bash",wsl_path(verify/"stage5_native_enter.sh"),
             wsl_path(study),wsl_path(verify),"sha256sum",*library_paths]
    measured=subprocess.check_output(command,text=True)
    actual={line.split(maxsplit=1)[1].strip():line.split()[0] for line in measured.splitlines()}
    require(actual==environment["native_files_sha256"],"native FIREFLY, FLARE, MPI or python binary changed")
    fixed["bridge.json"]=sha(study/"bridge.json")
    for case in campaign.CASES:
        key=f"native/cases/{case}/ready.json"
        fixed[key]=sha(study/key)
    registered=read(study/"wall-registration.json")
    require(registered["protocol_sha256"]==sha(campaign.PAPER/"protocol.ko.md") and
            sha(study/"wall-registration.json")==sha(campaign.PAPER/"wall-registration.json"),
            "preregistration differs from published copy")
    panel=read(inputs/"panel.json")
    campaign.verify_panel(registered["roles"],panel)
    require(registered["predictions"]=={c["case_id"]:c["predicted_family"] for c in panel["cases"]},
            "changed original prediction")
    jobs=campaign.jobs(registered["roles"])
    require(registered["jobs"]==[list(j) for j in jobs],"run order differs")
    return dict(code=snapshots,files=fixed,registration_sha256=sha(study/"wall-registration.json"),
                protocol_sha256=registered["protocol_sha256"],jobs=registered["jobs"],
                precision_rule_pp=1.,seeds=registered["seeds"],histories=registered["histories_per_seed"],
                interpreter=dict(path=sys.executable,sha256=sha(sys.executable)))


def execution(study,create=False):
    plan=study/"execution.json"
    expected=bindings(study)
    if not plan.exists():
        require(create and not (study/"native/results").exists(),"cannot create retrospective execution seal")
        for name,digest in expected["code"].items():
            target=study/"source-code"/name
            target.parent.mkdir(parents=True,exist_ok=True)
            if not target.exists():
                shutil.copyfile(campaign.ROOT/name,target)
            require(sha(target)==digest,"source snapshot differs")
        head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=campaign.ROOT,text=True).strip()
        public_save(plan,dict(schema="nf-stage5-native-execution.v1",created_utc=utc(),revision=head,**expected))
        (study/"execution.sha256").write_text(sha(plan)+"\n",encoding="ascii")
    record=read(plan)
    require(sha(plan)==(study/"execution.sha256").read_text().strip(),"changed execution seal")
    require(all(record[k]==v for k,v in expected.items()),"sealed code, inputs or ready dossiers changed")
    for name,digest in record["code"].items():
        require(sha(study/"source-code"/name)==digest,"changed source snapshot")
    return record,sha(plan)


def execute(study,job,seal):
    dest=receipt_path(study,job)
    if dest.exists():
        validate_receipt(study,job,seal,read(dest))
        return
    case,wall,seed=job
    attempts=study/"dispatch-attempts"/case/wall/str(seed)
    attempts.mkdir(parents=True,exist_ok=True)
    path=attempts/f"{1+max((int(d.name) for d in attempts.iterdir() if d.name.isdigit()),default=0):03d}"
    path.mkdir()
    cmd=[sys.executable,"-B","-u","-m","verification.stage5_native","heat","--study",str(study),
         "--case",case,"--wall",wall,"--seed",str(seed)]
    launch=dict(job=list(job),execution_sha256=seal,command=cmd,start_utc=utc())
    public_save(path/"launch.json",launch)
    start=perf_counter()
    print("START",case,wall[:10],seed,launch["start_utc"],flush=True)
    with (path/"stdout.log").open("xb") as log:
        status=subprocess.run(cmd,cwd=campaign.ROOT,env={**os.environ,"PYTHONDONTWRITEBYTECODE":"1"},
                              stdout=log,stderr=subprocess.STDOUT,check=False)
    completion=dict(elapsed_s=perf_counter()-start,ended_utc=utc(),returncode=status.returncode)
    public_save(path/"completion.json",completion)
    require(status.returncode==0,"native solver stopped; inspect "+str(path/"stdout.log"))
    result=dict(job=list(job),execution_sha256=seal,attempt=path.relative_to(study).as_posix(),
                launch_sha256=sha(path/"launch.json"),completion_sha256=sha(path/"completion.json"),
                **completion,**fingerprints(study,job))
    public_save(dest,result)
    print("DONE",case,wall[:10],seed,result["ended_utc"],flush=True)


def run(study,workers):
    with process_lock(study/"workflow.lock"),keep_awake():
        plan,seal=execution(study)
        jobs=[tuple(j) for j in plan["jobs"]]
        pending=[]
        for j in jobs:
            dest=receipt_path(study,j)
            if dest.exists():
                validate_receipt(study,j,seal,read(dest))
            else:
                pending.append(j)
        completed=len(jobs)-len(pending)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            it=iter(pending)
            futures={}
            def submit():
                j=next(it,None)
                if j is not None:
                    futures[pool.submit(execute,study,j,seal)]=j
            for _ in range(workers):
                submit()
            while futures:
                done,_=wait(futures,return_when=FIRST_COMPLETED)
                for task in done:
                    last=futures.pop(task)
                    task.result()
                    completed+=1
                    public_save(study/"progress.json",dict(completed=completed,total=len(jobs),last=last,utc=utc()))
                    submit()
        print("COMPLETE",completed,"native seed jobs",utc(),flush=True)


def check(study):
    plan,seal=execution(study)
    for j in plan["jobs"]:
        record=read(receipt_path(study,j))
        validate_receipt(study,j,seal,record)
        path=study/record["attempt"]
        require(sha(path/"launch.json")==record["launch_sha256"] and
                sha(path/"completion.json")==record["completion_sha256"] and record["returncode"]==0,
                "native process provenance differs")
    print("Verified",len(plan["jobs"]),"native dispatch receipts and original netCDFs")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("seal","run","check"))
    parser.add_argument("--study",type=Path,default=campaign.STUDY)
    parser.add_argument("--workers",type=int,choices=(1,2),default=2)
    args=parser.parse_args()
    study=args.study.resolve()
    if args.action=="seal":
        print("Execution seal",execution(study,create=True)[1])
    elif args.action=="run":
        run(study,args.workers)
    else:
        check(study)
