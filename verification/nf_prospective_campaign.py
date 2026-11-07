"""Register and supervise development and held-out divertor searches on one device."""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter_ns

from d3d_optimization.broadening_study import digest, load_json
from d3d_optimization.campaign_state import keep_awake, process_lock, public_save
from d3d_optimization.reference_transport import ReferencePoint
from verification import nf_anytime_campaign as prior
from verification.nf_anytime_final_audit import read_events
from verification.nf_anytime_reference import EvidenceStore, bind_evidence
from verification.nf_unused_live import INPUTS, load_cases
from verification.validate_boundary_transport import LEVELS

ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results/nf-prospective-search-20261107"
OLD_RESULT = ROOT / "results/nf-anytime-fullspace-20261106"
OLD_LIVE = Path("D:/DivertorOptimization/studies/nf-anytime-fullspace-20261106")
DEFAULT_DEV = Path("D:/DivertorOptimization/studies/unused-equilibrium-20261106")
DEFAULT_HOLDOUT = Path("D:/DivertorOptimization/studies/nf-prospective-inputs-20261107")
DEV_CASES = ("deb842cf1e", "fad1624516", "4d14a35fcb")
SEEDS = (20711, 20729, 20747)
CONFIGS = ((24, 8), (32, 12))
FILES = {"development":"development-plan.json", "holdout":"frozen.json"}
AUDITS = {"development":"development-audit.json", "holdout":"audit.json"}


def require(value, message):
    if not value:
        raise ValueError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def job_row(case, seed, policy, population, projected, budget, phase, directory):
    return dict(case_id=case["case_id"], input_root=case["root"], seed=seed, policy=policy,
                population=population, projected_seeds=projected, budget_s=budget,
                phase=phase, relative_dir=directory)


def development_jobs(cases):
    jobs = []
    for index, case in enumerate(cases):
        for size, projected in (CONFIGS if index % 2 == 0 else reversed(CONFIGS)):
            directory = f"development/{case['case_id']}/20611/p{size}-s{projected}"
            jobs.append(job_row(case, 20611, "fine_adaptive", size, projected, 600, "development", directory))
    return jobs


def holdout_jobs(cases, selected):
    fine = selected["population"], selected["projected_seeds"]
    require(selected["policy"] == "fine_adaptive" and fine in CONFIGS, "unregistered fine configuration")
    methods = (("coarse_adaptive", 64, 12), ("fine_adaptive", *fine))
    jobs = []
    for index, case in enumerate(cases):
        for j, seed in enumerate(SEEDS):
            budget = 1800 if seed == SEEDS[0] else 600
            for policy, count, projected in (methods if (index+j) % 2 == 0 else reversed(methods)):
                directory = f"live/{case['case_id']}/{seed}/{policy}"
                jobs.append(job_row(case, seed, policy, count, projected, budget, "holdout", directory))
    return jobs


def check_disjoint(development, holdout):
    require(len(development) == len(holdout) == 3, "three distinct inputs required in each stage")
    for key in (lambda case: case["case_id"], lambda case: case["hashes"][INPUTS[0]]):
        old, new = {key(c) for c in development}, {key(c) for c in holdout}
        require(len(old) == len(new) == 3 and not old.intersection(new), "reused case ID or equilibrium bytes")


def capture_files(store, paths):
    hashes = {}
    for path in sorted(paths):
        name = path.relative_to(ROOT).as_posix()
        hashes[name] = digest(path)
        dest = "code-snapshot/" + name
        if not store.exists(dest):
            store.write_bytes(dest, path.read_bytes())
        require(store.sha256(dest) == hashes[name], "snapshot differs from source: " + name)
    return hashes


def source_records(source, phase):
    names = set(prior.SOURCE_RECORDS)
    if phase == "holdout":
        names.update(p.name for p in source.glob("*.json"))
    return {name:digest(source/name) for name in sorted(names)}


def checked_development(source):
    archive = load_json(OLD_LIVE/"frozen.json")
    published = load_json(OLD_RESULT/"summary.json")["input_hashes"]["frozen.json"]
    require(digest(OLD_LIVE/"frozen.json") == published, "historical development inputs differ from published record")
    cases = load_cases(source)
    require(source_records(source, "development") == archive["source_records"], "development source metadata changed")
    require(cases == archive["cases"] and tuple(c["case_id"] for c in cases) == DEV_CASES,
            "historical development cases or hashes changed")
    return cases


def checked_holdout(source):
    cases = load_cases(source)
    receipt = load_json(source/"verification.json")
    registered = load_json(source/"registration.json")
    require(receipt["cases"] == [c["case_id"] for c in cases] and
            receipt["source_records"] == {name:digest(source/name) for name in prior.SOURCE_RECORDS},
            "held-out case order or verified input records changed")
    require(registered["staged_or_downloaded_candidates"] == 0 and
            registered["selection_sha256"] == receipt["source_records"]["selection.json"],
            "candidate names were not registered before acquisition")
    require(receipt["archived_sources"] >= 18 and receipt["archived_frames"] >= 3903 and
            receipt["valid_archived_frames"] == receipt["archived_frames"], "incomplete archival frame check")
    return cases


def new_or_matching(store, name, value):
    registered = store.read(name) if store.exists(name) else {}
    value = dict(value, created_utc=registered.get("created_utc") or utc())
    bind_evidence(store, name, value)
    return value


def register(study, source, phase):
    store = EvidenceStore(study)
    protocol = (RESULT_DIR/"protocol.ko.md").read_bytes()
    if not store.exists("protocol.ko.md"):
        store.write_bytes("protocol.ko.md", protocol)
    require(store.read_bytes("protocol.ko.md") == protocol, "written study rule changed")
    paths = prior.source_files() + [Path(__file__).resolve(), ROOT/"verification/nf_prospective_worker.py"]
    paths = set(prior.local_source_closure(paths))
    source_hashes = capture_files(store, paths)
    originals = capture_files(store, [ROOT/name for name in INPUTS])
    extra = {}
    if phase == "development":
        cases = checked_development(source)
        jobs, seeds = development_jobs(cases), [20611]
    else:
        dev = store.read(FILES["development"])
        verify_registration(store, dev)
        choice = selection(study, dev)
        require(choice["qualified"] and choice["selected"], "development comparator has not qualified")
        cases = checked_holdout(source)
        check_disjoint(dev["cases"], cases)
        previous_paths = load_json(Path(dev["source_study"])/"selection.json")["selected"]
        previous_ids = {Path(path).stem.removeprefix("d3d_shot_") for path in previous_paths} | {"1a3155370c"}
        require(not previous_ids.intersection(c["case_id"] for c in cases), "historical input in held-out selection")
        jobs, seeds = holdout_jobs(cases, choice["selected"]), list(SEEDS)
        extra = dict(development_choice=choice["selected"], development_choice_sha256=store.sha256("development-choice.json"),
                     development_plan_sha256=store.sha256(FILES["development"]))
    plan = dict(schema="nf-prospective-registration.v1", phase=phase, cases=cases, source_study=str(source),
                code_root=str(ROOT), source_records=source_records(source,phase), jobs=jobs, seeds=seeds,
                checkpoints_s=[200,400,600] if phase=="development" else [200,400,600,1800],
                source="OMP", radius_mm=20, reference=asdict(ReferencePoint()), levels={"coarse":list(LEVELS[0]),"fine":list(LEVELS[2])},
                protocol_sha256=store.sha256("protocol.ko.md"), source_hashes=source_hashes,
                original_reference_hashes=originals, python_version=platform.python_version(),
                packages=prior.versions(), hardware=prior.hardware(), **extra)
    return new_or_matching(store, FILES[phase], plan)


def verify_registration(store, plan):
    require(store.sha256("protocol.ko.md") == plan["protocol_sha256"], "protocol hash changed")
    for name, expected in (plan["source_hashes"] | plan["original_reference_hashes"]).items():
        require(digest(ROOT/name) == expected and store.sha256("code-snapshot/"+name) == expected, "code or snapshot changed: "+name)
    source = Path(plan["source_study"])
    for name, expected in plan["source_records"].items():
        require(digest(source/name) == expected, "registered input metadata changed: "+name)
    if plan["phase"] == "holdout":
        require(store.sha256("development-choice.json") == plan["development_choice_sha256"] and
                store.sha256(FILES["development"]) == plan["development_plan_sha256"], "development choice link changed")
    for case in plan["cases"]:
        prior.check_inputs(case)


def file_sha(path):
    return digest(path) if path.is_file() else None


def verify_process(out, expected_job):
    arm = EvidenceStore(out)
    require(arm.read("job.json") == expected_job, "job differs from phase registration")
    p = arm.read("process.json")
    require(p["job_sha256"] == arm.sha256("job.json") and
            p["registration_sha256"] == expected_job["registration_sha256"], "process/job binding differs")
    start, end = p["launch_ns"], p["exit_observed_ns"]
    require(isinstance(start, int) and isinstance(end, int) and 0 < start <= end and
            math.isclose(p["elapsed_s"],(end-start)/1e9,abs_tol=1e-9,rel_tol=0), "invalid process clock")
    for key, file in (("events_sha256", "events.jsonl"),("result_sha256","result.json"),("stdout_sha256","stdout.log")):
        require(p[key] == file_sha(Path(out)/file), "process evidence changed: "+file)
    require(p["stdout_sha256"] is not None, "process log missing")
    if p["returncode"] == 0:
        require(p["error"] is None and p["events_sha256"] and p["result_sha256"], "successful process lacks evidence")
    return p


def launch(out, job):
    out.mkdir(parents=True, exist_ok=False)
    arm = EvidenceStore(out)
    arm.write("job.json", job)
    job_hash = arm.sha256("job.json")
    caches = {"CUPY_CACHE_DIR":"cupy", "CUDA_CACHE_PATH":"cuda", "TEMP":"tmp", "TMP":"tmp",
              "MPLCONFIGDIR":"matplotlib","HF_HOME":"huggingface"}
    environment = {key:str(out/"cache"/leaf) for key,leaf in caches.items()}
    for directory in set(environment.values()):
        Path(directory).mkdir(parents=True, exist_ok=False)
    environment.update(NF_PROSPECTIVE_JOB_SHA256=job_hash, PYTHONDONTWRITEBYTECODE="1", PYTHONIOENCODING="utf8",
                       OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    command = [sys.executable, "-B", "-u", "-m","verification.nf_prospective_worker",
               "--policy",job["policy"],"--population",str(job["population"]),"--projected-seeds",str(job["projected_seeds"]),
               "--input-root",job["input_root"],"--seed",str(job["seed"]),"--out",str(out),"--budget-s",str(job["budget_s"]),"--started-ns"]
    print("START",job["phase"],job["case_id"],job["seed"],job["policy"],job["population"],job["budget_s"],utc(),flush=True)
    with (out/"stdout.log").open("x",encoding="utf8") as log:
        start_utc, start, returncode, error = utc(), perf_counter_ns(), None, None
        try:
            proc = subprocess.run([*command,str(start)],cwd=ROOT,env=os.environ | environment,
                                  stdout=log,stderr=subprocess.STDOUT,check=False,timeout=job["budget_s"]+600)
            returncode = proc.returncode
        except (OSError, subprocess.TimeoutExpired) as exc:
            error = f"{type(exc).__name__}: {exc}"
        finally:
            end = perf_counter_ns()
    receipt = dict(launch_ns=start, exit_observed_ns=end, elapsed_s=(end-start)/1e9, start_utc=start_utc, end_utc=utc(),
                   returncode=returncode, error=error, job_sha256=job_hash, registration_sha256=job["registration_sha256"],
                   environment=environment, events_sha256=file_sha(out/"events.jsonl"), result_sha256=file_sha(out/"result.json"),
                   stdout_sha256=file_sha(out/"stdout.log"))
    arm.write("process.json",receipt)
    print("FINISH",job["case_id"],job["seed"],job["policy"],"rc",returncode,"seconds",round(receipt["elapsed_s"],2),flush=True)


def assess(out, job, reference, levels):
    from verification.nf_prospective_worker import audit_arm
    p = verify_process(out, job)
    base = {key:job[key] for key in ("case_id","seed","policy","population","projected_seeds","budget_s","relative_dir")}
    base.update(process_sha256=digest(out/"process.json"), job_sha256=p["job_sha256"], events_sha256=p["events_sha256"],
                result_sha256=p["result_sha256"], parent_elapsed_s=p["elapsed_s"], launch_ns=p["launch_ns"], exit_observed_ns=p["exit_observed_ns"])
    if p["returncode"] not in (0,1) or not (out/"result.json").is_file():
        return dict(base, normal=False, startup_full_rank=False, eligible_thermal_de_trials=0, de_requests=0,
                    de_evidence="no verified worker result", peak_ratio=None, elapsed_s=None,
                    status="process_failed", error=p["error"] or f"child return code {p['returncode']}")
    record, (raw, _) = load_json(out/"result.json"), read_events(out/"events.jsonl")
    for key, expected in dict(case_id=job["case_id"], input_root=job["input_root"], launch_ns=p["launch_ns"],
                              registration_sha256=job["registration_sha256"], job_sha256=p["job_sha256"]).items():
        require(record[key]==expected, "child/parent identity differs: "+key)
    require(record["levels"]==levels, "registered calculation levels changed")
    if "reference" in record or any(e["stage"]=="prepare_fine" for e in raw):
        require(record.get("reference")==reference, "registered physics changed")
    require(record["elapsed_s"] <= p["elapsed_s"]+1e-8, "child exceeds parent duration")
    audit = audit_arm(record,raw, **{k:job[k] for k in ("policy","population","projected_seeds","budget_s","seed")})
    require(p["returncode"] == int(not audit["normal"]), "worker exit differs from audited completion")
    return dict(base | audit, status="completed" if audit["normal"] else "arm_failed", error=record.get("error"))


def audit_stage(study, plan):
    store = EvidenceStore(study)
    phase = plan["phase"]
    registration_hash = store.sha256(FILES[phase])
    rows, previous, baselines = [], 0, {}
    for registered in plan["jobs"]:
        job = dict(registered, registration_sha256=registration_hash)
        out = store.path(job["relative_dir"])
        require((out/"process.json").is_file(), "unclosed or unrun registered arm: "+job["relative_dir"])
        row = assess(out,job,plan["reference"],plan["levels"])
        require(row["launch_ns"] >= previous, "overlap or altered job ordering")
        previous = row["exit_observed_ns"]
        baseline = row.get("baseline_peak")
        if baseline is not None:
            original = baselines.setdefault(job["case_id"],baseline)
            require(math.isclose(original,baseline,rel_tol=1e-10,abs_tol=0), "original-wall peak changed between arms")
        rows.append(row)
    result = dict(schema="nf-prospective-audit.v1",phase=phase, registration_sha256=registration_hash,
                  rows=rows, planned_arm_count=len(plan["jobs"]), completed_arm_count=sum(r["normal"] for r in rows))
    bind_evidence(store,AUDITS[phase],result)
    return result


def selection(study, plan):
    from verification.nf_prospective_worker import choose_fine
    store = EvidenceStore(study)
    summary = audit_stage(study,plan)
    choice = choose_fine(summary["rows"])
    return new_or_matching(store,"development-choice.json",dict(schema="nf-prospective-choice.v1",
                           development_plan_sha256=store.sha256(FILES["development"]),
                           development_audit_sha256=store.sha256(AUDITS["development"]), **choice))


def run_stage(study, plan):
    store = EvidenceStore(study)
    verify_registration(store,plan)
    cases = {c["case_id"]:c for c in plan["cases"]}
    registration_hash = store.sha256(FILES[plan["phase"]])
    for index, registered in enumerate(plan["jobs"]):
        job = dict(registered, registration_sha256=registration_hash)
        case = cases[job["case_id"]]
        prior.check_inputs(case)
        out = store.path(job["relative_dir"])
        if not out.exists():
            launch(out, job)
        else:
            require((out/"process.json").is_file(), "preserving incomplete attempt; inspect before resuming: "+str(out))
        result = assess(out,job,plan["reference"],plan["levels"])
        prior.check_inputs(case)
        public_save(study/"progress.json",dict(phase=plan["phase"],finished=index+1,planned=len(plan["jobs"]),
                    latest=job["relative_dir"],status=result["status"],updated_utc=utc()))
    verify_registration(store,plan)
    return audit_stage(study,plan)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,required=True)
    parser.add_argument("--action",choices=("freeze-development","development","freeze-holdout","holdout","audit"),required=True)
    parser.add_argument("--development-source",type=Path,default=DEFAULT_DEV)
    parser.add_argument("--holdout-source",type=Path,default=DEFAULT_HOLDOUT)
    args=parser.parse_args()
    study=args.study.resolve()
    require(study.drive.lower()=="d:", "keep prospective study on D")
    study.mkdir(parents=True,exist_ok=True)
    with process_lock(study/"campaign.lock"), keep_awake():
        if args.action=="audit":
            store=EvidenceStore(study)
            for phase in FILES:
                plan=store.read(FILES[phase])
                verify_registration(store,plan)
                summary=audit_stage(study,plan)
                print("AUDIT",phase,summary["completed_arm_count"],"/",summary["planned_arm_count"],flush=True)
            return
        phase="development" if args.action.endswith("development") else "holdout"
        source=(args.development_source if phase=="development" else args.holdout_source).resolve()
        plan=register(study,source,phase)
        print("REGISTERED",phase,len(plan["jobs"]),"arms",flush=True)
        if args.action.startswith("freeze"):
            return
        summary=run_stage(study,plan)
        print("AUDIT",phase,summary["completed_arm_count"],"/",len(plan["jobs"]),flush=True)
        if phase=="development":
            result=selection(study,plan)
            print("DEVELOPMENT CHOICE",result["qualified"],result["selected"],flush=True)


if __name__=="__main__":
    main()

