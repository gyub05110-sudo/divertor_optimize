"""Seal and run the prospective DIII-D benchmark, one independent CLI process per arm."""
import argparse
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import runpy
import shutil
import subprocess
import sys
from time import perf_counter, sleep

ROOT = Path(__file__).resolve().parents[1]
SEEDS = (20711, 20729, 20747)
ARMS = ("floor-screened", "full-screened", "full-direct")
SETTINGS = dict(budget_seconds=600, radius_mm=20, power_mw=1., density_m3=1e19,
                chi_m2_s=.3, lambda_q_mm=3., gamma=7., kappa0=2000., backend="gpu")
PROTOCOL = "results/nf-stage3-20261108/protocol.ko.md"
INPUT_REGISTRATION_SHA = "35f77695555964d5bf6e6a8da71c94d6442624d25477802bc70f8fcb48709a83"
ARTIFACTS = {"result_sha256": "result.json", "events_sha256": "events.jsonl", "stdout_sha256": "stdout.log"}


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def canonical_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def require(value, message):
    if not value:
        raise ValueError(message)


def atomic_bytes(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("wb") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    for retry in range(9):
        try:
            os.replace(temp, path)
            return
        except PermissionError:
            if retry == 8:
                raise
            sleep(.01 * 2 ** retry)


def public_save(path, value):
    atomic_bytes(path, json.dumps(value, indent=2, allow_nan=False).encode("utf-8"))


@contextmanager
def resource_lock(path):
    """OS-released lock in the science process; a resumed worker waits for its predecessor."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as stream:
        if stream.tell() == 0:
            stream.write(b"0")
            stream.flush()
        stream.seek(0)
        if os.name == "nt":
            import errno
            import msvcrt
            while True:
                try:
                    msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError as exc:
                    if exc.errno not in (errno.EACCES, errno.EAGAIN, errno.EDEADLK):
                        raise
                    sleep(.25)
        else:
            import fcntl
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        yield


def build_tasks(cases):
    require(len(cases) == 12, "expected twelve registered slots")
    ids = [c["id"] for c in cases]
    require(len(set(ids)) == 12 and all(re.fullmatch(r"[A-Za-z0-9_-]+", i) for i in ids), "invalid case IDs")
    require([c["rank"] for c in cases] == list(range(1, 13)), "expected registered rank order")
    tasks = []
    for i, c in enumerate(cases):
        order = ARMS if i % 2 == 0 else ARMS[::-1]
        for j, seed in enumerate(SEEDS):
            offset = (i+j) % 3
            for arm in order[offset:] + order[:offset]:
                tasks.append(dict(case_id=c["id"], input_rank=c["rank"], registration_status=c["status"],
                                  registration_reason=c.get("reason"), equilibrium=c["path"], sha256=c["sha256"],
                                  time_ms=c["time_ms"], seed=seed, arm=arm, budget_seconds=600,
                                  family="floor-only" if arm == "floor-screened" else "full",
                                  method="direct" if arm == "full-direct" else "screened",
                                  relative_output=f"runs/{c['id']}/{seed}/{arm}"))
    return tasks


def study_files():
    paths = set(ROOT.glob("*.py")) | set((ROOT/"d3d_optimization").glob("*.py"))
    paths.update((ROOT/"verification").glob("nf_*.py"))
    paths.update(ROOT/p for p in [PROTOCOL, "verification/stage3_campaign.py", "verification/stage3_inputs.py",
                                 "verification/stage3_environment.py", "verification/stage3_report.py",
                                 "verification/stage3_replay.py", "verification/user_shape_replay.py",
                                 "verification/user_shape_results.py",
                                 "d3d_optimization/data/DIIID_geom.json", "baseline_review/optimization_region_policy.json",
                                 "baseline_review/baseline_review.json", "requirements-repro.txt"])
    paths.update((ROOT / PROTOCOL).with_name(n) for n in
                 ("method-freeze.ko.md", "execution-environment.json", "commissioning.json"))
    return sorted(paths)


def make_plan(study, *, root=ROOT, files=None, revision=None):
    source = study/"inputs"/"registration.json"
    cases = read(source)["cases"]
    if files is None:
        require(sha(source) == INPUT_REGISTRATION_SHA, "unexpected Stage 3 input registration")
        files = study_files()
    return dict(schema_version=1, registered_utc=utc(), registration_sha256=sha(source),
                revision=revision or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                files={p.relative_to(root).as_posix(): sha(p) for p in files},
                settings=SETTINGS, seeds=list(SEEDS), arms=list(ARMS), tasks=build_tasks(cases),
                python=dict(executable=str(Path(sys.executable).resolve()), sha256=sha(sys.executable),
                            version=sys.version), worker_environment=worker_env_settings(study))


def check_plan(study, plan, root=ROOT):
    source = study/"inputs"/"registration.json"
    require(sha(source) == plan["registration_sha256"], "input registration changed")
    if plan["registration_sha256"] == INPUT_REGISTRATION_SHA:
        require(set(plan["files"]) == {p.relative_to(root).as_posix() for p in study_files()},
                "production source inventory changed")
    require(plan["tasks"] == build_tasks(read(source)["cases"]), "schedule differs from registration")
    require(plan["settings"] == SETTINGS, "settings changed")
    require(plan["seeds"] == list(SEEDS) and plan["arms"] == list(ARMS), "method or seed declaration changed")
    require(plan["python"] == dict(executable=str(Path(sys.executable).resolve()), sha256=sha(sys.executable),
                                   version=sys.version) and plan["worker_environment"] == worker_env_settings(study),
            "worker interpreter or execution environment changed")
    for name, expected in plan["files"].items():
        require(sha(root/name) == expected, f"source digest mismatch: {name}")


def seal(study):
    path = study/"expected.json"
    fingerprint = study/"expected.sha256"
    with resource_lock(study/"schedule.lock"):
        if fingerprint.exists():
            return load_seal(study)
        require(not any((study/n).exists() for n in ("attempts", "runs", "observations")),
                "a missing plan seal cannot be repaired after run activity")
        plan = read(path) if path.is_file() else make_plan(study)
        check_plan(study, plan)
        for task in plan["tasks"][::9]:
            require(sha(task["equilibrium"]) == task["sha256"], "input bytes changed before seal")
        for name, expected in plan["files"].items():
            target = study/"source-code"/name
            if not target.exists():
                atomic_bytes(target, (ROOT/name).read_bytes())
            require(sha(target) == expected, f"source snapshot mismatch: {name}")
        if not path.exists():
            public_save(path, plan)
        atomic_bytes(fingerprint, (sha(path)+"\n").encode("ascii"))
    return plan


def load_seal(study):
    path = study/"expected.json"
    require((study/"expected.sha256").read_text().strip() == sha(path), "run plan checksum mismatch")
    plan = read(path)
    check_plan(study, plan)
    for name, expected in plan["files"].items():
        require(sha(study/"source-code"/name) == expected, f"source snapshot changed: {name}")
    return plan


def make_receipt(folder, task, plan, *, status, returncode=None, elapsed=None, start_utc=None,
                 end_utc=None, detail=None):
    return dict(status=status, case_id=task["case_id"], arm=task["arm"], seed=task["seed"],
                job_sha256=canonical_hash(task), registration_sha256=plan["registration_sha256"],
                plan_sha256=canonical_hash(plan),
                returncode=returncode, elapsed_s=elapsed, start_utc=start_utc, end_utc=end_utc, detail=detail,
                **{k: sha(folder/name) if folder and (folder/name).is_file() else None for k, name in ARTIFACTS.items()})


def excluded_receipt(task, plan):
    require(task["registration_status"] != "ready", "cannot exclude ready input")
    return make_receipt(None, task, plan, status="registration_excluded",
                        detail=task["registration_reason"] or task["registration_status"])


def verified_receipt(folder, task, plan):
    receipt = read(folder/"receipt.json")
    require(receipt["job_sha256"] == canonical_hash(task) and
            receipt["registration_sha256"] == plan["registration_sha256"] and
            receipt["plan_sha256"] == canonical_hash(plan), "receipt provenance mismatch")
    require(all(receipt[k] == task[k] for k in ("case_id", "seed", "arm")), "incorrect receipt identity")
    require(receipt["status"] in ("registration_excluded","execution_failure","completed"), "unknown receipt status")
    require((receipt["status"] == "registration_excluded") == (task["registration_status"] != "ready"),
            "receipt and registration status disagree")
    for key, name in ARTIFACTS.items():
        path = folder/name
        require(receipt[key] == (sha(path) if path.is_file() else None), f"changed {name} in {folder}")
    if receipt["status"] == "completed":
        require(all(receipt[k] for k in ARTIFACTS), "completed process lacks evidence")
    if task["registration_status"] == "ready":
        study = folder.parents[5] if folder.name == "output" else folder.parents[3]
        source = Path(receipt["source_attempt"])
        expected_parent = Path("attempts")/task["case_id"]/str(task["seed"])/task["arm"]
        require(source.parent == expected_parent and source.name.isdigit() and
                receipt["published_output"] == task["relative_output"], "attempt path mismatch")
        attempt = study/source
        require(sha(attempt/"launch.json") == receipt["launch_sha256"] and
                read(attempt/"command.json") == command(task, attempt/"output"), "changed launch evidence")
        launch = read(attempt/"launch.json")
        require(launch["job_sha256"] == canonical_hash(task) and launch["plan_sha256"] == canonical_hash(plan)
                and launch["started_utc"] == receipt["start_utc"], "changed launch identity")
        if receipt["timing_source"] == "cli_invocation":
            require(sha(attempt/"completion.json") == receipt["completion_sha256"], "changed completion evidence")
            measured = read(attempt/"completion.json")
            require((measured["elapsed_s"],measured["returncode"],measured["ended_utc"]) ==
                    (receipt["elapsed_s"],receipt["returncode"],receipt["end_utc"]), "changed completion timing")
        else:
            require(receipt["timing_source"] == "recovered_cli_result" and receipt["status"] == "completed" and
                    not (attempt/"completion.json").exists() and all(receipt[k] is None for k in
                    ("completion_sha256","elapsed_s","end_utc","returncode")), "invalid recovered timing")
        if receipt["status"] == "completed":
            result = validate_result(folder, task, attempt/"output")
            if receipt["returncode"] is not None:
                require(receipt["returncode"] == {"calculation_error": 1,"partial_result": 1,"input_error": 2}
                        .get(result["status"], 0), "return code differs from result")
    return receipt


def command(task, out):
    cmd = [sys.executable, "-B", "-u", str(ROOT/"optimize_divertor.py"),
           "--equilibrium", str(task["equilibrium"]), "--time-ms", str(task["time_ms"]),
           "--seed", str(task["seed"]), "--family", task["family"], "--method", task["method"],
           "--output", str(out)]
    for k, v in SETTINGS.items():
        cmd.extend(["--"+k.replace("_", "-"), str(v)])
    return cmd


def worker_env_settings(study):
    env = {}
    cache = study/"runtime"
    for key, leaf in (("TEMP","temp"), ("TMP","temp"), ("MPLCONFIGDIR","matplotlib"),
                      ("CUPY_CACHE_DIR","cupy"), ("CUDA_CACHE_PATH","cuda")):
        d = cache/leaf
        env[key] = str(d)
    deps = ROOT/".local_dependencies"/"d3d-fixed-equilibrium"
    env["PYTHONPATH"] = str(deps)
    env["DIVERTOR_GPU_RUNTIME"] = str(Path(os.environ.get("DIVERTOR_GPU_RUNTIME",
                                          "D:/DivertorOptimization/runtime/cupy")).resolve())
    env.update(OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1", PYTHONIOENCODING="utf-8")
    return env


def environment(study):
    settings = worker_env_settings(study)
    for key in ("TEMP", "MPLCONFIGDIR", "CUPY_CACHE_DIR", "CUDA_CACHE_PATH"):
        Path(settings[key]).mkdir(parents=True, exist_ok=True)
    return {**os.environ, **settings}


def validate_result(staging, task, attempted_output):
    result = read(staging/"result.json")
    flags = result["arguments"]
    require(all(flags[k] == v for k,v in SETTINGS.items()) and
            all(flags[k] == task[k] for k in ("family", "seed", "method", "time_ms")), "CLI arguments mismatch")
    require(Path(flags["equilibrium"]).resolve() == Path(task["equilibrium"]).resolve(), "input path mismatch")
    require(Path(flags["output"]).resolve() == attempted_output.resolve(), "output path mismatch")
    if result.get("inputs"):
        require(result["inputs"]["equilibrium_sha256"] == task["sha256"], "CLI used a different equilibrium")
    require((staging/"events.jsonl").is_file(), "completed CLI result lacks event log")
    for name in ("baseline", "selected"):
        value = result.get(name) or {}
        if "directory" in value:
            require(value["directory"] == name and value["sha256"], "unexpected exported evidence directory")
            for filename, expected in value["sha256"].items():
                require(Path(filename).name == filename and sha(staging/name/filename) == expected,
                        "exported file checksum differs")
    return result


def finish_attempt(study, attempt, task, plan):
    staging = attempt/"output"
    launch = read(attempt/"launch.json")
    require(launch["job_sha256"] == canonical_hash(task) and launch["plan_sha256"] == canonical_hash(plan)
            and read(attempt/"command.json") == command(task, staging), "attempt is not bound to scheduled CLI")
    staging.mkdir(exist_ok=True)
    shutil.copyfile(attempt/"stdout.log", staging/"stdout.log")
    complete = (attempt/"completion.json").is_file()
    observed = read(attempt/"completion.json") if complete else {}
    error = observed.get("error")
    rc = observed.get("returncode")
    result = None
    if (staging/"result.json").is_file():
        result = validate_result(staging, task, staging)
        expected_code = {"calculation_error": 1, "partial_result": 1, "input_error": 2}.get(result["status"], 0)
        if complete and rc != expected_code:
            error = f"exit code {rc} disagrees with status {result['status']}"
    status = "completed" if result is not None and not error else "execution_failure"
    receipt = make_receipt(staging, task, plan, status=status, returncode=rc, elapsed=observed.get("elapsed_s"),
                           start_utc=launch["started_utc"], end_utc=observed.get("ended_utc"),
                           detail=error or (result or {}).get("reason") or "no CLI result")
    receipt.update(source_attempt=attempt.relative_to(study).as_posix(), published_output=task["relative_output"],
                   timing_source="cli_invocation" if complete else "recovered_cli_result",
                   completion_sha256=sha(attempt/"completion.json") if complete else None,
                   launch_sha256=sha(attempt/"launch.json"))
    public_save(staging/"receipt.json", receipt)
    return receipt


def execute_cli(attempt, task, plan):
    staging = attempt/"output"
    cmd = command(task, staging)
    public_save(attempt/"command.json", cmd)
    public_save(attempt/"launch.json", dict(started_utc=utc(), pid=os.getpid(), job_sha256=canonical_hash(task),
                                           plan_sha256=canonical_hash(plan)))
    rc, error = None, None
    previous = sys.argv
    with (attempt/"stdout.log").open("x", encoding="utf8", buffering=1) as log:
        with redirect_stdout(log), redirect_stderr(log):
            start = perf_counter()
            try:
                sys.argv = cmd[3:]
                runpy.run_path(cmd[3], run_name="__main__")
                rc = 0
            except SystemExit as exc:
                rc = 0 if exc.code is None else exc.code
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
                import traceback
                traceback.print_exc()
            finally:
                sys.argv = previous
            elapsed = perf_counter()-start
            ended = utc()
            log.flush()
            os.fsync(log.fileno())
    public_save(attempt/"completion.json", dict(returncode=rc, elapsed_s=elapsed, ended_utc=ended, error=error))


def worker(study, index):
    # A detached coordinator may restart while the old worker is still writing its
    # receipt. Serialize workers as well as coordinators to avoid overlapping GPU use.
    with resource_lock(study/"worker.lock"):
        plan = load_seal(study)
        task = plan["tasks"][index]
        dest = study/task["relative_output"]
        if dest.exists():
            return verified_receipt(dest, task, plan)
        dest.parent.mkdir(parents=True, exist_ok=True)
        attempt_root = study/"attempts"/task["case_id"]/str(task["seed"])/task["arm"]
        attempt_root.mkdir(parents=True, exist_ok=True)
        for item in sorted(attempt_root.iterdir()):
            staging = item/"output"
            if (staging/"receipt.json").exists():
                receipt = verified_receipt(staging, task, plan)
                staging.rename(dest)
                return receipt
            if (staging/"result.json").is_file() or (item/"completion.json").is_file():
                receipt = finish_attempt(study, item, task, plan)
                staging.rename(dest)
                return receipt
        attempt = attempt_root/f"{max((int(i.name) for i in attempt_root.iterdir()), default=0)+1:03d}"
        attempt.mkdir()
        staging = attempt/"output"
        if task["registration_status"] != "ready":
            staging.mkdir()
            receipt = excluded_receipt(task, plan)
            public_save(staging/"receipt.json", receipt)
        else:
            require(sha(task["equilibrium"]) == task["sha256"], "equilibrium changed before launch")
            execute_cli(attempt, task, plan)
            receipt = finish_attempt(study, attempt, task, plan)
        staging.rename(dest)
        return receipt


def run(study):
    from d3d_optimization.campaign_state import keep_awake, process_lock
    with process_lock(study/"campaign.lock"), keep_awake():
        plan = load_seal(study)
        counts = {}
        for index, task in enumerate(plan["tasks"]):
            dest = study/task["relative_output"]
            if not dest.exists():
                print(f"START {index+1}/{len(plan['tasks'])} {task['relative_output']} {utc()}", flush=True)
                child = [sys.executable, "-B", "-u", "-m", "verification.stage3_campaign", "worker",
                         "--study", str(study), "--index", str(index)]
                # An orphaned *science worker* keeps the OS lock after coordinator loss.
                # Wait before starting the watchdog so an earlier live run owns its full budget.
                with resource_lock(study/"worker.lock"):
                    pass
                if not dest.exists():
                    observe_worker(study, task, child)
            row = verified_receipt(dest, task, plan)
            label = read(dest/"result.json")["status"] if row["status"] == "completed" else row["status"]
            counts[label] = counts.get(label, 0)+1
            public_save(study/"progress.json", dict(completed=index+1, total=len(plan["tasks"]),
                        last=task["relative_output"], counts=counts, utc=utc()))
            print(f"FINISH {index+1}/{len(plan['tasks'])} {task['relative_output']} {label} {utc()}", flush=True)
            if row["status"] == "execution_failure":
                print("Infrastructure failure retained for review:", row["detail"], flush=True)


def observe_worker(study, task, command_line):
    folder = study/"observations"/task["case_id"]/str(task["seed"])/task["arm"]
    folder.mkdir(parents=True, exist_ok=True)
    number = len(list(folder.glob("*.started.json"))) + 1
    prefix = f"{number:03d}"
    identity = dict(job_sha256=canonical_hash(task), number=number, started_utc=utc(),
                    clock_scope="parent process creation through worker exit", command=command_line)
    public_save(folder/(prefix+".started.json"), identity)
    start = perf_counter()
    try:
        result = subprocess.run(command_line, cwd=ROOT, env=environment(study),
                                timeout=SETTINGS["budget_seconds"]+600, check=False)
    except (subprocess.TimeoutExpired, OSError) as exc:
        public_save(folder/(prefix+".finished.json"), dict(**identity, ended_utc=utc(), elapsed_s=perf_counter()-start,
                                                           error=f"{type(exc).__name__}: {exc}", returncode=None))
        raise
    public_save(folder/(prefix+".finished.json"), dict(**identity, ended_utc=utc(), elapsed_s=perf_counter()-start,
                                                       returncode=result.returncode))
    result.check_returncode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("seal", "run", "worker", "check"))
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    study = args.study.resolve()
    if args.action == "seal":
        print("Sealed", len(seal(study)["tasks"]), "tasks")
    elif args.action == "run":
        run(study)
    elif args.action == "worker":
        if args.index is None or not 0 <= args.index < 108:
            parser.error("worker requires a planned --index from 0 to 107")
        worker(study, args.index)
    else:
        plan = load_seal(study)
        for task in plan["tasks"]:
            verified_receipt(study/task["relative_output"], task, plan)
        print("Verified", len(plan["tasks"]), "process or exclusion receipts")


if __name__ == "__main__":
    main()
