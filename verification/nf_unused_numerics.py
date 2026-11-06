"""Post-trial sensitivity of every admitted design to a fixed alternative numerical setting."""
import argparse
from pathlib import Path
from time import perf_counter
from dataclasses import asdict

from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.boundary_evaluator import prepare_boundary
from d3d_optimization.campaign_geometry import Candidate, evaluator
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_base_budget import score
from verification.nf_search_benchmark import measurement
from verification.validate_boundary_transport import LEVELS


def chosen_candidates(study, case_id, seeds):
    chosen={}
    source_hashes={}
    for seed in seeds:
        for policy in ("direct","screen"):
            path=study/"live"/case_id/str(seed)/f"{policy}.json"
            record=load_json(path)
            source_hashes[str(path.relative_to(study))]=digest(path)
            incumbent=record["incumbent"]
            if incumbent and incumbent["kind"]=="new_design":
                matches=[row for row in record["generated"] if row["id"]==incumbent["id"]]
                if len(matches)!=1:
                    raise ValueError("final design missing or duplicated")
                candidate=Candidate.from_record(matches[0]["candidate"])
                if candidate.identity!=incumbent["id"]:
                    raise ValueError("final design identity mismatch")
                if candidate.identity not in chosen:
                    chosen[candidate.identity]=dict(candidate=candidate.record(),native=incumbent,uses=[])
                chosen[candidate.identity]["uses"].append(dict(seed=seed,policy=policy))
    return chosen,source_hashes


def run(study):
    campaign=load_json(study/"live/protocol.json")
    summary=load_json(study/"summary.json")
    if not all(row["status"]=="completed" for row in summary["rows"]):
        raise ValueError("complete and audited paired trial records required")
    output=study/"numerics"
    output.mkdir(exist_ok=False)
    reference=ReferencePoint()
    _,count,step,nodes,wall_step=LEVELS[1]
    for case in campaign["cases"]:
        started=perf_counter()
        chosen,hashes=chosen_candidates(study,case["case_id"],campaign["seeds"])
        original=prepare_boundary(Path(case["root"]),"OMP",reference,count=count,step=step,nodes=nodes,
                                  wall_step=wall_step,method="quintic")
        baseline=measurement(original.baseline,reference)
        result=dict(case_id=case["case_id"],method="quintic",count=count,step=step,nodes=nodes,
                    wall_step=wall_step,reference=asdict(reference),baseline=baseline,
                    native_summary_sha256=digest(study/"summary.json"),native_arm_sha256=hashes,
                    driver_sha256=digest(__file__),protocol_sha256=digest(study/"numerics-protocol.ko.md"),
                    candidates={})
        gpu=evaluator(original,20,"gpu")
        cpu=evaluator(original,20,"cpu")
        for identity,row in chosen.items():
            candidate=Candidate.from_record(row["candidate"])
            alt=score(gpu,candidate)
            check=score(cpu,candidate)
            if alt["eligible"]!=check["eligible"]:
                raise ValueError("CPU/GPU eligibility mismatch")
            parity={}
            if alt["eligible"]:
                for key in ("peak","A2"):
                    error=abs(alt[key]/check[key]-1)
                    if error>1e-10:
                        raise ValueError("CPU/GPU "+key+" mismatch")
                    parity[key+"_relative_difference"]=error
            entry=dict(row,alternative=alt,cpu_check=check,parity=parity)
            if alt["eligible"] and baseline["eligible"]:
                entry["peak_reduction_pct"]=100*(1-alt["peak"]/baseline["peak"])
                entry["A2_change_pct"]=100*(alt["A2"]/baseline["A2"]-1)
            result["candidates"][identity]=entry
            print(case["case_id"],identity[:8],alt["eligible"],entry.get("peak_reduction_pct"),
                  entry.get("A2_change_pct"),flush=True)
        result["elapsed_s"]=perf_counter()-started
        save(output/(case["case_id"]+".json"),result,exclusive=True)
        del gpu,cpu,original


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",type=Path,required=True)
    args=p.parse_args()
    run(args.study.resolve())

