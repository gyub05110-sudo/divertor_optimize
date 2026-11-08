"""Reconstruct the four public CLI checks from their actual geometry and power exports."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from verification.user_shape_results import audit_measurement, fine_receipt, verify_geometry
from d3d_optimization.campaign_geometry import Candidate

EXPECTED = {"entry-demo": "qualified", "entry-release": "qualified",
            "entry-332a6d97fc": "constraint_excluded", "entry-51d83022ca": "constraint_excluded"}


def verify(folder):
    summary = {}
    for name, status in EXPECTED.items():
        run = folder / name
        row = json.loads((run/"result.json").read_text(encoding="utf-8"))
        if row["status"] != status:
            raise ValueError(f"unexpected return state: {name}")
        baseline = row["baseline"]
        ref = audit_measurement(run / baseline["directory"], baseline, row["reference_SI"]["power_W"])
        if not baseline["physics_measurable"] or not baseline["reference_only"] or baseline["candidate_eligible"]:
            raise ValueError(f"incorrect original wall state: {name}")
        shelf = row["fixed_shelf"]
        stored = dict(status=status, sha256=row["inputs"]["equilibrium_sha256"],
                      baseline_MW_m2=ref.metrics["peak_MW_m2"], A2_m2=ref.metrics["A2_m2"],
                      baseline_ready_s=row["timings"]["baseline_ready_s"],
                      fixed_actual_deg=shelf["actual_deg"], pitch_max_deg=shelf["max_deg"],
                      floor_deg=shelf["floor_deg"], unhit_W=ref.metrics["unhit_W"],
                      direct_subfloor_W=row["reference_SI"]["power_W"]*baseline["shelf_exposure"]["direct"]["subfloor_power_fraction"],
                      received_subfloor_midpoint_W=row["reference_SI"]["power_W"]*baseline["shelf_exposure"]["receiving"]["midpoint_subfloor_cell_power_fraction"])
        if status == "qualified":
            candidate = Candidate.from_record(row["candidate"])
            measured = audit_measurement(run/"selected", row["selected"], row["reference_SI"]["power_W"])
            verify_geometry(ref, measured, candidate)
            receipt = fine_receipt(run, row)
            if (candidate.identity != row["candidate"]["id"] or not row["selected"]["candidate_eligible"]
                    or not measured.metrics["peak_MW_m2"] < ref.metrics["peak_MW_m2"]):
                raise ValueError(f"incorrect selected geometry: {name}")
            stored.update(peak_MW_m2=measured.metrics["peak_MW_m2"], id=candidate.identity,
                          first_feasible_s=row["search"]["first_feasible_elapsed_s"], fine=receipt,
                          budget_s=row["arguments"]["budget_seconds"])
        else:
            if row["selected"] is not None or row["candidate"] is not None or not shelf["actual_deg"] < shelf["floor_deg"]:
                raise ValueError(f"incorrect excluded state: {name}")
            if name.endswith("51d83022ca") != (shelf["max_deg"] < shelf["floor_deg"]):
                raise ValueError("incorrect distinction between tangent and pitch limits")
        summary[name] = stored
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.study), indent=2, allow_nan=False))
