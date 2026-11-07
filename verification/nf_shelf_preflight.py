"""Check fixed shelf incidence from staged DIII-D inputs before starting a search."""
import argparse
import json
from pathlib import Path

from verification.nf_shelf_corner_check import check_case
from verification.nf_unused_live import SHOT, input_record


def inspect_input(root):
    root = Path(root).resolve()
    record = input_record(root)
    key = 'frame-' + record['hashes'][SHOT][:16]
    case = check_case(root, dict(record,case_id=key))
    angles = {method:values['end_60_deg'] for method,values in case['comparison'].items()}
    below = [m for m,x in angles.items() if x < case['angle_floor_deg']]
    statuses = ('at_or_above_floor','interpolants_disagree','below_floor_if_powered')
    return dict(input_root=str(root), shot_sha256=record['hashes'][SHOT], time_ms=record['time_ms'],
                shelf_corner_RZ_m=case['fixed_end_RZ_m'], minimum_angle_deg=case['angle_floor_deg'],
                incidence_deg=angles, interpolants_below_floor=below,
                corner_status=statuses[len(below)],
                geometry_or_candidate_feasibility='requires first-hit and receiving-power checks', heat_evaluations=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root',type=Path,required=True,nargs='+',help='staged physical input directories')
    args = parser.parse_args()
    print(json.dumps([inspect_input(root) for root in args.input_root],indent=2,allow_nan=False))


if __name__=='__main__':
    main()
