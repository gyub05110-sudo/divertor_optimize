"""Test the necessary incidence condition at the fixed DIII-D OMP shelf end."""
import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from d3d_optimization.broadening_study import digest, load_json
from d3d_optimization.case import load_case
from d3d_optimization.independent_sources import split_design_wall
from d3d_optimization.tangent_geometry import Curve
from verification.nf_anytime_reference import EvidenceStore, bind_evidence, require
from verification.nf_prospective_campaign import verify_registration
from verification.nf_unused_live import SHOT, REVIEW


def anchored_shelf(a, b):
    ends = [0., 1.]
    base = Curve(a,b,np.zeros(4))
    positions, tangents = [], []
    for vector in .1*np.eye(4):
        curve = Curve(a,b,vector)
        positions.append(float(np.abs(curve.point(ends)-base.point(ends)).max()))
        tangents.append(float(np.abs(curve.tangent(ends)-base.tangent(ends)).max()))
    return dict(basis_count=4, test_amplitude_m=.1, max_position_difference_m=max(positions),
                max_tangent_difference=max(tangents), linear_in_four_shelf_coefficients=True)


def first_passing_distance(incidence, floor, length):
    distances = np.linspace(0., length, 2001)
    values = np.asarray(incidence(distances)) - floor
    passing = np.flatnonzero(values >= 0.)
    if not len(passing):
        return None
    i = passing[0]
    if not i:
        return 0.
    return float(brentq(lambda x: float(incidence(x))-floor, distances[i-1], distances[i], xtol=1e-13))


def check_case(study, case):
    root = Path(case['root'])
    diag_path = study/'baseline-diagnostic/results'/case['case_id']/'prepare_fine.json'
    native = None
    if diag_path.exists():
        diag = EvidenceStore(study/'baseline-diagnostic').read(diag_path.relative_to(study/'baseline-diagnostic').as_posix())
        require(diag['all_primary_native_errors_reproduced'] and diag['input_sha256']==case['hashes'], 'diagnostic inputs or outcome differ')
        native = diag['native_observations']['continuous'][-1]['prescribed_direct']['surfaces']['shelf']
    review = load_json(root/REVIEW)
    shelf = next(s for s in split_design_wall(review)[1] if s['name']=='shelf')
    a,b = np.asarray(shelf['endpoints_RZ_m'])
    require(shelf['original_segment']==59 and np.array_equal(b,review['wall'][60]), 'different shelf or corner')
    bounds = anchored_shelf(a,b)
    require(bounds['max_position_difference_m'] < 1e-14 and bounds['max_tangent_difference'] < 1e-14,
            'shelf endpoints or tangents respond to design variables')
    length = np.linalg.norm(a-b)
    unit = (a-b)/length
    normal = np.array([-unit[1],unit[0]])
    floor = shelf['angle_floor_deg']
    by_method = {}
    for method in ('cubic','quintic'):
        field, _, _ = load_case(root/SHOT,case['time_ms'],method=method)

        def incidence(distance):
            ds = np.asarray(distance)
            points = b + np.atleast_1d(ds)[:,None]*unit
            values = field.incidence(points,np.broadcast_to(normal,points.shape))
            return float(values[0]) if ds.ndim==0 else values

        at_end = incidence(0.)
        first = first_passing_distance(incidence,floor,length)
        row = dict(end_60_deg=at_end, end_59_deg=incidence(length), margin_60_deg=at_end-floor,
                   baseline_distance_to_first_passing_angle_from_60_mm=None if first is None else 1000*first)
        if native:
            point = np.asarray(native['minimum_hit']['point_RZ_m'])
            calculated = float(field.incidence([point],[normal])[0])
            row['incidence_at_native_minimum_RZ_deg'] = calculated
            if method=='cubic':
                require(abs(calculated-native['minimum_angle_deg']) < 1e-7, 'native shelf incidence differs')
        by_method[method] = row
    return dict(case_id=case['case_id'],input_hashes=case['hashes'],fixed_end_RZ_m=b.tolist(),angle_floor_deg=floor,
                basis=bounds, original_shelf_length_m=length, original_wall_shelf_3mm=native,
                comparison=by_method, diagnostic_sha256=digest(diag_path) if native else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,required=True)
    study = parser.parse_args().study.resolve()
    require(study.drive.lower()=='d:', 'store the check on D')
    store = EvidenceStore(study)
    frozen = store.read('frozen.json')
    verify_registration(store,frozen)
    raw = EvidenceStore(study/'baseline-diagnostic').read('summary.json')
    require(raw is not None, 'finish primary baseline diagnostics first')
    result = dict(frozen_sha256=store.sha256('frozen.json'), diagnostic_sha256=digest(study/'baseline-diagnostic/summary.json'),
                  source_sha256=digest(Path(__file__)), clipped_or_changed_constraints=False,
                  scope='fixed-point incidence conditional on that shelf neighbourhood remaining powered; the distance is measured on the original straight shelf',
                  cases=[check_case(study,c) for c in frozen['cases']])
    bind_evidence(EvidenceStore(study/'fixed-shelf-corner'),'summary.json',result)
    for row in result['cases']:
        print(row['case_id'], {k: round(v['end_60_deg'],6) for k,v in row['comparison'].items()})


if __name__=='__main__':
    main()
