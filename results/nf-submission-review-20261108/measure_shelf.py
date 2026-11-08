"""Measure direct first-hit shelf power on the frozen failed walls, before broadening.

From the worktree root with its research Python environment:
python -m results.nf-submission-review-20261108.measure_shelf 332a6d97fc 51d83022ca
"""
import argparse
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from d3d_optimization.continuous_boundary import BoundaryAwareEvaluator
from d3d_optimization.case import load_case
from d3d_optimization.independent_sources import split_design_wall
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.tangent_geometry import Curve
from verification.nf_equilibrium_times import REVIEW, SHOT
from verification.nf_search_benchmark import prep


REPO = Path(__file__).resolve().parents[2]
BUNDLE = REPO / "results/nf-prospective-search-20261107"
LEVELS = {"cubic320": 0, "cubic640": 1, "cubic1280": 2}
FAILED = ("332a6d97fc", "51d83022ca")


def checked_root(case, input_study):
    root = input_study / "cases" / case["case_id"] / "input-root" if input_study else Path(case["input_root"])
    for name, expected in case["hashes"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"input digest mismatch: {name}")
    return root


def endpoint_probe(root, case):
    review = json.loads((root / REVIEW).read_text(encoding="utf8"))
    shelf = next(s for s in split_design_wall(review)[1] if s["name"] == "shelf")
    curve = Curve(*np.asarray(shelf["endpoints_RZ_m"]), [0, 0, 0, .039])
    field, _, _ = load_case(root / SHOT, case["time_ms"])

    def max_axisymmetric_incidence(method):
        magnetic_field = field if method == "cubic" else load_case(root / SHOT, case["time_ms"], method)[0]
        b = magnetic_field.magnetic([curve.b])[0]
        return float(np.degrees(np.arctan2(np.linalg.norm(b[:2]), abs(b[2]))))

    def angle(tangent):
        return float(field.incidence([curve.b], [[-tangent[1], tangent[0]]])[0])

    return dict(case=case["case_id"], coefficients_m=[0, 0, 0, .039],
                min_radius_mm=1000 / curve.max_curvature(), length_ratio=curve.length_ratio(),
                analytic_deg=angle(curve.tangent(1)),
                max_any_axisymmetric_normal_deg={m: max_axisymmetric_incidence(m) for m in ("cubic", "quintic")},
                last_chord_deg={n: angle(curve.b - curve.point(1 - 1/(n - 1))) for n in (161, 321, 641, 1281)})


def measure(root, case, level):
    observed = []
    native_evaluate = BoundaryAwareEvaluator.evaluate

    def observe(evaluator, coefficients):
        result = native_evaluate(evaluator, coefficients)
        if np.any(coefficients):
            raise ValueError("expected the unmodified reference wall")
        shelf = result["surface_segments"]["shelf"]
        source, mapping = next((s, m) for s, m in zip(evaluator.sources, result["maps"], strict=True)
                               if np.isclose(s.width, 0.003))
        power = np.asarray(mapping.hit_power)
        shelf_hits = np.isin(mapping.segment, shelf) & (power > 0)
        if not np.isfinite(mapping.angles[shelf_hits]).all():
            raise ValueError("nonfinite angle on illuminated shelf")
        below = shelf_hits & (mapping.angles < 1.5)
        if not np.isclose(source.weights.sum(), 1.0) or not np.isclose(power.sum(), 1.0):
            raise ValueError("expected the full unit-power direct map")
        observed.append(dict(case=case["case_id"], level=level, donors=len(power),
                             shelf_hits=int(shelf_hits.sum()), below_hits=int(below.sum()),
                             shelf_kW_at_1MW=1000 * float(power[shelf_hits].sum()),
                             below_kW_at_1MW=1000 * float(power[below].sum()),
                             min_shelf_deg=float(mapping.angles[shelf_hits].min())))
        return result

    with patch.object(BoundaryAwareEvaluator, "evaluate", observe):
        try:
            prep(root, "OMP", ReferencePoint(), LEVELS[level])
        except ValueError as error:
            if str(error) != "continuous-boundary baseline is infeasible":
                raise
        else:
            raise RuntimeError("expected the held-out baseline preparation failure")
    if len(observed) != 1:
        raise RuntimeError(f"expected a single reference-wall observation; got {len(observed)}")
    return observed[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", choices=FAILED)
    parser.add_argument("--levels", nargs="+", choices=list(LEVELS), default=list(LEVELS))
    parser.add_argument("--input-study", type=Path, help="alternate extracted prospective input study")
    parser.add_argument("--endpoint-only", action="store_true", help="only compare shelf endpoint normals")
    args = parser.parse_args()
    cases = json.loads((BUNDLE / "summary.json").read_text(encoding="utf8"))["cases"]
    for name in args.cases:
        case = next(c for c in cases if c["case_id"] == name)
        root = checked_root(case, args.input_study)
        if args.endpoint_only:
            print(json.dumps(endpoint_probe(root, case), allow_nan=False), flush=True)
            continue
        for level in args.levels:
            print(json.dumps(measure(root, case, level), allow_nan=False), flush=True)
