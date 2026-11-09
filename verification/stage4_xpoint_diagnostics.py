"""Development-only cubic stationary-point census for the three archived Stage 3 lower-X failures."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from shapely.geometry import Point, Polygon

from d3d_optimization.user_input import _parquet
from d3d_optimization.user_topology import hessian, upper_saddle_evidence


DEVELOPMENT_IDS = ("8d1575d8c4", "a7e51f17f7", "8bb901c2a1")
SUBDIVISIONS = 12


def _possible_boxes(field, boxes, precondition=False):
    # On a tensor cubic knot rectangle each gradient component has degree <= 3
    # in each variable. Recover its 4x4 Bernstein tensor on this subrectangle.
    t = np.linspace(0., 1., 4)
    basis = np.c_[(1-t)**3, 3*t*(1-t)**2, 3*t*t*(1-t), t**3]
    inverse = np.linalg.inv(basis)
    uv = np.stack(np.meshgrid(t, t, indexing="ij"), axis=-1)
    points = boxes[:, None, None, :2]*(1-uv) + boxes[:, None, None, 2:]*uv
    samples = field.gradient(points.reshape(-1, 2)).reshape(-1, 4, 4, 2)
    scale = np.maximum(1., np.abs(samples).max(axis=(1, 2)))
    if precondition:
        # Two nearly parallel gradient-component zero sets can each cross an
        # adjacent box. A fixed nonsingular combination per box shares every
        # common zero and is still bicubic. Respect the transform in the padding.
        centres = (boxes[:, :2]+boxes[:, 2:])/2
        matrices = np.array([hessian(field, c) for c in centres])
        good = np.linalg.cond(matrices) < 1e8
        transform = np.tile(np.eye(2), (len(boxes), 1, 1))
        transform[good] = np.linalg.inv(matrices[good])
        samples = np.einsum("nab,nijb->nija", transform, samples)
        scale = np.maximum(1., np.einsum("nab,nb->na", np.abs(transform), scale))
    coefficients = np.einsum("ai,nijc,bj->nabc", inverse, samples, inverse)
    # The convex-hull property bounds each component on the WHOLE box. Pad
    # recovered coefficients for floating point evaluation and the transform.
    pad = 128*np.finfo(float).eps*np.linalg.norm(inverse, ord=np.inf)**2*scale
    possible = (coefficients.min(axis=(1, 2)) <= pad) & (coefficients.max(axis=(1, 2)) >= -pad)
    return boxes[np.all(possible, axis=1)]


def _bisect(boxes):
    lo, hi = boxes[:, :2], boxes[:, 2:]
    mid = (lo+hi)/2
    children = []
    for right, up in ((0, 0), (0, 1), (1, 0), (1, 1)):
        choose = np.array([right, up], dtype=bool)
        children.append(np.c_[np.where(choose, mid, lo), np.where(choose, hi, mid)])
    return np.concatenate(children)


def _box_distances(boxes, point):
    excess = np.maximum(np.maximum(boxes[:, :2]-point, point-boxes[:, 2:]), 0.)
    return np.linalg.norm(excess, axis=1)


def stationary_enclosures(field):
    """Enclose possible cubic stationary points in physical knot rectangles, then polish roots in-grid."""
    if field.method != "cubic":
        raise ValueError("whole-grid Bernstein enclosure requires the cubic EFIT field")
    z, r = (np.unique(k) for k in field.spline.get_knots())
    rlo, zlo = np.meshgrid(r[:-1], z[:-1])
    rhi, zhi = np.meshgrid(r[1:], z[1:])
    boxes = np.c_[rlo.ravel(), zlo.ravel(), rhi.ravel(), zhi.ravel()]
    initial = len(boxes)
    boxes = _possible_boxes(field, boxes)
    for _ in range(SUBDIVISIONS):
        if not len(boxes):
            break
        if len(boxes) > 20000:
            raise ValueError("unresolved degenerate field exceeds the development census budget")
        boxes = _possible_boxes(field, _bisect(boxes))
    if len(boxes):
        boxes = _possible_boxes(field, boxes, precondition=True)
    bounds = ([r[0], z[0]], [r[-1], z[-1]])
    roots = []
    for seed in (boxes[:, :2]+boxes[:, 2:])/2:
        fit = least_squares(lambda p: field.gradient([p])[0], seed, jac=lambda p: hessian(field, p),
                            bounds=bounds, ftol=1e-14, xtol=1e-14, gtol=1e-14)
        if (fit.success and np.linalg.norm(fit.fun) <= 1e-7 and
                all(np.linalg.norm(root-fit.x) > 1e-7 for root in roots)):
            roots.append(fit.x)
    # Neighboring boxes can survive componentwise range tests even if the zero
    # curves do not intersect within that box. Report their spatial resolution.
    widths = np.linalg.norm(boxes[:, 2:]-boxes[:, :2], axis=1)
    nearest = np.full(len(boxes), np.inf)
    for root in roots:
        nearest = np.minimum(nearest, _box_distances(boxes, root))
    stats = dict(knot_rectangles=initial, subdivisions=SUBDIVISIONS, terminal_boxes=len(boxes),
                 max_terminal_box_diagonal_m=float(max(widths, default=0.)),
                 matched_root_distance_limit="one terminal box diagonal",
                 unmatched_terminal_boxes=int(np.count_nonzero(nearest > widths)),
                 distinct_polished_stationary_points=len(roots),
                 terminal_boxes_Rmin_Zmin_Rmax_Zmax_m=boxes.tolist())
    return boxes, roots, stats


def diagnose(field, lcfs, axis):
    boxes, roots, census = stationary_enclosures(field)
    polygon = Polygon(lcfs)
    bottom = lcfs[np.argmin(lcfs[:, 1])]
    boundary = field.values(lcfs)
    median = float(np.median(boundary))
    stationary = []
    for point in sorted(roots, key=lambda p: (p[1], p[0])):
        eigenvalues = np.linalg.eigvalsh(hessian(field, point))
        flux = float(field.values([point])[0])
        stationary.append(dict(RZ_m=point.tolist(), hessian_eigenvalues=eigenvalues.tolist(),
            saddle=bool(eigenvalues[0] < 0 < eigenvalues[1]), below_axis=bool(point[1] < axis[1]),
            distance_from_lcfs_lowest_m=float(np.linalg.norm(point-bottom)),
            distance_from_lcfs_boundary_m=float(polygon.boundary.distance(Point(point))),
            gradient_norm=float(np.linalg.norm(field.gradient([point]))), flux_Wb_per_rad=flux,
            flux_minus_lcfs_median_Wb_per_rad=flux-median))
    lower = [s for s in stationary if s["saddle"] and s["below_axis"]]
    census["possible_stationary_boxes_within_lower_search_radius"] = int(np.sum(_box_distances(boxes, bottom) <= .15))
    return dict(grid_bounds_RZ_m=[[field.r[0], field.z[0]], [field.r[-1], field.z[-1]]],
                grid_shape_RZ=[len(field.r), len(field.z)], axis_RZ_m=axis.tolist(),
                vertex_flux_median_Wb_per_rad=median, lcfs_vertex_flux_span_Wb_per_rad=float(np.ptp(boundary)),
                lcfs_lowest_RZ_m=bottom.tolist(), lower_search_radius_m=.15,
                census=census, stationary_points=stationary,
                nearest_lower_saddle=min(lower, key=lambda s: s["distance_from_lcfs_lowest_m"], default=None),
                upper_lcfs_evidence=upper_saddle_evidence(field, lcfs, axis))


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(study):
    plan = study/"expected.json"
    registration = study/"inputs"/"registration.json"
    if (digest(plan) != (study/"expected.sha256").read_text().strip() or
            digest(registration) != json.loads(plan.read_bytes())["registration_sha256"]):
        raise ValueError("Stage 3 expected plan or registration changed")
    registered = {case["id"]: case for case in json.loads(registration.read_bytes())["cases"]}
    cases = []
    for ident in DEVELOPMENT_IDS:
        case = registered[ident]
        path = study/"inputs"/"raw"/f"d3d_shot_{ident}.parquet"
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != case["sha256"] or case["time_ms"] != 1400.:
            raise ValueError("Stage 3 development input identity or time changed")
        field, lcfs, axis, _ = _parquet(raw, case["time_ms"])
        cases.append(dict(case_id=ident, equilibrium_sha256=case["sha256"], time_ms=case["time_ms"],
                          **diagnose(field, lcfs, axis)))
    return dict(schema="nf-stage4-development-lower-x.v1", stage3_expected_sha256=digest(plan),
                command="python -B -m verification.stage4_xpoint_diagnostics --stage3 STUDY --output OUTPUT",
                diagnostic_source_sha256=digest(Path(__file__)),
                scope="Three Stage 3 failures only. All physical cubic knot rectangles are screened by the "
                      "Bernstein coefficient ranges of both gradient components, padded for floating point. "
                      "Remaining boxes are quartered 12 times, checked again with a nonsingular preconditioned "
                      "gradient on well-conditioned boxes, and polished by in-grid least squares. These are "
                      "finite-precision enclosures and root counts, not exact arithmetic multiplicity proofs. "
                      "Upper LCFS evidence is descriptive. The lower-divertor acceptance criteria are unchanged.",
                cases=cases)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage3", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.stage3.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf8")
    for case in result["cases"]:
        counts = case["census"]
        print(case["case_id"], counts["distinct_polished_stationary_points"], "stationary points,",
              counts["possible_stationary_boxes_within_lower_search_radius"], "boxes in lower search radius")
