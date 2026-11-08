"""Lower-OMP, single-operating-point adapter for a supplied DIII-D equilibrium.

The original material is a thermal reference. Only a separately evaluated
replacement can satisfy the spline, wall, field-incidence and power constraints.
"""

from time import perf_counter

import numpy as np
import shapely
from scipy.optimize import brentq
from shapely.geometry import LineString, Polygon

from .boundary_evaluator import BoundaryEvaluator
from .campaign_geometry import evaluate_candidate
from .case import separatrix
from .continuous_boundary import BoundaryAwareEvaluator, BoundaryGuard
from .frozen_field import exponential_source
from .frozen_mapping import target_angle_status, trace_paths
from .independent_sources import split_design_wall
from .multisurface import MultiSurfaceEvaluator
from .tangent_geometry import Curve


SUPPORT = ("A positive-width periodic Gaussian from any positive first-hit donor "
           "has nonzero density on every open interval of the connected receiving wall.")


def _shelf(surfaces):
    return next(surface for surface in surfaces if surface["name"] == "shelf")


def shelf_limits(model_input, method="cubic"):
    """Necessary condition at fixed endpoints under the active receiving policy."""
    field, _, _ = model_input.field(method)
    shelf = _shelf(split_design_wall(model_input.review)[1])
    curve = Curve(*shelf["endpoints_RZ_m"], np.zeros(4))
    floor = shelf["angle_floor_deg"]
    points = curve.point([0., 1.])
    magnetic = field.magnetic(points)
    norm = np.linalg.norm(magnetic, axis=1)
    if not np.all(np.isfinite(norm) & (norm > 0)):
        raise ValueError("undefined magnetic field at fixed shelf endpoint")
    actual = field.incidence(points, np.broadcast_to(curve.n, points.shape))
    maximum = np.degrees(np.arcsin(np.clip(np.linalg.norm(magnetic[:, :2], axis=1) / norm, 0, 1)))
    ends = [dict(u=i, RZ_m=p.tolist(), normal_RZ=curve.n.tolist(), B_RZphi_T=b.tolist(),
                 actual_deg=float(a), max_deg=float(m), passes_fixed_normal=bool(a >= floor))
            for i, (p, b, a, m) in enumerate(zip(points, magnetic, actual, maximum))]
    worst = ends[int(np.argmin(actual))]
    return dict(admissible=all(e["passes_fixed_normal"] for e in ends), method=method,
                actual_deg=worst["actual_deg"], max_deg=worst["max_deg"], floor_deg=floor,
                limiting_u=worst["u"], endpoints=ends, support=SUPPORT,
                assumptions="Fixed endpoint positions and tangent of every allowed shelf spline; "
                            "continuous field and one-sided shelf normal; strictly positive receiving support. "
                            "An endpoint strictly below the floor implies a violating open neighbourhood.",
                scope="Both endpoints passing is a necessary condition only.")


class _SingleSourceMap(MultiSurfaceEvaluator):
    """Initialize the archived intersection kernel with exactly one imposed source."""

    def __init__(self, field, lcfs, wall, surfaces, paths, inside, source, nodes, solids):
        self.field, self.lcfs, self.wall = field, lcfs, wall
        self.paths, self.inside = paths, inside
        self.vessel_solids = solids
        self.surfaces, self.sources = surfaces, [source]
        self.dimension, self.bound = len(surfaces) * 4, .1
        self.targets = {s["name"]: np.linspace(*s["endpoints_RZ_m"], nodes) for s in surfaces}
        self.lengths = {name: float(np.linalg.norm(t[-1] - t[0])) for name, t in self.targets.items()}
        adjustable = {s["segment_index"] for s in surfaces}
        self.fixed_indices = np.array([i for i in range(len(wall)-1) if i not in adjustable])
        self.fixed_geometry = shapely.union_all([LineString(wall[i:i+2]) for i in self.fixed_indices])
        self.fixed_hits = self._fixed_hits()


class _UserDirect(BoundaryAwareEvaluator):
    """Attach the separatrix witness without imposing a design gate on the original."""

    def __init__(self, base, xp, psi, rz, step):
        self.base = base
        self.boundary_guard = BoundaryGuard(base, xp, psi, "OMP", rz, min(.0005, step / 4))


def _engineering_regions(review, design_field, design_psi):
    """Keep OMP portions immutable; label the two fixed inner portions for reporting."""
    leg = review["legs"]["inner"]
    if leg["segment"] != 69:
        return review, "supplied_fixed_regions"
    wall = np.asarray(review["wall"])
    a, b = wall[69:71]
    v = b-a
    hit = np.asarray(leg["strike_RZ_m"], dtype=float)
    fraction = float((hit-a) @ v / (v @ v))
    if (not np.isfinite(fraction) or not -1e-10 <= fraction <= 1+1e-10 or
            np.linalg.norm(hit-(a+fraction*v)) > 1e-7):
        raise ValueError("resolved inner first hit does not lie on its stated material segment")
    fraction = float(np.clip(fraction, 0, 1))
    grid = np.linspace(0, 1, 257)
    values = design_field.values(a+grid[:,None]*v) - design_psi
    bracket = np.flatnonzero(values[:-1]*values[1:] < 0)
    source = "resolved_inner_first_hit_on_69"
    if len(bracket) == 1:
        j = bracket[0]
        # This is the historical cubic root on the commissioned unique-root case.
        root = brentq(lambda u: float(design_field.values([a+u*v])[0]-design_psi),
                      grid[j], grid[j+1], xtol=1e-14)
        if np.linalg.norm(a+root*v-hit) <= 1e-4:
            fraction = root
            source = "unique_inner_flux_root_on_69"
    if fraction >= 1-1e-12:
        # A zero-length inactive named portion cannot be meshed. Keep the two
        # supplied reporting portions when the strike coincides with this corner.
        return review, "supplied_fixed_regions_at_inner_corner"
    va, vb = wall[70:72]
    cutoff = review["optimization_policy"]["imp_maximum_Z_m"]
    upper = va+(cutoff-va[1])/(vb[1]-va[1])*(vb-va)
    portions = {69: [(a+fraction*v).tolist(), b.tolist()], 70: [va.tolist(), upper.tolist()]}
    rows = [dict(row, points=portions[row["segment"]]) if row["segment"] in portions else row
            for row in review["optimization_footprints"]]
    return dict(review, optimization_footprints=rows), source


def prepare_user(model_input, reference, level):
    method, count, step, nodes, wall_step = level
    if (method not in ("cubic", "quintic") or type(count) is not int or count < 8
            or type(nodes) is not int or nodes < 9 or not np.isfinite([step, wall_step]).all()
            or min(step, wall_step) <= 0):
        raise ValueError("finite positive resolution, cubic/quintic field, >=8 donors and >=9 nodes required")
    field, lcfs, axis = model_input.field(method)
    rsep, xp, psi = separatrix(field, lcfs, axis)
    rz = [float(rsep), float(axis[1])]
    review = model_input.review
    # Both IMP segments remain fixed. Reproduce the historical split when
    # applicable, including when using quintic fields for sensitivity.
    design_field, design_psi = field, psi
    if method != "cubic" and review["legs"]["inner"]["segment"] == 69:
        design_field, design_lcfs, design_axis = model_input.field("cubic")
        design_psi = separatrix(design_field, design_lcfs, design_axis)[2]
    review, inner_region_basis = _engineering_regions(review, design_field, design_psi)
    wall, surfaces, original_segments = split_design_wall(review)
    crossing = LineString(wall).intersection(LineString([rz, [field.r[-1], axis[1]]]))
    points = list(crossing.geoms) if hasattr(crossing, "geoms") else [crossing]
    outward = [p.x-rsep for p in points if p.geom_type == "Point" and not p.is_empty and p.x > rsep]
    if not outward:
        raise ValueError("no outboard midplane intersection with supplied material wall")
    limit = min(.1, min(outward))
    edges = limit * np.linspace(0, 1, count+1)**2
    source = exponential_source(field, *rz, edges, reference.lambda_q_m)
    starts = np.column_stack([rsep+source.x, np.full(count, axis[1])])
    inside = shapely.contains(Polygon(wall), shapely.points(starts))
    if not np.all(inside):
        raise ValueError("OMP source starts outside supplied material boundary")
    paths = [LineString(p) for p in trace_paths(field, starts, step=step, max_length=8)]
    solids = shapely.union_all([Polygon(v["pts"]) for v in model_input.geometry["vv"]])
    base = _SingleSourceMap(field, lcfs, wall, surfaces, paths, inside, source, nodes, solids)
    direct = _UserDirect(base, xp, psi, rz, step)
    metadata = dict(source="OMP", source_RZ_m=rz, shot_sha256=model_input.equilibrium_sha256,
                    time_ms=model_input.time_ms, geometry_sha256=model_input.provenance["geometry_sha256"],
                    source_power_W=reference.power_W, widths_mm=[reference.lambda_q_m * 1000],
                    support_mm=limit * 1000, tubes=count, step_m=step, target_nodes=nodes,
                    dimension=base.dimension, interpolation=method, surfaces=surfaces,
                    split_to_original_segment=original_segments.tolist(), policy=review["optimization_policy"],
                    fixed_inner_region_basis=inner_region_basis)
    e = BoundaryEvaluator(direct, metadata, reference, wall_step)
    r = e.baseline
    if r["geometry_valid"] and "maps" in r:
        base.baseline = r
        base.baseline_peaks = np.array([float(np.max(m.q)) for m in r["maps"]])
        direct.baseline = r
        direct.baseline_peaks = np.array([r.get("direct_peak_q", np.nan)])
    shelf = _shelf(surfaces)
    curve = Curve(*shelf["endpoints_RZ_m"], np.zeros(4))
    e.baseline = _with_exposure(e, r, curve, reference_only=True)
    e.cache[tuple(np.zeros(base.dimension))] = e.baseline
    return e


def _parameters(curve, points):
    # Normal displacement is perpendicular to the reference chord. Projection
    # gives the exact longitudinal u on both the transport chord and the spline.
    return np.clip((np.asarray(points) - curve.a) @ curve.v / (curve.v @ curve.v), 0, 1)


def _signed_sine(field, curve, u):
    tangent = curve.tangent(np.atleast_1d(u))
    normal = np.column_stack([-tangent[:, 1], tangent[:, 0]])
    magnetic = field.magnetic(curve.point(np.atleast_1d(u)))
    value = np.sum(magnetic[:, :2] * normal, axis=1) / (
        np.linalg.norm(magnetic, axis=1) * np.linalg.norm(normal, axis=1))
    if not np.isfinite(value).all():
        raise ValueError("undefined analytic shelf incidence")
    return np.clip(value, -1, 1)


def _degrees(sine):
    return np.degrees(np.arcsin(np.abs(sine)))


def _refine_shelf_intervals(field, curve, u):
    """Bounded numerical minimization inside every sampled parameter interval."""
    left, right = u[:-1].copy(), u[1:].copy()
    golden = (5**.5 - 1) / 2
    low = right - golden * (right-left)
    high = left + golden * (right-left)
    f_low = np.abs(_signed_sine(field, curve, low))
    f_high = np.abs(_signed_sine(field, curve, high))
    # Maximum initial width is 1/160. After 32 reductions each remaining
    # bracket is <1.3e-9 in u. Each iteration evaluates one point per bracket.
    for _ in range(32):
        keep_left = f_low <= f_high
        left = np.where(keep_left, left, low)
        right = np.where(keep_left, high, right)
        probe = np.where(keep_left, right-golden*(right-left), left+golden*(right-left))
        measured = np.abs(_signed_sine(field, curve, probe))
        low, high = np.where(keep_left, probe, high), np.where(keep_left, low, probe)
        f_low, f_high = np.where(keep_left, measured, f_high), np.where(keep_left, f_low, measured)
    j = int(np.argmin(np.minimum(f_low, f_high)))
    if f_low[j] <= f_high[j]:
        return float(f_low[j]), float(low[j])
    return float(f_high[j]), float(high[j])


def _shelf_exposure(field, curve, floor, r):
    ids = r["surface_segments"]["shelf"]
    heat = r["maps"][0]
    direct_mask = np.isin(heat.segment, ids)
    direct_angles = np.full(len(direct_mask), np.nan)
    donor_u = _parameters(curve, heat.points[direct_mask])
    direct_angles[direct_mask] = _degrees(_signed_sine(field, curve, donor_u))
    direct = target_angle_status(direct_angles, direct_mask, heat.hit_power, floor)
    direct.update(angle_floor_deg=floor, power_fraction=float(heat.hit_power[direct_mask].sum()),
                  subfloor_power_fraction=float(heat.hit_power[direct_mask & (direct_angles < floor)].sum()))

    dep = r["deposition"]
    cells = np.flatnonzero(np.isin(dep["cell_parent"], ids))
    if not len(cells):
        raise ValueError("no receiving cells on shelf")
    samples = np.isin(dep["sample_parent"], ids)
    # The exposure test includes the entire shelf even when floating-point
    # integration rounds the deposited power on a cell to exactly zero.
    u = np.unique(np.r_[curve.spline.t, np.linspace(0, 1, 161), _parameters(curve, dep["points"][samples])])
    sine = _signed_sine(field, curve, u)
    magnitude = np.abs(sine)
    i = int(np.argmin(magnitude))
    minimum, location = float(magnitude[i]), float(u[i])
    crossings = np.flatnonzero(np.signbit(sine[1:]) != np.signbit(sine[:-1]))
    if len(crossings):
        j = crossings[0]
        location = float(brentq(lambda x: _signed_sine(field, curve, x)[0], u[j], u[j+1]))
        minimum = 0.
    elif minimum > 0:
        refined, u_refined = _refine_shelf_intervals(field, curve, u)
        if refined < minimum:
            minimum, location = refined, u_refined
    angle = float(_degrees(minimum))
    midpoint = (dep["wall"][cells] + dep["wall"][cells+1]) / 2
    midpoint_angles = _degrees(_signed_sine(field, curve, _parameters(curve, midpoint)))
    power = dep["cell_power"][cells]
    receiving = dict(minimum_angle_deg=angle, minimum_u=location, minimum_RZ_m=curve.point(location).tolist(),
                     endpoint_angles_deg=_degrees(sine[[0,-1]]).tolist(), passed=bool(angle >= floor),
                     checked_cells=len(cells), checked_zero_numeric_power_cells=int(np.count_nonzero(power == 0)),
                     power_fraction=float(power.sum()),
                     midpoint_subfloor_cell_power_fraction=float(power[midpoint_angles < floor].sum()),
                     refined_parameter_intervals=(len(u)-1 if not len(crossings) and magnitude[i] > 0 else 0),
                     maximum_parameter_interval=float(np.diff(u).max()),
                     minimum_method="analytic normals at endpoints, spline knots and receiving samples; "
                                    "sign-change roots or bounded minimization on every intervening interval",
                     minimum_scope="numerical minimum; multiple unresolved extrema in one interval are possible",
                     power_diagnostic="sum of represented cell powers with subfloor analytic incidence at midpoint")
    exposure = dict(floor_deg=floor, support=SUPPORT, normal="analytic spline tangent at projected longitudinal u",
                    direct=direct, receiving=receiving)
    return exposure, direct_angles


def _with_exposure(e, native, curve, reference_only=False):
    r = dict(native, reference_only=reference_only, physics_measurable=False, feasible=False)
    if "deposition" not in r:
        return r
    exposure, direct_angles = _shelf_exposure(e.field, curve, _shelf(e.surfaces)["angle_floor_deg"], r)
    r["shelf_exposure"] = exposure
    r["shelf_direct_incidence_deg"] = direct_angles
    r["angles"] = dict(r["angles"], shelf=exposure["direct"])
    r["receiving_angles"] = dict(r["receiving_angles"], shelf=exposure["receiving"])
    # Rebuild each archived non-angle predicate explicitly. The archived
    # aggregate includes chord normals and rounded positive-power masks.
    errors = [item for item in r["errors"] if item != "receiving_angle_constraint_failed"]
    ledger = r["deposition"]["ledger"]
    checks = {
        "unhit_power": r["missing"] <= .001 and ledger["unhit"] <= .001,
        "power_balance": r["power_balance_error"] <= 1e-4 and ledger["balance_error"] <= 2e-12,
        "unresolved_heat_peak": np.isfinite(r["peak_q"]) and r["peak_q"] > 0,
    }
    r["physics_measurable"] = bool(all(checks.values()))
    witness = r["boundary_witness"]
    checks.update(
        no_direct_design_power=any(any(p > 0 for p in row) for row in r["surface_power"].values()),
        unresolved_separatrix_witness=witness["resolved"],
        separatrix_angle_constraint_failed=witness.get("angle_passed", True),
        direct_shelf_angle_constraint_failed=exposure["direct"]["angle_constraint_passed"] is not False,
        receiving_angle_constraint_failed=exposure["receiving"]["passed"],
    )
    errors.extend(reason for reason, passed in checks.items() if not passed and reason not in errors)
    r["errors"] = errors
    r["feasible"] = bool(not reference_only and not errors and r["geometry_valid"]
                         and r["shape_constraints"]["passed"] and r["geometry_screen"]["passed"])
    r["penalty"] = 0. if not errors else 1000. * len(errors)
    return r


def evaluate_user_candidate(e, candidate):
    start = perf_counter()
    native = evaluate_candidate(e, candidate)
    if "deposition" not in native:
        return dict(native, physics_measurable=False, reference_only=False)
    curve = e.curves(candidate.coefficients)["shelf"]
    result = _with_exposure(e, native, curve)
    result["elapsed_seconds"] = perf_counter() - start
    return result
