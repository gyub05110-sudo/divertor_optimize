"""First-surface-hit accounting for prescribed direct magnetic heat flux."""
from dataclasses import dataclass

import numpy as np
import shapely
from shapely.geometry import LineString

from .shape import NormalDisplacementBSpline, ReferencePolyline

ANGLE_SCOPE = "positive_power_first_hits_on_explicit_target"


def target_angle_status(angles, target_hit, source_power=None, floor=1.5):
    """One angle policy: explicit target AND positive first-hit source power.

    Use source power, not q or an arbitrary heat cutoff: a grazing, powered
    target hit must still fail. Other surfaces and zero-power hits never enter
    the minimum, including nonfinite angles on excluded surfaces.
    """
    angles = np.asarray(angles, dtype=float)
    target_hit = np.asarray(target_hit, dtype=bool)
    power = np.ones_like(angles) if source_power is None else np.asarray(source_power, dtype=float)
    if (angles.ndim != 1 or target_hit.shape != angles.shape or power.shape != angles.shape or
            np.any(power < 0) or not np.all(np.isfinite(power)) or not np.isfinite(floor) or floor < 0):
        raise ValueError("matching angle/mask/nonnegative power arrays and finite floor required")
    selected = target_hit & (power > 0)
    summary = {"scope": ANGLE_SCOPE, "target_hit_count": int(np.count_nonzero(target_hit)),
               "applied_hit_count": int(np.count_nonzero(selected)),
               "excluded_zero_power_target_hit_count": int(np.count_nonzero(target_hit & (power == 0)))}
    if not np.any(selected):
        return dict(summary, minimum_angle_deg=None, angle_constraint_passed=None,
                    status="target_unilluminated")
    if not np.all(np.isfinite(angles[selected])):
        return dict(summary, minimum_angle_deg=None, angle_constraint_passed=False,
                    status="invalid_target_angle")
    minimum = float(np.min(angles[selected]))
    passed = bool(minimum >= floor)
    return dict(summary, minimum_angle_deg=minimum, angle_constraint_passed=passed,
                status="angle_satisfied" if passed else "target_angle_violation")


@dataclass
class Hits:
    segment: np.ndarray
    points: np.ndarray
    distance: np.ndarray


def _intersection_points(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Point":
        return [np.asarray(geometry.coords[0])]
    if geometry.geom_type == "LineString":
        return [np.asarray(geometry.coords[0]), np.asarray(geometry.coords[-1])]
    return [p for part in geometry.geoms for p in _intersection_points(part)]


def first_hits(paths, wall):
    """Intersect complete paths against the ordered wall, keeping first hits."""
    wall = np.asarray(wall, dtype=float)
    if wall.ndim != 2 or wall.shape[1] != 2 or len(wall) < 2:
        raise ValueError("wall must contain at least two R/Z points")
    lines = [p if isinstance(p, LineString) else LineString(p) for p in paths]
    intersections = shapely.intersection(np.asarray(lines, dtype=object), LineString(wall))
    segments = np.full(len(lines), -1, dtype=int)
    points = np.full((len(lines), 2), np.nan)
    distances = np.full(len(lines), np.inf)
    direction = np.diff(wall, axis=0)
    lengths2 = np.sum(direction**2, axis=1)
    if np.any(lengths2 <= 0):
        raise ValueError("duplicate wall vertices forbidden")
    for i, (line, intersection) in enumerate(zip(lines, intersections)):
        candidates = _intersection_points(intersection)
        if not candidates:
            continue
        ds = [line.project(shapely.Point(p)) for p in candidates]
        k = int(np.argmin(ds))
        point = candidates[k]
        projection = np.clip(np.sum((point-wall[:-1])*direction, axis=1)/lengths2, 0, 1)
        residual = np.linalg.norm(wall[:-1]+projection[:, None]*direction-point, axis=1)
        segments[i] = int(np.argmin(residual))
        points[i], distances[i] = point, ds[k]
    return Hits(segments, points, distances)


@dataclass
class HeatMap:
    points: np.ndarray
    segment: np.ndarray
    q: np.ndarray
    angles: np.ndarray
    segment_power: np.ndarray
    areas: np.ndarray
    mapped_power: float
    unmapped_power: float
    balance_error: float
    hit_power: np.ndarray

    def angle_mask(self, target_segments, source_power=None):
        targets = np.asarray(target_segments)
        if targets.ndim != 1 or not len(targets) or np.any(targets < 0) or not np.all(targets == targets.astype(int)):
            raise ValueError("explicit nonnegative target segment indices required")
        power = self.hit_power if source_power is None else np.asarray(source_power, dtype=float)
        if power.shape != self.segment.shape or not np.all(np.isfinite(power)) or np.any(power < 0):
            raise ValueError("finite nonnegative first-hit source power must match segments")
        return np.isin(self.segment, targets) & (power > 0)

    def angle_application(self, floor, target_segments, source_power=None):
        # Validate exactly the mask used by the audit and optimization path.
        self.angle_mask(target_segments, source_power)
        power = self.hit_power if source_power is None else source_power
        return target_angle_status(self.angles, np.isin(self.segment, target_segments), power, floor)

    def feasible_angle(self, floor, target_segments):
        return self.angle_application(floor, target_segments)["angle_constraint_passed"] is True


def map_heat(field, hits, weights, h, wall):
    wall = np.asarray(wall, dtype=float)
    weights, h = np.asarray(weights, dtype=float), np.asarray(h, dtype=float)
    if (weights.shape != hits.segment.shape or h.shape != weights.shape or
            np.any(weights < 0) or not np.all(np.isfinite(weights)) or
            np.any(h < 0) or not np.all(np.isfinite(h))):
        raise ValueError("finite nonnegative source arrays must match hit count")
    valid = hits.segment >= 0
    direction = np.diff(wall, axis=0)
    lengths = np.linalg.norm(direction, axis=1)
    areas = 2*np.pi*(wall[:-1, 0]+wall[1:, 0])/2*lengths
    if np.any(areas <= 0):
        raise ValueError("positive finite surface areas required")
    q, angles = np.zeros(len(weights)), np.full(len(weights), np.nan)
    normals = np.column_stack([-direction[:, 1], direction[:, 0]])/lengths[:, None]
    if np.any(valid):
        b = field.magnetic(hits.points[valid])
        q[valid] = h[valid]*np.abs(np.sum(b[:, :2]*normals[hits.segment[valid]], axis=1))
        angles[valid] = field.incidence(hits.points[valid], normals[hits.segment[valid]])
    power = np.bincount(hits.segment[valid], weights=weights[valid], minlength=len(areas))
    mapped, unmapped = float(np.sum(power)), float(np.sum(weights[~valid]))
    total = float(np.sum(weights))
    if total <= 0:
        raise ValueError("positive source power required")
    return HeatMap(hits.points, hits.segment, q, angles, power, areas, mapped, unmapped,
                   abs(mapped+unmapped-total)/total, np.where(valid, weights, 0.0))


def deform_target(target, coefficients):
    reference = ReferencePolyline(target)
    return NormalDisplacementBSpline(coefficients).deform(reference, reference.node_parameters)


def trace_paths(field, starts, step=0.002, max_length=6.0, extend_below=-1.50):
    """Vectorized RK4 poloidal arclength paths, directed toward lower divertor.

    Paths continue beyond the material wall to permit candidate first-hit tests.
    No field sample is taken outside the rectangular EFIT grid.
    """
    starts = np.asarray(starts, dtype=float)
    field.checked(starts)
    count = len(starts)
    pos = starts.copy()
    active = np.ones(count, dtype=bool)
    b0 = field.magnetic(starts)[:, :2]
    sign = np.where(b0[:, 1] < 0, 1.0, -1.0)
    history = [[p.copy()] for p in starts]

    def velocity(p, sg):
        g = field.gradient(p)
        v = np.column_stack([-g[:, 1], g[:, 0]])
        norm = np.linalg.norm(v, axis=1)
        return sg[:, None]*v/np.maximum(norm[:, None], 1e-15)

    margin = 2*step
    for _ in range(int(max_length/step)):
        active &= ((pos[:, 0] > field.r[0]+margin) & (pos[:, 0] < field.r[-1]-margin) &
                   (pos[:, 1] > max(field.z[0]+margin, extend_below)) &
                   (pos[:, 1] < field.z[-1]-margin))
        ids = np.flatnonzero(active)
        if len(ids) == 0:
            break
        p, sg = pos[ids], sign[ids]
        k1 = velocity(p, sg)
        k2 = velocity(p+step*k1/2, sg)
        k3 = velocity(p+step*k2/2, sg)
        k4 = velocity(p+step*k3, sg)
        pos[ids] = p+step*(k1+2*k2+2*k3+k4)/6
        for j in ids:
            history[j].append(pos[j].copy())
    return [np.asarray(path) for path in history]
