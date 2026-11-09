"""Identify positive evidence for an upper LCFS saddle after lower-divertor reconstruction fails."""
import numpy as np
from scipy.optimize import least_squares
from shapely.geometry import Point, Polygon


class LowerDivertorInputError(ValueError):
    def __init__(self, original, evidence):
        if evidence is None:
            label = "lower-divertor topology unresolved"
        else:
            r, z = evidence["saddle_RZ_m"]
            label = f"upper LCFS saddle at R={r:.7f} m, Z={z:.7f} m, required lower divertor unresolved"
        super().__init__(f"{label} ({original})")
        self.evidence = evidence


def hessian(field, point):
    r, z = field.checked([point])[0]
    spline = field.spline
    cross = spline.ev(z, r, dx=1, dy=1)
    return np.array([[spline.ev(z, r, dy=2), cross], [cross, spline.ev(z, r, dx=2)]])


def _local_stationary(field, point, tolerance):
    # Resolve only a root coincident with the supplied location, inside its physical grid.
    low = np.maximum([field.r[0], field.z[0]], point - tolerance)
    high = np.minimum([field.r[-1], field.z[-1]], point + tolerance)
    scale = np.linalg.norm(hessian(field, point), ord=2)
    if scale == 0.:
        return None
    def physical(delta):
        return np.clip(point+tolerance*delta, low, high)
    fit = least_squares(lambda d: field.gradient([physical(d)])[0]/(scale*tolerance), np.zeros(2),
                        jac=lambda d: hessian(field, physical(d))/scale,
                        bounds=((low-point)/tolerance, (high-point)/tolerance),
                        ftol=1e-14, xtol=1e-14, gtol=1e-14)
    root = physical(fit.x)
    gradient = field.gradient([root])[0]
    if not fit.success or np.linalg.norm(gradient) > 1e-7:
        return None
    try:
        correction_m = np.linalg.norm(np.linalg.solve(hessian(field, root), gradient))
    except np.linalg.LinAlgError:
        return None
    if correction_m > 1e-8 or np.linalg.norm(root-point)+correction_m > tolerance:
        return None
    return root


def upper_saddle_evidence(field, lcfs, axis):
    """Conservative descriptive match only. No saddle returned here can enter a lower-OMP calculation."""
    tolerance = 1e-5  # m, same displacement scale as the loader's two-step strike comparison.
    ring = lcfs[:-1] if np.array_equal(lcfs[0], lcfs[-1]) else lcfs
    top = int(np.argmax(ring[:, 1]))
    vertex = ring[top]
    polygon = Polygon(ring)
    if vertex[1] <= axis[1] or not polygon.is_valid or not polygon.contains(Point(axis)):
        return None
    saddle = _local_stationary(field, vertex, tolerance)
    if saddle is None:
        return None
    values, basis = np.linalg.eigh(hessian(field, saddle))
    if not values[0] < 0 < values[1]:
        return None
    centre = _local_stationary(field, axis, tolerance)
    if centre is None or not polygon.contains(Point(centre)):
        return None
    axis_hessian = hessian(field, centre)
    level, axis_level = map(float, field.values([saddle, centre]))
    boundary = field.values(ring)
    low, high = float(boundary.min()), float(boundary.max())
    roundoff = 64*np.finfo(float).eps*max(1., abs(level), abs(low), abs(high))
    if (np.linalg.det(axis_hessian) <= 0 or (level-axis_level)*np.trace(axis_hessian) <= 0 or
            low-roundoff <= axis_level <= high+roundoff or not low-roundoff <= level <= high+roundoff):
        return None
    # At a nondegenerate saddle, d^T H d = 0 gives two lines. Choose the two
    # rays pointing downwards and compare with both outgoing polygon chords.
    directions = np.array([np.sqrt(-values[0])*basis[:, 1]+sign*np.sqrt(values[1])*basis[:, 0]
                           for sign in (-1, 1)])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    directions[directions[:, 1] > 0] *= -1
    if np.any(directions[:, 1] >= 0):
        return None
    neighbours = ring[[(top-1) % len(ring), (top+1) % len(ring)]]
    chords = neighbours-saddle
    lengths = np.linalg.norm(chords, axis=1)
    if np.any(chords[:, 1] >= 0) or min(lengths) <= tolerance:
        return None
    angles = np.degrees(np.arccos(np.clip((chords/lengths[:, None]) @ directions.T, -1., 1.)))
    branches = np.argmin(angles, axis=1)
    errors = angles[np.arange(2), branches]
    half_separation = .5*np.degrees(np.arccos(np.clip(abs(directions[0] @ directions[1]), 0., 1.)))
    if branches[0] == branches[1] or max(errors) >= half_separation:
        return None
    return dict(classification="upper_lcfs_saddle", saddle_RZ_m=saddle.tolist(), lcfs_top_RZ_m=vertex.tolist(),
                saddle_distance_to_lcfs_top_m=float(np.linalg.norm(saddle-vertex)),
                position_limit_m=tolerance, saddle_gradient_norm=float(np.linalg.norm(field.gradient([saddle]))),
                saddle_hessian_determinant=float(np.prod(values)), saddle_flux_Wb_per_rad=level,
                lcfs_vertex_flux_range_Wb_per_rad=[low, high], axis_RZ_m=axis.tolist(),
                elliptic_root_RZ_m=centre.tolist(), elliptic_root_flux_Wb_per_rad=axis_level,
                axis_hessian_determinant=float(np.linalg.det(axis_hessian)), neighbours_RZ_m=neighbours.tolist(),
                downward_rays_RZ=directions.tolist(), neighbour_branch_indices=branches.tolist(),
                neighbour_branch_angles_deg=errors.tolist(), angular_assignment_limit_deg=float(half_separation))
