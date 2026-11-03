"""Fixed-endpoint CPD divertor face deformation, independent of transport.

Each selected source line or circular arc has its own normal-displacement
cubic spline. Original node corners remain corners: this is not a globally
C1 family. Collision checks use sampled curves; curvature extrema are refined
numerically within each smooth face and require resolution qualification.
"""
from copy import deepcopy
import math

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize_scalar
from shapely.geometry import Polygon
from shapely.validation import explain_validity
from .geometry import sample_three_point_arc


KNOTS = np.r_[np.zeros(4), np.arange(1, 5)/5, np.ones(4)]
BASIS = BSpline(KNOTS, np.eye(8), 3)


class FaceCurve:
    """Analytic reference line/circle plus an exact spline normal offset."""
    def __init__(self, face):
        self.face = face
        self.start = np.array(face.points_m[0], float, copy=True)
        self.end = np.array(face.points_m[-1], float, copy=True)
        self.arc = face.operation == 'threePointArc'
        if self.arc:
            mid = np.asarray(face.midpoint_m, float)
            u, v = mid-self.start, self.end-self.start
            det = np.linalg.det(np.vstack([u, v]))
            if abs(det) <= 1e-13*max(np.linalg.norm(u), np.linalg.norm(v))**2:
                raise ValueError('degenerate reference circle')
            self.center = self.start+np.linalg.solve(2*np.vstack([u, v]), [u@u, v@v])
            self.radius = float(np.linalg.norm(self.start-self.center))
            theta = np.arctan2(np.array([self.start, mid, self.end])[:, 1]-self.center[1],
                               np.array([self.start, mid, self.end])[:, 0]-self.center[0])
            self.theta0 = float(theta[0])
            self.direction = 1 if (theta[1]-theta[0]) % (2*np.pi) < (theta[2]-theta[0]) % (2*np.pi) else -1
            self.sweep = self.direction*((self.direction*(theta[2]-theta[0])) % (2*np.pi))
            self.mid_t = float(((self.direction*(theta[1]-theta[0])) % (2*np.pi))/abs(self.sweep))
            self.reference_length_m = self.radius*abs(self.sweep)
        else:
            self.vector = self.end-self.start
            self.reference_length_m = float(np.linalg.norm(self.vector))
        if self.reference_length_m <= 0:
            raise ValueError('positive reference face length required')
        self.rules = []
        for order in (16, 32):
            x, w = np.polynomial.legendre.leggauss(order)
            t = np.concatenate([left+.1*(x+1) for left in np.arange(5)/5])
            self.rules.append((t, np.tile(.1*w, 5)))

    def evaluate(self, t, coefficients):
        """Return positions and first/second derivatives with respect to t."""
        t = np.asarray(t, float).reshape(-1)
        c = np.asarray(coefficients, float)
        if c.shape != (4,) or np.any(~np.isfinite(c)) or np.any(~np.isfinite(t)) or np.any((t < 0) | (t > 1)):
            raise ValueError('four finite coefficients and parameters in [0,1] required')
        full = np.r_[0., 0., c, 0., 0.]
        d = BASIS(t)@full; d1 = BASIS(t, 1)@full; d2 = BASIS(t, 2)@full
        if self.arc:
            theta = self.theta0+self.sweep*t
            radial = np.column_stack([np.cos(theta), np.sin(theta)])
            angular = np.column_stack([-np.sin(theta), np.cos(theta)])
            r = self.center+self.radius*radial
            r1 = self.radius*self.sweep*angular
            r2 = -self.radius*self.sweep**2*radial
            n = -self.direction*radial
            n1 = -self.direction*self.sweep*angular
            n2 = self.direction*self.sweep**2*radial
        else:
            r = self.start+t[:, None]*self.vector
            r1 = np.broadcast_to(self.vector, r.shape)
            r2 = np.zeros_like(r)
            n = np.broadcast_to(np.array([-self.vector[1], self.vector[0]])/self.reference_length_m, r.shape)
            n1 = n2 = np.zeros_like(r)
        p = r+d[:, None]*n
        # Preserve engineering endpoints bit-for-bit.
        p[t == 0] = self.start; p[t == 1] = self.end
        return p, r1+d1[:, None]*n+d[:, None]*n1, r2+d2[:, None]*n+2*d1[:, None]*n1+d[:, None]*n2

    def curvature(self, t, coefficients):
        _, v, a = self.evaluate(t, coefficients)
        speed = np.linalg.norm(v, axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            result = abs(v[:, 0]*a[:, 1]-v[:, 1]*a[:, 0])/speed**3
        return np.where(speed > 1e-12*self.reference_length_m, result, np.inf)

    def metrics(self, coefficients):
        values = [float(w@np.linalg.norm(self.evaluate(t, coefficients)[1], axis=1)) for t, w in self.rules]
        error = abs(values[1]-values[0])
        # Include every knot and both endpoints, then refine all sampled maxima.
        t = np.linspace(0, 1, 161)
        kappa = self.curvature(t, coefficients)
        maximum = float(np.max(kappa))
        if np.isfinite(maximum) and maximum > 0:
            ids = np.flatnonzero((kappa[1:-1] >= kappa[:-2]) & (kappa[1:-1] >= kappa[2:]) &
                                ((kappa[1:-1] > kappa[:-2]) | (kappa[1:-1] > kappa[2:])))+1
            # A reference circle has constant curvature up to roundoff.
            if np.ptp(kappa) > 1e-10*maximum:
                for i in ids:
                    opt = minimize_scalar(lambda x: -float(self.curvature([x], coefficients)[0]),
                                          bounds=(t[i-1], t[i+1]), method='bounded', options={'xatol': 1e-11})
                    maximum = max(maximum, -float(opt.fun))
        return dict(reference_length_m=self.reference_length_m, length_m=values[1],
                    length_error_m=error, length_ratio=(max(values)+error)/self.reference_length_m,
                    min_curvature_radius_m=(1/maximum if maximum > 0 else math.inf),
                    curvature_method='161 samples plus bounded local maxima refinement')

    def sample(self, coefficients, max_chord_m):
        c = np.asarray(coefficients, float)
        if not np.any(c):
            if self.arc and np.max(np.linalg.norm(np.diff(self.face.points_m, axis=0), axis=1)) > max_chord_m*(1+1e-12):
                return sample_three_point_arc(self.start, self.face.midpoint_m, self.end, max_chord_m)
            return np.array(self.face.points_m, copy=True)
        spline = BSpline(KNOTS, np.r_[0., 0., c, 0., 0.], 3)
        # Convex-hull bounds on d and d' bound speed everywhere, giving a
        # conservative chord bound without missing an excursion between nodes.
        speed_bound = self.reference_length_m+np.max(abs(spline.derivative().c))
        if self.arc:
            speed_bound += abs(self.sweep)*np.max(abs(c))
        count = max(2, int(math.ceil(speed_bound/max_chord_m)))
        t = np.unique(np.r_[np.linspace(0, 1, count+1), KNOTS, self.mid_t if self.arc else []])
        points=self.evaluate(t,c)[0]
        # Distinct floating-point parameters (e.g. a knot and a linspace sample)
        # can round to the identical R/Z point. Removing only exact adjacent
        # duplicates changes no curve, length, material boundary or input.
        keep=np.r_[True,np.any(points[1:]!=points[:-1],axis=1)]
        return points[keep]


class DivertorDeformation:
    """40-DOF family on faces 8->9 through 17->18; all other faces fixed.

    Dome face 17->18 remains a target: the source anchors constrain its nodes,
    not every interior arc point. Its node positions/tangents stay fixed while
    its interior may deform. The narrow coupled slot faces 6->7->8 stay fixed.
    """
    def __init__(self, geometry, *, max_chord_m=.003, coefficient_bound_m=.1,
                 max_length_ratio=1.2, min_curvature_radius_m=.020, lcfs=None):
        if not all(np.isfinite(x) and x > 0 for x in
                   (max_chord_m, coefficient_bound_m, max_length_ratio, min_curvature_radius_m)):
            raise ValueError('positive finite deformation settings required')
        self.geometry = deepcopy(geometry)
        self.max_chord_m = float(max_chord_m)
        self.bound = float(coefficient_bound_m)
        self.max_length_ratio = float(max_length_ratio)
        self.min_curvature_radius_m = float(min_curvature_radius_m)
        self.lcfs = lcfs
        self.faces = self.geometry.contours['Divertor', 'outer'].faces
        selected = [f for f in self.faces if 8 <= f.start_node_id <= 17 and f.end_node_id == f.start_node_id+1]
        if len(selected) != 10:
            raise ValueError('expected ten CPD faces 8->9 through 17->18')
        self.slices = {f.face_id: slice(4*i, 4*i+4) for i, f in enumerate(selected)}
        self.curves = {f.face_id: FaceCurve(f) for f in selected}
        self.dimension = 4*len(selected)
        self.cavity = Polygon(self.geometry.contours['VV', 'inner'].points_m)
        self.basis_description = dict(type='cubic B-spline normal displacement per original face',
            knots=KNOTS.tolist(), active_coefficients_per_face=4, fixed_coefficients=[0, 1, 6, 7],
            coefficient_bound_m=self.bound, original_endpoints_and_tangents_fixed=True,
            globally_C1=False, fixed_faces=[f.face_id for f in self.faces if f.face_id not in self.slices],
            selected_faces=list(self.slices), max_length_ratio=self.max_length_ratio,
            dome_policy='17->18 interior movable; anchor node positions and tangents fixed',
            min_curvature_radius_m=self.min_curvature_radius_m)

    def deform(self, coefficients, *, lcfs=None, max_chord_m=None):
        c = np.asarray(coefficients, float)
        if c.shape != (self.dimension,) or np.any(~np.isfinite(c)) or np.any(abs(c) > self.bound+1e-12):
            raise ValueError('finite bounded coefficient vector required')
        step = self.max_chord_m if max_chord_m is None else float(max_chord_m)
        if not np.isfinite(step) or step <= 0:
            raise ValueError('positive finite max_chord_m required')
        points, provenance, node_ids, face_points, metrics, errors = [], [], [], {}, {}, []
        for face in self.faces:
            name = face.face_id
            if name in self.slices:
                curve = self.curves[name]; cc = c[self.slices[name]]
                part = curve.sample(cc, step)
                metrics[name] = curve.metrics(cc)
                if metrics[name]['length_ratio'] > self.max_length_ratio+1e-12:
                    errors.append('material_length_exceeds_limit:'+name)
                if metrics[name]['min_curvature_radius_m'] < self.min_curvature_radius_m-1e-12:
                    errors.append('curvature_radius_below_limit:'+name)
            else:
                part = np.array(face.points_m, copy=True)
            face_points[name] = part
            points.extend(part if not points else part[1:])
            provenance.extend([name]*(len(part)-1))
            node_ids.extend([face.start_node_id]*(len(part)-1))
        wall = np.asarray(points)
        solid = Polygon(wall)
        if np.any(wall[:, 0] <= 0):
            errors.append('nonpositive_radius')
        if not solid.is_valid or solid.area <= 0:
            errors.append('invalid_polygon:'+explain_validity(solid))
        else:
            if not self.cavity.covers(solid):
                errors.append('outside_vessel_cavity')
            for name in ('Blanket', 'Shield'):
                if solid.intersects(self.geometry.solids[name]):
                    errors.append('collision_'+name)
            plasma = self.lcfs if lcfs is None else lcfs
            if plasma is not None:
                plasma = plasma if hasattr(plasma, 'geom_type') else Polygon(plasma)
                if not plasma.is_valid:
                    raise ValueError('valid LCFS polygon required')
                if solid.intersects(plasma):
                    errors.append('lcfs_intrusion')
        return dict(wall=wall, solid=solid, face_points=face_points, segment_faces=tuple(provenance),
                    segment_face_ids=np.asarray(node_ids, int), feasible=not errors,
                    metrics=metrics, length_ratios={k:v['length_ratio'] for k,v in metrics.items()},
                    min_curvature_radius_m=min(v['min_curvature_radius_m'] for v in metrics.values()),
                    errors=errors, geometry_valid=not errors, max_chord_m=step,
                    max_coefficient_m=float(np.max(abs(c))),
                    basis_description=self.basis_description)

    def evaluate(self, coefficients, **kwargs):
        """Evaluator spelling used by the transport and search wrappers."""
        return self.deform(coefficients, **kwargs)
