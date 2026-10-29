"""Frozen EFIT field and prescribed, unit-power outer-branch heat source."""
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator


class FrozenField:
    """Psi is in Wb/rad, grid order is [Z,R], TF current is in kA."""

    def __init__(self, r, z, psi, tf_current_ka, method="cubic"):
        self.r, self.z, self.psi = (np.asarray(a, dtype=float) for a in (r, z, psi))
        if (self.psi.shape != (len(self.z), len(self.r)) or
                not all(np.all(np.isfinite(a)) for a in (self.r, self.z, self.psi)) or
                np.any(np.diff(self.r) <= 0) or np.any(np.diff(self.z) <= 0) or
                np.any(self.r <= 0)):
            raise ValueError("finite increasing R/Z and psi[Z,R] required")
        if method not in ("cubic", "linear", "quintic") or not np.isfinite(tf_current_ka):
            raise ValueError("valid interpolation and finite TF current required")
        self.method = method
        self.tf_current_ka = float(tf_current_ka)
        degree = 5 if method == "quintic" else 3
        self.spline = RectBivariateSpline(self.z, self.r, self.psi, kx=degree, ky=degree, s=0)
        if method == "linear":
            self.linear = RegularGridInterpolator((self.z, self.r), self.psi)

    def checked(self, points):
        p = np.asarray(points, dtype=float).reshape(-1, 2)
        if (not np.all(np.isfinite(p)) or np.any(p[:, 0] < self.r[0]) or
                np.any(p[:, 0] > self.r[-1]) or np.any(p[:, 1] < self.z[0]) or
                np.any(p[:, 1] > self.z[-1])):
            raise ValueError("sampling outside EFIT grid is forbidden")
        return p

    def values(self, points):
        p = self.checked(points)
        if self.method == "linear":
            return self.linear(p[:, ::-1])
        return self.spline.ev(p[:, 1], p[:, 0])

    def gradient(self, points):
        p = self.checked(points)
        if self.method == "linear":
            # Differentiate the same bilinear psi interpolant, rather than
            # interpolating an unrelated finite-difference gradient field.
            ir = np.clip(np.searchsorted(self.r, p[:, 0], side="right")-1, 0, len(self.r)-2)
            iz = np.clip(np.searchsorted(self.z, p[:, 1], side="right")-1, 0, len(self.z)-2)
            dr, dz = self.r[ir+1]-self.r[ir], self.z[iz+1]-self.z[iz]
            u, v = (p[:, 0]-self.r[ir])/dr, (p[:, 1]-self.z[iz])/dz
            f00, f10 = self.psi[iz, ir], self.psi[iz, ir+1]
            f01, f11 = self.psi[iz+1, ir], self.psi[iz+1, ir+1]
            return np.column_stack([((1-v)*(f10-f00)+v*(f11-f01))/dr,
                                    ((1-u)*(f01-f00)+u*(f11-f10))/dz])
        return np.column_stack([self.spline.ev(p[:, 1], p[:, 0], dy=1),
                                self.spline.ev(p[:, 1], p[:, 0], dx=1)])

    def magnetic(self, points):
        p = self.checked(points)
        g = self.gradient(p)
        return np.column_stack([-g[:, 1] / p[:, 0], g[:, 0] / p[:, 0],
                                2e-7 * 144 * self.tf_current_ka * 1000 / p[:, 0]])

    def incidence(self, points, normals):
        b = self.magnetic(points)
        n = np.asarray(normals, dtype=float)
        if n.shape != b[:, :2].shape or np.any(np.linalg.norm(n, axis=1) == 0):
            raise ValueError("finite nonzero R/Z normals required")
        n = n / np.linalg.norm(n, axis=1)[:, None]
        ratio = np.abs(np.sum(b[:, :2] * n, axis=1)) / np.linalg.norm(b, axis=1)
        return np.degrees(np.arcsin(np.clip(ratio, 0, 1)))


@dataclass
class HeatSource:
    x: np.ndarray
    weights: np.ndarray
    h: np.ndarray
    edges: np.ndarray
    width: float
    normalization: float


def exponential_source(field, r_sep, z_up, edges, width):
    """Integrate q_parallel/B * 2*pi*|dpsi| per tube, then normalize once.

    The finite prescribed source interval is common to every geometry. Missing
    surface hits do not change this normalization.
    """
    edges = np.asarray(edges, dtype=float)
    if width <= 0 or not np.isfinite(width) or edges[0] < 0 or np.any(np.diff(edges) <= 0):
        raise ValueError("positive width and increasing nonnegative source edges required")
    nodes, weights = np.polynomial.legendre.leggauss(4)
    centers = (edges[:-1] + edges[1:]) / 2
    half = np.diff(edges) / 2
    x = centers[:, None] + half[:, None] * nodes
    p = np.column_stack([r_sep + x.ravel(), np.full(x.size, z_up)])
    b = field.magnetic(p)
    dr = field.gradient(p)[:, 0]
    density = np.exp(-x.ravel()/width) * 2*np.pi*np.abs(dr)/np.linalg.norm(b, axis=1)
    intervals = half * (density.reshape(x.shape) @ weights)
    norm = float(np.sum(intervals))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("source must have positive finite integrated power")
    cp = np.column_stack([r_sep+centers, np.full(len(centers), z_up)])
    h = np.exp(-centers/width)/np.linalg.norm(field.magnetic(cp), axis=1)/norm
    return HeatSource(centers, intervals/norm, h, edges, width, norm)
