"""Read a user equilibrium on the supplied DIII-D material wall and OMP design."""
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from shapely.geometry import Point, Polygon

from .case import separatrix
from .continuous_boundary import extended_leg
from .frozen_field import FrozenField
from .frozen_mapping import first_hits
from .independent_sources import split_design_wall


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "d3d_optimization/data/DIIID_geom.json"
POLICY = ROOT / "baseline_review/optimization_region_policy.json"
REFERENCE = ROOT / "baseline_review/baseline_review.json"
# The material and engineering definitions used to commission the OMP family.
GEOMETRY_SHA256 = "d8374eb3512e5e7b1b291cfb1dabaea6583c866e04f815ca2fe7c23a7add34b1"
REFERENCE_SHA256 = "f0827adffce3f161f709ded564f8645f8d8798f9f8fa940b0486d832b5fc436b"
OMP_SURFACES = ("omp_entry", "shelf", "omp_connector", "floor")


@dataclass
class UserInput:
    time_ms: float
    equilibrium_sha256: str
    review: dict
    geometry: dict
    provenance: dict
    _field: FrozenField
    _lcfs: np.ndarray
    _axis: np.ndarray

    def field(self, method="cubic"):
        field = self._field
        if method == "quintic" and min(len(field.r), len(field.z)) < 6:
            raise ValueError("quintic EFIT grid requires at least six points per direction")
        if method != field.method:
            field = FrozenField(field.r, field.z, field.psi, field.tf_current_ka, method)
        return field, self._lcfs, self._axis


def _hash(raw):
    return sha256(raw).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _fixed_input(path, default, description):
    path = Path(path) if path is not None else default
    raw = path.read_bytes()
    value = json.loads(raw)
    default_raw = raw if path == default else default.read_bytes()
    if default == GEOMETRY and _hash(default_raw) != GEOMETRY_SHA256:
        raise ValueError("supplied DIII-D material geometry reference checksum differs")
    if path != default and _canonical(value) != _canonical(json.loads(default_raw)):
        raise ValueError(f"unsupported {description}; this optimizer uses the supplied DIII-D reference")
    return value, str(path.resolve()), _hash(raw)


def _parquet(raw, time_ms):
    import pyarrow as pa
    import pyarrow.parquet as pq

    frame_keys = ("efit_psirz", "efit_lcfs_n", "efit_lcfs_r", "efit_lcfs_z", "efit_r_axis", "efit_z_axis")
    columns = ("source", "efit_times", "efit_grid_R", "efit_grid_Z",
               "magnetics_time", "magnetics_bcoil", *frame_keys)
    parquet = pq.ParquetFile(pa.BufferReader(raw))
    missing = [name for name in columns if name not in parquet.schema_arrow.names]
    if missing:
        raise ValueError("missing DIII-D parquet columns: " + ", ".join(missing))
    if parquet.metadata.num_rows != 1:
        raise ValueError("file must contain exactly one DIII-D discharge")
    row = parquet.read(columns=list(columns))
    if row["source"][0].as_py() != "DIII-D":
        raise ValueError("expected a DIII-D discharge")
    arrays = {}
    for name in columns[1:]:
        scalar = row[name][0]
        dtype = scalar.type
        depth = 0
        while pa.types.is_list(dtype) or pa.types.is_large_list(dtype) or pa.types.is_fixed_size_list(dtype):
            dtype = dtype.value_type
            depth += 1
        expected = {"efit_psirz": 3, "efit_lcfs_r": 2, "efit_lcfs_z": 2}.get(name, 1)
        if (not scalar.is_valid or depth != expected or
                not (pa.types.is_integer(dtype) or pa.types.is_floating(dtype))):
            raise ValueError(f"expected a numeric {expected}-dimensional array in {name}")
        arrays[name] = scalar.values
    times = np.asarray(arrays["efit_times"], dtype=float)
    if not np.isfinite(times).all():
        raise ValueError("finite EFIT timestamps required")
    matches = np.flatnonzero(times == time_ms)
    if len(matches) != 1:
        raise ValueError("requested fixed EFIT time must exist exactly once")
    i = int(matches[0])
    if any(len(arrays[name]) != len(times) for name in frame_keys):
        raise ValueError("inconsistent number of EFIT frames")
    mag_times, mag_ka = (np.asarray(arrays[name], dtype=float) for name in ("magnetics_time", "magnetics_bcoil"))
    if (not len(mag_times) or mag_times.shape != mag_ka.shape or
            not np.isfinite(mag_times).all()):
        raise ValueError("finite, aligned magnetics timestamps required")
    if np.any(np.diff(mag_times) <= 0):
        raise ValueError("magnetics time must be strictly increasing")
    if not mag_times[0] <= time_ms <= mag_times[-1]:
        raise ValueError("magnetics time must cover the requested EFIT frame")
    upper = int(np.searchsorted(mag_times, time_ms))
    lower = upper if mag_times[upper] == time_ms else upper - 1
    bracket_times = mag_times[lower:upper+1]
    bracket_ka = mag_ka[lower:upper+1]
    if not np.isfinite(bracket_ka).all():
        raise ValueError("finite magnetics current required at both samples bracketing the requested time")
    current = float(np.interp(time_ms, bracket_times, bracket_ka))
    frame = {name: arrays[name][i].as_py() for name in frame_keys}
    n = frame["efit_lcfs_n"]
    if (n is None or not np.isfinite(n) or n != int(n) or n < 3 or
            any(frame[k] is None or len(frame[k]) < n for k in ("efit_lcfs_r", "efit_lcfs_z"))):
        raise ValueError("invalid EFIT LCFS point count")
    lcfs = np.asarray([frame["efit_lcfs_r"][:int(n)], frame["efit_lcfs_z"][:int(n)]], dtype=float).T
    axis = np.array([frame["efit_r_axis"], frame["efit_z_axis"]], dtype=float)
    if not np.isfinite(axis).all():
        raise ValueError("finite magnetic axis required")
    r, z = (np.asarray(arrays[name], dtype=float) for name in ("efit_grid_R", "efit_grid_Z"))
    if min(len(r), len(z)) < 4:
        raise ValueError("cubic EFIT grid requires at least four points per direction")
    field = FrozenField(r, z, frame["efit_psirz"], current)
    metadata = dict(format="DIII-D EFIT parquet", frame_index=i, time_source="exact efit_times match",
                    tf_current_kA=current, tf_interpolation_times_ms=bracket_times.tolist(),
                    tf_interpolation_currents_kA=bracket_ka.tolist(), toroidal_field_model="ideal 144-turn coil with interpolated magnetics_bcoil (kA)")
    return field, lcfs, axis, metadata


def _strike_legs(field, xp, psi_sep, wall):
    legs = {}
    for source, name in (("OMP", "outer"), ("IMP", "inner")):
        hits = []
        for spacing in (0.0005, 0.00025):
            hit = first_hits([extended_leg(field, xp, psi_sep, source, spacing)], wall)
            if hit.segment[0] < 0:
                raise ValueError(f"{name} separatrix leg has no resolved first material hit")
            hits.append((int(hit.segment[0]), hit.points[0]))
        delta = float(np.linalg.norm(hits[0][1] - hits[1][1]))
        if hits[0][0] != hits[1][0] or delta > 1e-5:
            raise ValueError(f"{name} separatrix first hit is unresolved at the two tracing steps")
        legs[name] = dict(segment=hits[1][0], strike_RZ_m=hits[1][1].tolist(), refinement_m=delta)
    return legs


def _outer_strike_location(strike, endpoints):
    """Describe the first hit against the four fixed material portions."""
    strike = np.asarray(strike)
    matches = []
    at_end = False
    for name, segment in endpoints.items():
        a, b = np.asarray(segment)
        direction = b - a
        u = np.clip((strike - a) @ direction / (direction @ direction), 0., 1.)
        if np.linalg.norm(strike - (a + u * direction)) <= 1e-7:
            matches.append(name)
            at_end |= min(np.linalg.norm(strike - p) for p in (a, b)) <= 1e-7
    if not matches:
        return dict(location="fixed_wall", surfaces=[])
    return dict(location="adjustable_endpoint" if at_end else "adjustable_interior", surfaces=matches)


def load_user_input(equilibrium, time_ms, geometry=None, design_policy=None):
    """Read one stored EFIT frame from a single-discharge DIII-D parquet.

    The cubic lower X-point, outboard midplane crossing and two lower leg first
    hits on the supplied wall must resolve. The four OMP portions stay fixed
    if the outer leg strikes elsewhere. Shelf exposure, first-hit incidence
    and the SOL-limit witness are evaluated for each geometry by the solver.
    """
    if not np.isfinite(time_ms):
        raise ValueError("finite equilibrium time in ms required")
    equilibrium = Path(equilibrium)
    material, geom_path, geom_hash = _fixed_input(geometry, GEOMETRY, "material geometry")
    policy, policy_path, policy_hash = _fixed_input(design_policy, POLICY, "design policy")
    reference_raw = REFERENCE.read_bytes()
    if _hash(reference_raw) != REFERENCE_SHA256:
        raise ValueError("DIII-D engineering design reference checksum differs")
    reference = json.loads(reference_raw)
    wall = np.asarray(reference["wall"], dtype=float)
    limiter = np.asarray(material["limiter"], dtype=float)
    if not np.array_equal(wall, np.vstack((limiter, limiter[0]))) or policy != reference["optimization_policy"]:
        raise ValueError("design reference disagrees with supplied material geometry or design policy")
    raw = equilibrium.read_bytes()
    if raw[:4] != b"PAR1":
        raise ValueError("expected a DIII-D EFIT parquet equilibrium")
    field, lcfs, axis, metadata = _parquet(raw, time_ms)
    if (lcfs.ndim != 2 or lcfs.shape[1] != 2 or not np.isfinite(lcfs).all() or
            not Polygon(lcfs).is_valid or not Polygon(lcfs).contains(Point(axis))):
        raise ValueError("a finite, simple LCFS enclosing the magnetic axis is required")
    field.checked(lcfs)
    field.checked(wall)
    try:
        rsep, xp, psi_sep = separatrix(field, lcfs, axis)
    except ValueError as error:
        if str(error) not in ("sampling outside EFIT grid is forbidden", "lower X-point cannot be validated near LCFS",
                              "stationary point is not a magnetic saddle"):
            raise
        from .user_topology import LowerDivertorInputError, upper_saddle_evidence
        try:
            evidence = upper_saddle_evidence(field, lcfs, axis)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError, OverflowError):
            # A secondary numerical diagnosis must preserve the actual reason for rejection.
            evidence = None
        raise LowerDivertorInputError(error, evidence) from error
    # The archived review supplies material segment limits only. No older field,
    # launch point, strike location, transport result or optimized design is used.
    portions = [{k: row[k] for k in ("source", "segment", "points")}
                for row in reference["optimization_footprints"]]
    review = dict(wall=wall.tolist(), optimization_footprints=portions, optimization_policy=policy,
                  legs=_strike_legs(field, xp, psi_sep, wall), time_ms=float(time_ms),
                  xpoint_RZ_m=xp.tolist(), axis_RZ_m=axis.tolist(),
                  sources={"OMP": {"midplane_RZ_m": [float(rsep), float(axis[1])]}, "IMP": {}})
    endpoints = {s["name"]: s["endpoints_RZ_m"] for s in split_design_wall(review)[1]
                 if s["name"] in OMP_SURFACES}
    outer_location = _outer_strike_location(review["legs"]["outer"]["strike_RZ_m"], endpoints)
    review["outer_strike_design"] = outer_location
    provenance = dict(baseline_lower_strike_legs=review["legs"], baseline_outer_strike_design=outer_location,
                      equilibrium_path=str(equilibrium.resolve()), equilibrium_sha256=_hash(raw),
                      geometry_path=geom_path, geometry_sha256=geom_hash,
                      design_policy_path=policy_path, design_policy_sha256=policy_hash,
                      wall_sha256=_hash(wall.astype("<f8").tobytes()),
                      omp_endpoints_RZ_m=endpoints, omp_endpoints_sha256=_hash(_canonical(endpoints)),
                      design_reference_sha256=_hash(reference_raw), time_ms=float(time_ms),
                      length_unit="m", poloidal_flux_unit="Wb/rad", magnetic_field_unit="T", **metadata)
    return UserInput(float(time_ms), provenance["equilibrium_sha256"], review, material, provenance, field, lcfs, axis)
