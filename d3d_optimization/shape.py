"""Normal-displacement shape parameterization for divertor polylines."""

from __future__ import annotations

import math
from numbers import Integral, Real
from typing import Dict, Mapping, Sequence, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]
NodeKey = Tuple[str, str, int]
NodeOffset = Tuple[float, float]


def _de_boor_value(
    knots: Sequence[float], coefficients: Sequence[float], degree: int, x: float
) -> float:
    """Evaluate one B-spline value without relying on SciPy."""

    knot_array = np.asarray(knots, dtype=np.float64)
    coefficient_array = np.asarray(coefficients, dtype=np.float64)
    count = len(coefficient_array)
    if degree < 0 or len(knot_array) != count + degree + 1:
        raise ValueError("invalid B-spline representation")
    if x == float(knot_array[count]):
        return float(coefficient_array[-1])

    span = int(np.searchsorted(knot_array, x, side="right") - 1)
    span = min(max(span, degree), count - 1)
    work = np.array(coefficient_array[span - degree : span + 1], copy=True)
    for recursion in range(1, degree + 1):
        for local in range(degree, recursion - 1, -1):
            index = span - degree + local
            denominator = knot_array[index + degree - recursion + 1] - knot_array[index]
            alpha = 0.0 if denominator == 0.0 else (x - knot_array[index]) / denominator
            work[local] = (1.0 - alpha) * work[local - 1] + alpha * work[local]
    return float(work[degree])


def _derivative_representation(
    knots: Sequence[float], coefficients: Sequence[float], degree: int
) -> Tuple[Tuple[float, ...], Tuple[float, ...], int]:
    """Return knots, coefficients, and degree for the derivative spline."""

    if degree <= 0:
        raise ValueError("cannot differentiate a degree-zero spline")
    knot_array = np.asarray(knots, dtype=np.float64)
    coefficient_array = np.asarray(coefficients, dtype=np.float64)
    denominators = (
        knot_array[degree + 1 : degree + len(coefficient_array)]
        - knot_array[1 : len(coefficient_array)]
    )
    if np.any(denominators <= 0.0):
        raise ValueError("invalid knot vector for spline differentiation")
    derived = degree * np.diff(coefficient_array) / denominators
    return (
        tuple(float(value) for value in knot_array[1:-1]),
        tuple(float(value) for value in derived),
        degree - 1,
    )


def _parameter_array(values: ArrayLike) -> Tuple[FloatArray, bool]:
    try:
        parameters = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("s must be a finite scalar or 1-D array") from exc
    scalar = parameters.ndim == 0
    if not scalar and parameters.ndim != 1:
        raise ValueError("s must be a finite scalar or 1-D array")
    if not np.all(np.isfinite(parameters)) or np.any(parameters < 0.0) or np.any(parameters > 1.0):
        raise ValueError("s must contain finite values in [0, 1]")
    return parameters, scalar


def _return_scalar_if_needed(values: FloatArray, scalar: bool) -> Union[float, FloatArray]:
    if scalar:
        return float(np.asarray(values))
    return np.asarray(values, dtype=np.float64)


class NormalDisplacementBSpline:
    """Cubic open-uniform spline with exactly the supplied active scalar DOF.

    Two fixed zero coefficients are prepended and appended.  Clamping the knot
    vector therefore enforces zero displacement and zero displacement slope at
    both endpoints without consuming active design variables.
    """

    degree = 3

    def __init__(self, coefficients: Sequence[Real]) -> None:
        try:
            active = np.asarray(coefficients, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError("coefficients must be a nonempty 1-D finite sequence") from exc
        if active.ndim != 1 or active.size < 1:
            raise ValueError("coefficients must contain at least one active value")
        if not np.all(np.isfinite(active)):
            raise ValueError("coefficients must contain only finite values")

        full = np.concatenate((np.zeros(2), active, np.zeros(2)))
        interior = np.linspace(0.0, 1.0, active.size + 2)[1:-1]
        knots = np.concatenate((np.zeros(self.degree + 1), interior, np.ones(self.degree + 1)))
        self._active = tuple(float(value) for value in active)
        self._full = tuple(float(value) for value in full)
        self._knots = tuple(float(value) for value in knots)

    @property
    def active_dof(self) -> int:
        return len(self._active)

    @property
    def coefficients(self) -> Tuple[float, ...]:
        return self._active

    @property
    def full_coefficients(self) -> FloatArray:
        return np.asarray(self._full, dtype=np.float64)

    def _evaluate(self, s: ArrayLike, derivative: int) -> Union[float, FloatArray]:
        parameters, scalar = _parameter_array(s)
        if derivative < 0 or derivative > self.degree:
            raise ValueError("unsupported spline derivative order")
        knots = self._knots
        coefficients = self._full
        degree = self.degree
        for _ in range(derivative):
            knots, coefficients, degree = _derivative_representation(knots, coefficients, degree)
        values = np.asarray(
            [
                _de_boor_value(knots, coefficients, degree, float(parameter))
                for parameter in np.atleast_1d(parameters)
            ],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("spline evaluation must be finite")
        return _return_scalar_if_needed(values, scalar)

    def displacement(self, s: ArrayLike) -> Union[float, FloatArray]:
        return self._evaluate(s, 0)

    def first_derivative(self, s: ArrayLike) -> Union[float, FloatArray]:
        return self._evaluate(s, 1)

    def second_derivative(self, s: ArrayLike) -> Union[float, FloatArray]:
        return self._evaluate(s, 2)

    def bending_second_derivative(self, s: ArrayLike) -> Union[float, FloatArray]:
        """Return ``d''(s)`` for the bending objective ``f3``."""

        return self.second_derivative(s)

    def deform(self, reference: "ReferencePolyline", s: ArrayLike) -> FloatArray:
        parameters, scalar = _parameter_array(s)
        samples = reference.sample(parameters)
        normals = reference.unit_normal(parameters)
        displacement = np.asarray(self.displacement(parameters), dtype=np.float64)
        result = samples + displacement[..., np.newaxis] * normals
        return np.asarray(result, dtype=np.float64).reshape(2) if scalar else result

    def deformed_tangent(self, reference: "ReferencePolyline", s: ArrayLike) -> FloatArray:
        parameters, scalar = _parameter_array(s)
        displacement = np.asarray(self.displacement(parameters), dtype=np.float64)
        slope = np.asarray(self.first_derivative(parameters), dtype=np.float64)
        result = (
            reference.tangent(parameters)
            + slope[..., np.newaxis] * reference.unit_normal(parameters)
            + displacement[..., np.newaxis] * reference.normal_derivative(parameters)
        )
        return np.asarray(result, dtype=np.float64).reshape(2) if scalar else result


class ReferencePolyline:
    """A finite positive-radius polyline parameterized by normalized arclength."""

    def __init__(self, points: ArrayLike) -> None:
        try:
            array = np.asarray(points, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError("base curve must be an N-by-2 numeric array") from exc
        if array.ndim != 2 or array.shape[1:] != (2,) or array.shape[0] < 2:
            raise ValueError("base curve must be an N-by-2 array with N >= 2")
        if not np.all(np.isfinite(array)):
            raise ValueError("base curve must contain only finite values")
        if np.any(array[:, 0] <= 0.0):
            raise ValueError("base curve radial coordinates R must be positive")
        segments = np.diff(array, axis=0)
        lengths = np.linalg.norm(segments, axis=1)
        if np.any(lengths <= 0.0):
            raise ValueError("base curve must not contain duplicate consecutive points")

        cumulative = np.concatenate(([0.0], np.cumsum(lengths)))
        parameters = cumulative / cumulative[-1]
        segment_units = segments / lengths[:, np.newaxis]
        node_units = np.empty_like(array)
        node_units[0] = segment_units[0]
        node_units[-1] = segment_units[-1]
        if len(array) > 2:
            combined = segment_units[:-1] + segment_units[1:]
            combined_norm = np.linalg.norm(combined, axis=1)
            if np.any(combined_norm <= np.finfo(np.float64).eps):
                raise ValueError("base curve has an undefined reversing vertex tangent")
            node_units[1:-1] = combined / combined_norm[:, np.newaxis]
        node_normals = np.column_stack((-node_units[:, 1], node_units[:, 0]))

        self._points = np.array(array, copy=True)
        self._segments = segments
        self._segment_lengths = lengths
        self._parameters = parameters
        self._node_normals = node_normals
        self._total_length = float(cumulative[-1])

    @property
    def points(self) -> FloatArray:
        return np.array(self._points, copy=True)

    @property
    def node_parameters(self) -> FloatArray:
        return np.array(self._parameters, copy=True)

    @property
    def total_length(self) -> float:
        return self._total_length

    def sample(self, s: ArrayLike) -> FloatArray:
        parameters, scalar = _parameter_array(s)
        flat = np.atleast_1d(parameters)
        result = np.column_stack(
            [np.interp(flat, self._parameters, self._points[:, axis]) for axis in range(2)]
        )
        return result[0] if scalar else result

    def _segment_indices(self, parameters: FloatArray) -> NDArray[np.int64]:
        indices = np.searchsorted(self._parameters, np.atleast_1d(parameters), side="right") - 1
        return np.clip(indices, 0, len(self._segments) - 1).astype(np.int64)

    def tangent(self, s: ArrayLike) -> FloatArray:
        parameters, scalar = _parameter_array(s)
        indices = self._segment_indices(parameters)
        result = self._segments[indices] / np.diff(self._parameters)[indices, np.newaxis]
        return result[0] if scalar else result

    def _raw_normal_and_derivative(self, parameters: FloatArray) -> Tuple[FloatArray, FloatArray]:
        flat = np.atleast_1d(parameters)
        indices = self._segment_indices(flat)
        left_s = self._parameters[indices]
        right_s = self._parameters[indices + 1]
        weight = (flat - left_s) / (right_s - left_s)
        left = self._node_normals[indices]
        right = self._node_normals[indices + 1]
        raw = left + weight[:, np.newaxis] * (right - left)
        derivative = (right - left) / (right_s - left_s)[:, np.newaxis]
        return raw, derivative

    def unit_normal(self, s: ArrayLike) -> FloatArray:
        parameters, scalar = _parameter_array(s)
        raw, _ = self._raw_normal_and_derivative(parameters)
        norms = np.linalg.norm(raw, axis=1)
        if np.any(norms <= np.finfo(np.float64).eps):
            raise ValueError("interpolated reference normal is undefined")
        result = raw / norms[:, np.newaxis]
        return result[0] if scalar else result

    def normal_derivative(self, s: ArrayLike) -> FloatArray:
        parameters, scalar = _parameter_array(s)
        raw, raw_derivative = self._raw_normal_and_derivative(parameters)
        norms = np.linalg.norm(raw, axis=1)
        unit = raw / norms[:, np.newaxis]
        projection = np.sum(unit * raw_derivative, axis=1)
        result = (raw_derivative - projection[:, np.newaxis] * unit) / norms[:, np.newaxis]
        return result[0] if scalar else result


def _integer_node_id(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError("member node IDs must be integers")
    return int(value)


def existing_node_overrides(
    component: str,
    contour: str,
    member_ids: Sequence[int],
    base_node_positions: Mapping[int, Sequence[Real]],
    coefficients: Sequence[Real],
) -> Dict[NodeKey, NodeOffset]:
    """Map a normal-displacement design onto an existing ordered node set."""

    if not isinstance(component, str) or not component:
        raise ValueError("component must be a nonempty string")
    if not isinstance(contour, str) or not contour:
        raise ValueError("contour must be a nonempty string")
    members = tuple(_integer_node_id(value) for value in member_ids)
    if len(set(members)) != len(members):
        raise ValueError("member node IDs must be unique")

    normalized_positions: Dict[int, Sequence[Real]] = {}
    for raw_node_id, position in base_node_positions.items():
        node_id = _integer_node_id(raw_node_id)
        if node_id in normalized_positions:
            raise ValueError("base node position IDs must be unique")
        normalized_positions[node_id] = position
    missing = [node_id for node_id in members if node_id not in normalized_positions]
    if missing:
        raise ValueError(f"base node positions are missing members: {missing}")

    reference = ReferencePolyline([normalized_positions[node_id] for node_id in members])
    spline = NormalDisplacementBSpline(coefficients)
    displacements = spline.displacement(reference.node_parameters)
    normals = reference.unit_normal(reference.node_parameters)
    offsets = np.asarray(displacements)[:, np.newaxis] * normals
    return {
        (component, contour, node_id): (float(offset[0]), float(offset[1]))
        for node_id, offset in zip(members, offsets)
    }
