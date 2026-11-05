"""Prepare another operating point on the same fixed field, wall and source sampling."""
from copy import copy
import numpy as np
from .boundary_evaluator import BoundaryEvaluator
from .reference_transport import ReferencePoint


def reuse_thermal_paths(original, reference):
    if not isinstance(reference,ReferencePoint):
        raise TypeError('a validated ReferencePoint is required')
    if reference.lambda_q_m != original.reference.lambda_q_m:
        raise ValueError('source width changed; build a fresh boundary evaluator')
    result=copy(original)
    result.reference=reference
    source=original.source
    result.qu=reference.power_W*np.exp(-np.asarray(source.x)/source.width)/source.normalization
    return result


def reuse_boundary(original, reference):
    """Share geometric integrals only. Every heat state and deposition is recomputed.

    Pass the unmodified baseline returned by prepare_boundary. Construct the usual
    constrained evaluator from the returned baseline for this new operating point.
    For a new field, wall, interpolation, source, sampling or width use prepare_boundary.
    """
    if type(original) is not BoundaryEvaluator:
        raise TypeError('unmodified baseline BoundaryEvaluator required')
    thermal=reuse_thermal_paths(original.thermal_paths,reference)
    # The wrapper delegates missing attributes to base; use its owned dict rather than
    # probing uninitialized delegated attributes during copy construction.
    d=object.__new__(type(original.direct))
    d.__dict__.update(vars(original.direct))
    d.base=copy(original.direct.base)
    d.boundary_guard=copy(original.direct.boundary_guard)
    d.boundary_guard.e=d.base
    meta=dict(original.metadata,source_power_W=reference.power_W)
    return BoundaryEvaluator(d,meta,reference,original.wall_step,thermal_paths=thermal)

