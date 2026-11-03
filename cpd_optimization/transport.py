"""Actual disconnected vacuum boundaries and conservative per-loop footprints.

This preserves the D3D phenomenological area-measure surface kernel separately
on each vacuum boundary component. It is not a cross-field plasma transport
solver and cannot redistribute heat through material or across a vacuum gap.
"""
import numpy as np
from shapely.geometry import Point
from shapely.ops import unary_union

from d3d_optimization.frozen_mapping import Hits, first_hits
from d3d_optimization.boundary_transport import deposit_boundary


def vacuum_boundaries(cavity, solids, seed):
    """Select the vacuum component containing seed, retaining every true hole."""
    if not cavity.is_valid or cavity.is_empty:
        raise ValueError('valid nonempty vessel cavity required')
    for solid in solids:
        if not solid.is_valid:
            raise ValueError('invalid material solid')
    vacuum = cavity.difference(unary_union(solids))
    parts = list(vacuum.geoms) if hasattr(vacuum, 'geoms') else [vacuum]
    selected = [p for p in parts if p.geom_type == 'Polygon' and p.contains(Point(seed))]
    if len(selected) != 1:
        raise ValueError('axis must lie in one connected vacuum region')
    region = selected[0]
    loops = [np.asarray(ring.coords, dtype=float)
             for ring in [region.exterior, *region.interiors]]
    return region, loops


def first_loop_hits(paths, loops):
    """Keep the earliest hit of all physical loops; segment IDs are loop-local."""
    count = len(paths)
    hits = Hits(np.full(count, -1, int), np.full((count, 2), np.nan), np.full(count, np.inf))
    assignment = np.full(count, -1, int)
    for k, loop in enumerate(loops):
        candidate = first_hits(paths, loop)
        earlier = candidate.distance < hits.distance
        hits.segment[earlier] = candidate.segment[earlier]
        hits.points[earlier] = candidate.points[earlier]
        hits.distance[earlier] = candidate.distance[earlier]
        assignment[earlier] = k
    return hits, assignment


def deposit_loops(field, hits, assignment, loops, weights, sigma_flux, grad_sep, step,
                  region_by_loop=None):
    """Return global power accounting without renormalizing loop subsets."""
    weights = np.asarray(weights, dtype=float)
    assignment = np.asarray(assignment, dtype=int)
    sigma_flux = np.asarray(sigma_flux, dtype=float)
    if (weights.ndim != 1 or assignment.shape != weights.shape or
            sigma_flux.shape != weights.shape or hits.segment.shape != weights.shape or
            not np.all(np.isfinite(weights)) or np.any(weights < 0) or
            not np.isfinite(weights.sum()) or weights.sum() <= 0 or
            np.any((assignment < -1) | (assignment >= len(loops))) or
            np.any((assignment < 0) != (hits.segment < 0))):
        raise ValueError('consistent hit assignments and finite nonnegative donor weights required')
    depositions = []
    for k, loop in enumerate(loops):
        ids = np.flatnonzero((assignment == k) & (weights > 0))
        if not len(ids):
            continue
        local = Hits(hits.segment[ids], hits.points[ids], hits.distance[ids])
        regions = None if region_by_loop is None else region_by_loop[k]
        dep = deposit_boundary(field, local, loop, weights[ids], sigma_flux[ids], grad_sep, step,
                               region_by_segment=regions)
        dep['loop_index'] = k
        dep['donor_ids'] = ids
        depositions.append(dep)
    prescribed = float(weights.sum())
    unhit = float(weights[assignment < 0].sum())
    deposited = float(sum(d['ledger']['deposited'] for d in depositions))
    balance = abs(deposited + unhit - prescribed) / prescribed
    if not np.isfinite(balance) or balance > 2e-12:
        raise FloatingPointError('global power ledger failed closure')
    peak = max((float(np.max(d['q'])) for d in depositions), default=0.)
    area90 = 0.
    if depositions:
        power = np.concatenate([d['cell_power'] for d in depositions])
        areas = np.concatenate([d['areas'] for d in depositions])
        density = power / areas
        order = np.argsort(-density)
        cumulative = np.cumsum(power[order])
        selected = min(int(np.searchsorted(cumulative, .9 * deposited)), len(order)-1)
        area90 = float(np.sum(areas[order[:selected]]))
        previous = cumulative[selected-1] if selected else 0.
        if power[order[selected]] > 0:
            area90 += float(areas[order[selected]] * (.9*deposited-previous)/power[order[selected]])
    return dict(depositions=depositions, peak_q=peak, area_90_m2=area90,
                effective_area_m2=deposited/peak if peak > 0 else 0.,
                ledger=dict(prescribed=prescribed, deposited=deposited, unhit=unhit,
                            unresolved=unhit, balance_error=balance),
                model='periodic area-measure footprint on each actual vacuum loop',
                missing_power_renormalized=False)
