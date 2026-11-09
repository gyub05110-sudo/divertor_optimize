"""Prove a complete-wall LCFS intersection from material common to every OMP candidate."""
from hashlib import sha256

import numpy as np
from shapely.geometry import LineString, Point, Polygon, mapping

from .campaign_geometry import CHORD_TOLERANCE_M, GEOMETRY_POLICY
from .independent_sources import split_design_wall
from .tangent_geometry import Curve
from .user_input import OMP_SURFACES
from .user_physics import _engineering_regions


def _samples(overlap):
    if overlap.geom_type == "Point":
        yield overlap
    elif overlap.geom_type == "LineString":
        yield Point(overlap.coords[0])
        yield overlap.interpolate(.5, normalized=True)
        yield Point(overlap.coords[-1])
    else:
        for part in overlap.geoms:
            yield from _samples(part)


def _contact(material, plasma, segment):
    if not material.intersects(plasma):
        return None
    overlap = material.intersection(plasma)
    boundary = plasma.boundary
    witness = max(_samples(overlap), key=boundary.distance)
    return dict(segment, intersection=mapping(overlap), overlap_length_m=float(overlap.length),
                witness_RZ_m=list(witness.coords[0]), witness_strictly_interior=bool(plasma.contains(witness)),
                witness_distance_to_LCFS_boundary_m=float(boundary.distance(witness)))


def fixed_wall_precheck(user):
    """Screen only immutable portions and both fixed ends of every movable face, before thermal work."""
    field, lcfs, _ = user.field("cubic")
    plasma = Polygon(field.checked(lcfs))
    xp = user.review["xpoint_RZ_m"]
    review, basis = _engineering_regions(user.review, field, float(field.values([xp])[0]))
    wall, surfaces, original = split_design_wall(review)
    field.checked(wall)
    by_index = {s["segment_index"]: s for s in surfaces}
    if [s["name"] for s in surfaces] != [*OMP_SURFACES, "imp_sloped", "imp_vertical"]:
        raise ValueError("fixed-wall check requires the supplied six-face OMP/IMP engineering layout")
    fixed_hits, endpoint_hits = [], []
    for index, (a, b) in enumerate(zip(wall[:-1], wall[1:])):
        surface = by_index.get(index)
        name = surface["name"] if surface else "fixed"
        segment = dict(original_material_segment=int(original[index]), split_segment=index, region=name,
                       segment_RZ_m=[a.tolist(), b.tolist()])
        if name in OMP_SURFACES:
            for side, end in enumerate((a, b)):
                hit = _contact(Point(end), plasma, dict(segment, endpoint=side))
                if hit is not None:
                    endpoint_hits.append(hit)
        else:
            if surface:
                # OMP candidates append eight zero IMP coefficients. At zero second
                # derivative assembled_wall uses its canonical 161-point chord mesh.
                curve = Curve(*surface["endpoints_RZ_m"], np.zeros(4))
                coords = curve.point(np.linspace(0., 1., 161))
                coords[0], coords[-1] = a, b
            else:
                coords = [a, b]
            hit = _contact(LineString(coords), plasma, segment)
            if hit is not None:
                fixed_hits.append(hit)
    excluded = bool(fixed_hits or endpoint_hits)
    return dict(excluded=excluded,
                classification="fixed_wall_lcfs_intersection" if excluded else "no_fixed_wall_lcfs_intersection",
                geometry_policy=GEOMETRY_POLICY, chord_tolerance_m=CHORD_TOLERANCE_M,
                predicate="LineString(assembled_candidate_wall).intersects(Polygon(supplied_LCFS))",
                scope="Any immutable contact excludes every OMP candidate. A clear check leaves feasibility undecided.",
                index_convention="Zero-based original and split material segment IDs and endpoint indices",
                engineering_region_basis=basis, split_segments=len(original),
                fixed_segments_tested=len(original)-len(OMP_SURFACES), adjustable_endpoints_tested=2*len(OMP_SURFACES),
                LCFS_vertices=len(lcfs), LCFS_RZ_little_endian_float64_sha256=sha256(lcfs.astype("<f8").tobytes()).hexdigest(),
                fixed_segment_intersections=fixed_hits, fixed_adjustable_endpoint_intersections=endpoint_hits)
