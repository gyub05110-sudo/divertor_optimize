"""Independent CPD material contours, sampled circular arcs and material hits.

Coordinates are always (R, Z) in metres after loading. Each input contour is a
closed material outline, including its own last-to-first line. No lines are
added between components. VV is its outer polygon minus its inner polygon;
Blanket, Shield and Divertor each use their supplied outer outline as a solid.
These outlines describe axisymmetric sections, not the separately specified
three-dimensional ports. Circular faces are exact three-point circles sampled
with a configurable maximum chord; intersection results use those chords.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
import math

import numpy as np
import yaml
from shapely.geometry import LineString, Point, Polygon
from shapely.validation import explain_validity

MATERIAL_COMPONENTS = ('Blanket', 'Divertor', 'Shield', 'VV')


def _point(value):
    point = np.asarray(value, dtype=float)
    if point.shape != (2,) or not np.isfinite(point).all():
        raise ValueError('Expected a finite (R, Z) point')
    return point


def sample_three_point_arc(start, midpoint, end, max_chord_m=0.005):
    """Sample the circular sweep from start through midpoint to end.

    Both clockwise and major arcs are supported. All three supplied points
    appear exactly in the returned array. Degenerate arcs fail explicitly.
    """
    a, b, c = map(_point, (start, midpoint, end))
    if not math.isfinite(max_chord_m) or max_chord_m <= 0:
        raise ValueError('max_chord_m must be finite and positive')
    u, v = b-a, c-a
    scale = max(np.linalg.norm(u), np.linalg.norm(v))
    determinant = u[0]*v[1]-u[1]*v[0]
    if scale == 0 or abs(determinant) <= 1e-13*scale*scale:
        raise ValueError('threePointArc has collinear or degenerate points')
    center = a + np.linalg.solve(2*np.vstack([u, v]), [u @ u, v @ v])
    radius = np.linalg.norm(a-center)
    angles = np.arctan2(np.array([a,b,c])[:,1]-center[1],
                        np.array([a,b,c])[:,0]-center[0])
    ccw_end = (angles[2]-angles[0]) % (2*np.pi)
    ccw_mid = (angles[1]-angles[0]) % (2*np.pi)
    direction = 1 if ccw_mid < ccw_end else -1
    step = 2*math.asin(min(max_chord_m/(2*radius), 1.))
    parts = []
    for i, (p, q) in enumerate(((a,b),(b,c))):
        sweep = direction*((direction*(angles[i+1]-angles[i])) % (2*np.pi))
        count = max(1, math.ceil(abs(sweep)/step))
        theta = np.linspace(angles[i], angles[i]+sweep, count+1)
        part = center+radius*np.column_stack((np.cos(theta), np.sin(theta)))
        part[0], part[-1] = p, q
        parts.append(part if i == 0 else part[1:])
    return np.vstack(parts)


@dataclass(frozen=True)
class SourceNode:
    node_id: int
    point_m: tuple[float, float]
    operation: str
    midpoint_m: tuple[float, float] | None = None


@dataclass(frozen=True)
class SourceFace:
    component: str
    contour: str
    start_node_id: int
    end_node_id: int
    operation: str
    points_m: np.ndarray
    midpoint_m: tuple[float, float] | None = None
    is_closing: bool = False

    @property
    def face_id(self):
        return f'{self.component}.{self.contour}:{self.start_node_id}->{self.end_node_id}'

    @property
    def line(self):
        return LineString(self.points_m)


@dataclass(frozen=True)
class SampledContour:
    component: str
    name: str
    nodes: tuple[SourceNode, ...]
    faces: tuple[SourceFace, ...]
    points_m: np.ndarray
    # For each points_m[i] -> points_m[i+1], the original source face.
    segment_faces: tuple[SourceFace, ...]


@dataclass(frozen=True)
class MaterialHit:
    point_m: tuple[float, float]
    distance_m: float
    component: str
    face: SourceFace | None


def _sample_contour(component, name, raw, scale, max_chord_m):
    nodes = []
    for item in raw['nodes']:
        point = tuple(_point((item['x'], item['z']))*scale)
        midpoint = item.get('arc_midpoint')
        midpoint = None if midpoint is None else tuple(_point((midpoint['x'],midpoint['z']))*scale)
        nodes.append(SourceNode(int(item['node_id']), point, item['operation'], midpoint))
    if len(nodes) < 3 or nodes[0].operation != 'moveTo':
        raise ValueError(f'{component}.{name} must start with moveTo and contain at least 3 nodes')
    if len({node.node_id for node in nodes}) != len(nodes):
        raise ValueError(f'{component}.{name} contains duplicate node IDs')
    faces, points, provenance = [], [np.array(nodes[0].point_m)], []
    pairs = [(nodes[i-1], nodes[i], False) for i in range(1,len(nodes))]
    if nodes[-1].point_m != nodes[0].point_m:
        pairs.append((nodes[-1], nodes[0], True))
    for start, end, closing in pairs:
        operation = 'lineTo' if closing else end.operation
        midpoint = None if closing else end.midpoint_m
        if operation == 'threePointArc':
            if midpoint is None:
                raise ValueError(f'{component}.{name} node {end.node_id} arc has no midpoint')
            sampled = sample_three_point_arc(start.point_m, midpoint, end.point_m, max_chord_m)
        elif operation == 'lineTo':
            # Straight faces need no subdivision for exact line intersection.
            sampled = np.array([start.point_m, end.point_m])
        else:
            raise ValueError(f'{component}.{name} unsupported operation {operation}')
        if np.linalg.norm(sampled[-1]-sampled[0]) == 0:
            raise ValueError(f'{component}.{name} zero length face ending at {end.node_id}')
        face = SourceFace(component, name, start.node_id, end.node_id, operation,
                          sampled, midpoint, closing)
        faces.append(face)
        points.extend(sampled[1:])
        provenance.extend([face]*(len(sampled)-1))
    return SampledContour(component, name, tuple(nodes), tuple(faces),
                          np.asarray(points), tuple(provenance))


def _intersection_points(geometry):
    if geometry.is_empty:
        return
    if geometry.geom_type == 'Point':
        yield geometry
    elif geometry.geom_type in ('LineString', 'LinearRing'):
        yield Point(geometry.coords[0])
        yield Point(geometry.coords[-1])
    elif hasattr(geometry, 'geoms'):
        for part in geometry.geoms:
            yield from _intersection_points(part)


@dataclass
class MaterialGeometry:
    contours: dict[tuple[str, str], SampledContour]
    solids: dict[str, Polygon]
    source_data: dict
    max_chord_m: float

    @property
    def faces(self):
        return tuple(face for contour in self.contours.values() for face in contour.faces)

    def first_hit(self, start_rz, end_rz):
        """First actual material contact along a finite directed segment.

        Tangency and overlap count as contact. Starting in material returns a
        zero-distance hit with face=None unless the start is on a boundary.
        Component ordering deterministically resolves shared-boundary ties.
        """
        start, end = map(_point, (start_rz, end_rz))
        origin = Point(start)
        segment = LineString([start, end])
        candidates = []
        for name, solid in self.solids.items():
            if solid.covers(origin):
                candidates.append((0., name, origin))
            elif segment.length:
                candidates.extend((float(segment.project(point)), name, point)
                                  for point in _intersection_points(segment.intersection(solid)))
        if not candidates:
            return None
        distance, name, point = min(candidates, key=lambda item: item[0])
        component_faces = [face for face in self.faces if face.component == name]
        nearest = min(component_faces, key=lambda face: face.line.distance(point))
        face = nearest if nearest.line.distance(point) <= 1e-8 else None
        return MaterialHit(tuple(point.coords[0]), distance, name, face)

    def ray_first_hit(self, origin, direction, max_distance_m=100.):
        """First contact along a normalized ray truncated at max_distance_m."""
        direction = _point(direction)
        norm = np.linalg.norm(direction)
        if norm == 0 or not math.isfinite(max_distance_m) or max_distance_m <= 0:
            raise ValueError('Ray direction must be nonzero and max_distance_m positive')
        return self.first_hit(origin, _point(origin)+direction/norm*max_distance_m)

    def with_divertor_nodes(self, positions_m: Mapping[int, tuple], *, arc_midpoints_m=None):
        """Return independent rebuilt geometry from SI node/midpoint updates.

        Unspecified arc midpoints remain fixed. Callers deform them explicitly
        when moving curved faces. This low-level operation validates contours;
        engineering anchors, collisions and clearance are caller constraints.
        """
        data = deepcopy(self.source_data)
        scale = {'m': 1., 'mm': .001}[data['metadata']['units']]
        midpoint_updates = arc_midpoints_m or {}
        nodes = data['components']['Divertor']['contours']['outer']['nodes']
        ids = {node['node_id'] for node in nodes}
        if (set(positions_m) | set(midpoint_updates))-ids:
            raise ValueError('Unknown divertor node ID')
        for node in nodes:
            node_id = node['node_id']
            if node_id in positions_m:
                node['x'], node['z'] = map(float, _point(positions_m[node_id])/scale)
            if node_id in midpoint_updates:
                if node['operation'] != 'threePointArc':
                    raise ValueError(f'Node {node_id} is not an arc endpoint')
                x, z = map(float, _point(midpoint_updates[node_id])/scale)
                node['arc_midpoint'] = {'x': x, 'z': z}
        return build_geometry(data, max_chord_m=self.max_chord_m)


def build_geometry(data, *, max_chord_m=0.005):
    """Build independent closed material solids; reject invalid input, never repair."""
    if not math.isfinite(max_chord_m) or max_chord_m <= 0:
        raise ValueError('max_chord_m must be finite and positive')
    units = data.get('metadata', {}).get('units')
    if units not in ('mm', 'm'):
        raise ValueError(f'Unsupported or missing geometry units: {units!r}')
    scale = .001 if units == 'mm' else 1.
    contours, solids = {}, {}
    for component in MATERIAL_COMPONENTS:
        raw_contours = data['components'][component]['contours']
        for name in (('outer', 'inner') if component == 'VV' else ('outer',)):
            contour = _sample_contour(component, name, raw_contours[name], scale, max_chord_m)
            contours[(component,name)] = contour
            polygon = Polygon(contour.points_m)
            if not polygon.is_valid or polygon.area <= 0:
                raise ValueError(f'{component}.{name} invalid polygon: {explain_validity(polygon)}')
        outer = Polygon(contours[(component,'outer')].points_m)
        if component == 'VV':
            inner = Polygon(contours[(component,'inner')].points_m)
            if not outer.contains(inner):
                raise ValueError('VV inner contour must lie strictly inside its outer contour')
            outer = Polygon(outer.exterior.coords, [inner.exterior.coords])
        solids[component] = outer
    return MaterialGeometry(contours, solids, deepcopy(data), max_chord_m)


def load_geometry(path: str | Path, *, max_chord_m=0.005):
    with Path(path).open(encoding='utf-8-sig') as stream:
        return build_geometry(yaml.safe_load(stream), max_chord_m=max_chord_m)
