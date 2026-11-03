"""Render CPD material contours and write reproducible geometry diagnostics."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cpd_optimization.geometry import load_geometry
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from shapely.geometry import Point, Polygon
from shapely.geometry.polygon import orient
from shapely.ops import nearest_points, unary_union

COLORS = {'VV': '#9dabb9', 'Shield': '#688e83', 'Blanket': '#dba448', 'Divertor': '#5d84b2'}


def material_patch(solid, color):
    polygon = orient(solid)
    rings = [polygon.exterior, *polygon.interiors]
    vertices, codes = [], []
    for ring in rings:
        coords = list(ring.coords)
        vertices.extend(coords)
        codes.extend([MplPath.MOVETO]+[MplPath.LINETO]*(len(coords)-2)+[MplPath.CLOSEPOLY])
    return PathPatch(MplPath(vertices, codes), facecolor=color, alpha=.6, edgecolor=color, lw=.7)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('input/CPD_5/CPD_5_nodes.yaml'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    geometry = load_geometry(args.source)
    figure, axes = plt.subplots(1, 2, figsize=(16, 10), layout='constrained')
    for ax in axes:
        for name in ('VV', 'Shield', 'Blanket', 'Divertor'):
            ax.add_patch(material_patch(geometry.solids[name], COLORS[name]))
        for (component, contour_name), contour in geometry.contours.items():
            if component == 'VV' and contour_name == 'outer':
                continue
            for node in contour.nodes:
                x, z = node.point_m
                ax.plot(x, z, '.', color=COLORS[component], ms=3)
                if ax is axes[0] and component == 'Divertor':
                    continue
                if ax is axes[1] and (x < 2 or x > 5.1 or z < -3.8 or z > -1.3):
                    continue
                ax.annotate(str(node.node_id), (x,z), xytext=(3,3), textcoords='offset points',
                            fontsize=7 if ax is axes[0] else 9, color=COLORS[component],
                            bbox={'facecolor':'white','alpha':.7,'edgecolor':'none','pad':.4})
        ax.set_aspect('equal')
        ax.set_xlabel('R [m]')
        ax.set_ylabel('Z [m]')
        ax.grid(alpha=.15)
    axes[0].set(xlim=(1.7,7.7), ylim=(-3.9,3.9), title='CPD_5: independent material solids')
    axes[1].set(xlim=(2,5.05), ylim=(-3.75,-1.3), title='Detached divertor: source node IDs')
    for face in geometry.contours[('Divertor','outer')].faces:
        if 6 <= face.start_node_id <= 17:
            axes[1].plot(face.points_m[:,0], face.points_m[:,1], color='#154573', lw=2)
    axes[1].annotate('20 mm open inboard slot', (2.6,-1.523), (2.1,-1.39), fontsize=9,
                     arrowprops={'arrowstyle':'->'})
    axes[1].annotate('Slanted open slot', (4.3,-2.24), (4.5,-1.8), fontsize=9,
                     arrowprops={'arrowstyle':'->'})
    figure.suptitle('CPD_5 supplied geometry | mm converted to m | no component bridges', fontsize=15)
    handles = [plt.Line2D([],[],color=color,lw=7,label=name) for name,color in COLORS.items()]
    axes[0].legend(handles=handles, loc='upper right')
    image_path = args.output/'cpd5_geometry_nodes.png'
    figure.savefig(image_path, dpi=180)
    plt.close(figure)

    gaps = {}
    for other in ('Blanket','Shield','VV'):
        a, b = nearest_points(geometry.solids['Divertor'],geometry.solids[other])
        gaps[f'Divertor--{other}'] = {'distance_m': a.distance(b),
                                    'nearest_points_m': [list(a.coords[0]),list(b.coords[0])]}
    inner = Polygon(geometry.contours[('VV','inner')].points_m)
    free = inner.difference(unary_union([geometry.solids[n] for n in ('Blanket','Divertor','Shield')]))
    fine = load_geometry(args.source, max_chord_m=.0025)
    report = {
        'source': str(args.source.resolve()),
        'source_sha256': hashlib.sha256(args.source.read_bytes()).hexdigest(),
        'units': 'm', 'arc_max_chord_m': geometry.max_chord_m,
        'solids': {name: {'valid': solid.is_valid, 'area_m2': solid.area,
                         'bounds_m': list(solid.bounds),
                         'area_change_at_half_chord_m2': abs(solid.area-fine.solids[name].area),
                         'boundary_hausdorff_at_half_chord_m': solid.boundary.hausdorff_distance(fine.solids[name].boundary)}
                   for name,solid in geometry.solids.items()},
        'divertor_gaps': gaps,
        'vacuum': {'type': free.geom_type, 'valid': free.is_valid,
                   'connected_chamber_contains_4_0': free.covers(Point(4,0)),
                   'holes': len(free.interiors) if free.geom_type == 'Polygon' else None},
        'divertor_faces': [{'face_id':face.face_id, 'operation':face.operation,
                            'sampled_length_m':face.line.length, 'closing_face':face.is_closing}
                           for face in geometry.contours[('Divertor','outer')].faces],
        'interpretation': [
            'All contours close only within their source component, using the implicit last-to-first edge.',
            'Anchors identify slanted Divertor 6->7->8, interior 8->...->17, and dome 17->18 as the likely plasma-facing chain. This is geometric identification, not a strike prediction.',
            'Divertor 1->2->3->4->5 faces the vessel; 5->6 and 18->1 are material end faces and remain collision surfaces.',
            'Inboard 20 mm slot and slanted Blanket/Divertor slot remain open; detached backside vacuum communicates with the chamber.',
            'The supplied VV.inner contour is closed. No physical vacuum closure across Blanket/Divertor endpoints is invented.',
            'Three-dimensional ports are not represented by the supplied axisymmetric contours.',
            'Anchor displacement groups define follow-motion, not an independent admissible optimization range; optimizer must enforce anchors and explicit design bounds.',
        ],
    }
    (args.output/'geometry_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'preview':str(image_path),'gaps':gaps,'vacuum':report['vacuum']},indent=2))


if __name__ == '__main__':
    main()
