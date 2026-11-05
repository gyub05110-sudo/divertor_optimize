"""Attribute OMP 640 cubic/quintic differences on a fixed physical wall.

Run from the repo root: python -m verification.nf_local_attribution --help.
Uses archived shot splines, source quadrature and deposits; no path tracing.
"""
import argparse
import csv
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq
from shapely.geometry import LineString

from d3d_optimization.boundary_transport import AreaBoundary, wrapped_density, wrapped_mass
from d3d_optimization.campaign_geometry import Candidate
from d3d_optimization.case import load_case, separatrix
from d3d_optimization.independent_sources import midplane_source
from verification.preimpact_routing import Wall

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'viewers/d3d-progression-20261104'))
from distribution_policy_v3_final import distribution_compare, from_npz

METHODS = ('cubic', 'quintic')
FACTORS = ('center', 'physical_width', 'weight')
SHOT = ROOT / 'input/D3D/d3d/d3d_shot_1a3155370c.parquet'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def shapley(cube):
    """Average signed increments over all factor orderings, donor by donor."""
    count = len(next(iter(cube)))
    paths = list(permutations(range(count)))
    effects = np.zeros((count, *cube[(0,) * count].shape))
    for order in paths:
        bits = [0] * count
        for index in order:
            before = cube[tuple(bits)]
            bits[index] = 1
            effects[index] += (cube[tuple(bits)] - before) / len(paths)
    np.testing.assert_allclose(effects.sum(axis=0), cube[(1,) * count] - cube[(0,) * count],
                               rtol=2e-12, atol=2e-15)
    return effects


def source_inputs(reference, count, review):
    records = []
    for method in METHODS:
        field, lcfs, axis = load_case(SHOT, review['time_ms'], method)
        rsep, xpoint, psi_sep = separatrix(field, lcfs, axis)
        line = LineString([[field.r[0], axis[1]], [field.r[-1], axis[1]]])
        crossing = line.intersection(LineString(review['wall']))
        points = crossing.geoms if hasattr(crossing, 'geoms') else [crossing]
        limit = min(.1, min(point.x - rsep for point in points if point.x > rsep))
        edges = limit * np.linspace(0., 1., count + 1) ** 2
        source = midplane_source(field, rsep, axis[1], edges, reference['lambda_q_m'], 'OMP')
        starts = np.column_stack([rsep + source.x, np.full(count, axis[1])])
        records.append(dict(field=field, source=source, starts=starts, rsep=rsep,
                            xpoint=xpoint, psi_sep=psi_sep))
    return records


def replay(z, source):
    wall = z['wall']
    period = z['area_edges'][-1]
    mesh = AreaBoundary(wall, max(np.linalg.norm(np.diff(wall, axis=0), axis=1)) * 2)
    segments = z['hit_segment']
    if np.any(segments < 0):
        raise ValueError('assumes complete first hits')
    centers = mesh.hit_area(z['hit_points'], segments)
    np.testing.assert_allclose(centers, z['donor_area_m2'], atol=2e-12, rtol=0)
    width = 2 * np.pi * z['hit_points'][:, 0] * z['sigma_physical_m']
    np.testing.assert_allclose(width, z['sigma_area_m2'], rtol=1e-14)
    np.testing.assert_allclose(source.x, z['source_x_m'], rtol=0, atol=1e-16)
    weights = source.weights
    collapsed = z['donor_region_power'].sum(axis=1)
    np.testing.assert_allclose(weights, collapsed, atol=2e-18, rtol=2e-13)
    edges = z['area_edges']
    power = np.zeros(len(edges) - 1)
    for c, sigma, weight in zip(centers, width, weights):
        power += weight * wrapped_mass(edges[:-1], edges[1:], c, sigma, period)
    scale = float(z['assigned_power_W']) / 1e6
    err = power - z['cell_power']
    max_q = np.max(abs(err) * scale / np.diff(edges))
    if max_q > 1e-9:
        raise ValueError(f'failed native cell replay: {max_q} MW/m²')
    return dict(cells=len(power), max_cell_power_error_fraction=max(abs(err)),
                max_cell_q_error_MW_m2=max_q, max_weight_error=max(abs(weights-collapsed)),
                source_weight_sum=weights.sum(), total_cell_power_fraction=power.sum())


def wall_location(profile, edges, u):
    cell = int(np.clip(np.searchsorted(profile['u'], u, side='right') - 1, 0, len(edges)-2))
    t = (u - profile['u'][cell]) / np.diff(profile['u'])[cell]
    a, b = profile['wall'][cell:cell+2]
    ds = np.linalg.norm(b-a)
    rz = a + t*(b-a)
    area = edges[cell] + np.pi * ds * t * (a[0] + rz[0])
    return dict(s_m=u*profile['length'], RZ_m=rz.tolist(), area_m2=area, cell=cell)


def cube_at(zs, inputs, pos, cell, edges):
    """Each corner returns [point density, cell average] for each donor."""
    cube = {}
    for ci, si in product(range(2), repeat=2):
        widths = 2*np.pi*zs[ci]['hit_points'][:, 0]*zs[si]['sigma_physical_m']
        centers = zs[ci]['donor_area_m2']
        density = []
        for mu, width in zip(centers, widths):
            mass = wrapped_mass(edges[cell], edges[cell+1], mu, width, edges[-1])
            density.append((wrapped_density(pos, mu, width, edges[-1]), mass/(edges[cell+1]-edges[cell])))
        for wi, inp in enumerate(inputs):
            weights = inp['source'].weights * float(zs[wi]['assigned_power_W'])/1e6
            cube[(ci, si, wi)] = weights[:, None]*np.asarray(density)
    return cube


def flux_root_check(k, inputs, walls, hit_s, zs):
    """Check the same observed wall branch independently using constant psi."""
    roots = {}
    for fi, wall in enumerate(walls):
        segment = int(wall.arc_id(hit_s[fi][k]))
        a, b = wall.turn[segment:segment+2]
        for xi, inp in enumerate(inputs):
            psi = wall.flux.values(inp['starts'][k:k+1])[0]
            roots[(fi, xi)] = np.array(brentq(lambda s: wall.flux.values(wall.points_at([s]))[0] - psi,
                                               a, b, xtol=2e-14))
    movement = shapley(roots)
    slope = [float(wall.flux.gradient(z['hit_points'][k:k+1])[0] @ wall.t[z['hit_segment'][k]])
             for z, wall in zip(zs, walls)]
    return dict(donor=k, psi_derivative_along_wall=slope, cubic_hit_vs_flux_root_um=1e6*(hit_s[0][k]-roots[(0, 0)]),
                quintic_hit_vs_flux_root_um=1e6*(hit_s[1][k]-roots[(1, 1)]),
                predicted_hit_shift_um=1e6*(roots[(1, 1)]-roots[(0, 0)]),
                interpolant_at_fixed_launch_um=1e6*movement[0], launch_RZ_shift_um=1e6*movement[1])


def archived_identity(grid_dir, stage, paths):
    design = read(grid_dir/'design.json')
    grid = read(grid_dir/'OMP_chi0.3/result.json')
    selected, = grid['selected']
    finite = Candidate('OMP', design['radius_mm'], tuple(design['families']['OMP'][selected]))
    assert finite.identity == grid['fine'][selected]['id']
    validation = read(stage/'validation.json')['payload']
    old = Candidate.from_record(validation['candidate'])
    assert old.identity == validation['identity'] == stage.name
    assert validation['reference'] == grid['reference']
    assert read(grid_dir/'OMP_chi0.6/result.json')['selected'] == [selected]
    diffs = {}
    for method, p in zip(METHODS, paths):
        archive = stage/f'{method}_640.npz'
        with np.load(archive) as a, np.load(p) as b:
            np.testing.assert_array_equal(a['hit_segment'], b['hit_segment'])
            diffs[method] = dict(archive=str(archive), sha256=digest(archive),
                max_wall_vertex_distance_m=max(np.linalg.norm(a['wall']-b['wall'], axis=1)),
                max_cell_power_diff_fraction=max(abs(a['cell_power']-b['cell_power'])),
                max_saved_q_diff_MW_m2=max(abs(a['q_MW_m2']-b['q_MW_m2'])))
    return dict(finite_grid_id=finite.identity, historical_stage18_id=old.identity,
                vector_max_abs_difference=max(abs(np.array(finite.vector)-old.vector)),
                reference=grid['reference'], comparisons=diffs)


def pair(paths, expected, inputs, key):
    with np.load(paths[0]) as f, np.load(paths[1]) as g:
        zs = [dict(f), dict(g)]
    for name in ('wall', 'receiving_wall', 'area_edges', 'cell_parent', 'assigned_power_W'):
        np.testing.assert_array_equal(zs[0][name], zs[1][name])
    replayed = {m:replay(z, inp['source']) for m, z, inp in zip(METHODS, zs, inputs)}
    prof = [from_npz(p, float(z['assigned_power_W'])/1e6) for p, z in zip(paths, zs)]
    policy = distribution_compare(*prof)
    for metric in ('linf_relative', 'l1_relative', 'a2_relative', 'power_relative'):
        np.testing.assert_allclose(policy[metric], expected['comparison'][metric], rtol=3e-11, atol=1e-14)
    assert policy['failures'] == expected['comparison']['failures']
    edges = zs[0]['area_edges']
    supports = [p['continuous'] for p in prof]
    u = np.unique(np.r_[supports[0][0], supports[1][0]])
    curves = [np.interp(u, *v) for v in supports]
    j = int(np.argmax(abs(curves[1]-curves[0])))
    loc = wall_location(prof[0], edges, u[j])
    cell = loc['cell']
    cube = cube_at(zs, inputs, loc['area_m2'], cell, edges)
    effects = shapley(cube)
    delta = cube[(1,1,1)]-cube[(0,0,0)]
    endpoints = [cube[(i,i,i)].sum(axis=0) for i in (0,1)]
    for p, (_, qring) in zip(prof, endpoints):
        np.testing.assert_allclose(qring, p['q'][cell], atol=1e-9)
    point_residual = [endpoints[k][0]-curves[k][j] for k in (0,1)]
    if max(abs(np.array(point_residual))) > 1e-4:
        raise ValueError('point replay fails to approximate saved piecewise-linear q')
    peaks = []
    for method, p, (grid, values) in zip(METHODS, prof, supports):
        i = np.argmax(values)
        peak = wall_location(p, edges, grid[i])
        peak.update(q_MW_m2=values[i], A2_cell_m2=np.dot(p['q'],p['area'])**2 / np.dot(p['q']**2,p['area']))
        np.testing.assert_allclose([peak['q_MW_m2'], peak['A2_cell_m2']],
                                   [expected[method]['peak'], expected[method]['A2']], rtol=2e-11)
        peaks.append(peak)

    # Report the two local maxima separately when their locations exchange order.
    # Use equal neighborhoods around each global maximum on the common physical wall.
    neighborhoods = []
    if abs(peaks[0]['s_m'] - peaks[1]['s_m']) > .004:
        for reference_peak in peaks:
            center = reference_peak['s_m']
            record = {}
            for method, p, (grid, q) in zip(METHODS, prof, supports):
                s = grid*p['length']
                ids = np.flatnonzero(abs(s - center) <= .002)
                best = ids[np.argmax(q[ids])]
                assert ids[0] < best < ids[-1]
                record[method] = dict(s_m=float(s[best]), q_MW_m2=float(q[best]))
            neighborhoods.append(dict(center_m=center, radius_mm=2., maxima=record))

    walls = [Wall(z['wall'], inp['field']) for z, inp in zip(zs, inputs)]
    hit_s = [wall.hit_s(z['hit_points'], z['hit_segment']) for z, wall in zip(zs, walls)]
    shift = (hit_s[1]-hit_s[0] + walls[0].s[-1]/2) % walls[0].s[-1] - walls[0].s[-1]/2
    changes = zs[1]['hit_segment']-zs[0]['hit_segment']
    arcs = [int(np.sum(w.arc_id(hit_s[0]) != w.arc_id(hit_s[1]))) for w in walls]
    order = np.argsort(-abs(delta[:, 0]), kind='stable')
    ranked = order[:10]
    top = int(ranked[0])
    donors = []
    for i in range(len(shift)):
        row = dict(case=key, donor=i, point_delta=delta[i,0], ring_delta=delta[i,1],
                   hit_shift_um=1e6*shift[i])
        for k, (z, method) in enumerate(zip(zs, METHODS)):
            row.update({f'{method}_segment': int(z['hit_segment'][i]), f'{method}_hit_s_m': hit_s[k][i],
                f'{method}_R_m': z['hit_points'][i,0], f'{method}_Z_m': z['hit_points'][i,1],
                f'{method}_source_x_mm': 1e3*z['source_x_m'][i],
                f'{method}_width_mm': 1e3*z['sigma_physical_m'][i],
                f'{method}_target_gradient': z['target_gradient'][i],
                f'{method}_Lparallel_m': z['thermal_length_m'][i],
                f'{method}_weight_fraction': inputs[k]['source'].weights[i]})
        for f, effect in zip(FACTORS, effects):
            row[f'point_{f}'], row[f'ring_{f}'] = effect[i]
        donors.append(row)
    cells = [int(np.argmax(abs(prof[1]['q']-prof[0]['q']))), cell]
    cell_details = [dict(cell=i, parent=int(zs[0]['cell_parent'][i]),
                         arc_edges_m=(prof[0]['u'][i:i+2]*prof[0]['length']).tolist(),
                         q_MW_m2=[p['q'][i] for p in prof]) for i in cells]
    launch_limit = []
    for i, (inp, w) in enumerate(zip(inputs, walls)):
        turn = int(w.arc_id(hit_s[i][0]))
        segment = w.turn[turn:turn+2]
        root = brentq(lambda s: w.flux.values(w.points_at([s]))[0]-inp['psi_sep'], *segment,
                      xtol=2e-14)
        launch_limit.append(dict(s_m=root, RZ_m=w.points_at([root])[0].tolist(),
                                 policy_location_minus_flux_limit_mm=1e3*(loc['s_m']-root)))
    tags = ('point', 'ring')
    return dict(files={m:{'path':str(p), 'sha256':digest(p)} for m,p in zip(METHODS,paths)},
        policy=policy, replay=replayed, peaks=dict(zip(METHODS, peaks)),
        local_peak_neighborhoods=neighborhoods, near_separatrix_floor_flux_root=launch_limit,
        peak_shift_mm=1e3*(peaks[1]['s_m']-peaks[0]['s_m']),
        peak_change_relative=peaks[1]['q_MW_m2']/peaks[0]['q_MW_m2']-1,
        A2_change_relative=peaks[1]['A2_cell_m2']/peaks[0]['A2_cell_m2']-1,
        location=loc, policy_point_q_MW_m2=[curve[j] for curve in curves],
        exact_kernel_minus_saved_interpolation_MW_m2=point_residual,
        max_native_difference=cell_details[0], policy_location_ring=cell_details[1],
        factorial_q_MW_m2={''.join(map(str,bits)):dict(zip(tags,arr.sum(axis=0))) for bits,arr in cube.items()},
        shapley_MW_m2={factor:dict(zip(tags,arr.sum(axis=0))) for factor,arr in zip(FACTORS,effects)},
        top10_point_donors_by_absolute_delta=ranked.tolist(),
        top10_point_signed_delta_MW_m2=delta[ranked,0].sum(),
        donors_8_to_25_point_signed_delta_MW_m2=delta[8:26,0].sum(),
        top_donor=donors[top], top_donor_local_flux_test=flux_root_check(top,inputs,walls,hit_s,zs),
        donor16_local_flux_test=flux_root_check(16,inputs,walls,hit_s,zs),
        median_abs_hit_shift_um=np.median(abs(shift))*1e6, max_abs_hit_shift_um=max(abs(shift))*1e6,
        donor_with_max_abs_hit_shift=int(np.argmax(abs(shift))), segment_changed_ids=np.flatnonzero(changes).tolist(),
        max_abs_segment_increment=int(max(abs(changes))),
        monotone_arc_crossings_when_measured_on_each_field=dict(zip(METHODS,arcs)),
        max_source_x_offset_difference_nm=1e9*max(abs(zs[0]['source_x_m']-zs[1]['source_x_m'])),
        max_source_weight_difference=max(abs(inputs[0]['source'].weights-inputs[1]['source'].weights)),
        nearest_hit_edge_offset_mm=[1e3*(loc['s_m']-s[0]) for s in hit_s]), donors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True, help='parent of grid/ and sensitivity/')
    parser.add_argument('--historical', type=Path, required=True, help='historical stage 18 validation directory')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.study/'sensitivity'
    design = read(args.study/'grid/design.json')
    if digest(SHOT) != design['inputs']['input/D3D/d3d/d3d_shot_1a3155370c.parquet']:
        raise ValueError('equilibrium hash differs from search inputs')
    inputs = source_inputs(design['references'][0], 640, read(ROOT/'baseline_review/baseline_review.json'))
    report = dict(evidence_sha256={'shot':digest(SHOT), 'grid_design':digest(args.study/'grid/design.json'),
        'driver':digest(Path(__file__)), 'v3_policy':digest(ROOT/'viewers/d3d-progression-20261104/distribution_policy_v3_final.py')},
        interpolation_sources={m:dict(separatrix_R_m=x['rsep'], Xpoint_RZ_m=x['xpoint'].tolist(),
                                      psi_sep=x['psi_sep']) for m,x in zip(METHODS,inputs)},
        cases={})
    donors = []
    for chi in (.3,.6):
        label = f'OMP_chi{chi}'
        summary = read(root/f'{label}.json')
        for case in ('candidate', 'baseline'):
            paths = [root/f'{label}_{method}_{case}.npz' for method in METHODS]
            if (chi,case) == (.3,'candidate'):
                report['historical_identity'] = archived_identity(args.study/'grid',args.historical,paths)
            key = f'{label}_{case}'
            report['cases'][key], rows = pair(paths, summary[case], inputs, key)
            donors.extend(rows)
            print(key, 'replayed', flush=True)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'local-attribution.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    with (args.output/'local-donors.csv').open('w',newline='',encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=donors[0])
        writer.writeheader()
        writer.writerows(donors)


if __name__ == '__main__':
    main()
