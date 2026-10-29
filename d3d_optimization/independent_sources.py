"""Actual IMP/OMP adapters around the preserved direct magnetic mapping kernel."""
import json
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Polygon, LineString
from .case import load_case, WIDTHS_MM, sha256, separatrix
from scipy.optimize import brentq
from .frozen_field import exponential_source
from .frozen_mapping import trace_paths
from .multisurface import MultiSurfaceEvaluator


def midplane_source(field, r_sep, z_up, edges, width, source):
    if source not in ('IMP', 'OMP'):
        raise ValueError('source must be IMP or OMP')
    if source == 'OMP':
        return exponential_source(field, r_sep, z_up, edges, width)
    class ReflectedSampling:
        def points(self, p):
            p = np.array(p, dtype=float, copy=True)
            p[:, 0] = 2*r_sep-p[:, 0]
            return p
        def magnetic(self, p): return field.magnetic(self.points(p))
        def gradient(self, p): return field.gradient(self.points(p))
    return exponential_source(ReflectedSampling(), r_sep, z_up, edges, width)


def split_design_wall(review):
    wall = np.asarray(review['wall'], dtype=float)
    selected = {}
    for row in review['optimization_footprints']:
        i = row['segment']
        a, b = wall[i:i+2]; direction = b-a
        alpha = np.clip((np.asarray(row['points'])-a) @ direction / (direction @ direction), 0, 1)
        if np.ptp(alpha) > 1e-12:
            selected[i] = (float(min(alpha)), float(max(alpha)))
    for i in review['optimization_policy']['existing_design_segments']: selected[i] = (0., 1.)
    names = {58:'omp_entry', 59:'shelf', 60:'omp_connector', 68:'floor', 69:'imp_sloped', 70:'imp_vertical'}
    vertices = [wall[0]]; original = []; surfaces = []
    for i, (a, b) in enumerate(zip(wall[:-1], wall[1:])):
        limits = selected.get(i)
        cuts = sorted(set([0., 1.] + (list(limits) if limits else [])))
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            if hi-lo < 1e-12: continue
            index = len(original); endpoints = np.array([a+lo*(b-a), a+hi*(b-a)])
            vertices.append(endpoints[1]); original.append(i)
            if limits and abs(lo-limits[0]) < 1e-12 and abs(hi-limits[1]) < 1e-12:
                surfaces.append(dict(name=names[i], segment_index=index, original_segment=i,
                                     endpoints_RZ_m=endpoints.tolist(), angle_floor_deg=1.5 if i==59 else None))
    return np.asarray(vertices), surfaces, np.asarray(original)


def prepare_evaluator(root, source, count=320, step=.002, nodes=161, method='cubic', bound=.1,
                      boundary_mode='legacy', imp_buffer_m=0.):
    root = Path(root)
    if source not in ('IMP','OMP') or count < 8 or step <= 0 or nodes < 9:
        raise ValueError('invalid source or numerical resolution')
    if (boundary_mode not in ('legacy','continuous') or not np.isfinite(imp_buffer_m)
            or not 0<=imp_buffer_m<=.02 or (boundary_mode=='legacy' and imp_buffer_m!=0)):
        raise ValueError('invalid boundary policy or IMP buffer (maximum 20 mm)')
    review = json.loads((root/'baseline_review/baseline_review.json').read_text(encoding='utf-8'))
    shot = root/'input/D3D/d3d/d3d_shot_1a3155370c.parquet'
    reference = root/'input/D3D/optimized/d3d-divertor-floor-shelf'
    old = json.loads((reference/'manifest.json').read_text(encoding='utf-8'))
    if sha256(shot) != old['shot_sha256']: raise ValueError('restored shot hash mismatch')
    field, lcfs, axis = load_case(shot, old['time_ms'], method=method)
    info = review['sources'][source]
    r_omp, xp, psi_sep = separatrix(field, lcfs, axis)
    if source == 'OMP': r_sep = r_omp
    else:
        rr = np.linspace(field.r[0]+1e-5, axis[0]-.02, 500)
        values = field.values(np.c_[rr,np.full(len(rr),axis[1])])-psi_sep
        roots = [brentq(lambda r: float(field.values([[r,axis[1]]])[0]-psi_sep),a,b,xtol=1e-13)
                 for a,b,u,v in zip(rr[:-1],rr[1:],values[:-1],values[1:]) if u*v<0]
        if len(roots)!=1: raise ValueError('ambiguous IMP separatrix')
        r_sep = roots[0]
    sg = -1 if source=='IMP' else 1
    crossing = shapely.intersection(LineString(review['wall']),LineString([[field.r[0],axis[1]],[field.r[-1],axis[1]]]))
    points = list(crossing.geoms) if hasattr(crossing,'geoms') else [crossing]
    limit = min(.1,min((p.x-r_sep)*sg for p in points if (p.x-r_sep)*sg>0))
    edges = limit*np.linspace(0,1,count+1)**2
    starts = np.c_[r_sep+(-1 if source=='IMP' else 1)*(edges[:-1]+edges[1:])/2, np.full(count,axis[1])]
    if boundary_mode=='continuous':
        from .continuous_boundary import continuous_design, BoundaryAwareEvaluator
        # Engineering design endpoints are fixed by the cubic reference case,
        # independent of the interpolant used for a physics sensitivity check.
        if method=='cubic': design_field, design_psi=field,psi_sep
        else:
            design_field,design_lcfs,design_axis=load_case(shot,old['time_ms'],method='cubic')
            design_psi=separatrix(design_field,design_lcfs,design_axis)[2]
        review,strike=continuous_design(review,design_field,design_psi,imp_buffer_m)
    wall, surfaces, original = split_design_wall(review)
    inside = shapely.contains(Polygon(wall), shapely.points(starts))
    if not inside.all(): raise ValueError('source starts outside material boundary')
    paths = [LineString(p) for p in trace_paths(field, starts, step=step, max_length=8)]
    geometry = json.loads((reference/'DIIID_geom.json').read_text(encoding='utf-8'))
    vessel = shapely.union_all([Polygon(v['pts']) for v in geometry['vv']])
    # Initialize preserved geometry/cached-hit machinery with its historical OMP input.
    e = MultiSurfaceEvaluator(field,lcfs,wall,paths,inside,edges,
                             r_omp,axis[1],surfaces,
                             target_nodes=nodes,coefficient_bound=bound,vessel_solids=vessel)
    e.dimension = 4*len(surfaces)
    e.sources = [midplane_source(field,r_sep,axis[1],edges,mm/1000,source) for mm in WIDTHS_MM]
    e.baseline = e.evaluate(np.zeros(e.dimension))
    e.baseline_peaks = np.array([float(max(m.q)) for m in e.baseline['maps']])
    if boundary_mode=='continuous':
        e=BoundaryAwareEvaluator(e,xp,psi_sep,source,[r_sep,float(axis[1])],min(.0005,step/4))
    meta = dict(source=source,source_power=1.,source_split=None,shot_path=str(shot),shot_sha256=old['shot_sha256'],
                time_ms=old['time_ms'],geometry_sha256=sha256(reference/'DIIID_geom.json'),
                source_RZ_m=[r_sep,float(axis[1])],support_mm=limit*1000,tubes=count,step_m=step,
                target_nodes=nodes,interpolation=method,dimension=e.dimension,surfaces=surfaces,
                split_to_original_segment=original.tolist(),policy=review['optimization_policy'],widths_mm=WIDTHS_MM.tolist())
    meta.update(boundary_mode=boundary_mode,imp_buffer_m=imp_buffer_m,
                objective_version='sol-limit-guard.v1' if boundary_mode=='continuous' else 'sampled-peak.v1')
    return e, meta
