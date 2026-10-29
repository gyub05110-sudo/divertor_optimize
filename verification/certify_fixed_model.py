"""Replayable whole-box geometry certificate; NOT a continuous physics proof."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from verification import global_certificate_geometry as g

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'input/D3D/optimized/global-certificate-20261024'

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':')).encode()).hexdigest()

def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')

def geometry(data, witnesses=None):
    """Recompute enclosures, never trust stored intervals or success flags."""
    if not data['surfaces'] or any(len(s['parameters'])<2 for s in data['surfaces']):
        raise ValueError('nonempty surfaces with at least two nodes required')
    try:hit=g.first_fixed_intersection(data['path'],data['fixed'])
    except ValueError as exc:
        if str(exc)!='stored path has no fixed-wall intersection':raise
        if witnesses not in (None,[]):raise ValueError('unexpected witnesses without fixed hit')
        return dict(all_moving_edges_excluded=False,reason=str(exc),fixed_index=None,
                    edge_count=0,failed_edges=[],conservative_separation_lower_bound_m=None),[]
    prefix=g.hull(hit['prefix'])
    bound=float(data['coefficient_bound_with_tolerance_m'])
    if not 0<bound<=.100000000002:raise ValueError('unsupported coefficient box')
    records=[];failed=[];minimum=None
    expected=sum(len(s['parameters'])-1 for s in data['surfaces'])
    if witnesses is not None and len(witnesses)!=expected:raise ValueError('incomplete edge witnesses')
    for surface in data['surfaces']:
        ss,ps,ns=surface['parameters'],surface['positions'],surface['normals']
        if len(ss)!=len(ps) or len(ps)!=len(ns):raise ValueError('incomplete node data')
        nodes=[g.node_envelope(p,n,g.displacement_interval(s,[(-bound,bound)]*4)) for s,p,n in zip(ss,ps,ns)]
        for index in range(len(nodes)-1):
            edge=g.hull(nodes[index]+nodes[index+1])
            if witnesses is None:axis=g.separator(prefix,edge)
            else:
                row=witnesses[len(records)]
                if (row['surface'],row['edge'])!=(surface['name'],index):raise ValueError('edge witness identity mismatch')
                axis=row['axis']
                if axis is not None and not g.check_separator(prefix,edge,axis):raise ValueError('invalid separating witness')
            row=dict(surface=surface['name'],edge=index,axis=axis)
            if axis is None:failed.append(dict(surface=surface['name'],edge=index))
            else:
                # ||axis||_2 <= ||axis||_1, so this is a distance LOWER bound.
                margin=g.gap(prefix,edge,axis)/sum(abs(v) for v in g.point(axis))
                minimum=margin if minimum is None else min(minimum,margin)
            records.append(row)
    lower=math.nextafter(float(minimum),-math.inf) if minimum is not None else None
    return dict(all_moving_edges_excluded=not failed,edge_count=len(records),failed_edges=failed,
                fixed_index=hit['fixed_index'],prefix_edges=len(hit['prefix'])-1,
                first_fixed_point_m=[float(v) for v in hit['point']],
                conservative_separation_lower_bound_m=lower),records

def replay(path):
    report=json.loads(Path(path).read_text(encoding='utf-8'))
    if report.get('schema')!='fixed-polyline-geometry-certificate.v1':raise ValueError('unknown certificate schema')
    if digest(report['frozen'])!=report['frozen_sha256']:raise ValueError('frozen data digest mismatch')
    result,_=geometry(report['frozen'],report['separating_witnesses'])
    if result!=report['geometry']:raise ValueError('replayed geometry result mismatch')
    # This is a geometry replay only. Numerical heat and library assumptions are
    # explicitly outside this exact-rational verifier.
    return result

def build(source,count,step,nodes,method):
    from d3d_optimization.independent_sources import prepare_evaluator
    from d3d_optimization.shape import ReferencePolyline
    from d3d_optimization.run_provenance import REQUIRED,digest as file_digest
    start=time.monotonic()
    fingerprints={p:file_digest(ROOT/p) for p in REQUIRED}
    proof_files=['verification/global_certificate_geometry.py','verification/certify_fixed_model.py']
    proof_fingerprints={p:file_digest(ROOT/p) for p in proof_files}
    e,meta=prepare_evaluator(ROOT,source,count,step,nodes,method,boundary_mode='continuous')
    frozen=dict(path=np.asarray(e.boundary_guard.path.coords).tolist(),
                fixed=[(int(i),e.wall[i:i+2].tolist()) for i in e.fixed_indices],
                coefficient_bound_with_tolerance_m=float(e.bound+1e-12),surfaces=[])
    for surface in e.surfaces:
        ref=ReferencePolyline(e.targets[surface['name']]);s=ref.node_parameters
        index=surface['segment_index']
        if not np.array_equal(ref.sample(s)[[0,-1]],e.wall[index:index+2]):
            raise ValueError('stitched target endpoints do not match frozen wall')
        frozen['surfaces'].append(dict(name=surface['name'],parameters=s.tolist(),
                                      positions=ref.sample(s).tolist(),normals=ref.unit_normal(s).tolist()))
    if len(frozen['surfaces'])!=6 or e.dimension!=24:raise ValueError('expected all six surfaces and 24 controls')
    print(f'{source} {method} {count}: inputs frozen; checking whole-box geometry',flush=True)
    result,witnesses=geometry(frozen)
    baseline=e.evaluate(np.zeros(e.dimension));w=baseline['boundary_witness']
    equality=[bool(q==b) for q,b in zip(w.get('q',[]),e.baseline_peaks)]
    fixed_match=bool(w.get('resolved') and w.get('surface')=='fixed' and w.get('split_segment')==result['fixed_index'])
    roots=[]
    if fixed_match:
        from scipy.optimize import brentq
        a,b=e.wall[result['fixed_index']:result['fixed_index']+2]
        residual=lambda s:float(e.field.values([a+s*(b-a)])[0]-e.boundary_guard.psi_sep)
        grid=np.linspace(0.,1.,13);values=[residual(s) for s in grid]
        roots=[float(s) for s,v in zip(grid,values) if abs(v)<1e-14]
        roots += [brentq(residual,l,h,xtol=1e-14) for l,h,u,v in zip(grid[:-1],grid[1:],values[:-1],values[1:]) if u*v<0]
    unique_root=len(roots)==1
    conditional=bool(result['all_moving_edges_excluded'] and fixed_match and unique_root and any(equality) and baseline['feasible'])
    for p,expected in {**fingerprints,**proof_fingerprints}.items():
        if file_digest(ROOT/p)!=expected:raise ValueError('inputs changed during certificate construction')
    report=dict(schema='fixed-polyline-geometry-certificate.v1',source=source,
        configuration=dict(interpolation=method,tubes=count,step_m=step,nodes=nodes,imp_buffer_m=0.),
        frozen=frozen,frozen_sha256=digest(frozen),geometry=result,separating_witnesses=witnesses,
        provenance=fingerprints,verifier_sha256=proof_fingerprints,
        baseline=dict(feasible=bool(baseline['feasible']),peak_q=e.baseline_peaks.tolist(),witness=w,
                      invariant_peak_equals_baseline=equality,fixed_selection_agrees=fixed_match,
                      fixed_refiner_root_candidates=roots,unique_operational_refiner_root=unique_root),
        conclusion=dict(exact_stored_geometry_exclusion=result['all_moving_edges_excluded'],
                        conditional_frozen_model_global_minimum=1. if conditional else None,
                        global_optimum_proven=False,continuous_model_proven=False,
                        experimental_validation=False),
        assumptions=[
            'Frozen source normalization and field/root evaluations are treated as fixed numerical constants.',
            'Operational GEOS first-hit selection agrees with exact rational stored-polyline geometry for every admissible candidate.',
            'IEEE binary64 round-to-nearest basic operations, no unmodelled arithmetic transformations.',
            'Baseline feasibility is the numerical production predicate, not a validated continuous-geometry proof.'
        ],
        scope='24 coefficient box; stored polygonal targets and sampled separatrix only. No validated ODE or field interpolation error bound.',
        elapsed_seconds=time.monotonic()-start)
    destination=OUT/f'{source}-{method}-{count}.json';save(destination,report)
    check=replay(destination)
    print(json.dumps(dict(file=str(destination),geometry=check,conclusion=report['conclusion'],
                          baseline=report['baseline'],elapsed_seconds=report['elapsed_seconds']),ensure_ascii=False),flush=True)

def challenge(count,step,nodes,method):
    from d3d_optimization.independent_sources import prepare_evaluator
    from verification.validate_continuous_run import audit_level
    from d3d_optimization.run_provenance import verify
    run=ROOT/'input/D3D/optimized/continuous-boundary-incumbent-20261024/OMP'
    provenance=verify(ROOT,run)
    m=json.loads((run/'manifest.json').read_text(encoding='utf-8'))
    c=np.array(m['verified_coefficients_m'])
    e,meta=prepare_evaluator(ROOT,'OMP',count,step,nodes,method,boundary_mode='continuous')
    rows=[]
    for delta in [0.,-.00001,-.0001,.00001,.0001]:
        x=c.copy();x[15]+=delta
        row,_=audit_level(e,meta,x)
        row.update(coefficient_index=15,delta_m=delta,coefficients_m=x.tolist())
        rows.append(row)
    original=rows[0]['worst_all_wall_peak_ratio']
    better=[r for r in rows[1:] if r['feasible'] and r['audit_passed'] and r['worst_all_wall_peak_ratio']<original]
    report=dict(source='OMP',interpolation=method,tubes=count,step_m=step,nodes=nodes,rows=rows,
                incumbent_has_feasible_better_neighbor=bool(better),
                best_tested_ratio=min([original]+[r['worst_all_wall_peak_ratio'] for r in better]),
                global_optimum_proven=False,scope='Counterexample for the stated floating-point discrete evaluator only.',
                incumbent_provenance=provenance)
    destination=OUT/f'OMP-challenge-{method}-{count}.json';save(destination,report)
    print(json.dumps(dict(file=str(destination),old=original,best=report['best_tested_ratio'],
                          counterexample=bool(better))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',choices=['IMP','OMP'],default='IMP')
    p.add_argument('--level',choices=['coarse','fine','quintic'],default='coarse')
    p.add_argument('--replay',type=Path)
    p.add_argument('--challenge',action='store_true')
    a=p.parse_args()
    if a.replay:print(json.dumps(replay(a.replay)));sys.exit()
    method,count,step,nodes={'coarse':('cubic',320,.002,161),'fine':('cubic',2560,.00025,1281),
                             'quintic':('quintic',1280,.0005,641)}[a.level]
    if a.challenge:challenge(count,step,nodes,method)
    else:build(a.source,count,step,nodes,method)
