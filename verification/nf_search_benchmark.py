"""Matched finite design-space benchmark of the production coarse search and refined peak.
Run from the repository root with python -m verification.nf_search_benchmark --output NEW_DIR.
All output paths are write-once. The experiment takes no FIREFLY data as input.
"""
import argparse
import gc
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import numpy as np
from d3d_optimization.broadening_study import digest, load_json, save
from d3d_optimization.boundary_evaluator import prepare_boundary
from d3d_optimization.campaign_geometry import Candidate, evaluator, evaluate_candidate
from d3d_optimization.reference_transport import ReferencePoint
from verification.validate_boundary_transport import LEVELS

ROOT = Path(__file__).resolve().parents[1]
ANCHORS = ROOT/'verification/fixtures/nf_search_anchors.json'


def triangle(anchors, subdivisions):
    a = np.asarray(anchors, float)
    if (a.ndim != 2 or len(a) != 3 or not np.isfinite(a).all()
            or np.linalg.matrix_rank(a[1:] - a[0]) != 2 or not isinstance(subdivisions, int) or subdivisions < 1):
        raise ValueError('three finite noncollinear design vectors and a positive integer required')
    return [(a[0]*(subdivisions-i-j) + i*a[1] + j*a[2])/subdivisions
            for i in range(subdivisions+1) for j in range(subdivisions+1-i)]


def shortlist(rows, count):
    eligible = [r for r in rows if r['eligible']]
    if count < 1 or any(not np.isfinite(r['peak']) or r['peak'] <= 0 for r in eligible):
        raise ValueError('positive shortlist size and positive finite eligible peaks required')
    return [r['i'] for r in sorted(eligible, key=lambda r: (r['peak'], r['i']))[:count]]


def measurement(result, reference):
    ledger = result.get('deposition', {}).get('ledger', {})
    row = dict(eligible=bool(result['feasible']), errors=result['errors'], ledger=ledger)
    if result.get('deposition') is not None:
        d = result['deposition']
        p = d['cell_power']
        row.update(peak=float(result['peak_q']*reference.power_W/1e6),
                   A2=float(p.sum()**2/np.sum(p*p/d['areas'])))
    return row


def prep(root, src, ref, level):
    method, count, step, nodes, wall_step = LEVELS[level]
    t = perf_counter()
    original = prepare_boundary(root, src, ref, count=count, step=step, nodes=nodes, wall_step=wall_step, method=method)
    if not original.baseline['feasible']:
        raise ValueError('supplied reference baseline failed: '+str(original.baseline['errors']))
    return original, perf_counter()-t


def evaluate(original, vectors, ids, backend):
    start = perf_counter()
    e = evaluator(original, 20, backend)
    rows = []
    for i in ids:
        c = Candidate(original.metadata['source'], 20, tuple(vectors[i]))
        out = evaluate_candidate(e, c)
        rows.append(dict(i=i, id=c.identity, **measurement(out, original.reference)))
    return rows, perf_counter()-start


def run_case(root, folder, src, ref, backend):
    vectors = load_json(folder.parent/'design.json')['families'][src]
    ids = list(range(len(vectors)))
    low, prep320 = prep(root, src, ref, 0)
    coarse, coarse_s = evaluate(low, vectors, ids, backend)
    chosen = shortlist(coarse, 3)
    save(folder/'shortlist.json', dict(reference=asdict(ref), source=src, ids=chosen, coarse=coarse), exclusive=True)
    del low
    gc.collect()
    high, prep1280 = prep(root, src, ref, 2)
    base = measurement(high.baseline, ref)
    finalists, final_s = evaluate(high, vectors, chosen, backend)
    # New evaluator with an empty geometry/thermal cache for the exhaustive comparator.
    all_high, exhaust_s = evaluate(high, vectors, ids, backend)
    for row in finalists:
        replay = all_high[row['i']]
        if row['eligible'] != replay['eligible'] or row.get('peak') != replay.get('peak'):
            raise RuntimeError('inconsistent precise replay')
    best = shortlist(all_high, 1)
    selected = shortlist(finalists, 1)
    fast_s = coarse_s + final_s
    delta = None
    if best and selected:
        delta = all_high[selected[0]]['peak']/all_high[best[0]]['peak']-1
    report = dict(source=src, reference=asdict(ref), baseline=base, coarse=coarse, fine=all_high,
                  shortlisted=chosen, selected=selected, finite_best=best, excess=delta,
                  fine_feasible=sum(r['eligible'] for r in all_high),
                  warm_fast_s=fast_s, warm_exhaustive_s=exhaust_s,
                  cold_fast_s=prep320+prep1280+fast_s, cold_exhaustive_s=prep1280+exhaust_s,
                  phase_seconds=dict(prep320=prep320, prep1280=prep1280, coarse=coarse_s,
                                     finalists=final_s, exhaustive=exhaust_s))
    save(folder/'result.json', report, exclusive=True)
    print(src, ref.chi_perp_m2_s, 'selected / finite best / excess', selected,best,delta,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--input-root', default=ROOT, type=Path)
    p.add_argument('--backend', choices=('gpu','cpu'), default='gpu')
    args=p.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    records=load_json(ANCHORS)
    families={}
    for src, anchors in records.items():
        candidates=[Candidate.from_record(a['candidate']) for a in anchors]
        if any(c.source!=src or c.identity!=a['original_id'] for c,a in zip(candidates,anchors)):
            raise ValueError('invalid anchor provenance')
        families[src]=[v.tolist() for v in triangle([c.vector for c in candidates],6)]
    refs=[ReferencePoint(chi_perp_m2_s=chi) for chi in (.3,.6)]
    save(out/'design.json', dict(families=families, references=[asdict(r) for r in refs], radius_mm=20,
         anchors_sha256=digest(ANCHORS), driver_sha256=digest(Path(__file__)), backend=args.backend,
         inputs={k:digest(args.input_root/k) for k in ('input/D3D/d3d/d3d_shot_1a3155370c.parquet',
                    'baseline_review/baseline_review.json','input/D3D/optimized/d3d-divertor-floor-shelf/manifest.json')},
         levels=[LEVELS[i] for i in (0,2)], k=3, finite_grid_excess_tolerance=.01), exclusive=True)
    for ref in refs:
        for src in ('OMP','IMP'):
            folder=out/f'{src}_chi{ref.chi_perp_m2_s:g}'
            folder.mkdir()
            run_case(args.input_root, folder, src, ref, args.backend)
            gc.collect()

if __name__=='__main__':main()

