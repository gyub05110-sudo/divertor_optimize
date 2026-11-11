"""Inspect the frozen source, wall, and receiving domain of the Stage 5 comparison."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from verification.stage5_native import audit_source, audit_tally


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def geometry(folder):
    with (folder/'geometry.csv').open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    return np.array([[float(r['R_m']), float(r['Z_m'])] for r in rows]), np.array([r['outgoing_region'] for r in rows[:-1]])


def selected_vertices(wall, names, surface):
    return np.array([wall[i:i+2] for i in np.flatnonzero(names == surface)])


def run(study):
    panel = read(study/'inputs/panel.json')
    execution = read(study/'execution.json')
    assert sha(study/'execution.json') == '3151ee0a78584acd70714847b3957a16730c13347053d3ef2d2add070770e023'
    assert sha(study/'inputs/panel.json') == execution['files']['inputs/panel.json']
    cases = []
    for case in panel['cases']:
        cid = case['case_id']
        inp = study/'inputs'/cid
        assert sha(inp/'case.json') == panel['case_record_sha256'][cid]
        source = audit_source(inp)
        donors = np.genfromtxt(inp/'OMP.csv', names=True, delimiter=',')
        rsep, zcut = case['topology']['omp_separatrix_RZ_m']
        assert np.all(donors['R_m'] > rsep) and np.all(donors['Z_m'] == zcut)
        baseline = next(r for r in case['roles'] if r['role'] == 'baseline')
        ref, names = geometry(Path(baseline['stage4_directory']))
        geometries = {}
        inner_roundoff = 0.
        for role in case['roles']:
            folder = Path(role['stage4_directory'])
            assert sha(folder/'geometry.csv') == role['files_sha256']['geometry.csv']
            wall, labels = geometry(folder)
            for surface in ('imp_sloped', 'imp_vertical'):
                before = selected_vertices(ref, names, surface)
                after = selected_vertices(wall, labels, surface)
                assert before.shape == after.shape
                delta = float(np.max(abs(before-after)))
                inner_roundoff = max(inner_roundoff, delta)
                # Reconstructing a zero-coefficient spline introduces a few ulps.
                assert delta <= 32*np.finfo(float).eps*float(np.max(abs(before)))
            ident = role['wall_id']
            if ident in geometries:
                continue
            assert np.array_equal(wall, np.loadtxt(inp/role['wall_file']))
            pfc = study/'native/cases'/cid/(ident+'.pfc.nc')
            assert sha(pfc) == read(study/'native/cases'/cid/'ready.json')['walls'][ident]['pfc_sha256']
            with h5py.File(pfc) as file:
                assert set(file) == {'wall', 'reflect_midplane'}
                plane = file['reflect_midplane/rz'][:]
                assert np.all(plane[:,1] == zcut)
                assert plane[:,0].min() < wall[:,0].min() and plane[:,0].max() > wall[:,0].max()
                assert np.array_equal(file['wall/rz'][:], wall)
            geometries[ident] = []
            for seed in execution['seeds']:
                raw = study/'native/results'/cid/ident/(str(seed)+'.nc')
                receipt = read(raw.with_suffix('.json'))
                assert receipt['raw_sha256'] == sha(raw)
                assert receipt['source_sha256'] == sha(inp/'OMP.source')
                assert receipt['pfc_sha256'] == sha(pfc)
                ledger = audit_tally(raw, wall)
                assert ledger['wall_above_midplane_W'] == 0.
                assert ledger['cap_W'] == ledger['domain_loss_W'] == ledger['censored_W'] == 0.
                geometries[ident].append(dict(seed=seed, raw_sha256=receipt['raw_sha256'], **ledger))
        folder = Path(baseline['stage4_directory'])
        assert sha(folder/'boundary.npz') == baseline['files_sha256']['boundary.npz']
        with np.load(folder/'boundary.npz') as raw:
            assert np.array_equal(raw['source_x_m'], donors['x_mid_m'])
            assert np.allclose(donors['R_m'], rsep+raw['source_x_m'], rtol=0., atol=1e-14)
            direct = np.loadtxt(inp/baseline['hit_file'])
            assert np.array_equal(raw['hit_points'], direct[:,:2])
            hit_names = names[raw['hit_segment']]
            hit_watts = {n: float(donors['power_W'][hit_names == n].sum()) for n in np.unique(hit_names)}
            receiver = raw['receiving_wall']
            above = np.minimum(receiver[:-1,1], receiver[1:,1]) > zcut
            upper = float((raw['cell_power'][above]*1e6).sum())
            temperature = {k:dict(min=float(raw[k].min()), max=float(raw[k].max()),
                                 power_weighted_mean=float(np.average(raw[k], weights=donors['power_W'])))
                           for k in ('thermal_Tu_eV', 'thermal_Tt_eV', 'thermal_nt_m3')}
        cases.append(dict(case_id=cid, source='OMP', source_sha256=sha(inp/'OMP.source'),
                          source_csv_sha256=sha(inp/'OMP.csv'), **source,
                          source_anchor_RZ_m=[rsep,zcut], all_source_R_above_separatrix=True,
                          independently_fixed_inner_surfaces=True, inner_vertex_roundoff_m=inner_roundoff,
                          first_hit_power_W=hit_watts,
                          reduced_wholly_above_midplane_W=upper, reduced_thermal=temperature,
                          native_surfaces=['wall', 'reflect_midplane'], native_reflecting_plane_Z_m=zcut,
                          native_fixed_background=dict(n0=1e19, T0=100., chi=.45), geometries=geometries))
    assert len(cases) == 7 and sum(len(c['geometries']) for c in cases) == 46
    assert sum(len(rows) for c in cases for rows in c['geometries'].values()) == 276
    return dict(schema='independent-source-domain-audit.v1', execution_sha256=sha(study/'execution.json'),
                driver_sha256=sha(Path(__file__)), cases=cases)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = (json.dumps(run(args.study), indent=2, allow_nan=False)+'\n').encode()
    args.output.parent.mkdir(exist_ok=True, parents=True)
    if args.output.exists() and args.output.read_bytes() != data:
        raise ValueError('audit output already exists with different bytes')
    args.output.write_bytes(data)
    print('Verified 7 OMP source experiments, 46 fixed geometries, and 276 native output ledgers.')
