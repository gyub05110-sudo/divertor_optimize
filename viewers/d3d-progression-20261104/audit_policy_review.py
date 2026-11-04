"""Read-only audit of sealed OMP studies; writes a new evidence directory only."""
import argparse
import json
import shutil
import time
from pathlib import Path

import numpy as np
from scipy.special import ndtr

from qualify_omp_stages import read, save, sha, checked_read, check_preservation
from distribution_policy_v3_final import from_npz, distribution_compare, qualify, POLICY, LIMITS

HERE = Path(__file__).resolve().parent
ROOT = Path('D:/DivertorOptimization/studies')
FINAL = ROOT / 'omp-distribution-final-20261104'
OLD = ROOT / 'omp-stage-qualification-20261104'
NAMES = ['cubic_320', 'cubic_640', 'cubic_1280', 'quintic_640', 'wall_refinement']
PAIRS = [(1, 2), (1, 3), (2, 3), (0, 4)]


def verify_seal(folder):
    marker = read(folder / 'COMPLETED.json')
    assert marker['summary_sha256'] == sha(folder / 'summary.json')
    assert marker['artifact_index_sha256'] == sha(folder / 'artifact_sha256.json')
    files = read(folder / 'artifact_sha256.json')
    for rel, digest in files.items():
        assert sha(folder / rel) == digest, str(folder / rel)
    return len(files)


def cell_a2(path):
    with np.load(path, allow_pickle=False) as z:
        p = z['cell_power']
        area = np.diff(z['area_edges'])
        assert np.all(p >= 0) and np.all(area > 0)
        assert abs(float(p.sum()) - 1) < 2e-12
        return float(p.sum() ** 2 / np.sum(p * p / area))


def assess(paths):
    profiles = [from_npz(p) for p in paths]
    checks = [dict(pair=[NAMES[i], NAMES[j]],
                   **distribution_compare(profiles[i], profiles[j], profiles[2]))
              for i, j in PAIRS]
    return dict(passed=all(c['passed'] for c in checks), checks=checks), profiles


def direct_density(path, profile, u):
    """Independent real-space normal sum, including broad periodic kernels."""
    with np.load(path, allow_pickle=False) as z:
        wall = z['receiving_wall']
        i = np.clip(np.searchsorted(profile['u'], u, side='right') - 1, 0, len(wall)-2)
        t = (u - profile['u'][i]) / (profile['u'][i+1] - profile['u'][i])
        ds = np.linalg.norm(wall[i+1] - wall[i], axis=1)
        r, dr = wall[i, 0], wall[i+1, 0] - wall[i, 0]
        area = z['area_edges'][i] + 2*np.pi*ds*(r*t + .5*dr*t*t)
        period = float(z['area_edges'][-1])
        weights = z['donor_region_power'].sum(axis=1)
        assert abs(float(weights.sum()) - 1) < 2e-12
        q = np.zeros_like(u)
        for weight, center, sigma in zip(weights, z['donor_area_m2'], z['sigma_area_m2']):
            if weight == 0:
                continue
            assert sigma > 0
            # At least 10 sigma tail coverage; no normalization of saved weights.
            images = int(np.ceil(10*sigma/period)) + 1
            for k in range(-images, images+1):
                q += weight*np.exp(-.5*((area-center+k*period)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        return q


def run(out):
    started = time.perf_counter()
    if out.exists():
        raise FileExistsError('Use a new audit directory; do not overwrite evidence')
    config = read(FINAL / 'config.json')
    seal_counts = {str(p): verify_seal(p) for p in (OLD, FINAL)}
    for key in ('protected_sha256', 'code_sha256'):
        assert check_preservation(config[key])['unchanged'], key
    # Preserve exactly the mutable presentation files before correction.
    out.mkdir(parents=True)
    for name in ('dist/app.mjs', 'dist/index.html', 'README.md'):
        target = out / 'before' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(HERE / name, target)
    shutil.copyfile(__file__, out / 'audit_policy_review.py')
    protected = dict(config['protected_sha256'])
    protected.update(config['code_sha256'])
    for p in FINAL.rglob('*'):
        if p.is_file():
            protected[str(p)] = sha(p)
    save(out / 'preserved_sha256.json', protected)
    saved = read(FINAL / 'summary.json')
    base_paths = [OLD / 'baseline' / n / 'baseline.npz' if n != 'wall_refinement'
                  else OLD / 'baseline/wall_refinement.npz' for n in NAMES]
    base, _ = assess(base_paths)
    assert base['passed'] == saved['baseline']['passed']
    base_a2 = [cell_a2(p) for p in base_paths]
    rows = []
    direct = None
    for item in sorted(saved['candidates'], key=lambda r: (r['radius_mm'], r['stage'])):
        v3 = read(FINAL / 'validations' / item['identity'] / 'validation.json')
        prior_path = Path(v3['prior_validation_path'])
        assert sha(prior_path) == v3['prior_validation_sha256']
        prior = checked_read(prior_path)
        paths = [prior_path.parent / (n+'.npz') for n in NAMES]
        for rel, digest in prior['artifacts'].items():
            assert sha(prior_path.parent / rel) == digest
        selected, profiles = assess(paths)
        result = qualify(prior, dict(policy=POLICY, passed=base['passed'] and selected['passed'],
                                     baseline=base, selected=selected))
        assert result['verdict']['improvement_confirmed'] == item['passed']
        maximum = max(c['linf_relative'] for c in selected['checks'])
        assert abs(maximum-item['max_linf_relative']) < 1e-12
        a2 = [cell_a2(p) for p in paths]
        peaks = [l['rows']['selected']['peak_MW_m2'] for l in prior['levels']]
        baseline_peaks = [l['rows']['baseline']['peak_MW_m2'] for l in prior['levels']]
        rows.append(dict(identity=item['identity'], radius_mm=item['radius_mm'], stage=item['stage'],
                         v3_passed=item['passed'], levels=NAMES, peak_MW_m2=peaks,
                         peak_reduction_percent=[100*(1-p/b) for p, b in zip(peaks, baseline_peaks)],
                         A2_m2=a2, A2_gain_percent=[100*(a/b-1) for a, b in zip(a2, base_a2)],
                         observed_fine_peak_range_MW_m2=max(peaks[1:])-min(peaks[1:]),
                         max_linf_relative=maximum,
                         max_l1_relative=max(c['l1_relative'] for c in selected['checks']),
                         checks=selected['checks']))
        if item['radius_mm'] == 20 and item['stage'] == 18:
            a, b = profiles[1], profiles[3]
            u = np.unique(np.r_[a['continuous'][0], b['continuous'][0]])
            delta = abs(np.interp(u, *a['continuous']) - np.interp(u, *b['continuous']))
            pos = u[delta.argmax()]
            sample = np.linspace(pos-.002, pos+.002, 10001)
            qa = direct_density(paths[1], a, sample)
            qb = direct_density(paths[3], b, sample)
            scale = max(a['continuous'][1].max(), b['continuous'][1].max())
            direct = dict(stage=18, identity=item['identity'], samples=len(sample),
                          RZ_m=[float(np.interp(pos, a['u'], a['wall'][:,j])) for j in (0, 1)],
                          direct_cubic640_MW_m2=float(qa[5000]), direct_quintic640_MW_m2=float(qb[5000]),
                          local_linf_relative=float(max(abs(qa-qb))/scale),
                          cubic_interpolation_error_MW_m2=float(max(abs(qa-np.interp(sample, *a['continuous'])))),
                          quintic_interpolation_error_MW_m2=float(max(abs(qb-np.interp(sample, *b['continuous'])))),
                          scope='Direct periodic Gaussian sum over local interval; not a global continuum certificate')
        print('RECHECK', item['radius_mm'], item['stage'], item['passed'], flush=True)
    assert len(rows) == 20 and direct is not None
    preservation = check_preservation(protected)
    assert preservation['unchanged']
    result = dict(complete=True, policy=POLICY, limits_unchanged=LIMITS, seal_file_counts=seal_counts,
                  reference=saved['reference'], baseline_A2_m2=base_a2, candidates=rows,
                  direct_kernel_check=direct, preservation=preservation,
                  elapsed_seconds=time.perf_counter()-started,
                  limitations=['A2 is computed from native cell averages; not a continuum integral certificate',
                               'Observed resolution/interpolation spread is not a rigorous error bound',
                               'No new optimization, changed tolerance, changed pass flag, or promoted final'])
    save(out / 'audit.json', result)
    print('AUDIT COMPLETE', result['elapsed_seconds'], direct, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
