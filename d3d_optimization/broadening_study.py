"""Evidence, persistence and acceptance rules shared by the new study runner."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from .run_provenance import REQUIRED

PRESERVATION = 'docs/superpowers/plans/2026-10-27-broadening-preservation.json'
STUDY_FILES = (
    'd3d_optimization/broadening.py', 'd3d_optimization/broadened_evaluator.py',
    'd3d_optimization/broadening_study.py', 'compare_broadening.py',
    'verification/validate_broadening.py', 'verification/plot_broadening.py',
    PRESERVATION,
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, data, exclusive=False):
    text = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if exclusive:
        with path.open('x', encoding='utf-8') as stream:
            stream.write(text)
    else:
        temporary = path.with_suffix(path.suffix + '.tmp')
        temporary.write_text(text, encoding='utf-8')
        temporary.replace(path)


def preserve_check(root):
    root = Path(root)
    original = load_json(root / PRESERVATION)
    changed = [name for name, value in original.items()
               if not (root / name).is_file() or digest(root / name) != value]
    if changed:
        raise ValueError(f'preserved production files changed: {changed}')
    return dict(unchanged=True, files=len(original))


def capture_provenance(root, destination, extra=()):
    root, destination = Path(root), Path(destination)
    preserve_check(root)
    names = sorted(set(REQUIRED) | set(STUDY_FILES) | set(extra))
    record = dict(schema='connection-broadening-provenance.v1',
                  captured_at_utc=datetime.now(timezone.utc).isoformat(),
                  files={name: digest(root/name) for name in names}, extra_files=list(extra))
    save(destination/'provenance.json', record, exclusive=True)
    return record


def verify_provenance(root, destination):
    root = Path(root)
    record = load_json(Path(destination)/'provenance.json')
    expected = set(REQUIRED) | set(STUDY_FILES) | set(record.get('extra_files', []))
    if record.get('schema') != 'connection-broadening-provenance.v1' or set(record.get('files', {})) != expected:
        raise ValueError('incomplete study provenance')
    for name, value in record['files'].items():
        if not (root/name).is_file() or digest(root/name) != value:
            raise ValueError(f'study code/input changed: {name}')
    preserve_check(root)
    return record


def acceptance(levels):
    """Four levels: coarse, fine, finer, and fine with quintic field.

    This empirical screen measures ratios AND absolute peak convergence. It
    provides neither a rigorous error bound nor physical validation.
    """
    failed = dict(numerically_verified_feasible=False, resolution_screen_passed=False,
                  improvement_confirmed=False, verification_status='incomplete_or_invalid_evidence',
                  physical_validation=False, global_optimum_proven=False)
    if len(levels) != 4 or not all(r.get('feasible') and r.get('audit_passed') for r in levels):
        return failed
    try:
        ratios = np.array([r['peak_ratios'] for r in levels], dtype=float)
        peaks = np.array([r['peak_q'] for r in levels], dtype=float)
        baselines = np.array([r['baseline_peak_q'] for r in levels], dtype=float)
    except (KeyError, ValueError, TypeError):
        return failed
    if (ratios.shape != (4, 5) or peaks.shape != (4, 5) or baselines.shape != (4, 5)
            or not all(np.all(np.isfinite(x)) for x in (ratios, peaks, baselines))
            or np.any(peaks <= 0) or np.any(baselines <= 0)):
        return failed
    ratio_delta = float(max(np.max(abs(ratios[2]-ratios[1])), np.max(abs(ratios[3]-ratios[1]))))
    absolute_delta = float(max(np.max(abs(values[j]-values[1])/values[1])
                               for values in (peaks, baselines) for j in (2, 3)))
    screen = ratio_delta <= .002 and absolute_delta <= .005
    nonworsening = bool(np.all(ratios <= 1+1e-10))
    reduction = 1-float(np.max(ratios[2]))
    return dict(numerically_verified_feasible=True, resolution_screen_passed=bool(screen),
                improvement_confirmed=bool(screen and nonworsening and reduction > 2*ratio_delta+.0001),
                nonworsening=nonworsening, finest_worst_peak_ratio=float(np.max(ratios[2])),
                finest_reduction_percent=100*reduction, max_ratio_change=ratio_delta,
                max_absolute_peak_relative_change=absolute_delta,
                verification_status='passed' if screen else 'resolution_failed',
                physical_validation=False, global_optimum_proven=False)


def gate_original_screen(verdict, original_baseline, original_candidate):
    result = dict(verdict)
    passed = bool(original_baseline.get('baseline_admissible') and original_candidate.get('feasible'))
    result['original_screen_passed'] = passed
    result['refined_verification_status'] = result['verification_status']
    if not passed:
        result.update(numerically_verified_feasible=False, improvement_confirmed=False,
                      verification_status='original_screen_failed')
    return result


def summarize(evaluator, result):
    row = dict(geometry_valid=bool(result['geometry_valid']), feasible=bool(result['feasible']),
               baseline_admissible=bool(evaluator.baseline_admissible))
    if not result['geometry_valid'] or 'peak_q' not in result:
        row['model_errors'] = result.get('model_errors', [])
        return row
    peak = np.asarray(result['peak_q'])
    base = evaluator.baseline_peaks
    if not np.all(np.isfinite(peak)) or base is None or np.any(base <= 0):
        row.update(feasible=False, model_errors=['invalid_peak_or_baseline'])
        return row
    ratios = peak/base
    row.update(peak_q=peak.tolist(), baseline_peak_q=base.tolist(), peak_ratios=ratios.tolist(),
               worst_peak_ratio=float(max(ratios)), missing_power=float(result['missing']),
               balance_error=float(result['power_balance_error']),
               surface_power=result['surface_power'])
    if 'deposition' in result:
        dep = result['deposition']
        row['power_ledger'] = {key: np.asarray(value).tolist() for key, value in dep['ledger'].items()}
        row['model_errors'] = result.get('model_errors', [])
        row['receiving_angles'] = result['receiving_angles']
        points, q = dep['points'], dep['q']
        row['peak_points_RZ_m'] = points[np.argmax(q, axis=0)].tolist()
        length = np.asarray(dep['lengths_m'])
        weights = np.column_stack([s.weights for s in evaluator.sources])
        finite = np.isfinite(length)
        row['connection_length_m'] = None
        if np.any(finite):
            denominator = weights[finite].sum(axis=0)
            average = np.divide(np.sum(length[finite, None]*weights[finite], axis=0), denominator,
                                out=np.zeros_like(denominator), where=denominator > 0)
            row['connection_length_m'] = dict(min=float(np.min(length[finite])),
                max=float(np.max(length[finite])), power_weighted_mean=average.tolist())
    return row
