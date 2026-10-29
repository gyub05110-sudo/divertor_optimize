"""Immutable reference-study evidence and scalar nominal-width acceptance."""
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from .broadening_study import load_json,save,digest
from .run_provenance import REQUIRED

PRESERVATION='docs/superpowers/plans/2026-10-28-reference-preservation.json'
NEW_FILES=('d3d_optimization/reference_transport.py','d3d_optimization/transport_deposition.py',
           'd3d_optimization/transport_evaluator.py','d3d_optimization/reference_study.py',
           'compare_reference_transport.py','verification/validate_reference_transport.py',
           'verification/plot_reference_transport.py',PRESERVATION)


def verify_preserved(root):
    root=Path(root);original=load_json(root/PRESERVATION)
    changed=[name for name,value in original.items() if not (root/name).is_file() or digest(root/name)!=value]
    if changed:raise ValueError(f'preserved files changed: {changed}')
    return dict(unchanged=True,files=len(original))


def capture(root,folder,extra=()):
    root=Path(root);verify_preserved(root)
    names=sorted(set(REQUIRED)|set(NEW_FILES)|set(load_json(root/PRESERVATION))|set(extra))
    record=dict(schema='reference-transport-provenance.v1',extra=list(extra),
                captured_at_utc=datetime.now(timezone.utc).isoformat(),files={n:digest(root/n) for n in names})
    save(Path(folder)/'provenance.json',record,exclusive=True)
    return record


def verify(root,folder):
    root=Path(root);record=load_json(Path(folder)/'provenance.json')
    expected=set(REQUIRED)|set(NEW_FILES)|set(load_json(root/PRESERVATION))|set(record.get('extra',[]))
    if record.get('schema')!='reference-transport-provenance.v1' or set(record['files'])!=expected:
        raise ValueError('incomplete reference provenance')
    if any(not (root/n).is_file() or digest(root/n)!=v for n,v in record['files'].items()):
        raise ValueError('reference-study code/input changed')
    verify_preserved(root)


def acceptance(levels,original_passed):
    result=dict(numerically_verified_feasible=False,improvement_confirmed=False,
                physical_validation=False,global_optimum_proven=False,
                original_screen_passed=bool(original_passed),verification_status='incomplete_or_invalid_evidence')
    if not original_passed:
        result['verification_status']='original_screen_failed'
        return result
    if len(levels)!=4 or not all(r.get('feasible') and r.get('audit_passed') for r in levels):return result
    try:
        ratios=np.array([r['peak_ratios'] for r in levels],float)
        peaks=np.array([r['peak_q'] for r in levels],float)
        bases=np.array([r['baseline_peak_q'] for r in levels],float)
    except (KeyError,ValueError,TypeError):return result
    if any(x.shape!=(4,1) or not np.all(np.isfinite(x)) or np.any(x<=0) for x in [ratios,peaks,bases]):return result
    ratio_change=float(max(abs(ratios[2,0]-ratios[1,0]),abs(ratios[3,0]-ratios[1,0])))
    absolute_change=float(max(np.max(abs(x[j]-x[1])/x[1]) for x in [peaks,bases] for j in [2,3]))
    passed=ratio_change<=.002 and absolute_change<=.005
    gain=1-float(ratios[2,0])
    result.update(numerically_verified_feasible=True,resolution_screen_passed=passed,
        finest_peak_ratio=float(ratios[2,0]),finest_reduction_percent=100*gain,
        max_ratio_change=ratio_change,max_absolute_peak_relative_change=absolute_change,
        improvement_confirmed=bool(passed and np.all(ratios<=1+1e-10) and gain>2*ratio_change+1e-4),
        verification_status='passed' if passed else 'resolution_failed')
    return result


def export(path,e,result):
    if not result['geometry_valid'] or 'peak_q' not in result:return
    if 'deposition' in result:
        dep=result['deposition'];points,q=dep['points'],dep['q']
        extra={k:dep[k] for k in ['sigma_m','lengths_m','cell_power','cell_parent','cell_xi_m']}
        extra.update({'thermal_'+k:v for k,v in result['thermal'].items()})
    else:
        points=result['maps'][0].points;q=result['maps'][0].q[:,None];extra={}
        w=result.get('boundary_witness')
        if w and w['resolved']:
            points=np.vstack([points,w['point']]);q=np.vstack([q,w['q']])
    np.savez_compressed(path,wall=result['wall'],points=points,q_over_P=q,
                        q_MW_m2=q*e.reference.power_W/1e6,peak_q=result['peak_q'],
                        source_x_m=e.sources[0].x,**extra)
