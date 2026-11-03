"""Audit saved CPD results and render figures; never reruns transport/search.

Run from the repository with the original D: studies available. All exported
NPZs are byte-for-byte copies. The additional geometry check is postprocessing
and is not counted as optimization or transport qualification time.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad
from shapely.geometry import Polygon

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[2]
sys.path.insert(0, str(REPO))
from cpd_optimization.geometry import load_geometry
from cpd_optimization.deformation import DivertorDeformation

RAW = Path('D:/DivertorOptimization/studies/cpd5-20261103')
RUN = RAW / 'resume-remaining-v1'
SNAP = RAW / 'resume-execution-v2'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
study = json.loads((RUN / 'study.json').read_text())
parent = json.loads((RAW / 'campaign-3h/study.json').read_text())
manifest = json.loads((SNAP / 'execution-manifest.json').read_text())
assert study['status'] == 'completed'
assert study['deadline_overrun_seconds'] == 0
assert study['total_consumed_seconds'] <= 10800
assert abs(study['prior_consumed_seconds'] + study['campaign_elapsed_seconds'] - study['total_consumed_seconds']) < .01
assert sha(RAW / 'campaign-3h/study.json') == study['resume']['parent_study_sha256']
assert sha(RAW / 'campaign-3h/candidates.jsonl') == study['resume']['parent_candidate_log_sha256']
for name, digest in manifest['files'].items():
    assert sha(SNAP / name) == digest, name
for key, name in [('equilibrium_sha256', 'equilibrium/equilibrium.npz'), ('geometry_sha256', 'CPD_5_nodes.yaml')]:
    assert sha(REPO / 'input/CPD_5' / name) == study[key] == parent[key]
assert study['fixed_inputs'] == parent['fixed_inputs']
shutil.copy2(RUN / 'study.json', OUT / 'study.json')
shutil.copy2(RUN / 'status.json', OUT / 'status.json')

geometry = load_geometry(REPO / 'input/CPD_5/CPD_5_nodes.yaml', max_chord_m=.00025)
with np.load(REPO / 'input/CPD_5/equilibrium/equilibrium.npz') as f:
    lcfs = f['lcfs'].copy()
deformation = DivertorDeformation(geometry, max_chord_m=.00025, lcfs=Polygon(lcfs))
audit = dict(status='passed', execution_files_verified=len(manifest['files']),
    original_checkpoint_unchanged=True, fixed_inputs_unchanged=True,
    original_budget_seconds=10800, prior_consumed_seconds=study['prior_consumed_seconds'],
    resumed_seconds=study['campaign_elapsed_seconds'], total_consumed_seconds=study['total_consumed_seconds'],
    qualification_seconds=study['campaign_elapsed_seconds']-study['search_elapsed_seconds'],
    geometry_audit_chord_m=.00025, sources={}, artifacts={})
data = {}
for source, source_data in study['sources'].items():
    rows = source_data['qualification_rows']
    coefficients = np.asarray(source_data['best_coefficients'])
    geom = deformation.evaluate(coefficients)
    assert geom['feasible'], geom['errors']
    assert np.max(abs(coefficients)) <= .1
    dense_min_radius = np.inf
    quad_length_ratio = 0.
    endpoint_error = tangent_error = 0.
    for name, curve in deformation.curves.items():
        c = coefficients[deformation.slices[name]]
        p, v, _ = curve.evaluate([0., 1.], c)
        p0, v0, _ = curve.evaluate([0., 1.], np.zeros(4))
        endpoint_error = max(endpoint_error, float(np.max(abs(p-p0))))
        tangent_error = max(tangent_error, float(np.max(abs(v-v0))))
        dense_min_radius = min(dense_min_radius, 1/np.max(curve.curvature(np.linspace(0, 1, 10001), c)))
        length = quad(lambda t: float(np.linalg.norm(curve.evaluate([t], c)[1][0])), 0, 1,
                      points=[.2,.4,.6,.8], epsabs=1e-11, epsrel=1e-11)[0]
        quad_length_ratio = max(quad_length_ratio, length/curve.reference_length_m)
    assert endpoint_error == tangent_error == 0
    assert dense_min_radius >= .020 and quad_length_ratio <= 1.2
    ga = dict(feasible=True, endpoint_max_error_m=endpoint_error, tangent_max_error=tangent_error,
        max_abs_coefficient_m=float(np.max(abs(coefficients))),
        max_length_ratio_independent_quadrature=quad_length_ratio,
        min_radius_dense_10001_samples_m=float(dense_min_radius),
        min_radius_refined_m=geom['min_curvature_radius_m'],
        vessel_clearance_m=geom['solid'].distance(deformation.cavity.boundary),
        blanket_clearance_m=geom['solid'].distance(geometry.solids['Blanket']),
        shield_clearance_m=geom['solid'].distance(geometry.solids['Shield']),
        lcfs_clearance_m=geom['solid'].distance(Polygon(lcfs)))
    checks = []
    for count in (320,640):
        row = rows[str(count)]
        assert not row.get('errors')
        assert row['metadata']['reference'] == study['fixed_inputs']['reference']
        assert row['metadata']['count'] == count
        assert row['metadata']['source_RZ_m'] == study['fixed_inputs']['sources'][source]['source_RZ_m']
        for key in ('geometry_sha256','equilibrium_sha256'):
            assert row['metadata'][key] == study[key]
        for label in ('baseline','candidate'):
            path = RUN / f'{label}-{source}-{count}.npz'
            shutil.copy2(path, OUT / path.name)
            f = dict(np.load(path, allow_pickle=False))
            desc = row[label]
            assert desc['feasible'] and not desc['errors']
            assert len(f['source_weights']) == count
            assert float(f['source_power_W']) == 1e6
            for key in ('geometry_sha256','equilibrium_sha256'):
                assert str(f[key]) == study[key]
            assert np.array_equal(f['coefficients'], coefficients if label=='candidate' else np.array([]))
            qkeys = [k for k in f if k.endswith('_q')]
            peak = max(float(np.max(f[k])) for k in qkeys) * float(f['source_power_W'])/1e6
            power = np.concatenate([f[k] for k in f if k.endswith('_cell_power')])
            areas = np.concatenate([f[k] for k in f if k.endswith('_areas')])
            order = np.argsort(-power/areas)
            cumulative = np.cumsum(power[order]); n = np.searchsorted(cumulative, .9*power.sum())
            area90 = areas[order[:n]].sum() + areas[order[n]]*(.9*power.sum()-(cumulative[n-1] if n else 0))/power[order[n]]
            unhit = f['source_weights'][f['hit_loop']<0].sum()
            balance = abs(power.sum()+unhit-f['source_weights'].sum())
            assert abs(peak-desc['peak_MW_m2']) < 1e-12
            assert abs(area90-desc['area_90_m2']) < 1e-11
            assert unhit == 0 and balance < 2e-12
            assert abs(power.sum()-desc['ledger']['deposited']) < 1e-12
            assert abs(desc['divertor_front_power_fraction']+desc['divertor_rear_and_top_cap_power_fraction']-desc['divertor_first_hit_power_fraction']) < 1e-12
            if label=='candidate':
                assert np.allclose(f[f"loop_{int(f['divertor_loop_index'])}"][[0,-1]], geom['wall'][[0,-1]],atol=1e-14)
            checks.append(dict(count=count,case=label,peak_MW_m2=peak,area_90_m2=float(area90),
                deposited_W=float(power.sum()*1e6),unhit_W=float(unhit*1e6),balance_error=float(balance)))
            if count==640: data[source,label] = f
            audit['artifacts'][path.name] = sha(path)
    b0,c0 = [rows['320'][k]['peak_MW_m2'] for k in ('baseline','candidate')]
    b1,c1 = [rows['640'][k]['peak_MW_m2'] for k in ('baseline','candidate')]
    gain=b1-c1; uncertainty=abs(b1-b0)+abs(c1-c0)
    assert gain > uncertainty and source_data['qualification']['qualified_improvement']
    assert abs(uncertainty-source_data['qualification']['estimated_numerical_uncertainty_MW_m2']) < 1e-14
    assert rows['640']['candidate']['area_90_m2'] > rows['640']['baseline']['area_90_m2']
    audit['sources'][source] = dict(geometry=ga, array_checks=checks,qualification=source_data['qualification'])
    (OUT/f'geometry-{source}.json').write_text(json.dumps(dict(source=source,
        source_power_W=1e6,basis=source_data['deformation_basis'],coefficients_m=coefficients.tolist(),
        section_polyline_m=geom['wall'].tolist(), original_geometry_sha256=study['geometry_sha256']),indent=2),encoding='utf-8')

# Compact complete search history: only scalar outcomes, never raw donor data.
history=[]
for directory,offset in [('campaign-3h',0),('resume-remaining-v1',study['prior_consumed_seconds'])]:
    for line in (RAW/directory/'candidates.jsonl').read_text().splitlines():
        r=json.loads(line)
        if r.get('stage')=='qualification': continue
        desc=r['description']
        history.append(dict(source=r['source'],evaluation=r['evaluation'],
            elapsed_seconds=offset+r['campaign_elapsed_seconds'],accepted=r['accepted'],
            feasible=desc['feasible'],peak_MW_m2=desc.get('peak_MW_m2'),
            evaluation_seconds=desc.get('elapsed_seconds'),errors=desc.get('errors',[])))
assert len(history)==study['search_evaluations']
(OUT/'search-history.json').write_text(json.dumps(history,separators=(',',':')),encoding='utf-8')
for source,d in audit['sources'].items():
    rr=[r for r in history if r['source']==source]
    times=[r['evaluation_seconds'] for r in rr if r['feasible'] and r['evaluation_seconds'] is not None]
    final=study['sources'][source]['best_search']['peak_MW_m2']
    near=[r for r in rr if r['accepted'] and r['peak_MW_m2'] <= final*1.001]
    d['search']=dict(evaluations=len(rr),feasible_evaluation_seconds_median=float(np.median(times)),
        feasible_evaluation_seconds_p95=float(np.percentile(times,95)),
        first_within_0p1_percent_of_final_search_peak_seconds=near[0]['elapsed_seconds'],
        candidate_evaluation_errors=sum(bool(r['errors']) and any('Error:' in e for e in r['errors']) for r in rr))
(OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':160})
colors={'baseline':'#526879','candidate':'#d95732'}
def save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=170)
    fig.savefig(OUT/f'{name}.svg')
    plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained')
for ax,source in zip(axes,('IMP','OMP')):
    for label in ('baseline','candidate'):
        f=data[source,label]; wall=f[f"loop_{int(f['divertor_loop_index'])}"]
        ax.plot(*wall.T,color=colors[label],lw=1.5,label=label.title())
    ax.plot(*lcfs.T,color='#55916d',lw=1,label='Frozen LCFS')
    for name in ('Blanket','Shield'):
        ax.plot(*geometry.solids[name].exterior.xy,color='.75',lw=1)
    ax.set(xlim=(1.95,4.65),ylim=(-3.7,-1.4),xlabel='R (m)',ylabel='Z (m)',title=f'{source}: separate 1 MW optimization')
    ax.set_aspect('equal');ax.legend(fontsize=9,loc='lower right')
fig.suptitle('CPD_5 divertor section: original and qualified candidate')
save(fig,'geometry-comparison')

fig,axes=plt.subplots(2,2,figsize=(12,7.5),layout='constrained')
for i,source in enumerate(('IMP','OMP')):
    vmax=study['sources'][source]['qualification_rows']['640']['baseline']['peak_MW_m2']
    for j,label in enumerate(('baseline','candidate')):
        ax=axes[i,j];f=data[source,label];k=int(f['divertor_loop_index']);wall=f[f'loop_{k}']
        points=f[f'loop_{k}_points'];q=f[f'loop_{k}_q']*float(f['source_power_W'])/1e6
        ax.plot(*wall.T,color='.8',lw=.7)
        select=np.flatnonzero(q>vmax*.001)
        select=select[np.argsort(q[select])]
        sc=ax.scatter(*points[select].T,c=q[select],s=6,cmap='inferno',vmin=0,vmax=vmax,rasterized=True)
        desc=study['sources'][source]['qualification_rows']['640'][label]
        ax.plot(*desc['peak_RZ_m'],'x',color='#209ec2',ms=8,mew=1.7)
        limits=((2.86,3.05),(-2.23,-2.00)) if source=='IMP' else ((3.65,3.94),(-2.23,-2.01))
        ax.set(xlim=limits[0],ylim=limits[1],xlabel='R (m)',ylabel='Z (m)',
            title=f"{source} {label} | peak {desc['peak_MW_m2']:.4f} MW/m²")
        ax.set_aspect('equal')
    fig.colorbar(sc,ax=axes[i,:],label='Surface heat flux (MW/m²)',shrink=.8)
fig.suptitle('640 donors | strike-region zoom, identical axes and color scale within each source\nColored points: q > 0.1% of baseline peak; gray: divertor outline; ×: peak')
save(fig,'heat-load-comparison')

fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for j,source in enumerate(('IMP','OMP')):
    ax=axes[0,j]
    for label in ('baseline','candidate'):
        f=data[source,label];power=np.concatenate([f[k] for k in f if k.endswith('_cell_power')]);areas=np.concatenate([f[k] for k in f if k.endswith('_areas')])
        order=np.argsort(-power/areas);a=np.r_[0,np.cumsum(areas[order])];p=np.r_[0,np.cumsum(power[order])]/power.sum()
        ax.plot(a,p,color=colors[label],label=label.title())
    ax.axhline(.9,color='.5',ls='--',lw=1);ax.set(xlim=(0,1.5),ylim=(0,1.01),xlabel='Accumulated receiving area (m²), highest flux first',ylabel='Fraction of deposited power',title=f'{source}: wider footprint');ax.legend()
    ax=axes[1,j];d=study['sources'][source]
    rr=[r for r in history if r['source']==source]
    best=d['baseline']['peak_MW_m2'];xx=[0];yy=[best]
    for r in rr:
        if r['accepted']: best=min(best,r['peak_MW_m2'])
        xx.append(r['elapsed_seconds']/60);yy.append(best)
    ax.plot(xx,yy,color=colors['candidate'],label='Best 160-donor search peak')
    ax.axvline(study['prior_consumed_seconds']/60,color='.5',ls=':',label='Checkpoint resume')
    ax.plot(study['total_consumed_seconds']/60,d['qualification_rows']['640']['candidate']['peak_MW_m2'],'D',color='#2a8096',label='640-donor qualified peak')
    ax.set(xlim=(0,180),xlabel='Cumulative campaign time (min), downtime excluded',ylabel='Peak (MW/m²)',title=f'{source}: search progress');ax.legend(fontsize=8)
fig.suptitle('CPD_5 fixed-input comparison | IMP and OMP each receive 1 MW independently')
save(fig,'footprint-and-runtime')
print(json.dumps(audit,indent=2))
