"""Independent full-wall first-hit and finer-resolution audit of saved candidates."""
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import quad
from shapely.geometry import LineString, Polygon
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.frozen_mapping import first_hits, map_heat
from d3d_optimization.case import separatrix, sha256

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'input/D3D/optimized/independent-production'

def audit(source):
    manifest=json.loads((OUT/source/'manifest.json').read_text(encoding='utf-8'))
    coefficients=np.asarray(manifest['verified_coefficients_m'])
    levels=[]
    for count,step,nodes in [(320,.002,161),(640,.001,321),(1280,.0005,641),(2560,.00025,1281)]:
        start=time.monotonic()
        e,meta=prepare_evaluator(ROOT,source,count,step,nodes,bound=manifest['bound_m'])
        if any(meta[k]!=manifest[k] for k in ['shot_sha256','geometry_sha256','policy','surfaces']):
            raise ValueError('saved input mismatch')
        r=e.evaluate(coefficients)
        row=dict(tubes=count,step_m=step,target_nodes=nodes,geometry_valid=r['geometry_valid'],feasible=r['feasible'])
        if not r['geometry_valid']:levels.append(row);continue
        full=first_hits(e.paths,r['wall'])
        maps=[map_heat(e.field,full,s.weights,s.h,r['wall']) for s in e.sources]
        row.update(full_first_hit_segment_mismatches=int(np.count_nonzero(full.segment!=r['maps'][0].segment)),
                   full_first_hit_point_max_delta_m=float(np.nanmax(np.linalg.norm(full.points-r['maps'][0].points,axis=1))),
                   independent_q_max_delta=float(max(np.max(abs(a.q-b.q)) for a,b in zip(maps,r['maps']))),
                   wall_simple=bool(LineString(r['wall']).is_simple),wall_polygon_valid=bool(Polygon(r['wall']).is_valid),
                   lcfs_intersection=bool(LineString(r['wall']).intersects(Polygon(e.lcfs))),
                   vessel_intersection=bool(LineString(r['wall']).intersects(e.vessel_solids)),
                   unmapped_power=r['missing'],accounting_error=r['power_balance_error'],
                   baseline_peaks=e.baseline_peaks.tolist(),candidate_peaks=[float(max(m.q)) for m in maps],
                   ratios=[float(max(m.q)/b) for m,b in zip(maps,e.baseline_peaks)],angles=r['angles'])
        candidate_fixed=[]
        for i in e.fixed_indices:
            ids=np.flatnonzero(r['original_segment']==i)
            candidate_fixed.append(np.max(abs(r['wall'][ids[0]:ids[0]+2]-e.wall[i:i+2])))
        row['fixed_segment_max_coordinate_change_m']=float(max(candidate_fixed))
        peak_ids=[int(np.argmax(m.q)) for m in maps]
        row['peak_hits']=[dict(point=m.points[i].tolist(),original_segment=meta['split_to_original_segment'][r['original_segment'][m.segment[i]]],
                              movable=bool(any(m.segment[i] in ids for ids in r['surface_segments'].values()))) for m,i in zip(maps,peak_ids)]
        row['source_weight_sums']=[float(s.weights.sum()) for s in e.sources]
        rsep,zup=meta['source_RZ_m'];sg=-1 if source=='IMP' else 1
        errors=[]
        for s in e.sources:
            def density(x):
                p=[[rsep+sg*x,zup]]
                return float(np.exp(-x/s.width)*2*np.pi*abs(e.field.gradient(p)[0,0])/np.linalg.norm(e.field.magnetic(p)[0]))
            integral,_=quad(density,0,s.edges[-1],epsabs=1e-13,epsrel=1e-10,limit=100)
            errors.append(abs(integral-s.normalization)/integral)
        row['source_normalization_relative_errors_vs_adaptive_quad']=errors
        _,xp,psi_sep=separatrix(e.field,e.lcfs,np.asarray(json.loads((ROOT/'baseline_review/baseline_review.json').read_text())['axis_RZ_m']))
        review=json.loads((ROOT/'baseline_review/baseline_review.json').read_text())
        span=abs(psi_sep-float(e.field.values([review['axis_RZ_m']])[0]))
        drift=0.
        for path,hit,distance in zip(e.paths,full.points,full.distance):
            if not np.isfinite(distance):continue
            ps=np.asarray(path.coords);ss=np.r_[0,np.cumsum(np.linalg.norm(np.diff(ps,axis=0),axis=1))]
            sampled=ps[ss<=distance][::max(1,len(ps)//150)]
            sampled=np.vstack([sampled,hit])
            drift=max(drift,float(np.max(abs(e.field.values(sampled)-e.field.values([ps[0]])[0]))/span))
        row['max_psiN_drift_to_first_hit']=drift
        if count==640:
            base=max(row['ratios']);neighbors=[]
            for delta in [.0001,.001]:
                for j in range(e.dimension):
                    for sign in [-1,1]:
                        c=coefficients.copy();c[j]+=sign*delta
                        if max(abs(c))>e.bound:continue
                        result=e.evaluate(c)
                        if not result['feasible']:continue
                        score=max(float(max(m.q)/b) for m,b in zip(result['maps'],e.baseline_peaks))
                        if score<base-1e-8:neighbors.append(dict(coefficient=j,delta_m=sign*delta,ratio=score))
            row['local_probe_better_feasible_candidates']=len(neighbors)
            row['best_local_probe']=min(neighbors,key=lambda x:x['ratio']) if neighbors else None
        row['elapsed_seconds']=time.monotonic()-start
        levels.append(row)
        (OUT/source/'validity_audit.json').write_text(json.dumps(dict(source=source,levels=levels),indent=2,allow_nan=False),encoding='utf-8')
        print(source,count,'ratio',max(row['ratios']),'full-hit mismatch',row['full_first_hit_segment_mismatches'],'peak',row['peak_hits'][-1],flush=True)
    return dict(source=source,coefficients_m=coefficients.tolist(),levels=levels,
                audit_script_sha256=sha256(Path(__file__)),physical_scope='prescribed direct magnetic mapping; no transverse transport; no experimental validation')

if __name__=='__main__':
    summaries=[]
    for source in ['IMP','OMP']:
        report=audit(source)
        (OUT/source/'validity_audit.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
        finest=report['levels'][-1]
        ratio=max(finest['ratios'])
        note=('추가 2560경로 검증: IMP 최대 열유속 감소 0%. 이전 약 1% 개선은 해상도 의존.' if source=='IMP'
              else f'추가 2560경로 검증: 최대 열유속 감소 {100*(1-ratio):.2f}%.')
        summary=dict(source=source,finest_worst_peak_ratio=ratio,finest_reduction_percent=100*(1-ratio),
                     improvement_confirmed=ratio<1-1e-6,validity_note=note,feasible=finest['feasible'],
                     physical_scope=report['physical_scope'],global_optimum_proven=False)
        summaries.append(summary)
        for name in ['manifest.json','viewer_result.json']:
            path=OUT/source/name
            saved=json.loads(path.read_text(encoding='utf-8'));saved['validity_audit_summary']=summary;saved['validity_note']=note
            path.write_text(json.dumps(saved,indent=2,allow_nan=False),encoding='utf-8')
    (OUT/'validity_summary.json').write_text(json.dumps(summaries,indent=2,allow_nan=False),encoding='utf-8')
