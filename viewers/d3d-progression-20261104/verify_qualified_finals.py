"""Independent stored-array and finite-ranking checks. No new thermal calculation."""
from pathlib import Path
import json,math
import numpy as np
from shapely.geometry import LineString
from export_qualified_finals import HERE,STUDY,verified_study
from qualify_omp_stages import read,save,sha,checked_read,QualifiedStore,context,passes_improvement
from qualification_selection import ranked_improvements


def verify():
    report=verified_study();index=read(HERE/'dist/qualified-finals.json');data=read(HERE/'dist/data.json')
    assert index['data_sha256']==sha(HERE/'dist/data.json')==report['config']['data_sha256']
    rows=[]
    for entry in index['entries']:
        final=entry['final'];mm=final['candidate']['radius_mm'];identity=final['id']
        assert sha(HERE/'dist'/entry['profile']['profile'])==entry['profile']['sha256']
        assert sha(HERE/'dist'/entry['evidence']['path'])==entry['evidence']['sha256']
        profile=read(HERE/'dist'/entry['profile']['profile']);manifest=read(HERE/'dist'/entry['evidence']['path'])
        assert passes_improvement(manifest['validation'])
        ranked=ranked_improvements(report['ranking'][str(mm)]);rank=next(i for i,r in enumerate(ranked) if r['identity']==identity)
        attempts=report['attempts'][str(mm)]
        assert [a['identity'] for a in attempts]==[r['identity'] for r in ranked[:rank+1]]
        assert all(not a['passed'] for a in attempts[:-1]) and attempts[-1]['passed']
        for a in attempts:
            v=checked_read(STUDY/f'R{mm:03d}'/'validation'/a['identity']/'validation.json')
            assert passes_improvement(v)==a['passed']
        best=QualifiedStore(STUDY/f'R{mm:03d}',context(report['config'],'OMP',mm)).best()
        assert best['candidate']==manifest['candidate']
        npz=STUDY/f'R{mm:03d}'/manifest['artifact_bundle']/'selected_refined.npz'
        assert profile['npz_sha256']==sha(npz)
        with np.load(npz) as z:
            wall=z['receiving_wall'];area=np.diff(z['area_edges']);cell=z['cell_power'];q=z['q_MW_m2']
            geometric_area=np.pi*(wall[:-1,0]+wall[1:,0])*np.linalg.norm(np.diff(wall,axis=0),axis=1)
            np.testing.assert_allclose(area,geometric_area,rtol=3e-9,atol=3e-13)
            assert abs(cell.sum()-1)<2e-12
            assert math.isclose(float(q.max()),final['peak'],rel_tol=2e-11)
            hausdorff=LineString(final['wall']).hausdorff_distance(LineString(z['wall']))
            assert hausdorff<=2e-5
            seg=np.array(profile['segments']);shown_area=np.pi*(seg[:,0]+seg[:,2])*np.linalg.norm(seg[:,2:4]-seg[:,:2],axis=1)
            np.testing.assert_allclose(shown_area,profile['areas_m2'],rtol=3e-9,atol=3e-13)
            display_power=float(np.sum(shown_area*seg[:,4])*1e6/data['reference']['power_W'])
            assert abs(display_power-1)<2e-12
            # Re-aggregate exact stored cells into the independently reconstructed display mesh.
            cumulative=np.r_[0,np.cumsum(cell)];mesh_arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(wall,axis=0),axis=1))]
            cursor=0
            for s,a in zip(seg,profile['areas_m2']):
                assert np.linalg.norm(wall[cursor]-s[:2])<1e-10
                end=cursor+1
                while end<len(wall) and np.linalg.norm(wall[end]-s[2:4])>=1e-10:end+=1
                assert end<len(wall)
                expected=float(cell[cursor:end].sum()*data['reference']['power_W']/1e6/a)
                assert math.isclose(float(s[4]),expected,rel_tol=2e-11,abs_tol=1e-13)
                cursor=end
            assert cursor==len(cell)
        geom=manifest['validation']['levels'][2]['geometry_audit']
        rows.append(dict(radius_mm=mm,identity=identity,archived_stage=final['archived_stage'],
            ranked_candidates=len(report['ranking'][str(mm)]),improving_candidates=len(ranked),full_qualification_attempts=len(attempts),
            selected_fine_rank=rank+1,peak_MW_m2=final['peak'],baseline_MW_m2=final['baseline'],reduction_percent=100*(1-final['ratio']),
            raw_power_fraction=float(cell.sum()),display_power_fraction=display_power,
            profile_wall_hausdorff_m=hausdorff,display_cells=len(seg),geometry=geom,
            verdict=manifest['validation']['verdict'],cpu_gpu=manifest['cpu_gpu_distribution'],
            cpu_gpu_peak_relative_difference=manifest['cpu_gpu_peak_relative_difference'],
            levels=[dict(method=l['method'],tubes=l['tubes'],selected_peak=l['rows']['selected']['peak_MW_m2'],baseline_peak=l['rows']['baseline']['peak_MW_m2']) for l in manifest['validation']['levels']]))
    assert len(rows)==sum(r['status']=='qualified_improvement' for r in report['results'].values())
    assert all(s['final'] is None for s in data['series'] if s['run']=='campaign' and s['source']=='OMP')
    result=dict(passed=True,checks='Committed study, finite ranking, full qualification, candidate identity, raw/display toroidal area and power, geometry alignment',
        cases=rows,protected_files_unchanged=report['preservation']['unchanged'],global_optimum_proven=False,
        summary_sha256=sha(STUDY/'summary.json'),overlay_sha256=sha(HERE/'dist/qualified-finals.json'))
    save(HERE/'qualified-final-verification.json',result)
    print(json.dumps({**result,'cases':[{k:r[k] for k in ('radius_mm','archived_stage','selected_fine_rank','peak_MW_m2','reduction_percent','raw_power_fraction','display_power_fraction')} for r in rows]},indent=2))
    return result


if __name__=='__main__':verify()
