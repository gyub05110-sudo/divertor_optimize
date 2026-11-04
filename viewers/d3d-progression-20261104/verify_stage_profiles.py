"""Independent inventory, source-array and toroidal area checks for every displayed stage."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from shapely.geometry import LineString

HERE=Path(__file__).resolve().parent;DIST=HERE/'dist'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()
def run():
    started=time.perf_counter();data=read(DIST/'data.json');index=read(DIST/'stage-profiles.json')
    audit=read(HERE/'stage-evaluations/audit.json');saved=read(DIST/'full-saved-profiles.json')
    assert index['complete'] and audit['complete'] and audit['protected_unchanged']
    assert index['data_sha256']==saved['data_sha256']==sha(DIST/'data.json')
    assert all(audit['parity'][s]['passed'] for s in ['IMP','OMP'])
    rows={};stage_count=0;max_peak_error=0.;max_power_error=0.;cell_count=0;max_shape_error=0.
    for series in data['series']:
        for row in series['stages']:
            entry=index['cases'][row['id']];assert entry['source']==series['source'] and entry['donors']==320
            assert row['validation']=='geometry' and row['profile'] is None
            max_peak_error=max(max_peak_error,abs(entry['peak']/row['peak']-1));stage_count+=1
            rows[row['id']]=row
    assert len(rows)==len(index['cases']) and max_peak_error<1e-7
    sources=[]
    for identity,row in rows.items():
        entry=index['cases'][identity];meta=read(HERE/'stage-evaluations'/identity/'metadata.json')
        assert meta['identity']==identity and meta['candidate']==row['candidate']
        assert meta['fingerprint']==index['fingerprint']==audit['fingerprint']
        assert meta['description']['feasible'] and meta['description']['geometry_screen']['passed']
        path=HERE/'stage-evaluations'/identity/'selected.npz'
        assert sha(path)==meta['npz_sha256']
        sources.append((entry['profile'],entry['sha256'],path,row['wall'],row['peak'],identity))
    for series in data['series']:
        row=series['final']
        if row:
            e=saved['profiles'][row['profile']]
            assert e['identity']==row['id'] and e['source']==series['source']
            sources.append((e['profile'],e['sha256'],Path(e['input']),row['wall'],row['peak'],row['id']))
    for source in ['IMP','OMP']:
        e=saved['profiles'][data['baseline']['profiles'][source]]
        assert e['source']==source
        sources.append((e['profile'],e['sha256'],Path(e['input']),data['baseline']['wall'],data['baseline']['refined'][source],None))
    for rel,h,path,wall,peak,identity in sources:
        p=read(DIST/rel);assert sha(DIST/rel)==h and p['identity']==identity
        segments=np.array(p['segments']);assert np.all(np.isfinite(segments)) and np.min(segments[:,4])>=0
        # Independent surface-of-revolution formula, not exported area_edges.
        geom_area=np.pi*(segments[:,0]+segments[:,2])*np.linalg.norm(segments[:,2:4]-segments[:,:2],axis=1)
        np.testing.assert_allclose(geom_area,p['areas_m2'],rtol=2e-8,atol=1e-12)
        power=float(geom_area@segments[:,4])*1e6/data['reference']['power_W']
        max_power_error=max(max_power_error,abs(power-1));assert abs(power-1)<2e-10
        with np.load(path) as raw:
            z={key:raw[key] for key in ['wall','q_MW_m2','cell_power','cell_parent','receiving_wall','area_edges']}
            delta=LineString(wall).hausdorff_distance(LineString(z['wall']))
            max_shape_error=max(max_shape_error,delta);assert delta<=2e-5
            np.testing.assert_allclose(p['peak'],max(z['q_MW_m2']),rtol=1e-13)
            assert abs(p['peak']/peak-1)<1e-7
            cell_count+=len(segments)
            # Check each bin against the underlying original cell power independently.
            start=0;starts=[];ends=[]
            for segment in segments:
                end=start+1
                while end<min(start+4,len(z['cell_power'])) and z['cell_parent'][end]==z['cell_parent'][start]:end+=1
                starts.append(start);ends.append(end);start=end
            assert start==len(z['cell_power'])
            np.testing.assert_allclose(segments[:,:2],z['receiving_wall'][starts],rtol=0,atol=1e-13)
            np.testing.assert_allclose(segments[:,2:4],z['receiving_wall'][ends],rtol=0,atol=1e-13)
            expected=np.add.reduceat(z['cell_power'],starts)*data['reference']['power_W']/1e6/np.add.reduceat(np.diff(z['area_edges']),starts)
            np.testing.assert_allclose(segments[:,4],expected,rtol=1e-12,atol=1e-12)
    for path,h in audit['protected_sha256'].items():assert sha(Path(path))==h,path
    for e in saved['profiles'].values():assert sha(Path(e['input']))==e['input_sha256']
    result=dict(passed=True,stage_occurrences=stage_count,unique_stage_profiles=len(rows),saved_profiles=8,
        checked_display_cells=cell_count,max_recorded_peak_relative_error=max_peak_error,
        max_power_fraction_error=max_power_error,max_geometry_distance_m=max_shape_error,
        cpu_gpu_parity=audit['parity'],original_files_unchanged=True,final_qualifications_unchanged=True,
        data_sha256=index['data_sha256'],stage_index_sha256=sha(DIST/'stage-profiles.json'),
        full_saved_index_sha256=sha(DIST/'full-saved-profiles.json'),seconds=time.perf_counter()-started)
    (HERE/'stage-profiles-verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
