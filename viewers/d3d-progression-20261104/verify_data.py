"""Validate exported evidence without optimization or heat solves."""
from export_data import *
from collections import defaultdict
from urllib.request import urlopen

def main():
    d=read(HERE/'dist/data.json');audit=read(HERE/'data-audit.json');requested=defaultdict(dict)
    stages=0;finals=0
    for series in d['series']:
        best=1.;seen=set()
        for row in series['stages']:
            assert 0<row['ratio']<best-1e-9
            best=row['ratio'];assert row['id'] not in seen;seen.add(row['id'])
            assert row['profile'] is None and row['validation']=='geometry'
            requested[row['origin']][row['line']]=(row,series)
            stages+=1
        for row in series['stages']+([series['final']] if series['final'] else []):
            c=Candidate.from_record(row['candidate'])
            assert c.identity==row['id'] and c.source==series['source'] and c.radius_mm==series['rmin']
            assert row['constraints']['passed'] and row['geometry_screen']['passed']
            assert max(row['constraints']['length_ratios'].values())<=1.2+1e-10
            assert row['geometry_screen']['policy']=='full-wall-chord-1e-7.v2'
            assert math.isclose(row['peak']/row['baseline'],row['ratio'],abs_tol=1e-12)
            assert np.all(np.isfinite(row['wall']))
        if series['final']:
            row=series['final'];finals+=1
            assert 0<row['ratio']<1
            p=read(HERE/'dist'/row['profile'])
            assert math.isclose(p['peak'],row['peak'],abs_tol=1e-10)
            assert abs(p['power_fraction']-1)<1e-10
            if row['validation']=='current':
                validation=read(row['origin'])
                assert validation['identity']==row['id'] and validation['complete'] and validation['eligible']
                assert validation['verdict']['resolution_screen_passed'] and validation['verdict']['improvement_confirmed']
    for path,lines in requested.items():
        with Path(path).open('rb') as f:
            for i,line in enumerate(f,1):
                if i not in lines:
                    if i>max(lines):break
                    continue
                exported,series=lines[i];original=json.loads(line)
                assert original['feasible']
                c=Candidate.from_record(original['candidate']) if 'candidate' in original else Candidate(series['source'],series['rmin'],tuple(original['coefficients_m'][:16]))
                assert c.identity==exported['id']
                assert math.isclose(original.get('objective',original.get('peak_ratio')),exported['ratio'],abs_tol=1e-7)
    assert not audit['old_viewer_research_data_used'] and audit['new_thermal_solves']==audit['new_optimizations']==0
    assert audit['inputs_unchanged_verified']
    assert all('divertor-viewer-public' not in p for p in audit['input_sha256'])
    assert d['baseline']['wall']==read(REPO/'verification/fixtures/campaign_geometry_20261029.json')['wall']
    mime={}
    for filename,expected in [('app.mjs','text/javascript'),('vendor/three.module.js','text/javascript'),('data.json','application/json')]:
        with urlopen('http://127.0.0.1:8780/'+filename) as response:
            mime[filename]=response.headers.get_content_type()
            assert mime[filename]==expected
    save(HERE/'verification.json',{'data_checks':'passed','series':len(d['series']),'stages':stages,'finals':finals,
                                  'exact_source_records_checked':stages,'http_mime':mime,
                                  'baseline_preserved':True,'new_thermal_solves':0,'new_optimizations':0})
    print(json.dumps({'series':len(d['series']),'stages':stages,'finals':finals,'checks':'passed'}))

if __name__=='__main__':main()
