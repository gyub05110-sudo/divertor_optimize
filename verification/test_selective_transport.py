import copy
import importlib.util
import pytest


def test_module_exists():
    assert importlib.util.find_spec('d3d_optimization.selective_transport') is not None


def test_centered_spline_range_retains_constant_offset_correlation():
    from d3d_optimization.selective_transport import tight_spline_range
    from d3d_optimization.enclosure_math import I,D
    value=tight_spline_range([0,0,1,1],[I(1001),I(1002)],1,I('.4','.6'))
    assert value.lo<=D('1001.4') and value.hi>=D('1001.6')
    assert float(value.hi-value.lo)<.201


def test_tight_trace_contracts_analytic_curved_path_and_replays():
    from d3d_optimization.selective_transport import trace_tight,replay_tight
    from d3d_optimization.continuous_path_bounds import trace
    from verification.test_continuous_path_bounds import field_model
    m=field_model(True,3);guide=[[2,0],['1.99','.2'],['1.9','.6'],['1.6','1.2'],['1.2','1.6']]
    old=trace(m,guide,1);new=trace_tight(m,guide,1)
    assert new['complete'] and replay_tight(new)
    def width(p):
        a,b=map(float,p['cells'][-1]['length_end']);return b-a
    assert width(new)<width(old)
    bad=copy.deepcopy(new);bad['cells'][-1]['length_end']=['0','0']
    with pytest.raises(ValueError,match='replay'):replay_tight(bad)


def test_refinement_keeps_source_and_budget_stop():
    from d3d_optimization.selective_transport import refine_donor,select_cells
    from d3d_optimization.continuous_path_bounds import trace
    from verification.test_continuous_path_bounds import field_model
    p=trace(field_model(),[[2,1],[2,0],[2,-1]],1)
    chosen=select_cells(p,[0],limit=2)
    q=refine_donor(p,chosen,factor=2)
    assert q['inputs']['guide'][0]==p['inputs']['guide'][0]
    assert q['direction']==p['direction'] and q['complete']


def test_refinement_preserves_exact_long_decimal_terminal_plane():
    from d3d_optimization.selective_transport import refine_donor
    from d3d_optimization.continuous_path_bounds import trace
    from verification.test_continuous_path_bounds import field_model
    z='-1.000247237678733913668338573188520967960357666015625'
    p=trace(field_model(),[[2,1],[2,0],[2,z]],1)
    q=refine_donor(p,{1},factor=4)
    axis=p['cells'][-1]['axis']
    assert q['cells'][-1]['axis']==axis
    assert q['cells'][-1]['end'][axis]==p['cells'][-1]['end'][axis]


def test_replacement_replays_only_new_record_and_retains_model_binding():
    from d3d_optimization.continuous_path_bounds import trace
    from d3d_optimization.continuous_transport_certificate import ValidatedPaths,certify_box
    from verification.test_geometry_transport_bounds import slab_model
    from d3d_optimization.selective_transport import trace_tight
    m=slab_model();p=trace(m['field'],[[2,1],[2,'.4'],[2,0]],1)
    cache=ValidatedPaths(m,[p]);q=trace_tight(m['field'],[[2,1],[2,'.4'],[2,0]],1)
    new=cache.replacing(m,0,q)
    assert new._hashes!=cache._hashes
    assert certify_box(m,m['domains'],new,field_enclosure='centered-v1')['kernel']
    wrong=trace_tight(m['field'],[[2,1],[2,2]],1)
    with pytest.raises(ValueError,match='branch'):cache.replacing(m,0,wrong)
    tampered=copy.deepcopy(q);tampered['cells'][0]['integral_end']=['0','0']
    with pytest.raises(ValueError,match='replay'):cache.replacing(m,0,tampered)


def test_selection_targets_large_peak_bound_and_preserves_all_donors():
    from d3d_optimization.selective_transport import next_donor
    rows=[{'index':0,'status':'enclosed','weight_MW':['.9','.9'],'sigma_area_m2':['1','2']},
          {'index':1,'status':'enclosed','weight_MW':['.1','.1'],'sigma_area_m2':['.001','.1']}]
    original=copy.deepcopy(rows)
    assert next_donor(rows,{})==1
    assert next_donor(rows,{1:2},max_visits=2)==0
    assert next_donor(rows,{0:2,1:2},max_visits=2) is None
    assert rows==original


def test_previous_study_manifest_detects_preexisting_weight_tampering(tmp_path):
    import json
    from verification.run_selective_transport_study import verified_previous
    from verification.run_optimality_bounds import digest
    p=tmp_path/'model.json';p.write_text('{"weight":1}')
    settings=tmp_path/'settings.json';settings.write_text('{"protected_sha256":{}}')
    (tmp_path/'artifact_sha256.json').write_text(json.dumps({'model.json':digest(p),'settings.json':digest(settings)}))
    assert str(p) in verified_previous(tmp_path)
    p.write_text('{"weight":2}')
    with pytest.raises(ValueError,match='recorded'):
        verified_previous(tmp_path)
