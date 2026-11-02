"""Analytic adversarial tests; numerical samples are not proof oracles."""
import copy
import importlib.util
import math
import pytest


def field_model(circle=False,K=3):
    r=[1,1,1,3,3,3];z=[-2,-2,-2,2,2,2]
    rc=[1,3,9] if circle else [1,2,3]
    zc=[4,-4,4] if circle else [0,0,0]
    return dict(r_knots=r,z_knots=z,degree_r=2,degree_z=2,
                coefficients=[[x+y for x in rc] for y in zc],toroidal_K=K)


def test_continuous_module_exists():
    assert importlib.util.find_spec('d3d_optimization.continuous_path_bounds') is not None


def test_constant_helical_length_and_replay():
    from d3d_optimization.continuous_path_bounds import trace,replay
    r=trace(field_model(),[[2,1],[2,0],[2,-1]],grad_sep=1,max_cells=32)
    assert r['complete'] and r['continuous_prefix_proven']
    end=r['cells'][-1]
    assert float(end['end'][0][0])<=2<=float(end['end'][0][1])
    assert float(end['length_end'][0])<=2*math.sqrt(10)<=float(end['length_end'][1])
    assert replay(r)
    bad=copy.deepcopy(r);bad['cells'][-1]['length_end']=['0','0']
    with pytest.raises(ValueError,match='replay'):replay(bad)


def test_curved_contour_and_chart_switch():
    from d3d_optimization.continuous_path_bounds import trace
    r=trace(field_model(True,0),[[2,0],['1.99','.2'],['1.9','.6'],['1.6','1.2'],['1.2','1.6']],1,max_cells=128)
    assert r['complete']
    assert {c['axis'] for c in r['cells']}=={0,1}
    p=r['cells'][-1]['end']
    assert float(p[0][0])<=1.2<=float(p[0][1])
    assert float(p[1][0])<=1.6<=float(p[1][1])
    length=2*math.acos(.6)
    assert float(r['cells'][-1]['length_end'][0])<=length<=float(r['cells'][-1]['length_end'][1])


def test_wrong_branch_and_budget_are_explicit():
    from d3d_optimization.continuous_path_bounds import trace
    r=trace(field_model(),[[2,1],[2,0],[2,1]],1,max_cells=16,max_depth=2)
    assert not r['complete'] and r['unresolved']
    r=trace(field_model(),[[2,1],[2,0],[2,-1]],1,max_cells=1)
    assert not r['complete'] and len(r['cells'])==1


def test_non_C1_and_uncertain_field_rejected():
    from d3d_optimization.continuous_path_bounds import ContinuousField
    m=field_model();m['r_knots']=[1,1,1,2,2,3,3,3]
    m['coefficients']=[[1,2,3,4,5]]*3
    with pytest.raises(ValueError,match='C1'):ContinuousField(m)
    m=field_model();m['coefficients'][0][0]=[1,2]
    with pytest.raises(ValueError,match='exact'):ContinuousField(m)


def test_zero_gradient_refuses_trace():
    from d3d_optimization.continuous_path_bounds import trace
    m=field_model();m['coefficients']=[[0]*3]*3
    with pytest.raises(ValueError):trace(m,[[2,1],[2,0]],1)


def test_variable_terminal_time_and_switch_preserve_incoming_set():
    from d3d_optimization.continuous_path_bounds import ContinuousField,chart_step
    from d3d_optimization.enclosure_math import I
    # psi=R-Z on a larger rectangle. Every incoming state has psi in [1,2].
    m=dict(r_knots=[1,1,5,5],z_knots=[-1,-1,3,3],degree_r=1,degree_z=1,
           coefficients=[[2,6],[-2,2]],toroidal_K=0)
    f=ContinuousField(m)
    a=chart_step(f,[I(1,2),I(0)],I(3),0,1,I(1,2),I(1),I(1))
    assert float(a['time'][0])==1 and float(a['time'][1])==2
    assert float(a['end'][1][0])<=1 and float(a['end'][1][1])>=2
    b=chart_step(f,[I(2),I(0,1)],I(2),1,1,I(1,2),I(1),I(1))
    assert float(b['end'][0][0])<=3 and float(b['end'][0][1])>=4


def test_continuous_certificate_and_cache_tampering():
    assert importlib.util.find_spec('d3d_optimization.continuous_transport_certificate') is not None
    from verification.test_geometry_transport_bounds import slab_model
    from d3d_optimization.continuous_path_bounds import trace
    from d3d_optimization.continuous_transport_certificate import ValidatedPaths,certify_box,replay
    m=slab_model();guide=[[2,1],[2,'.4'],[2,'.1'],[2,0]]
    path=trace(m['field'],guide,m['grad_sep'])
    cache=ValidatedPaths(m,[path]);r=certify_box(m,m['domains'],cache,area_panels=2,thermal_panels=4)
    assert r['continuous_fieldline_enclosure_proven'] and r['transport_parameter_enclosure_proven']
    assert r['kernel'] and replay(r,cache)
    bad=copy.deepcopy(path);bad['cells'][0]['tube']=[['2','2'],['0','1']]
    with pytest.raises(ValueError):ValidatedPaths(m,[bad])
    other=copy.deepcopy(m);other['field']['toroidal_K']=99
    with pytest.raises(ValueError,match='cache'):certify_box(other,other['domains'],cache)


def test_complete_subdivision_no_dropped_unresolved_child():
    from verification.test_geometry_transport_bounds import slab_model
    from d3d_optimization.continuous_path_bounds import trace
    from d3d_optimization.continuous_transport_certificate import ValidatedPaths,cover,replay_cover
    m=slab_model();cache=ValidatedPaths(m,[trace(m['field'],[[2,1],[2,'.8']],1)])
    r=cover(m,cache,splits=2,area_panels=2,thermal_panels=2)
    assert len(r['leaves'])==3 and not r['all_boxes_resolved'] and replay_cover(r,cache)
    bad=copy.deepcopy(r);bad['leaves'].pop()
    with pytest.raises(ValueError):replay_cover(bad,cache)


def test_cache_rejects_opposite_physical_branch():
    from verification.test_geometry_transport_bounds import slab_model
    from d3d_optimization.continuous_path_bounds import trace
    from d3d_optimization.continuous_transport_certificate import ValidatedPaths
    m=slab_model();up=trace(m['field'],[[2,1],[2,'1.6'],[2,'1.9'],[2,2]],1)
    with pytest.raises(ValueError,match='branch'):ValidatedPaths(m,[up])
    down=trace(m['field'],[[2,1],[2,'.7'],[2,'.2'],[2,0]],1)
    assert ValidatedPaths(m,[down])


def test_continuous_tangency_retained_and_domain_exit_not_a_hit():
    from d3d_optimization.continuous_path_bounds import trace
    from d3d_optimization.continuous_transport_certificate import first_hits
    from d3d_optimization.geometry_transport_bounds import Surface
    from verification.test_geometry_transport_bounds import line
    r=trace(field_model(True,0),[['1.6','-1.2'],['1.9','-.6'],[2,0],['1.9','.6'],['1.6','1.2']],1)
    assert r['complete']
    hit=first_hits(r,[Surface(line([2,-.5],[2,.5]),[]),Surface(line([1,1],[3,1]),[])])
    assert hit['guaranteed_hit'] and {p['surface'] for p in hit['possible_hits']}=={0,1}
    missed=first_hits(r,[Surface(line([1.5,0],[1.7,0]),[])])
    assert not missed['guaranteed_hit']


def test_C1_knot_crossing():
    from d3d_optimization.continuous_path_bounds import trace
    rc=[1,1.5,2.5,3];zc=[-2,0,2]
    m=dict(r_knots=[1,1,1,2,3,3,3],z_knots=[-2,-2,-2,2,2,2],degree_r=2,degree_z=2,
           coefficients=[[r-z for r in rc] for z in zc],toroidal_K=1)
    r=trace(m,[[1.5,-.5],[2.5,.5]],1)
    assert r['complete']
    assert float(r['cells'][-1]['end'][1][0])<=.5<=float(r['cells'][-1]['end'][1][1])


def test_independent_path_audit():
    assert importlib.util.find_spec('verification.run_continuous_path_study') is not None
    from verification.run_continuous_path_study import numeric_audit
    from d3d_optimization.continuous_path_bounds import trace
    result=numeric_audit(trace(field_model(),[[2,1],[2,0],[2,-1]],1))
    assert result['passed'] and result['checks']>10 and not result['is_proof']


def test_proved_hit_before_unresolved_tail_and_ambient_precision():
    from verification.test_geometry_transport_bounds import slab_model
    from d3d_optimization.continuous_path_bounds import trace
    from d3d_optimization.continuous_transport_certificate import ValidatedPaths,certify_box
    from decimal import localcontext
    m=slab_model();guide=[[2,1],[2,'.4'],[2,'.1'],[2,0],[2,1]]
    for precision in [6,28,60]:
        with localcontext() as ctx:
            ctx.prec=precision
            p=trace(m['field'],guide,1,max_depth=1)
            assert not p['complete'] and p['unresolved']
            r=certify_box(m,m['domains'],ValidatedPaths(m,[p]),area_panels=2,thermal_panels=4)
            assert r['continuous_fieldline_enclosure_proven']


def test_continuous_mapping_integral_nonconstant_analytic_value():
    from d3d_optimization.continuous_transport_certificate import mapping_average
    from d3d_optimization.enclosure_math import I,record
    from fractions import Fraction as F
    cells=[]
    for j in range(64):
        a=I(j)/64;b=I(j+1)/64
        cells.append({'integral_start':record(a),'integral_end':record(b),
                      'mapping':record(I((1+2*a).lo,(1+2*b).hi))})
    r=mapping_average(cells,I(1),{'target_y':I(1),'energy':I(3)},128)
    assert F(r.lo)<=F(19,9)<=F(r.hi)
    assert float(r.hi-r.lo)<.06


def test_numeric_audit_supports_linear_fields():
    from verification.test_geometry_transport_bounds import slab_model
    from d3d_optimization.continuous_path_bounds import trace
    from verification.run_continuous_path_study import numeric_audit
    m=slab_model()
    assert numeric_audit(trace(m['field'],[[2,1],[2,'.5']],1))['passed']


def test_local_guide_refinement_keeps_source_and_recovers_small_surface():
    from verification.run_continuous_path_study import refined_guide
    from d3d_optimization.continuous_path_bounds import trace
    from d3d_optimization.continuous_transport_certificate import first_hits
    from d3d_optimization.geometry_transport_bounds import Surface
    from verification.test_geometry_transport_bounds import line
    p=trace(field_model(),[[2,1],[2,0]],1)
    wall=[Surface(line([1.99,.5],[2.01,.5]),[])]
    assert not first_hits(p,wall)['guaranteed_hit']
    guide=refined_guide(p,{0},32)
    assert guide[0]==p['inputs']['guide'][0]
    q=trace(field_model(),guide,1)
    assert first_hits(q,wall)['guaranteed_hit']
