"""Counterexamples and analytic checks for geometry-to-transport inclusion."""
from decimal import Decimal as D
import importlib.util
import math
import pytest


def test_enclosure_module_exists():
    assert importlib.util.find_spec('d3d_optimization.enclosure_math') is not None


def test_primitives():
    from d3d_optimization.enclosure_math import I, positive_root, sin_cos, spline_range
    for x in ['1e-60', '2', '1e80']:
        r=positive_root(I(x),7)
        from fractions import Fraction as F
        assert F(r.lo)**7 <= F(x) <= F(r.hi)**7
    s,c=sin_cos(I('1.56','1.58'))
    assert s.lo<=1<=s.hi
    assert c.lo<0<c.hi
    t=[0,0,0,0,1,1,1,1]
    for a,b in [(0,1),('.2','.3')]:
        value=spline_range(t,[I(0),I(1),I(2),I(3)],3,I(a,b))
        assert value.lo<=D(str(a))*3 and value.hi>=D(str(b))*3


def test_field_polynomial():
    from d3d_optimization.enclosure_math import I, SplineField
    # psi(Z,R)=R^2+3Z. Bernstein controls on R=[1,2], Z=[0,1].
    f=SplineField({'r_knots':[1,1,1,2,2,2], 'z_knots':[0,0,1,1],
        'degree_r':2,'degree_z':1,'coefficients':[[1,2,4],[4,5,7]],'toroidal_K':2})
    g=f.gradient([I('1.2','1.4'),I('.2','.3')])
    assert g[0].lo<=D('2.4') and g[0].hi>=D('2.8')
    assert g[1].lo<=3<=g[1].hi


def test_nonfinite_rejected():
    from d3d_optimization.enclosure_math import interval
    for value in [float('nan'),float('inf'),'-Infinity']:
        with pytest.raises(ValueError):interval(value)


def line(a,b,controls=None):
    dx=b[0]-a[0];dy=b[1]-a[1]
    return {'a':a,'b':b,'normal':[-dy,dx], 'knots':[0,0,1,1],
            'controls':controls or [0,0], 'degree':1}


def test_first_hit_switch_and_tangency_coverage():
    from d3d_optimization.geometry_transport_bounds import Surface, first_hits
    from d3d_optimization.enclosure_math import I
    surfaces=[Surface(line([1,0],[3,0],[{'variable':0}]*2),[I('-.2','.2')]),
              Surface(line([1,-.1],[3,-.1]),[])]
    r=first_hits([[2,1],[2,.3],[2,-.3],[2,-1]],surfaces,depth=7)
    assert r['guaranteed_hit']
    assert {h['surface'] for h in r['possible_hits']}=={0,1}
    # A tangent root before the guaranteed straight-wall hit cannot disappear.
    tangent={'a':[1,0],'b':[3,0],'normal':[0,1], 'knots':[0]*3+[1]*3,
             'controls':[1,-1,1],'degree':2}
    r=first_hits([[1,0],[3,0],[3,-1]], [Surface(tangent,[]),Surface(line([2,-.5],[4,-.5]),[])])
    assert any(h['surface']==0 for h in r['possible_hits'])


def test_sign_change_outside_material_domain_is_not_hit():
    from d3d_optimization.geometry_transport_bounds import Surface, first_hits
    r=first_hits([[1.5,1],[3,1],[3,-1],[1.5,-1]], [Surface(line([1,0],[2,0]),[])])
    assert not r['guaranteed_hit']
    assert not r['possible_hits']


def test_shared_endpoint_and_earlier_fixed_surface():
    from d3d_optimization.geometry_transport_bounds import Surface, first_hits
    surfaces=[Surface(line([1,0],[2,0]),[]),Surface(line([2,0],[3,0]),[]),
              Surface(line([1,-.5],[3,-.5]),[])]
    r=first_hits([[2,1],[2,-1]],surfaces)
    assert {h['surface'] for h in r['possible_hits']} >= {0,1}
    assert all(h['surface']!=2 for h in r['possible_hits'])


def test_multiple_controls_and_area():
    from d3d_optimization.geometry_transport_bounds import Surface
    from d3d_optimization.enclosure_math import I
    s=Surface(line([1,0],[1,1],[{'variable':0},{'variable':1}]),[I('-.1','.1'),I('-.2','.2')])
    p=s.point(I('.5'))
    assert p[0].lo<=D('.85') and p[0].hi>=D('1.15')
    q=Surface(line([1,0],[1,1]),[])
    area=q.area(I(1))
    assert area.lo < D(str(2*math.pi)) < area.hi or float(area.lo)==pytest.approx(2*math.pi)


def test_metric_prefix_toroidal_and_unused_invalid_tail():
    from d3d_optimization.thermal_box_bounds import PathThermal
    from d3d_optimization.enclosure_math import I
    m={'path':[[1,0],[1,1],[1,2],[1,3]],'metric':[5,5,5,'NaN'],
       'integrand':[5,5,5,'NaN'],'mapping':[2,2,2,'NaN']}
    p=PathThermal(m,1)
    L,J,avg=p.at(1,I('.2','.6'))
    assert L.lo<=6 and L.hi>=8 and J.lo<=6 and J.hi>=8
    assert avg.lo==avg.hi==2
    with pytest.raises(ValueError):PathThermal(m,2)


def test_thermal_root_bounds_and_sigma():
    from d3d_optimization.thermal_box_bounds import thermal_width, sheath_root
    from d3d_optimization.enclosure_math import I,power
    y=sheath_root(I(3),I(4)) # y^3+3y=4 -> y=1
    assert y.lo<=1<=y.hi
    reference={'density_m3':1e19,'chi_perp_m2_s':.3,'gamma':7,'kappa0':2000}
    result=thermal_width(I(1e7),I(1),I(5,6),I(2),reference)
    assert result['sigma_m'].lo>0
    # Independent production root is diagnostic, not the certificate oracle.
    from d3d_optimization.reference_transport import solve_two_point,ReferencePoint,E_CHARGE
    for j in [5,5.5,6]:
        q=solve_two_point(1e7,1,j,ReferencePoint())
        sigma=math.sqrt(3*E_CHARGE*1e19*q['Tu_eV']*.3/2000)*j/(q['Tu_eV']**1.75+q['Tt_eV']**1.75)*2
        assert float(result['sigma_m'].lo)<=sigma<=float(result['sigma_m'].hi)


def test_field_zero_gradient_blocks_width():
    from d3d_optimization.enclosure_math import SplineField,I
    f=SplineField({'r_knots':[1,1,2,2],'z_knots':[0,0,1,1],
        'degree_r':1,'degree_z':1,'coefficients':[[0,0],[0,0]],'toroidal_K':2})
    with pytest.raises(ValueError,match='gradient'):f.quantities([I('1.5'),I('.5')])


def slab_model():
    surfaces=[line([1,0],[3,0]),line([3,0],[3,2]),line([3,2],[1,2]),line([1,2],[1,0])]
    return {'schema':'geometry-transport-input.v1','domains':[['-.1','.1']],
        'surfaces':surfaces,'field':{'r_knots':[1,1,3,3],'z_knots':[0,0,2,2],
        'degree_r':1,'degree_z':1,'coefficients':[[1,3],[1,3]],'toroidal_K':2},
        'grad_sep':1,'reference':{'density_m3':1e19,'chi_perp_m2_s':.3,'gamma':7,'kappa0':2000},
        'donors':[{'path':[[2,1],[2,0]],'metric':[5,5],'integrand':[5,5],
        'mapping':[2,2],'b0':1,'qu':1e7,'weight_MW':1}]}


def test_certificate_bridge_and_replay():
    from d3d_optimization.geometry_transport_certificate import certify_box,replay
    from d3d_optimization.transport_peak_bounds import KernelBounds
    import copy
    r=certify_box(slab_model(),[['-.01','.01']],root_depth=4,area_panels=2)
    assert r['transport_parameter_enclosure_proven']
    assert not r['continuous_fieldline_enclosure_proven']
    assert not r['geometry_global_optimum_proven']
    KernelBounds(r['kernel'])
    assert replay(r)
    bad=copy.deepcopy(r);bad['donors'][0]['length_m'][0]='0'
    with pytest.raises(ValueError):replay(bad)
    missing=slab_model();missing['donors'][0]['path']=[[2,1],[2,.5]]
    r=certify_box(missing,[['-.01','.01']],root_depth=2,area_panels=2)
    assert not r['transport_parameter_enclosure_proven'] and r['kernel'] is None


def test_domain_cover_and_tamper():
    from d3d_optimization.geometry_transport_certificate import cover,replay_cover
    import copy
    r=cover(slab_model(),splits=1,root_depth=2,area_panels=2)
    assert len(r['leaves'])==2 and replay_cover(r)
    bad=copy.deepcopy(r);bad['leaves'].pop()
    with pytest.raises(ValueError):replay_cover(bad)


def test_junction_is_not_frozen_and_all_input_domains_checked():
    from d3d_optimization.geometry_transport_bounds import Surface
    from d3d_optimization.geometry_transport_certificate import certify_box
    from d3d_optimization.enclosure_math import I
    model={'a':[1,0],'b':[2,0],'normal':[0,1], 'degree':3,
        'knots':[0]*4+[1]*4,'controls':[0,0,{'junction':
            {'variable':8,'angle_low':0,'angle_high':1,'sign':-1}},0]}
    q=Surface(model,[I(0)]*8+[I('.2','.4')])
    assert float(q.c[2].lo)<=-math.tan(.4)/15
    assert float(q.c[2].hi)>=-math.tan(.2)/15
    with pytest.raises(ValueError,match='outside declared'):
        certify_box(slab_model(),[['-.2','.2']])


def test_actual_adapter_module_exists():
    assert importlib.util.find_spec('verification.run_geometry_transport_bounds') is not None


def test_trig_outward_negation_ignores_ambient_context():
    from d3d_optimization.enclosure_math import I,sin_cos
    from decimal import localcontext
    for precision in [6,28,60]:
        with localcontext() as c:
            c.prec=precision
            s,_=sin_cos(I('0','2.2222222222222222222222222222e-20'))
            assert s.lo<=0


def test_unclamped_wall_cannot_claim_closed_geometry():
    from d3d_optimization.geometry_transport_certificate import certify_box
    m=slab_model();s=m['surfaces'][0]
    s.update(degree=2,knots=[-1,-1,0,1,2,2],controls=[0,'.1',0])
    with pytest.raises(ValueError,match='clamped'):certify_box(m,m['domains'])


def test_uncertain_knots_are_rejected_not_collapsed():
    from d3d_optimization.enclosure_math import I,spline_range,SplineField
    with pytest.raises(ValueError,match='exact'):
        spline_range([0,0,[1,2],[1,2]],[I(0),I(1)],1,I('.5'))
    m=slab_model()['field'];m['r_knots']=[1,1,[3,4],[3,4]]
    with pytest.raises(ValueError,match='exact'):SplineField(m)


def test_invalid_cache_marker_is_unresolved_not_crash():
    from d3d_optimization.geometry_transport_certificate import certify_box
    m=slab_model();m['donors'][0]['metric'][-1]=None
    r=certify_box(m,m['domains'],root_depth=2,area_panels=2)
    assert r['unresolved'] and r['kernel'] is None


def test_universal_kernel_bounds_cover_broad_width_family():
    from d3d_optimization.geometry_transport_certificate import universal_peak_bounds
    m={'period_m2':['2','3'],'donors':[{'weight_MW':['1','1'],
        'center_fraction':['0','1'],'sigma_m2':['.01','1e10']}]}
    r=universal_peak_bounds(m)
    assert float(r[0])<=1/3 and float(r[1])>1/(math.sqrt(2*math.pi)*.01)


def test_transformed_integral_encloses_analytic_nonconstant_mapping():
    from d3d_optimization.thermal_box_bounds import PathThermal
    from d3d_optimization.enclosure_math import I
    p=PathThermal({'path':[[1,0],[1,1]],'metric':[1,1],
                  'integrand':[1,1],'mapping':[1,3]},0)
    # J=1, y=1, E=3 -> yu=2, d/S=1/3. r(s)=1+2s.
    # Integral r((1-u)(1+u/3)) du = 1+2*(1/2+1/18)=19/9.
    b=p.average(I(1),I(1),I(3),64)
    from fractions import Fraction as F
    assert F(b.lo)<=F(19,9)<=F(b.hi)
    assert b.hi-b.lo<D('.1')


def test_straight_area_closed_form_is_tight():
    from d3d_optimization.geometry_transport_bounds import Surface
    m=line([1,0],[2,1]);m['exact_straight_area']=True
    q=Surface(m,[],panels=1)
    assert q.total_area.hi-q.total_area.lo<D('1e-30')


def test_spline_enclosure_spans_interior_knot_peak():
    from d3d_optimization.enclosure_math import I,spline_range
    r=spline_range([0,0,'.5',1,1],[I(0),I(1),I(0)],1,I('.25','.75'))
    assert r.lo<=D('.5') and r.hi>=1


def test_sheath_extreme_scales_and_parameter_corners():
    from d3d_optimization.thermal_box_bounds import sheath_root
    from d3d_optimization.enclosure_math import I
    from fractions import Fraction as F
    for e in ['1e-80','3','1e80']:
        for rhs in ['1e-120','4','1e120']:
            r=sheath_root(I(e),I(rhs))
            assert F(r.lo)**3+F(e)*F(r.lo)<=F(rhs)
            assert F(r.hi)**3+F(e)*F(r.hi)>=F(rhs)
    r=sheath_root(I('1','10'),I('2','100'))
    assert F(r.lo)**3+10*F(r.lo)<=2
    assert F(r.hi)**3+F(r.hi)>=100
