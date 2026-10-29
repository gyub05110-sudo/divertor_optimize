from fractions import Fraction as F
import importlib.util
import numpy as np
import pytest
from d3d_optimization.shape import NormalDisplacementBSpline

def module():
    assert importlib.util.find_spec('verification.global_certificate_geometry'), 'whole-domain certificate primitives missing'
    from verification import global_certificate_geometry as g
    return g

def test_outward_arithmetic_contains_exact_rational_corner_values():
    g=module()
    for a,b in [((.1,.3),(-.7,.2)),((-1e-200,2e-200),(.2,.9)),((1.,1.),(.1,.1))]:
        for op,exact in [(g.iadd,lambda x,y:x+y),(g.imul,lambda x,y:x*y)]:
            lo,hi=op(a,b)
            for x in a:
                for y in b:assert F(lo)<=exact(F(x),F(y))<=F(hi)

def test_whole_coefficient_box_encloses_real_kernel_and_fixed_endpoints():
    g=module();rng=np.random.default_rng(31)
    lo=np.array([-.1,-.03,.01,-.05]);hi=np.array([.02,.1,.09,.03])
    s=np.linspace(0,1,81)
    bounds=np.array([g.displacement_interval(x,list(zip(lo,hi))) for x in s])
    for c in rng.uniform(lo,hi,(25,4)):
        actual=NormalDisplacementBSpline(c).displacement(s)
        assert np.all(bounds[:,0]<=actual) and np.all(actual<=bounds[:,1])
    np.testing.assert_array_equal(bounds[[0,-1]],0.)

def test_exact_separator_does_not_certify_touching_or_crossing():
    g=module()
    a=g.hull([(0,0),(1,0),(1,1),(0,1)])
    for b in [[(1,0),(2,0),(2,1),(1,1)],[(.5,-1),(.5,2)]]:
        assert g.separator(a,g.hull(b)) is None
    b=g.hull([(1+2**-40,0),(2,0),(2,1),(1+2**-40,1)])
    witness=g.separator(a,b)
    assert witness is not None and g.check_separator(a,b,witness)
    assert not g.check_separator(a,b,[0,1])

def test_exact_first_fixed_hit_selects_nearest_and_prefix_endpoint():
    g=module()
    result=g.first_fixed_intersection([[0,0],[2,0],[3,1]],[(7,[[1,-1],[1,1]]),(8,[[2,-1],[2,1]])])
    assert result['fixed_index']==7
    assert result['prefix'][-1]==(F(1),F(0))

def test_roundoff_envelope_contains_deformed_nodes():
    g=module()
    from d3d_optimization.shape import ReferencePolyline
    ref=ReferencePolyline(np.linspace([1.1,-1.3],[1.0,-1.2],17))
    parameters=ref.node_parameters;positions=ref.sample(parameters);normals=ref.unit_normal(parameters)
    c=[-.07,.09,-.01,.05]
    actual=NormalDisplacementBSpline(c).deform(ref,parameters)
    for i,s in enumerate(parameters):
        envelope=g.node_envelope(positions[i],normals[i],g.displacement_interval(s,[(-.1,.1)]*4))
        assert g.contains(envelope,tuple(F(float(v)) for v in actual[i]))

def test_first_fixed_hit_fails_closed_on_unresolved_collinear_overlap():
    g=module()
    with pytest.raises(ValueError,match='collinear'):
        g.first_fixed_intersection([[0,0],[2,0]],[(0,[[1,0],[3,0]])])

def frozen_example(y=2.):
    return dict(path=[[0,0],[2,0]],fixed=[(7,[[1,-1],[1,1]])],
                coefficient_bound_with_tolerance_m=.1,
                surfaces=[dict(name='test',parameters=[0.,.5,1.],
                               positions=[[0,y],[.5,y],[1,y]],normals=[[0,1]]*3)])

def test_whole_edge_envelopes_fail_closed_if_possible_path_collision():
    from verification.certify_fixed_model import geometry
    assert geometry(frozen_example())[0]['all_moving_edges_excluded']
    result,_=geometry(frozen_example(.01))
    assert not result['all_moving_edges_excluded'] and len(result['failed_edges'])==2

def test_replay_rejects_missing_forged_and_modified_proof(tmp_path):
    import json
    from verification.certify_fixed_model import geometry,digest,replay
    data=frozen_example();result,witnesses=geometry(data)
    report=dict(schema='fixed-polyline-geometry-certificate.v1',frozen=data,frozen_sha256=digest(data),
                geometry=result,separating_witnesses=witnesses)
    path=tmp_path/'certificate.json'
    path.write_text(json.dumps(report),encoding='utf-8')
    assert replay(path)==result
    saved=list(witnesses)
    report['separating_witnesses']=saved[:-1]
    path.write_text(json.dumps(report),encoding='utf-8')
    with pytest.raises(ValueError,match='incomplete'):replay(path)
    report['separating_witnesses']=saved
    saved[0]['axis']=[1.,0.]
    path.write_text(json.dumps(report),encoding='utf-8')
    with pytest.raises(ValueError,match='invalid'):replay(path)
    data['path'][0][1]=.5
    path.write_text(json.dumps(report),encoding='utf-8')
    with pytest.raises(ValueError,match='digest'):replay(path)

def test_empty_geometry_is_not_a_certificate():
    from verification.certify_fixed_model import geometry
    data=frozen_example();data['surfaces']=[]
    with pytest.raises(ValueError,match='surface'):geometry(data)

def test_no_fixed_hit_is_a_failed_obligation_not_a_certificate():
    from verification.certify_fixed_model import geometry
    data=frozen_example();data['fixed']=[(7,[[3,-1],[3,1]])]
    result,witnesses=geometry(data)
    assert not result['all_moving_edges_excluded']
    assert result['reason']=='stored path has no fixed-wall intersection'
    assert witnesses==[]
