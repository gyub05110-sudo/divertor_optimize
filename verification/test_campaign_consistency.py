"""Archived geometry regressions; no field tracing or transport solves."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from shapely.geometry import shape
from d3d_optimization import campaign_geometry as geometry
from d3d_optimization.tangent_direct import TangentMultiSurface
from d3d_optimization.tangent_evaluator import TangentEvaluator
from d3d_optimization.imp_tangent import IMPDirect, IMPEvaluator


def geometry_example(name='OMP/R020', nodes=161):
    data=json.loads((Path(__file__).parent/'fixtures/campaign_geometry_20261029.json').read_text())
    candidate=geometry.Candidate.from_record(data['candidates'][name])
    direct=object.__new__(IMPDirect if candidate.source=='IMP' else TangentMultiSurface)
    direct.wall=np.asarray(data['wall']);direct.surfaces=data['surfaces']
    direct.targets={s['name']:np.linspace(*direct.wall[s['segment_index']:s['segment_index']+2],nodes) for s in direct.surfaces}
    direct.junction_fraction=candidate.junction_fraction
    e=object.__new__(IMPEvaluator if candidate.source=='IMP' else TangentEvaluator)
    e.direct=SimpleNamespace(base=direct)
    e.wall=direct.wall;e.surfaces=direct.surfaces;e.targets=direct.targets
    e.lcfs=np.asarray(data['lcfs']);e.vessel_solids=shape(data['vessel'])
    e.metadata={'source':candidate.source};e.rmin=candidate.radius_mm/1000;e.dimension=24;e.bound=.1
    e.cache={};e.lip=np.array([1.37256,-1.3247])
    s=next(s for s in e.surfaces if s['name']=='floor');i=s['segment_index']
    e.left_direction=e.wall[i+2]-e.wall[i+1]
    from d3d_optimization.tangent_geometry import Curve
    e.gap_min=Curve(*e.targets['floor'][[0,-1]],np.zeros(4)).minimum_distance(e.lip)[0]
    return e,candidate


def test_archived_vessel_collision_is_rejected_before_heat_evaluation():
    e,candidate=geometry_example()
    c=candidate.bind(e)
    assert e.constraints(c)['passed']
    assert geometry.geometry_check(e,c,None)  # Historical coarse false acceptance.
    assert not geometry.geometry_check(e,c,1e-7)
    def forbidden(_):
        pytest.fail('thermal evaluation must not start for a fine-wall collision')
    e.evaluate=forbidden
    result=geometry.evaluate_candidate(e,candidate)
    assert not result['geometry_valid']
    assert 'vessel_intersection' in result['geometry_screen']['errors']


@pytest.mark.parametrize('name',['OMP/R100','IMP/R020','IMP/R050'])
def test_shared_screen_preserves_valid_archived_shapes_across_meshes(name):
    for nodes in [81,161,321]:
        e,candidate=geometry_example(name,nodes)
        assert geometry.geometry_screen(e,candidate.bind(e))['passed']


def test_search_and_independent_audit_use_same_geometry_policy():
    from d3d_optimization.campaign_validation import geometric_audit
    e,candidate=geometry_example()
    c=candidate.bind(e)
    targets={n:q.point(np.linspace(0,1,len(e.targets[n]))) for n,q in e.curves(c).items()}
    e.baseline={'targets':e.targets}
    report=geometric_audit(e,candidate,{'targets':targets})
    assert not report['passed_independent']
    assert report['geometry_screen']==geometry.geometry_screen(e,c)


def test_actual_search_score_retains_collision_reason_without_thermal_work():
    from collections import OrderedDict
    from d3d_optimization.campaign_search import SearchCase
    e,candidate=geometry_example()
    case=object.__new__(SearchCase);case.source='OMP';case.mm=20;case.e=e
    case.memo=OrderedDict();case.state={'parity':None}
    e.evaluate=lambda c:pytest.fail('search must reject before heat evaluation')
    row=case.score(candidate.vector,'regression')
    assert not row['feasible'] and not row['physics_evaluated']
    assert row['geometry_screen']['errors']==['vessel_intersection']
    assert case.score(candidate.vector,'again')['cached']
