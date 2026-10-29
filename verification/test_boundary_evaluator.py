import importlib.util,json
from pathlib import Path
import numpy as np
import pytest
from d3d_optimization.reference_transport import ReferencePoint

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def evaluator():
    assert importlib.util.find_spec('d3d_optimization.boundary_evaluator') is not None
    from d3d_optimization.boundary_evaluator import prepare_boundary
    return prepare_boundary(ROOT,'OMP',ReferencePoint(lambda_q_m=.0037),count=48,nodes=81,wall_step=.002)


def test_fixed_omp_baseline_conserves_at_arbitrary_source_width(evaluator):
    e=evaluator;r=e.baseline
    assert e.source.width==.0037
    assert r['geometry_valid'] and r['feasible']
    assert r['deposition']['ledger']['balance_error']<1e-12
    assert r['deposition']['ledger']['unresolved']==0
    np.testing.assert_array_equal(r['wall'],e.direct.baseline['wall'])


def test_refined_length_failure_is_caught_before_thermal_calculation(evaluator):
    e=evaluator
    m=json.loads((ROOT/'input/D3D/optimized/reference-comparison-20261028/OMP/manifest.json').read_text())
    r=e.evaluate(m['candidates']['direct_best'])
    assert not r['feasible']
    assert 'material_length_exceeds_limit' in r['errors']
    assert r['precise_length_ratios']['floor']>1.2


def test_repeated_candidates_reuse_cache_and_do_not_mutate_baseline(evaluator):
    e=evaluator;wall=e.baseline['wall'].copy();c=np.zeros(e.dimension);c[-1]=.0003
    first=e.evaluate(c);second=e.evaluate(c)
    assert first is second
    np.testing.assert_array_equal(e.baseline['wall'],wall)
    assert first['deposition']['diagnostics']['missing_power_renormalized'] is False
