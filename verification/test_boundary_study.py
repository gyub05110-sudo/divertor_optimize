import importlib.util
from types import SimpleNamespace
import numpy as np
import pytest


def api():
    assert importlib.util.find_spec('optimize_boundary_transport') is not None
    import optimize_boundary_transport as study
    return study


def test_cli_accepts_fixed_arbitrary_reference_and_rejects_bad_inputs():
    s=api();a=s.parse_args(['--output','unused','--lambda-mm','3.7','--gamma','6','--power-mw','2'])
    assert a.lambda_mm==3.7 and a.gamma==6 and a.power_mw==2
    for option,value in [('--chi','0'),('--power-mw','nan'),('--lambda-mm','-1'),('--tubes','2')]:
        with pytest.raises(SystemExit):s.parse_args(['--output','unused',option,value])


def test_search_keeps_baseline_and_deduplicates_proposals(tmp_path):
    s=api()
    class Fake:
        dimension=2;bound=.1;lengths={'plate':.1}
        baseline={'feasible':True,'peak_q':1.,'penalty':0.}
        baseline_peak=1.;metadata={'source':'OMP'}
        def __init__(self):self.calls=[]
        def evaluate(self,c):
            self.calls.append(tuple(c));return dict(feasible=True,peak_q=1+float(np.dot(c,c)),penalty=0.)
        def describe(self,r):return dict(feasible=r['feasible'],peak_q=r['peak_q'])
    e=Fake();a=SimpleNamespace(seed=17,sobol_power=2,population=5,generations=1)
    result=s.search(e,tmp_path,a,[np.zeros(2)])
    np.testing.assert_array_equal(result['best_coefficients_m'],[0,0])
    assert len(e.calls)==len(set(e.calls))==result['evaluated']
    assert result['global_optimum_proven'] is False


def test_refinement_verdict_rejects_missing_or_drifting_distribution():
    from verification.validate_boundary_transport import verdict
    assert not verdict([])['improvement_confirmed']
    row=dict(audit_passed=True,feasible=True,peak_ratio=.8,peak_MW_m2=2.,peak_RZ_m=[1.,0.],
             region_power={'floor':1.},cdf_wall_fraction=[0.,.5,1.])
    levels=[dict(rows={'baseline':dict(row,peak_ratio=1.,peak_MW_m2=2.5),'selected':dict(row)}) for _ in range(4)]
    wall=dict(max_peak_relative_change=0.,rows={'baseline':dict(row),'selected':dict(row)})
    assert verdict(levels,wall)['improvement_confirmed']
    wall['rows']['selected']['feasible']=False
    assert not verdict(levels,wall)['improvement_confirmed']
    wall['rows']['selected']['feasible']=True
    levels[2]['rows']['selected']['cdf_wall_fraction']=[0.,.6,1.]
    assert not verdict(levels,wall)['improvement_confirmed']
