from pathlib import Path
import numpy as np
import pytest
from shapely.geometry import LineString
from d3d_optimization.broadening import BroadeningConfig, PathLengths, build_deposition
from d3d_optimization.frozen_mapping import first_hits
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.transport_deposition import deposit_sigma
from d3d_optimization.transport_evaluator import TransportEvaluator

ROOT=Path(__file__).resolve().parents[1]


def test_explicit_width_kernel_reproduces_reviewed_original():
    from verification.test_broadening import LinearField
    field=LinearField();wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    paths=[LineString([[2.,1.],[2.,-.5]])];lengths=PathLengths(field,paths)
    hits=first_hits(paths,wall);config=BroadeningConfig(.1,length_ref_m=5,wall_step_m=.025)
    args=(field,paths,lengths,hits,wall,[[1.]],[1.,1.],1,None,[],config)
    old=build_deposition(*args);new=deposit_sigma(*args,sigma_m=[.1])
    for key in ['q','cell_power','sigma_m','lengths_m','quadrature_integrated_power']:
        np.testing.assert_array_equal(new[key],old[key])
    for key in old['ledger']:np.testing.assert_array_equal(new['ledger'][key],old['ledger'][key])
    with pytest.raises(ValueError):deposit_sigma(*args,sigma_m=[-1.])


@pytest.fixture(scope='module')
def direct():
    return prepare_evaluator(ROOT,'IMP',count=48,nodes=41,step=.0005,boundary_mode='continuous')


def test_direct_adapter_preserves_nominal_original_and_source(direct):
    d,meta=direct;original=d.baseline['peak_q'].copy()
    e=TransportEvaluator(d,meta,ReferencePoint(),mode='direct')
    assert e.baseline['peak_q'][0]==original[2]
    assert e.baseline['feasible']==d.baseline['feasible']
    assert len(e.sources)==1
    assert len(d.sources)==5
    np.testing.assert_array_equal(d.baseline['peak_q'],original)


def test_thermal_and_frozen_agree_at_baseline_and_preserve_input_power(direct):
    d,meta=direct;ref=ReferencePoint()
    t=TransportEvaluator(d,meta,ref,mode='thermal',wall_step_m=.003)
    f=TransportEvaluator(d,meta,ref,mode='frozen_transport',wall_step_m=.003,thermal_paths=t.thermal_paths)
    np.testing.assert_array_equal(t.baseline['peak_q'],f.baseline['peak_q'])
    ledger=t.baseline['deposition']['ledger']
    np.testing.assert_allclose(ledger['prescribed'],[1.],atol=1e-14)
    np.testing.assert_allclose(ledger['deposited']+ledger['unresolved'],[1.],atol=1e-12)
    c=np.zeros(d.dimension);c[-4:]=.001
    rt,rf=t.evaluate(c),f.evaluate(c)
    assert rt['geometry_valid'] and rf['geometry_valid']
    np.testing.assert_array_equal(rf['deposition']['sigma_m'],f.baseline['deposition']['sigma_m'])
    assert np.nanmax(abs(rt['deposition']['sigma_m']-t.baseline['deposition']['sigma_m']))>1e-12
    np.testing.assert_array_equal(t.sources[0].weights,d.sources[2].weights)


def test_unsupported_mode_rejected(direct):
    d,meta=direct
    with pytest.raises(ValueError):TransportEvaluator(d,meta,ReferencePoint(),mode='unknown')


def test_nominal_scenario_not_rejected_by_unused_source_width(direct,monkeypatch):
    d,meta=direct;original=d.evaluate
    def unused_width_failed(c):
        r=dict(original(c))
        if r['geometry_valid']:
            r.update(feasible=False,missing=.8,power_balance_error=.2,penalty=800)
        return r
    monkeypatch.setattr(d,'evaluate',unused_width_failed)
    e=TransportEvaluator(d,meta,ReferencePoint(),mode='direct')
    assert e.baseline['feasible']
    assert e.baseline['missing']==e.baseline['maps'][0].unmapped_power
    assert e.baseline['penalty']==0
