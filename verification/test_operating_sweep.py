from dataclasses import replace
import numpy as np
import pytest
from shapely.geometry import LineString
from d3d_optimization.reference_transport import ReferencePoint, ThermalPaths
from d3d_optimization.frozen_mapping import Hits
from d3d_optimization.operating_sweep import reuse_thermal_paths, reuse_boundary
from verification.test_reference_transport import SlabField, source

def test_field_path_reuse_recomputes_heat_and_matches_a_fresh_solver():
    field=SlabField();paths=[LineString([[1.,0.],[3.,0.]])]
    ref=ReferencePoint()
    initial=ThermalPaths(field,paths,[1.,0.],source(),ref)
    hits=Hits(np.array([0]),np.array([[2.5,0.]]),np.array([1.5]))
    before=initial.evaluate(hits)
    after=replace(ref,power_W=2e6,chi_perp_m2_s=.6,density_m3=1.7e19,gamma=6.)
    shared=reuse_thermal_paths(initial,after)
    assert shared.data is initial.data and shared.qu is not initial.qu
    fresh=ThermalPaths(field,paths,[1.,0.],source(),after)
    for key,value in shared.evaluate(hits).items():np.testing.assert_array_equal(value,fresh.evaluate(hits)[key])
    assert not np.array_equal(before['sigma_m'],shared.evaluate(hits)['sigma_m'])
    np.testing.assert_array_equal(before['sigma_m'],initial.evaluate(hits)['sigma_m'])

def test_refuses_to_reuse_a_different_source_width_or_a_candidate_evaluator():
    ref=ReferencePoint()
    initial=ThermalPaths(SlabField(),[LineString([[1.,0.],[3.,0.]])],[1.,0.],source(),ref)
    with pytest.raises(ValueError,match='width'):reuse_thermal_paths(initial,replace(ref,lambda_q_m=.004))
    with pytest.raises(TypeError):reuse_boundary(object(),ref)

