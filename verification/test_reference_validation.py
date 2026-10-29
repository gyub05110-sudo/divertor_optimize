from pathlib import Path
import numpy as np
import pytest
from shapely.geometry import LineString
from d3d_optimization.frozen_mapping import Hits
from d3d_optimization.reference_transport import ReferencePoint,ThermalPaths
from verification.test_reference_transport import SlabField,source
from verification.validate_reference_transport import independent_thermal,validate_source


def test_independent_quadrature_matches_constant_field():
    field=SlabField();paths=[LineString([[1.,0.],[2.,0.],[3.,0.]])]
    hit=Hits(np.array([0]),np.array([[2.2,0.]]),np.array([1.2]))
    ref=ReferencePoint();solver=ThermalPaths(field,paths,[1.,0.],source(),ref)
    state=solver.evaluate(hit)
    independent=independent_thermal(field,paths,hit,[1.,0.],solver.qu,ref)
    for key in ['sigma_m','length_m','integral_Bratio_m']:
        np.testing.assert_allclose(independent[key],state[key],rtol=1e-8)


def test_validation_refuses_partial_existing_evidence(tmp_path):
    (tmp_path/'validation.json').write_text('{"complete":false}')
    with pytest.raises(FileExistsError):validate_source(tmp_path)
    assert (tmp_path/'validation.json').read_text()=='{"complete":false}'
