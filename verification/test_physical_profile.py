import numpy as np
import pytest
from d3d_optimization.physical_profile import read_physical_profile

def artifact(path, recorded=None):
    wall=np.array([[1.,0.],[2.,0.],[2.,1.],[1.,0.]])
    area=np.pi*(wall[:-1,0]+wall[1:,0])*np.linalg.norm(np.diff(wall,axis=0),axis=1)
    p=area/area.sum()
    kw={} if recorded is None else {'assigned_power_W':recorded}
    np.savez(path,receiving_wall=wall,cell_power=p,area_edges=np.r_[0.,np.cumsum(area)],
             points=wall[:-1],q_MW_m2=np.full(3,2/area.sum()),**kw)

def test_recorded_two_megawatts_has_matching_cell_and_continuous_units(tmp_path):
    path=tmp_path/'two.npz';artifact(path,2e6)
    d=read_physical_profile(path)
    np.testing.assert_allclose(d['cell_q_MW_m2'],d['continuous_q_MW_m2'])
    assert np.isclose(sum(d['cell_power_MW']),2)
    with pytest.raises(ValueError,match='mismatch'):read_physical_profile(path,1e6)

def test_legacy_input_needs_its_actual_operating_power(tmp_path):
    path=tmp_path/'old.npz';artifact(path)
    with pytest.raises(ValueError,match='operating power'):read_physical_profile(path)
    assert read_physical_profile(path,2e6)['assigned_power_W']==2e6

