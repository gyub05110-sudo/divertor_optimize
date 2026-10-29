import numpy as np
import pytest
from types import SimpleNamespace


def test_regions_expand_without_changing_opposite_coefficients():
    from d3d_optimization.regional_transport import RegionalEvaluator
    names=['omp_entry','shelf','omp_connector','floor','imp_sloped','imp_vertical']
    fake=SimpleNamespace(surfaces=[{'name':n} for n in names],dimension=24,metadata={},lengths={n:1. for n in names})
    for source,indices in [('IMP',np.arange(16,24)),('OMP',np.arange(16))]:
        e=RegionalEvaluator(fake,source);c=np.arange(e.dimension)*.001;full=e.expand(c)
        np.testing.assert_array_equal(full[indices],c)
        np.testing.assert_array_equal(full[np.setdiff1d(np.arange(24),indices)],0.)
        np.testing.assert_array_equal(e.project(full),c)
        with pytest.raises(ValueError):e.expand(np.zeros(e.dimension+1))


def test_gpu_density_matches_cpu_narrow_broad_and_seam():
    from d3d_optimization.gpu_footprint import density_sum
    from d3d_optimization.boundary_transport import wrapped_density
    a=np.r_[np.linspace(0,1,1001),1e-10];mu=np.array([1e-10,.2,.7,.9]);sigma=np.array([1e-7,.01,.12,10.]);w=np.array([.1,.2,.3,.4])
    expected=sum(p*wrapped_density(a,c,s,1.) for c,s,p in zip(mu,sigma,w))
    actual=density_sum(a,mu,sigma,w,1.,backend='gpu')
    np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=1e-12)


def test_gpu_deposition_matches_preserved_cpu():
    from d3d_optimization.accelerated_transport import deposit_accelerated
    from d3d_optimization.boundary_transport import deposit_boundary
    from verification.test_broadening import LinearField
    from d3d_optimization.frozen_mapping import Hits
    wall=np.array([[1.,0],[3,0],[3,2],[1,2],[1,0]])
    h=Hits(np.array([0,1,-1]),np.array([[2.,0],[3,1],[np.nan,np.nan]]),np.array([1.,1.,np.inf]))
    args=(LinearField(),h,wall,[.4,.5,.1],[.03,.3,np.nan],1.,.05)
    a=deposit_boundary(*args);b=deposit_accelerated(*args)
    np.testing.assert_allclose(b['q'],a['q'],rtol=2e-12,atol=1e-12)
    np.testing.assert_array_equal(b['cell_power'],a['cell_power'])
    assert b['ledger']==a['ledger']
