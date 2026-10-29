import importlib.util
import numpy as np
import pytest
from scipy.integrate import quad


def api():
    assert importlib.util.find_spec('d3d_optimization.boundary_transport') is not None
    from d3d_optimization import boundary_transport
    return boundary_transport


@pytest.mark.parametrize('sigma',[.001,.1,.119,.12,.4,2.,1e6])
def test_wrapped_mass_density_and_uniform_limit(sigma):
    b=api();edges=np.linspace(0,1,83);mu=.002
    mass=b.wrapped_mass(edges[:-1],edges[1:],mu,sigma,1.)
    assert mass.min()>=0
    assert abs(mass.sum()-1)<2e-13
    # Independent adaptive quadrature of density, with narrow peak split points.
    value=quad(lambda a:float(b.wrapped_density(np.array([a]),mu,sigma,1)[0]),0,1,
               points=[.002,.003,.01,.5,.99,.999],epsabs=1e-10,limit=300)[0]
    assert abs(value-1)<1e-9
    if sigma>2:np.testing.assert_allclose(mass,np.diff(edges),atol=1e-14)


def test_exact_area_partial_hit_inverse_and_orientation():
    b=api();wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    mesh=b.AreaBoundary(wall,.13)
    assert np.all(np.diff(mesh.edges)>0)
    assert abs(mesh.total-32*np.pi)<1e-12
    ah=mesh.hit_area(np.array([[2.,0.]]),np.array([0]))[0]
    assert abs(ah-3*np.pi)<1e-12
    points,cells=mesh.points_at(np.array([ah,mesh.total-1e-8]))
    np.testing.assert_allclose(points[0],[2,0],atol=1e-13)
    target=np.array([[2.5,0.]])
    a=mesh.hit_area(target,[0])[0]
    for other,seg in [(np.roll(wall[:-1],2,axis=0),2),(wall[:-1][::-1],2)]:
        other=np.vstack([other,other[0]]);m=b.AreaBoundary(other,.09)
        h=m.hit_area([[2,0]],[seg])[0];x=m.hit_area(target,[seg])[0]
        np.testing.assert_allclose(b.wrapped_density([a],ah,.5,mesh.total),b.wrapped_density([x],h,.5,m.total),rtol=1e-12)


def test_boundary_rejects_open_crossed_and_disconnected_input():
    b=api()
    for wall in [[[1,0],[2,0],[2,1]],[[1,0],[2,1],[1,1],[2,0],[1,0]]]:
        with pytest.raises(ValueError):b.AreaBoundary(wall,.1)
    with pytest.raises(ValueError):b.wrapped_mass([0],[1],.5,0,1)


def test_concave_boundary_does_not_create_spatial_chord_shortcut():
    b=api();wall=np.array([[1,0],[4,0],[4,3],[3,3],[3,1],[2,1],[2,3],[1,3],[1,0]],float)
    mesh=b.AreaBoundary(wall,.1)
    a=mesh.hit_area([[3,2]],[3])[0];c=mesh.hit_area([[2,2]],[5])[0]
    # Opposite lips are 1m apart in space but 3m along the connected wall.
    assert min(abs(c-a),mesh.total-abs(c-a))>2*np.pi*2*2.9
    assert b.wrapped_density([c],a,.1,mesh.total)[0]<1e-20


def test_deposition_complete_and_unhit_never_renormalized():
    b=api();assert hasattr(b,'deposit_boundary')
    from verification.test_broadening import LinearField
    from d3d_optimization.frozen_mapping import Hits
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    hits=Hits(np.array([0,-1]),np.array([[2.,0.],[np.nan,np.nan]]),np.array([1.,np.inf]))
    r=b.deposit_boundary(LinearField(),hits,wall,[.7,.3],[.2,np.nan],1.,.01)
    assert abs(r['ledger']['deposited']-.7)<1e-12
    assert r['ledger']['unhit']==.3
    assert r['cell_power'].min()>=0
    assert abs(r['quadrature_integrated_power']-.7)<1e-8


def test_broad_limit_uniform_heat_and_no_flux_flat_wall_collapse():
    b=api();assert hasattr(b,'deposit_boundary')
    from verification.test_broadening import LinearField
    from d3d_optimization.frozen_mapping import Hits
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    # psi=R, so both vertical walls are exactly isoflux. They still receive area mass.
    hits=Hits(np.array([1]),np.array([[3.,1.]]),np.array([1.]))
    r=b.deposit_boundary(LinearField(),hits,wall,[1.],[1e8],1.,.1)
    np.testing.assert_allclose(r['q'],1/(32*np.pi),rtol=1e-14)
    np.testing.assert_allclose(r['cell_power'],r['areas']/(32*np.pi),rtol=1e-12)
    assert r['diagnostics']['broad_kernel_power_fraction']==1


@pytest.mark.parametrize('sigma',[1e-7,1e-9])
def test_narrow_cyclic_seam_closes_without_renormalizing(sigma):
    b=api();mass=b.wrapped_mass([0,.5],[.5,1],1e-10,sigma,1.)
    assert abs(mass.sum()-1)<2e-13


def test_continuous_mixture_peak_is_found_between_donors():
    b=api()
    from verification.test_broadening import LinearField
    from d3d_optimization.frozen_mapping import Hits
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    s=1e-5;hits=Hits(np.array([1,1]),np.array([[3.,1.],[3.,1+.5*s]]),np.ones(2))
    r=b.deposit_boundary(LinearField(),hits,wall,[.5,.5],[s,s],1.,.1)
    exact=np.exp(-.5*.25**2)/(np.sqrt(2*np.pi)*2*np.pi*3*s)
    assert abs(r['q'].max()/exact-1)<1e-8
    # Distance estimate must disclose its cell-size uncertainty; far-tail CDF is exact.
    assert r['diagnostics']['wall_distance_estimate_error_bound_m']>=.049
    assert r['diagnostics']['power_beyond_0p1m_wall_distance']<1e-12


def test_overflowing_density_rejected_even_when_mass_conserves():
    b=api()
    from verification.test_broadening import LinearField
    from d3d_optimization.frozen_mapping import Hits
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    hits=Hits(np.array([0]),np.array([[2.,0.]]),np.ones(1))
    with pytest.raises((ValueError,FloatingPointError)):
        b.deposit_boundary(LinearField(),hits,wall,[1.],[1e-320],1.,.1)
