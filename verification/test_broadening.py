"""Analytic evidence for the connection-length deposition surrogate."""
import numpy as np
import pytest
from scipy.special import ndtr
from shapely.geometry import LineString, Polygon
from d3d_optimization.broadening import BroadeningConfig, PathLengths, build_deposition
from d3d_optimization.frozen_mapping import first_hits


class LinearField:
    def values(self, p): return np.asarray(p)[:, 0]
    def gradient(self, p): return np.tile([1., 0.], (len(p), 1))
    def magnetic(self, p): return np.tile([0., 1., np.sqrt(24.)], (len(p), 1))


def case(sigma=.1, solid=None, unhit=False):
    field=LinearField()
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    paths=[LineString([[2.,1.],[2.,-.5]])]
    hits=first_hits(paths,wall)
    if unhit:
        hits.segment[:]=-1; hits.distance[:]=np.inf; hits.points[:]=np.nan
    return build_deposition(field,paths,PathLengths(field,paths),hits,wall,
        np.array([[1.,.7]]),[1.,1.],1,None,[] if solid is None else [solid],
        BroadeningConfig(sigma, length_ref_m=5., wall_step_m=.025))


def test_helical_length_and_fractional_hit():
    p=[LineString([[2,0],[2,.3],[2,1.]])]
    lengths=PathLengths(LinearField(),p)
    assert lengths.at([.65])[0] == pytest.approx(3.25)
    assert np.isnan(lengths.at([np.inf])[0])
    with pytest.raises(ValueError): lengths.at([1.1])


def test_length_uses_partial_trapezoid():
    class Variable(LinearField):
        def magnetic(self,p):
            p=np.asarray(p); metric=2+p[:,1]
            return np.column_stack([np.zeros(len(p)),np.ones(len(p)),np.sqrt(metric**2-1)])
    lengths=PathLengths(Variable(),[LineString([[2,0],[2,1]])])
    assert lengths.at([.5])[0] == pytest.approx(1.125)


def test_gaussian_cdf_mass_and_density_area_jacobian():
    r=case()
    bottom=r['cell_parent']==0
    wall=r['wall']; a=wall[:-1,0][bottom]; b=wall[1:,0][bottom]
    np.testing.assert_allclose(r['cell_power'][bottom,0],ndtr((b-2)/.1)-ndtr((a-2)/.1),atol=1e-14)
    k=np.argmin(np.linalg.norm(r['points']-[2,0],axis=1))
    assert r['q'][k,0] == pytest.approx(1/(np.sqrt(2*np.pi)*.1*2*np.pi*2))
    assert r['lengths_m'][0] == pytest.approx(5)
    assert r['sigma_m'][0] == pytest.approx(.1)
    np.testing.assert_allclose(r['ledger']['deposited']+r['ledger']['blocked']+r['ledger']['branch_tail']+r['ledger']['unhit'],[1,.7],atol=1e-12)
    assert np.all(r['cell_power'][r['cell_parent']==2]==0) # same psi, other branch
    assert r['valid']


def test_no_tail_renormalization():
    r=case(.8)
    assert r['ledger']['branch_tail'][0] == pytest.approx(2*ndtr(-1/.8))
    assert r['ledger']['deposited'][0] < .8
    assert not r['valid']


def test_blocked_mass_is_explicit_and_not_redistributed():
    blocker=Polygon([[1.8,.001],[2.2,.001],[2.2,.003],[1.8,.003]])
    r=case(solid=blocker)
    assert r['ledger']['blocked'][0]>.9
    assert not r['valid']
    assert all(np.all(r['ledger'][k]>=0) for k in ['deposited','blocked','branch_tail','unhit'])


def test_unhit_is_accounted_without_infinite_kernel():
    r=case(unhit=True)
    np.testing.assert_array_equal(r['ledger']['unhit'],[1,.7])
    assert np.sum(r['q'])==0
    assert not r['valid']


@pytest.mark.parametrize('kwargs',[{'sigma_ref_m':-1},{'sigma_ref_m':np.nan},{'sigma_ref_m':.1,'length_ref_m':0},{'sigma_ref_m':.1,'exponent':2},{'sigma_ref_m':.1,'wall_step_m':0},{'sigma_ref_m':.1,'anchor_m':np.inf}])
def test_invalid_config(kwargs):
    with pytest.raises(ValueError): BroadeningConfig(**kwargs)


def test_invalid_bp_fails_closed():
    class Bad(LinearField):
        def magnetic(self,p): return np.tile([0,0,1.],(len(p),1))
    with pytest.raises(ValueError): PathLengths(Bad(),[LineString([[2,0],[2,1]])]).at([.5])

def test_single_solid_geometry_supported():
    field=LinearField(); wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    paths=[LineString([[2.,1.],[2.,-.5]])]
    blocker=Polygon([[1.8,.001],[2.2,.001],[2.2,.003],[1.8,.003]])
    r=build_deposition(field,paths,PathLengths(field,paths),first_hits(paths,wall),wall,[[1.]],[1.,1.],1,None,blocker,BroadeningConfig(.1,length_ref_m=5.,wall_step_m=.025))
    assert r['ledger']['blocked'][0]>.9


def test_independent_surface_quadrature_matches_cdf_mass():
    r=case()
    np.testing.assert_allclose(r['quadrature_integrated_power'],r['ledger']['deposited'],atol=2e-7)


def test_genuine_branch_boundary_hit_rejected():
    field=LinearField(); wall=np.array([[1.,0.],[2.,0.],[2.,2.],[1.,2.],[1.,0.]])
    paths=[LineString([[1.5,1.],[2.5,-1.]])]
    with pytest.raises(ValueError,match='ambiguous.*branch'):
        build_deposition(field,paths,PathLengths(field,paths),first_hits(paths,wall),wall,[[1.]],[1.,1.],1,None,[],BroadeningConfig(.1,wall_step_m=.025))

def test_interior_flux_extremum_is_split_into_distinct_branches():
    from d3d_optimization.broadening import _wall_mesh
    class Bowl(LinearField):
        def values(self,p): return (np.asarray(p)[:,0]-2.013)**2
        def gradient(self,p): return np.column_stack([2*(np.asarray(p)[:,0]-2.013),np.zeros(len(p))])
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    mesh,parent,ds,tangent,xi,branch=_wall_mesh(Bowl(),wall,1.,0.,.15)
    bottom=np.flatnonzero(parent==0)
    assert np.min(abs(mesh[:,0]-2.013))<1e-12
    assert len(np.unique(branch[bottom]))==2
    assert np.all(ds>0)


def test_zero_weight_invalid_length_is_excluded():
    field=LinearField(); wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    paths=[LineString([[2.,1.],[2.,-.5]]),LineString([[1.5,1.],[1.5,-.5]])]
    lengths=PathLengths(field,paths); lengths.metric[1][:]=np.nan
    r=build_deposition(field,paths,lengths,first_hits(paths,wall),wall,[[1.],[0.]],[1.,1.],1,None,[],BroadeningConfig(.1,length_ref_m=5.,wall_step_m=.025))
    assert np.isnan(r['lengths_m'][1])
    assert r['valid']


def test_exponential_gaussian_convolution_matches_analytic_density():
    from scipy.stats import exponnorm
    field=LinearField(); wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    edges=np.linspace(0,.6,401); centers=(edges[1:]+edges[:-1])/2
    lam=.1; sigma=.05
    weights=(np.exp(-edges[:-1]/lam)-np.exp(-edges[1:]/lam))/(1-np.exp(-edges[-1]/lam))
    paths=[LineString([[1.5+x,1.],[1.5+x,-.5]]) for x in centers]
    r=build_deposition(field,paths,PathLengths(field,paths),first_hits(paths,wall),wall,weights[:,None],[1.5,1.],1,None,[],BroadeningConfig(sigma,length_ref_m=5.,wall_step_m=.025))
    selected=(r['sample_parent']==0)&(r['points'][:,0]>1.52)&(r['points'][:,0]<2.)
    x=r['points'][selected,0]-1.5
    expected=(exponnorm.pdf(x,lam/sigma,scale=sigma)-np.exp(-.6/lam)*exponnorm.pdf(x-.6,lam/sigma,scale=sigma))/(1-np.exp(-.6/lam))
    actual=r['q'][selected,0]*2*np.pi*r['points'][selected,0]
    np.testing.assert_allclose(actual,expected,rtol=3e-5)


def test_imp_xi_increases_outward_with_negative_source_sign():
    field=LinearField(); wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    paths=[LineString([[2.,1.],[2.,-.5]])]
    r=build_deposition(field,paths,PathLengths(field,paths),first_hits(paths,wall),wall,[[1.]],[2.5,1.],-1,None,[],BroadeningConfig(.1,length_ref_m=5.,wall_step_m=.025))
    assert r['source_xi_m'][0]==pytest.approx(.5)
    assert r['ledger']['deposited'][0]==pytest.approx(1.)

def test_solid_between_visibility_probe_chords_blocks_whole_cell():
    blocker=Polygon([[1.9529,.0024],[1.9533,.0024],[1.9533,.0026],[1.9529,.0026]])
    r=case(solid=blocker)
    selected=(r['sample_parent']==0)&np.isclose(r['points'][:,0],1.90625,atol=1e-12)
    assert np.count_nonzero(selected)==1
    cell=int(r['sample_cell'][selected][0])
    # The quarter-cell receiving ray intersects the material although the
    # previous endpoint/midpoint visibility probes missed this rectangle.
    assert LineString([[2.,.005],[1.90625,0.]]).intersects(blocker)
    assert r['q'][selected,0].item()==0
    assert r['cell_power'][cell,0]==0
    assert r['ledger']['blocked'][0]>.06
    assert not r['valid']

def test_paired_hidden_extrema_in_one_receiving_cell_are_all_found():
    from d3d_optimization.broadening import _wall_mesh
    class Paired(LinearField):
        def values(self,p):
            t=np.asarray(p)[:,0]-1
            return t**3/3-.01*t**2/2+.000024*t
        def gradient(self,p):
            t=np.asarray(p)[:,0]-1
            return np.column_stack([t*t-.01*t+.000024,np.zeros(len(p))])
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    mesh,parent,ds,tangent,xi,branch=_wall_mesh(Paired(),wall,1.,0.,.15)
    bottom=np.flatnonzero(parent==0)
    assert len(np.unique(branch[bottom]))==3
    assert np.min(abs(mesh[:,0]-1.004))<1e-10
    assert np.min(abs(mesh[:,0]-1.006))<1e-10


def test_exact_vertex_length_ignores_metric_beyond_hit():
    lengths=PathLengths(LinearField(),[LineString([[2.,0.],[2.,.5],[2.,1.]])])
    lengths.metric[0][-1]=np.nan
    assert lengths.at([.5])[0]==pytest.approx(2.5)
    assert lengths.at([0.])[0]==0


def test_nonpolynomial_generic_receiving_field_fails_closed():
    from d3d_optimization.broadening import _wall_mesh
    class Oscillatory(LinearField):
        def values(self,p): return np.sin(200*np.asarray(p)[:,0])
        def gradient(self,p): return np.column_stack([200*np.cos(200*np.asarray(p)[:,0]),np.zeros(len(p))])
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    with pytest.raises(ValueError,match='polynomial'):
        _wall_mesh(Oscillatory(),wall,1.,0.,.15)

@pytest.mark.parametrize('method',['cubic','quintic'])
def test_spline_patch_roots_and_knot_boundaries_are_resolved(method):
    from d3d_optimization.broadening import _wall_mesh
    from d3d_optimization.frozen_field import FrozenField
    r=np.linspace(1.,3.,13); z=np.linspace(0.,2.,13)
    t=r[None,:]-1
    psi=np.broadcast_to(t**3/3-.01*t**2/2+.000024*t,(len(z),len(r))).copy()
    field=FrozenField(r,z,psi,1.,method=method)
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    mesh,parent,ds,tangent,xi,branch=_wall_mesh(field,wall,1.,0.,.15)
    bottom=np.flatnonzero(parent==0)
    assert len(np.unique(branch[bottom]))==3
    assert np.min(abs(mesh[:,0]-1.004))<1e-9
    assert np.min(abs(mesh[:,0]-1.006))<1e-9
    for knot in np.unique(field.spline.get_knots()[1]):
        assert np.min(abs(mesh[np.r_[bottom,bottom[-1]+1],0]-knot))<1e-12


def test_resolved_nonzero_spread_approaches_direct_density_away_from_edge():
    field=LinearField()
    wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    edges=np.linspace(0,.6,401)
    centers=(edges[:-1]+edges[1:])/2
    lam=.1
    weights=(np.exp(-edges[:-1]/lam)-np.exp(-edges[1:]/lam))/(1-np.exp(-.6/lam))
    paths=[LineString([[1.5+x,1.],[1.5+x,-.5]]) for x in centers]
    lengths=PathLengths(field,paths)
    hits=first_hits(paths,wall)
    errors=[]
    for sigma in [.02,.01,.005]:
        r=build_deposition(field,paths,lengths,hits,wall,weights[:,None],[1.5,1.],1,None,[],
                           BroadeningConfig(sigma,length_ref_m=5.,wall_step_m=.025))
        ids=np.flatnonzero((r['sample_parent']==0)&(abs(r['points'][:,0]-1.7)<.005))
        j=ids[np.argmin(abs(r['points'][ids,0]-1.7))]
        x=r['points'][j,0]-1.5
        expected=np.exp(-x/lam)/(lam*(1-np.exp(-.6/lam)))
        actual=r['q'][j,0]*2*np.pi*r['points'][j,0]
        errors.append(abs(actual/expected-1))
    assert errors[2]<errors[1]<errors[0]
    assert errors[2]<.002
