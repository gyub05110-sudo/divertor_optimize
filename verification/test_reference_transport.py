from dataclasses import replace
import numpy as np
import pytest
from shapely.geometry import LineString
from d3d_optimization.reference_transport import ReferencePoint, ThermalPaths, solve_two_point, E_CHARGE, M_DEUTERIUM
from d3d_optimization.frozen_field import HeatSource
from d3d_optimization.frozen_mapping import Hits


class ConstantField:
    def magnetic(self, p):
        return np.tile([1., 0., np.sqrt(24.)], (len(p), 1))
    def gradient(self, p):
        return np.column_stack([np.zeros(len(p)), np.asarray(p)[:, 0]])
    def values(self, p):
        return np.zeros(len(p))


class SlabField(ConstantField):
    def gradient(self, p):
        return np.tile([1., 0.], (len(p), 1))


def source():
    return HeatSource(np.array([0.]), np.array([1.]), np.array([.2]), np.array([0., .003]), .003, 1.)


def test_two_point_satisfies_independent_energy_sheath_and_pressure():
    ref = ReferencePoint()
    out = solve_two_point(1e8, .8, 25., ref)
    tu,tt,nt=out['Tu_eV'],out['Tt_eV'],out['nt_m3']
    assert tu > tt > 0
    assert tu**3.5-tt**3.5 == pytest.approx(3.5*1e8*25/ref.kappa0, rel=1e-10)
    assert 2*nt*tt == pytest.approx(ref.density_m3*tu, rel=1e-10)
    assert ref.gamma*E_CHARGE*nt*tt*np.sqrt(2*E_CHARGE*tt/M_DEUTERIUM) == pytest.approx(.8e8,rel=1e-10)


def test_constant_field_matches_analytic_spreading_and_partial_hit():
    ref=ReferencePoint()
    # Gradient is independent of R: general metric g=1, B/Bp=5.
    paths=[LineString([[1.,0.],[2.,0.],[3.,0.]])]
    solver=ThermalPaths(SlabField(),paths,[1.,0.],source(),ref)
    hits=Hits(np.array([0]),np.array([[2.2,0.]]),np.array([1.2]))
    out=solver.evaluate(hits)
    assert out['length_m'][0] == pytest.approx(6.)
    tu,tt=out['Tu_eV'][0],out['Tt_eV'][0]
    expected=np.sqrt(1.5*E_CHARGE*ref.density_m3*tu*ref.chi_perp_m2_s/ref.kappa0)*12/(tu**1.75+tt**1.75)
    assert out['sigma_m'][0]*np.sqrt(2) == pytest.approx(expected,rel=1e-7)
    assert out['conduction_relative_residual'][0] < 1e-10
    assert out['sheath_relative_residual'][0] < 1e-10


def test_width_reacts_to_reference_power_and_sqrt_chi_without_clipping():
    field=SlabField();paths=[LineString([[1.,0.],[3.,0.]])]
    hits=Hits(np.array([0]),np.array([[3.,0.]]),np.array([2.]))
    ref=ReferencePoint()
    base=ThermalPaths(field,paths,[1.,0.],source(),ref).evaluate(hits)
    chi=ThermalPaths(field,paths,[1.,0.],source(),replace(ref,chi_perp_m2_s=4*ref.chi_perp_m2_s)).evaluate(hits)
    power=ThermalPaths(field,paths,[1.,0.],source(),replace(ref,power_W=2*ref.power_W)).evaluate(hits)
    assert chi['sigma_m'] == pytest.approx(2*base['sigma_m'])
    assert power['Tt_eV'][0] > base['Tt_eV'][0]
    assert power['sigma_m'][0] < base['sigma_m'][0]


@pytest.mark.parametrize('key,value',[('power_W',0),('density_m3',-1),('chi_perp_m2_s',float('nan')),('gamma',0),('kappa0',float('inf'))])
def test_invalid_reference_rejected(key,value):
    with pytest.raises(ValueError):
        replace(ReferencePoint(),**{key:value})


def test_cold_tail_root_no_temperature_floor():
    out=solve_two_point(1e-4,1.,30.,ReferencePoint())
    assert 0 < out['Tt_eV'] < 1
    assert out['sheath_relative_residual'] < 1e-9


def test_unhit_donor_keeps_nan_state_instead_of_fabricating_temperature():
    solver=ThermalPaths(SlabField(),[LineString([[1.,0.],[3.,0.]])],[1.,0.],source(),ReferencePoint())
    out=solver.evaluate(Hits(np.array([-1]),np.array([[np.nan,np.nan]]),np.array([np.inf])))
    assert np.isnan(out['sigma_m'][0])


def test_variable_field_conserves_q_over_B_and_matches_analytic_energy_integral():
    class Variable(SlabField):
        def magnetic(self,p):
            b=2+np.asarray(p)[:,1]
            return np.column_stack([np.zeros(len(p)),np.ones(len(p)),np.sqrt(b*b-1)])
    z=np.linspace(0,1,1001)
    path=LineString(np.column_stack([np.ones(len(z)),z]))
    solver=ThermalPaths(Variable(),[path],[1.,0.],source(),ReferencePoint())
    out=solver.evaluate(Hits(np.array([0]),np.array([[1.,.7]]),np.array([.7])))
    assert out['qt_W_m2'][0]/out['qu_W_m2'][0] == pytest.approx(2.7/2)
    assert out['integral_Bratio_m'][0] == pytest.approx(2*.7+.7**2+.7**3/6,rel=5e-8)
    assert out['length_m'][0] == pytest.approx(2*.7+.7**2/2,rel=1e-12)


def test_very_cold_tail_spreading_resolves_integrable_endpoint_without_floor():
    ref=ReferencePoint(power_W=1e-8)
    path=LineString([[1.,0.],[3.,0.]])
    out=ThermalPaths(SlabField(),[path],[1.,0.],source(),ref).evaluate(
        Hits(np.array([0]),np.array([[3.,0.]]),np.array([2.])))
    tu,tt=out['Tu_eV'][0],out['Tt_eV'][0]
    expected=np.sqrt(1.5*E_CHARGE*ref.density_m3*tu*ref.chi_perp_m2_s/ref.kappa0)*20/(tu**1.75+tt**1.75)/np.sqrt(2)
    assert out['sigma_m'][0]==pytest.approx(expected,rel=1e-10)
    from d3d_optimization.reference_transport import M_ELECTRON
    exact_bulk_target=ref.gamma/2*np.sqrt(M_ELECTRON/M_DEUTERIUM)
    assert out['max_free_streaming_ratio'][0]>=exact_bulk_target*(1-1e-10)
