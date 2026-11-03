import importlib
import unittest
import tempfile
from pathlib import Path
import numpy as np

class EquilibriumContract(unittest.TestCase):
    def module(self):
        try:
            return importlib.import_module('cpd_optimization.equilibrium')
        except ModuleNotFoundError:
            self.fail('CPD equilibrium implementation is required')

    def test_reference_units_and_beta(self):
        m = self.module()
        s = m.reference_inputs()
        self.assertAlmostEqual(s['Ip_A'], 1388124.0684837103)
        self.assertAlmostEqual(s['Bt_T'], 2e-7*144*123784.60159301758/1.715993046760559)
        self.assertAlmostEqual(m.normalized_beta(10000.,1.1,s),100*(2*m.MU0*10000./s['Bt_T']**2)*1.1*s['Bt_T']/1.3881240684837103)

    def test_frozen_field_grid_orientation_and_derivatives(self):
        m = self.module()
        r,z=np.linspace(2,6,17),np.linspace(-3,3,25)
        rr,zz=np.meshgrid(r,z,indexing='xy')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'f.npz'
            np.savez(p,R=r,Z=z,psi=3*rr**2+2*zz**2,fvac=8.,psi_axis=0.,psi_lcfs=1.)
            f=m.FrozenField(p)
            br,bz,bt=f.field(4.,1.)
            self.assertAlmostEqual(br,-1.,places=8)
            self.assertAlmostEqual(bz,6.,places=8)
            self.assertAlmostEqual(bt,2.,places=8)
            np.testing.assert_allclose(f.gradient([[4.,1.]]),[[24.,4.]])
            np.testing.assert_allclose(f.magnetic([[4.,1.]]),[[-1.,6.,2.]])
            np.testing.assert_allclose(f.values([[4.,1.]]),[50.])
            self.assertAlmostEqual(f.incidence([[4.,1.]],[[1.,0.]])[0],np.degrees(np.arcsin(1/np.sqrt(41))))
            np.testing.assert_array_equal(f.r,r)
            np.testing.assert_array_equal(f.z,z)
            with self.assertRaises(ValueError): f.field(1.,0.)
            with self.assertRaises(ValueError): f.psi[0,0]=4

    def test_distributed_solenoid_is_total_ampere_turns(self):
        m=self.module()
        from freegs.machine import Solenoid
        c=m.DistributedSolenoid(.9,-4.,4.,quadrature_points=51)
        legacy=Solenoid(.9,-4.,4.,51)
        self.assertAlmostEqual(c.controlPsi(4.,0.)*51,legacy.controlPsi(4.,0.))
        self.assertEqual(c.turns,1)
        c.current=1234.
        self.assertAlmostEqual(c.psi(4.,0.),1234.*c.controlPsi(4.,0.))

    def test_reference_is_verified_against_parquet(self):
        m=self.module()
        p=Path(__file__).resolve().parents[1]/'input/D3D/d3d/d3d_shot_1a3155370c.parquet'
        provenance=m.verify_reference_source(p)
        self.assertTrue(provenance['values_verified'])
        self.assertEqual(provenance['efit_index'],60)
        self.assertEqual(len(provenance['sha256']),64)

    def test_validation_rejects_topology_even_with_exact_operating_match(self):
        m=self.module()
        good=dict(Ip_relative_error=0.,betaN_relative_error=0.,lower_single_null=True,lcfs_valid_polygon=True,xpoint_target_distance_m=.0001,material_overlap_area_m2={'Blanket':0.},material_clearance_m={'Blanket':.1},frozen_field_parity={'p95_Bp_relative_error_away_from_null':1e-5})
        self.assertEqual(m.equilibrium_failures(good),[])
        self.assertIn('lower_single_null',m.equilibrium_failures(dict(good,lower_single_null=False)))
        self.assertIn('material_overlap',m.equilibrium_failures(dict(good,material_overlap_area_m2={'Blanket':1e-5})))
        self.assertIn('betaN',m.equilibrium_failures(dict(good,betaN_relative_error=.01)))
        self.assertIn('material_checks_missing',m.equilibrium_failures({k:v for k,v in good.items() if not k.startswith('material_')}))
if __name__=='__main__': unittest.main()
