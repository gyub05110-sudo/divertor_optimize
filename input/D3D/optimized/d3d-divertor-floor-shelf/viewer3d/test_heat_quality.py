import unittest
import numpy as np
import heat_quality as q

class QualityTests(unittest.TestCase):
    def test_refines_spatial_gap_and_surface_switch(self):
        edges=np.array([0.,.01,.02,.03])
        ids=q.refinement_cells(edges,[0.,.03,.031],['floor','floor','shelf'],np.array([[1.,1.,1.]]),gap=.002)
        self.assertTrue({0,1,2}.issubset(set(ids)))
    def test_no_zero_assigned_to_unreached_point(self):
        self.assertEqual(q.connection_status(False,False,False),'outside_source_scope')
        self.assertEqual(q.connection_status(True,False,False),'trace_unresolved')
        self.assertEqual(q.connection_status(True,True,False),'blocked_first_hit')
        self.assertEqual(q.connection_status(True,True,True),'direct_hit')
    def test_surface_integral_uses_area_not_assigned_weights(self):
        a=np.array([[1.,0.],[1.,1.]])
        self.assertAlmostEqual(q.surface_integral(a,np.array([[2.]]))[0],4*np.pi)
    def test_convergence_rejects_peak_location_jump(self):
        a={'peaks':[1.], 'locations':[[1.,0.]], 'ratios':[.6]}
        b={'peaks':[1.001], 'locations':[[1.01,0.]], 'ratios':[.6001]}
        self.assertFalse(q.compare_levels(a,b)['passed'])
    def test_convergence_accepts_small_value_position_changes(self):
        a={'peaks':[1.], 'locations':[[1.,0.]], 'ratios':[.6]}
        b={'peaks':[1.001], 'locations':[[1.001,0.]], 'ratios':[.6001]}
        self.assertTrue(q.compare_levels(a,b)['passed'])

if __name__=='__main__': unittest.main()
