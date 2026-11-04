import unittest
import numpy as np
from stage_profile_data import make_profile

class ProfileChecks(unittest.TestCase):
    def arrays(self):
        w=np.array([[1.,-1.2],[1.,-1.3],[1.,-1.4]])
        return dict(wall=w,receiving_wall=w,area_edges=np.array([0.,1.,2.]),
            cell_power=np.array([.4,.6]),cell_parent=np.array([0,1]),
            points=np.array([[1.,-1.35]]),q_MW_m2=np.array([.7]))
    def test_power_and_coordinates_are_preserved(self):
        z=self.arrays();p=make_profile(z,z['wall'],.7,1e6)
        self.assertEqual(p['segments'],[[1.,-1.2,1.,-1.3,.4],[1.,-1.3,1.,-1.4,.6]])
        self.assertAlmostEqual(p['power_fraction'],1.)
    def test_wrong_geometry_or_recorded_peak_rejected(self):
        z=self.arrays()
        with self.assertRaises(ValueError):make_profile(z,z['wall']+.01,.7,1e6)
        with self.assertRaises(ValueError):make_profile(z,z['wall'],.8,1e6)
    def test_missing_power_or_nonfinite_cell_rejected(self):
        z=self.arrays();z['cell_power'][0]=.3
        with self.assertRaises(ValueError):make_profile(z,z['wall'],.7,1e6)
        z=self.arrays();z['cell_power'][0]=np.nan
        with self.assertRaises(ValueError):make_profile(z,z['wall'],.7,1e6)

if __name__=='__main__':unittest.main()
