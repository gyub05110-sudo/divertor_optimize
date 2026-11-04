import unittest,copy
import numpy as np
from distribution_policy_v3_final import profile,distribution_compare,qualify
from qualify_omp_stages import checked_read
from pathlib import Path

class ReviewRegression(unittest.TestCase):
    def test_absolute_linear_error_integrated_exactly_across_zero(self):
        wall=np.array([[1.,0.],[1.,1.]])
        a=profile(wall,[2.5]);b=profile(wall,[2.8])
        a['continuous']=(np.array([0.,1.]),np.array([2.,3.]))
        b['continuous']=(np.array([0.,1.]),np.array([1.8,3.8]))
        self.assertAlmostEqual(distribution_compare(a,b)['l1_relative'],.34/2.65,places=13)
    def test_nonimprovement_cannot_hide_independent_audit_failure(self):
        path=Path('D:/DivertorOptimization/studies/omp-stage-qualification-20261104/R020/validation/0eb363d2b21c1d391948f380244813840964a06edab792fffca167be74755cc1/validation.json')
        v=checked_read(path);v['levels'][0]['rows']['selected']['audit_passed']=False
        v['levels'][2]['rows']['selected']['peak_ratio']=1.1
        r=qualify(v,{'passed':True})['verdict'];self.assertFalse(r['resolution_screen_passed']);self.assertFalse(r['improvement_confirmed'])
    def test_reconstruction_must_match_unmodified_cell_power(self):
        a=profile([[1,0],[1,1]],[1]);a['continuous']=(np.array([0,1]),np.array([2,2]))
        r=distribution_compare(a,a);self.assertFalse(r['passed']);self.assertIn('continuous_reconstruction_power',r['failures'])

if __name__=='__main__':unittest.main()
