"""Actual v3 primary-policy adversaries, separate from frozen run driver tests."""
import unittest,copy
import numpy as np
from distribution_policy_v3_final import profile,distribution_compare,qualify
from qualify_omp_stages import checked_read
from pathlib import Path

def smooth(n=600,swap=0,shift=0,spike=False,scale=1):
    def q(x):
        base=.1+(1+swap)*np.exp(-((x-.3-shift)/.04)**2)+(1-swap)*np.exp(-((x-.7-shift)/.04)**2)
        if spike:base+=.6*np.exp(-((x-.5)/.0001)**2)
        return scale*base
    edges=np.linspace(0,1,n+1);g,w=np.polynomial.legendre.leggauss(12)
    mid=(edges[:-1]+edges[1:])/2;h=np.diff(edges)/2
    means=(q(mid[:,None]+h[:,None]*g)*w).sum(axis=1)/2
    p=profile(np.column_stack((np.ones(len(edges)),edges)),means)
    u=np.unique(np.r_[np.linspace(0,1,20001),edges]);p['continuous']=(u,q(u));return p

class AdversarialPrimaryTests(unittest.TestCase):
    def test_rank_exchange_allowed(self):
        self.assertTrue(distribution_compare(smooth(swap=.0005),smooth(swap=-.0005))['passed'])
    def test_smooth_function_different_averaging_windows_allowed(self):
        self.assertTrue(distribution_compare(smooth(600),smooth(1300))['passed'])
    def test_material_translation_fails(self):
        self.assertFalse(distribution_compare(smooth(),smooth(shift=.01))['passed'])
    def test_secondary_narrow_peak_below_global_max_fails(self):
        self.assertFalse(distribution_compare(smooth(),smooth(spike=True))['passed'])
    def test_lost_power_fails(self):
        self.assertFalse(distribution_compare(smooth(),smooth(scale=.98))['passed'])
    def test_reversed_wall_fails(self):
        a=smooth();b=smooth();b['wall']=b['wall'][::-1].copy()
        self.assertFalse(distribution_compare(a,b)['passed'])
    def test_no_independent_audit_bypass(self):
        path=Path('D:/DivertorOptimization/studies/omp-stage-qualification-20261104/R020/validation/0eb363d2b21c1d391948f380244813840964a06edab792fffca167be74755cc1/validation.json')
        v=checked_read(path);d={'passed':True}
        self.assertTrue(qualify(v,d)['verdict']['improvement_confirmed'])
        for changed in ('geometry','thermal','wall'):
            bad=copy.deepcopy(v)
            if changed=='geometry':bad['levels'][0]['geometry_audit']['passed_independent']=False
            elif changed=='thermal':bad['levels'][1]['rows']['selected']['audit_passed']=False
            else:bad['wall_only_refinement']['rows']['baseline']['feasible']=False
            self.assertFalse(qualify(bad,d)['verdict']['improvement_confirmed'])
    def test_bad_existing_numeric_gate_remains_failure(self):
        path=Path('D:/DivertorOptimization/studies/omp-stage-qualification-20261104/R020/validation/0eb363d2b21c1d391948f380244813840964a06edab792fffca167be74755cc1/validation.json')
        v=checked_read(path);v['verdict']['max_wall_cdf_change']=.1
        self.assertFalse(qualify(v,{'passed':True})['verdict']['improvement_confirmed'])

if __name__=='__main__':unittest.main()
