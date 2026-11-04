import unittest
import numpy as np
from distribution_policy_v3 import compare, profile, continuous_compare

def fixture(n=2000,shift=0,swap=0):
    x=np.linspace(0,1,n+1);s=(x[:-1]+x[1:])/2
    q=.1+(1+swap)*np.exp(-((s-.3-shift)/.04)**2)+(1-swap)*np.exp(-((s-.7-shift)/.04)**2)
    return profile(np.column_stack((np.ones(n+1),x)),q)

class DistributionPolicyTests(unittest.TestCase):
    def test_secondary_subcell_peak_rejected_even_if_averages_match(self):
        a=fixture();b=fixture();u=np.array([0,.2,.2001,.2002,.5,1.])
        a['continuous']=(u,np.array([.1,.2,.9,.2,1.1,.1]))
        b['continuous']=(u,np.array([.1,.2,.2,.2,1.1,.1]))
        self.assertTrue(compare(a,b)['passed'])
        self.assertFalse(continuous_compare(a,b)['passed'])
    def test_near_equal_peaks_exchange_order(self):
        a,b=fixture(swap=.0005),fixture(swap=-.0005)
        self.assertTrue(compare(a,b)['passed'])
    def test_actual_translation_rejected(self):
        self.assertFalse(compare(fixture(),fixture(shift=.01))['passed'])
    def test_disappearing_peak_rejected(self):
        self.assertFalse(compare(fixture(),fixture(swap=1))['passed'])
    def test_one_cell_spike_cannot_hide_in_integral(self):
        a=fixture();q=a['q'].copy();q[1000]+=1
        b=profile(a['wall'],q)
        result=compare(a,b);self.assertFalse(result['passed']);self.assertGreater(result['linf_relative'],.01)
    def test_power_loss_rejected(self):
        a=fixture();self.assertFalse(compare(a,profile(a['wall'],a['q']*.98))['passed'])
    def test_identical_constant_regridding(self):
        a=profile(np.column_stack((np.ones(8),np.linspace(0,1,8))),np.ones(7))
        b=profile(np.column_stack((np.ones(14),np.linspace(0,1,14))),np.ones(13))
        self.assertTrue(compare(a,b)['passed'])
    def test_wall_shift_is_not_registered_away(self):
        a=fixture();w=a['wall'].copy();w[:,1]+=.01
        self.assertFalse(compare(a,profile(w,a['q']))['passed'])
    def test_nonfinite_negative_and_zero_cells_rejected(self):
        a=fixture()
        for bad in (np.nan,-1):
            q=a['q'].copy();q[4]=bad
            with self.assertRaises(ValueError):profile(a['wall'],q)
        w=a['wall'].copy();w[4]=w[3]
        with self.assertRaises(ValueError):profile(w,a['q'])

if __name__=='__main__':unittest.main()
