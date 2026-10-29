import unittest
import numpy as np
from shapely.geometry import LineString
from history_data import normalize_psi,truncate_path,select_history_records
class HistoryDataTests(unittest.TestCase):
 def test_infeasible_reference_keeps_improving_feasible_frontier(self):
  records=[{'case':1,'stage':'baseline','feasible':False,'worst_all_wall_peak_ratio':1.},
           {'case':2,'stage':'search','feasible':True,'worst_all_wall_peak_ratio':1.2},
           {'case':3,'stage':'search','feasible':True,'worst_all_wall_peak_ratio':1.15}]
  self.assertEqual([r['case'] for r in select_history_records(records,{1:True,2:True,3:True})],[1,2,3])
 def test_baseline_remains_reference_when_it_fails_angle_constraint(self):
  records=[{'case':1,'stage':'baseline','feasible':False,'worst_all_wall_peak_ratio':1.},
           {'case':2,'stage':'search','feasible':True,'worst_all_wall_peak_ratio':.95},
           {'case':3,'stage':'search','feasible':True,'worst_all_wall_peak_ratio':.97},
           {'case':4,'stage':'search','feasible':True,'worst_all_wall_peak_ratio':.9},
           {'case':5,'stage':'search','feasible':True,'worst_all_wall_peak_ratio':.94}]
  result=select_history_records(records,{1:True,2:True,3:True,4:False,5:True})
  self.assertEqual([r['case'] for r in result],[1,2,5]);self.assertFalse(result[0]['feasible'])
 def test_axis_separatrix_values_are_zero_and_one_for_either_sign(self):
  np.testing.assert_allclose(normalize_psi(np.array([5.,3.,1.]),5.,1.),[0.,.5,1.])
  with self.assertRaises(ValueError):normalize_psi([1.],1.,1.)
 def test_sol_curve_stops_at_first_hit_not_the_trace_extension(self):
  path=LineString([[2.,0.],[2.,-1.],[2.,-2.]])
  curve=truncate_path(path,[2.,-.8],stride=1)
  np.testing.assert_allclose(curve[-1],[2.,-.8]);self.assertLessEqual(len(curve),3)
if __name__=='__main__':unittest.main()
