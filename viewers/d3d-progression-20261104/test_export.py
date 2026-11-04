import importlib.util
from pathlib import Path
import unittest

class SelectionTests(unittest.TestCase):
    def setUp(self):
        spec=importlib.util.spec_from_file_location('export_data',Path(__file__).with_name('export_data.py'))
        self.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.m)

    def test_original_records_are_not_rewritten_after_geometry_rejection(self):
        rows=[{'id':'a','ratio':.8},{'id':'bad','ratio':.1},{'id':'b','ratio':.7},{'id':'worse','ratio':.75}]
        records=self.m.RunningBest()
        result=[r for r in rows if records.consider(r['ratio']) and r['id']!='bad']
        self.assertEqual([r['id'] for r in result],['a'])

    def test_heat_uses_cell_area_and_power_not_unordered_point_samples(self):
        import numpy as np
        z={'receiving_wall':np.array([[1,-1.1],[1,-1.2],[1,-1.3]]),
           'area_edges':np.array([0,1,3]),'cell_power':np.array([.2,.8]),'cell_parent':np.array([0,0]),
           'points':np.array([[99,99]]),'q_MW_m2':np.array([999])}
        cells=self.m.heat_cells(z,1e6)
        self.assertEqual(len(cells),1)
        self.assertAlmostEqual(cells[0][-1],1/3)
        self.assertEqual(cells[0][:4],[1,-1.1,1,-1.3])

    def test_no_worse_duplicate_or_nonfinite_candidate(self):
        records=self.m.RunningBest()
        ratios=[1.1,.9,.9,float('nan'),0,1,None,float('inf')]
        self.assertEqual([q for q in ratios if records.consider(q)],[.9])

    def test_saved_profile_cannot_be_attached_to_different_shape(self):
        self.assertFalse(self.m.same_shape([[1,0],[2,0]],[[1,0],[2,.1]]))
        self.assertTrue(self.m.same_shape([[1,0],[2,0]],[[1,0],[1.5,0],[2,0]]))

if __name__=='__main__':unittest.main()
