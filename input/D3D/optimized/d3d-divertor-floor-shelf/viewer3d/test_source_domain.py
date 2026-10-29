import unittest
import numpy as np
from source_domain import physical_source_limit,require_physical_extension,scope_reason

class SourceDomainTests(unittest.TestCase):
    def test_limit_stops_at_first_outboard_wall(self):
        wall=[[1,-1],[2,-1],[2,1],[1,1],[1,-1]]
        result=physical_source_limit(wall,1.4,0.,2.5)
        self.assertAlmostEqual(result['physical_max_offset_m'],.6)
    def test_grid_limit_is_respected(self):
        wall=[[1,-1],[2,-1],[2,1],[1,1],[1,-1]]
        self.assertAlmostEqual(physical_source_limit(wall,1.4,0.,1.8)['physical_max_offset_m'],.4)
    def test_cannot_extend_source_into_material(self):
        with self.assertRaisesRegex(ValueError,'wall'):
            require_physical_extension(.15,.0988)
    def test_scope_reasons_do_not_assume_field_sign(self):
        self.assertEqual(scope_reason(-.01,0.,.02),'separatrix_side')
        self.assertEqual(scope_reason(.03,0.,.02),'beyond_source_end')
        self.assertEqual(scope_reason(.03,.02,0.),'separatrix_side')

if __name__=='__main__':unittest.main()
