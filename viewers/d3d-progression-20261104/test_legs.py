import unittest
import numpy as np

class LegTests(unittest.TestCase):
    def test_clips_at_first_wall_crossing(self):
        from export_legs import clip_leg
        # A descending branch crosses a closed rectangle twice.
        leg=[[1,0],[1,-.5],[1,-1],[1,-2],[1,-3]]
        wall=[[.5,-1],[1.5,-1],[1.5,-2],[.5,-2],[.5,-1]]
        r=clip_leg(leg,wall)
        self.assertTrue(r['resolved'])
        np.testing.assert_allclose(r['strike_RZ_m'],[1,-1])
        self.assertAlmostEqual(r['poloidal_length_m'],1)

    def test_lower_target_extends_same_branch(self):
        from export_legs import clip_leg
        leg=[[1,0],[1.1,-1],[1.3,-2]]
        one=clip_leg(leg,[[.5,-.5],[2,-.5]])
        two=clip_leg(leg,[[.5,-1.5],[2,-1.5]])
        self.assertGreater(two['poloidal_length_m'],one['poloidal_length_m'])
        np.testing.assert_allclose(two['strike_RZ_m'],[1.2,-1.5])

    def test_no_intersection_does_not_fabricate_strike(self):
        from export_legs import clip_leg
        r=clip_leg([[1,0],[1,-1]],[[2,0],[2,-1]])
        self.assertFalse(r['resolved'])
        self.assertIsNone(r['strike_RZ_m'])

if __name__=='__main__':unittest.main()
