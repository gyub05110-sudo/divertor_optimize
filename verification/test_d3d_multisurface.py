import unittest
import numpy as np
from shapely.geometry import LineString
from d3d_optimization.frozen_field import FrozenField


class MultiSurfaceTests(unittest.TestCase):
    def fixture(self):
        from d3d_optimization.multisurface import MultiSurfaceEvaluator
        r = np.linspace(1., 3., 33); z = np.linspace(-2., 2., 33)
        rr, zz = np.meshgrid(r, z)
        field = FrozenField(r, z, rr**2+zz**2, 100.)
        wall = np.array([[1.1, 1.], [2.9, 1.], [2.9, -1.], [1.1, -1.], [1.1, 1.]])
        surfaces = [{"name": "floor", "segment_index": 1, "angle_floor_deg": None},
                    {"name": "shelf", "segment_index": 2, "angle_floor_deg": 1.5}]
        paths = [LineString([[2., 0.], [3., 0.]]), LineString([[2., 0.], [2., -1.8]])]
        return MultiSurfaceEvaluator(field, np.array([[1.2,.3],[1.3,.4],[1.2,.5]]),
                                     wall, paths, np.ones(2,dtype=bool), [0.,.001,.002], 1.6,0.,surfaces,41)

    def test_floor_is_design_surface_without_shelf_angle_floor(self):
        e = self.fixture(); result = e.baseline
        self.assertTrue(result["feasible"])
        self.assertLess(result["maps"][0].angles[0], 1.5)
        self.assertGreater(result["angles"]["shelf"]["minimum_angle_deg"], 1.5)
        self.assertIsNone(result["angles"]["floor"]["angle_floor_deg"])
        self.assertTrue(np.isin(result["maps"][0].segment[0], result["surface_segments"]["floor"]))
        self.assertAlmostEqual(sum(result["surface_power"][n][0] for n in ["floor","shelf"]), 1.)

    def test_each_surface_changes_with_its_own_four_coefficients(self):
        e = self.fixture()
        both = e.evaluate(np.full(8,.002))
        self.assertTrue(both["geometry_valid"])
        for name in ["floor","shelf"]:
            self.assertGreater(np.max(np.linalg.norm(both["targets"][name]-e.targets[name],axis=1)),0.)
            np.testing.assert_allclose(both["targets"][name][[0,-1]],e.targets[name][[0,-1]])
        floor_only=e.evaluate([.002]*4+[0.]*4)
        np.testing.assert_allclose(floor_only["targets"]["shelf"],e.targets["shelf"])

    def test_masks_distinguish_design_membership_and_angle_application(self):
        e=self.fixture(); result=e.baseline
        design,angle=e.hit_masks(result)
        np.testing.assert_array_equal(design,[True,True])
        np.testing.assert_array_equal(angle,[False,True])

    def test_coefficients_cannot_exceed_explicit_bound(self):
        e=self.fixture()
        with self.assertRaises(ValueError):
            e.evaluate([.101]*8)

    def test_design_surface_cannot_enter_vessel_solid(self):
        from shapely.geometry import Polygon
        from d3d_optimization.multisurface import MultiSurfaceEvaluator
        old=self.fixture()
        steel=Polygon([[2.92,-.5],[3.2,-.5],[3.2,.5],[2.92,.5]])
        e=MultiSurfaceEvaluator(old.field,old.lcfs,old.wall,old.paths,old.inside,
             [.0,.001,.002],1.6,0.,old.surfaces,41,vessel_solids=steel)
        self.assertTrue(e.baseline['geometry_valid'])
        self.assertFalse(e.evaluate([.1]*4+[0.]*4)['geometry_valid'])

    def test_cached_first_hits_match_full_changed_wall_intersection(self):
        from d3d_optimization.frozen_mapping import first_hits
        e=self.fixture(); r=e.evaluate([.003]*4+[-.003]*4)
        full=first_hits(e.paths,r["wall"])
        np.testing.assert_allclose(full.points,r["maps"][0].points)
        np.testing.assert_array_equal(full.segment,r["maps"][0].segment)

    def test_unheated_shelf_has_no_floor_violation_when_floor_is_heated(self):
        from d3d_optimization.multisurface import MultiSurfaceEvaluator
        old=self.fixture()
        e=MultiSurfaceEvaluator(old.field,old.lcfs,old.wall,
             [LineString([[2.,0.],[3.,0.]]),LineString([[2.,.1],[3.,.1]])],
             np.ones(2,dtype=bool),[0.,.001,.002],1.6,0.,old.surfaces,41)
        self.assertTrue(e.baseline["feasible"])
        self.assertEqual(e.baseline["angle_penalty"],0.)
        self.assertIsNone(e.baseline["angles"]["shelf"]["minimum_angle_deg"])
        np.testing.assert_array_equal(e.hit_masks(e.baseline)[1],[False,False])


if __name__ == "__main__":
    unittest.main()
