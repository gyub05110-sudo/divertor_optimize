"""Independent analytic checks for frozen-equilibrium direct heat mapping."""
import unittest

import numpy as np

from d3d_optimization.frozen_field import FrozenField, exponential_source
from d3d_optimization.frozen_mapping import first_hits, map_heat, deform_target, trace_paths


class FrozenHeatTests(unittest.TestCase):
    def setUp(self):
        self.r = np.linspace(1.0, 3.0, 33)
        self.z = np.linspace(-2.0, 2.0, 33)
        rr, zz = np.meshgrid(self.r, self.z)
        self.field = FrozenField(self.r, self.z, rr ** 2 + zz ** 2, 100.0)

    def test_analytic_field_signs_and_tf_units(self):
        b = self.field.magnetic(np.array([[2.0, 0.5]]))[0]
        np.testing.assert_allclose(b, [-0.5, 2.0, 1.44], atol=1e-10)

    def test_outside_grid_is_rejected(self):
        with self.assertRaises(ValueError):
            self.field.magnetic(np.array([[0.9, 0.0]]))

    def test_full_field_angle_uses_toroidal_component(self):
        angle = self.field.incidence(np.array([[2.0, 0.0]]), np.array([[0.0, 1.0]]))[0]
        self.assertAlmostEqual(angle, np.degrees(np.arcsin(2.0 / np.hypot(2.0, 1.44))))

    def test_every_width_has_same_input_power(self):
        edges = np.linspace(0, 0.1, 1001)
        for width in [0.001, 0.002, 0.003, 0.005, 0.01]:
            source = exponential_source(self.field, 1.5, 0.0, edges, width)
            self.assertAlmostEqual(np.sum(source.weights), 1.0, places=12)
            self.assertTrue(np.all(source.weights >= 0))
            self.assertTrue(np.all(source.h > 0))

    def test_first_surface_receives_heat_and_unmapped_is_preserved(self):
        paths = [np.array([[1.5, 0.0], [1.5, -1.8]]),
                 np.array([[2.5, 0.0], [2.5, -1.8]])]
        wall = np.array([[1.0, -1.0], [2.0, -1.0], [2.0, -1.5], [1.0, -1.5]])
        hit = first_hits(paths, wall)
        self.assertEqual(hit.segment[0], 0)
        self.assertEqual(hit.segment[1], -1)
        np.testing.assert_allclose(hit.points[0], [1.5, -1.0])
        mapped = map_heat(self.field, hit, np.array([0.3, 0.7]), np.ones(2), wall)
        self.assertAlmostEqual(mapped.mapped_power, 0.3)
        self.assertAlmostEqual(mapped.unmapped_power, 0.7)
        self.assertAlmostEqual(mapped.balance_error, 0.0)
        self.assertEqual(mapped.segment_power[2], 0.0)

    def test_surface_area_and_angle_floor_are_independent(self):
        paths = [np.array([[1.5, 0.0], [1.5, -1.8]])]
        wall = np.array([[1.0, -1.0], [2.0, -1.0]])
        hit = first_hits(paths, wall)
        mapped = map_heat(self.field, hit, np.array([1.0]), np.ones(1), wall)
        self.assertAlmostEqual(mapped.areas[0], 3.0 * np.pi)
        self.assertEqual(mapped.feasible_angle(1.5, [0]), bool(mapped.angles[0] >= 1.5))

    def test_spline_endpoints_are_fixed_and_zero_is_identity(self):
        target = np.column_stack([np.linspace(1.0, 2.0, 51), np.zeros(51)])
        np.testing.assert_allclose(deform_target(target, np.zeros(4)), target)
        changed = deform_target(target, np.full(4, 0.02))
        np.testing.assert_allclose(changed[[0, -1]], target[[0, -1]])
        self.assertLessEqual(np.max(np.linalg.norm(changed-target, axis=1)), 0.02 + 1e-12)

    def test_trace_preserves_analytic_flux_and_moves_to_lower_branch(self):
        path = trace_paths(self.field, [[2.0, 0.0]], step=0.01, max_length=0.5)[0]
        self.assertLess(path[-1, 1], 0)
        self.assertLess(np.max(np.abs(self.field.values(path)-4.0)), 1e-8)

    def test_case_contract_has_one_angle_and_five_widths(self):
        from d3d_optimization.case import ANGLE_FLOOR, WIDTHS_MM
        self.assertEqual(ANGLE_FLOOR, 1.5)
        np.testing.assert_array_equal(WIDTHS_MM, [1, 2, 3, 5, 10])

    def test_independent_surface_quadrature_matches_source_power(self):
        # psi=R^2: straight vertical tubes. Exact dmu=4*pi*R*dR;
        # q=h*2 and dA=2*pi*R*dR give the same power independently.
        rr, zz = np.meshgrid(self.r, self.z)
        field = FrozenField(self.r, self.z, rr**2, 100.0)
        edges = np.linspace(0, 0.1, 4001)
        source = exponential_source(field, 1.5, 0, edges, 0.003)
        radii = 1.5+source.x
        surface_power = np.sum(source.h * 2 * 2*np.pi*radii * np.diff(edges))
        self.assertAlmostEqual(surface_power, 1.0, places=5)

    def test_constraint_summary_rejects_subfloor_candidate(self):
        from d3d_optimization.case import candidate_feasible
        self.assertFalse(candidate_feasible(1.499, 0, 0))
        self.assertTrue(candidate_feasible(1.5, 0, 0))
        self.assertFalse(candidate_feasible(2.0, 0.002, 0))
        self.assertFalse(candidate_feasible(2.0, 0, 0.001))

    def test_linear_gradient_is_derivative_of_same_interpolated_psi(self):
        linear = FrozenField(self.r, self.z, self.field.psi, 100.0, "linear")
        point = np.array([[1.531, 0.281]])
        gradient = linear.gradient(point)[0]
        eps = 1e-6
        finite_difference = [(linear.values(point + np.eye(2)[i]*eps)[0] -
                              linear.values(point - np.eye(2)[i]*eps)[0])/(2*eps) for i in range(2)]
        np.testing.assert_allclose(gradient, finite_difference, atol=1e-8)

    def test_quintic_field_reproduces_analytic_field(self):
        quintic = FrozenField(self.r, self.z, self.field.psi, 100.0, "quintic")
        np.testing.assert_allclose(quintic.magnetic([[2, 0.5]]), [[-0.5, 2, 1.44]], atol=1e-10)

    def test_refined_infeasibility_invalidates_coarse_feasible_claim(self):
        from d3d_optimization.case import verified_feasibility
        self.assertFalse(verified_feasibility(True, [{"feasible": True}, {"feasible": False}]))
        self.assertFalse(verified_feasibility(True, []))
        self.assertTrue(verified_feasibility(True, [{"feasible": True}, {"feasible": True}]))

    def test_angle_constraint_ignores_non_target_hits(self):
        from d3d_optimization.case import target_angle_status
        status = target_angle_status(np.array([0.4, 2.0, 1.7]), np.array([False, True, True]))
        self.assertEqual(status["minimum_angle_deg"], 1.7)
        self.assertTrue(status["angle_constraint_passed"])
        self.assertEqual(status["target_hit_count"], 2)

    def test_unilluminated_target_is_not_zero_degree_violation(self):
        from d3d_optimization.case import target_angle_status
        status = target_angle_status(np.array([0.4, 0.8]), np.array([False, False]))
        self.assertIsNone(status["minimum_angle_deg"])
        self.assertIsNone(status["angle_constraint_passed"])
        self.assertEqual(status["status"], "target_unilluminated")

    def test_explicit_target_matches_geometry_hash_and_coordinates(self):
        from d3d_optimization.case import select_target
        wall = np.array([[1., 0.], [2., 0.], [2., -1.], [1., 0.]])
        config = {"name": "reference shelf", "geometry_sha256": "abc", "segment_index": 1,
                  "endpoints_RZ_m": [[2., 0.], [2., -1.]],
                  "physical_identity_validated": False}
        self.assertEqual(select_target(wall, "abc", config), 1)
        with self.assertRaises(ValueError):
            select_target(wall, "changed", config)
        with self.assertRaises(ValueError):
            select_target(wall, "abc", dict(config, endpoints_RZ_m=[[1., 0.], [2., 0.]]))
        with self.assertRaises(ValueError):
            select_target(wall, "abc", None)

    def test_evaluator_keeps_non_target_heat_without_constraining_its_angle(self):
        from shapely.geometry import LineString
        from d3d_optimization.case import ShapeEvaluator
        wall = np.array([[1.1, 1.], [2.9, 1.], [2.9, -1.], [1.1, -1.], [1.1, 1.]])
        paths = [LineString([[2., 0.], [2., -1.8]]), LineString([[2., 0.], [3., 0.]])]
        lcfs = np.array([[1.2, .3], [1.3, .4], [1.2, .5]])
        e = ShapeEvaluator(self.field, lcfs, wall, paths, np.ones(2, dtype=bool),
                           np.array([0., .001, .002]), 1.6, 0., 2)
        result = e.baseline
        self.assertTrue(result["feasible"])
        self.assertLess(result["maps"][0].angles[1], 1.5)
        self.assertGreater(result["angle"], 1.5)
        self.assertGreater(result["maps"][0].mapped_power-result["target_power"][0], 0)
        np.testing.assert_allclose(result["maps"][0].points, [[2., -1.], [2.9, 0.]])

    def test_report_does_not_infer_component_identity_or_peak_coincidence(self):
        import tempfile
        from pathlib import Path
        from d3d_optimization.case import write_report
        row = {"lambda_q_mm": 1, "minimum_target_angle_deg": 2., "q_over_P_peak_m2inv": 1.,
               "all_wall_peak_ratio": 1., "target_power": .1, "other_wall_power": .9,
               "unmapped_power": 0., "power_balance_error": 0.}
        manifest = {"results": [dict(row, geometry="baseline"), dict(row, geometry="candidate")],
                    "time_ms": 1400., "feasible_candidate": True, "preflight_only": True,
                    "angle_application": {"status": "angle_satisfied", "target_hit_count": 1},
                    "target_configuration": {"name": "arbitrary reference", "segment_index": 5,
                        "endpoints_RZ_m": [[1., 0.], [2., 0.]], "physical_identity_validated": False},
                    "angle_application_diagnostic": {"minimum_angle_hit_RZ_m": [1.5, 0.],
                        "minimum_angle_source_x_mm": 20., "minimum_angle_is_lambda1mm_peak_location": False,
                        "comparison_shelf_angle_deg": 2.}}
        with tempfile.TemporaryDirectory() as folder:
            write_report(Path(folder), manifest, [])
            report = (Path(folder)/"REPORT.md").read_text(encoding="utf-8")
        self.assertIn("세분화 검증 전", report)
        self.assertNotIn("자동 선택했습니다", report)
        self.assertNotIn("선택된 하부 outer target", report)
        self.assertNotIn("λq=1 mm 최대 열유속 위치와 같은 점입니다", report)
        self.assertNotIn("떨어진 첫 tube", report)

    def test_zero_power_target_hit_is_excluded_from_lower_bound(self):
        from d3d_optimization.frozen_mapping import Hits
        wall = np.array([[1., -1.], [3., -1.]])
        hits = Hits(np.array([0, 0]), np.array([[1.5, -1.], [2., -1.]]), np.ones(2))
        heat = map_heat(self.field, hits, [0., 1.], [0., 1.], wall)
        heat.angles[:] = [0.1, 2.]
        self.assertTrue(heat.feasible_angle(1.5, [0]))
        status = heat.angle_application(1.5, [0])
        self.assertEqual(status["minimum_angle_deg"], 2.)
        self.assertEqual(status["excluded_zero_power_target_hit_count"], 1)

    def test_no_target_hits_have_no_angle_penalty(self):
        from d3d_optimization.case import constraint_penalties
        status = {"status": "target_unilluminated", "minimum_angle_deg": None}
        penalties = constraint_penalties(status, 0.)
        self.assertEqual(penalties["angle_lower_bound"], 0.)
        self.assertGreater(penalties["target_unilluminated"], 0.)

    def test_low_angle_on_positive_power_target_still_fails(self):
        from d3d_optimization.frozen_mapping import Hits
        wall = np.array([[1., -1.], [3., -1.]])
        hits = Hits(np.array([0, 0]), np.array([[1.5, -1.], [2., -1.]]), np.ones(2))
        heat = map_heat(self.field, hits, [1e-40, 1.], [1e-40, 1.], wall)
        heat.angles[:] = [0.1, 2.]
        self.assertFalse(heat.feasible_angle(1.5, [0]))
        self.assertEqual(heat.angle_application(1.5, [0])["minimum_angle_deg"], .1)

    def test_scope_audit_excludes_other_surface_and_records_first_hit(self):
        from shapely.geometry import LineString
        from d3d_optimization.case import ShapeEvaluator, angle_scope_rows
        wall = np.array([[1.1, 1.], [2.9, 1.], [2.9, -1.], [1.1, -1.], [1.1, 1.]])
        paths = [LineString([[2., 0.], [2., -1.8]]), LineString([[2., 0.], [3., 0.]])]
        lcfs = np.array([[1.2, .3], [1.3, .4], [1.2, .5]])
        e = ShapeEvaluator(self.field, lcfs, wall, paths, np.ones(2, dtype=bool),
                           np.array([0., .001, .002]), 1.6, 0., 2)
        rows = angle_scope_rows(e, e.baseline, "baseline")
        self.assertTrue(rows[0]["angle_floor_applied"])
        self.assertFalse(rows[1]["angle_floor_applied"])
        self.assertEqual(rows[1]["excluded_reason"], "non_target_surface")
        self.assertEqual(rows[1]["original_wall_segment"], 1)
        self.assertEqual(rows[1]["first_hit_R_m"], 2.9)

    def test_shadowed_target_crossing_has_no_angle_lower_bound_penalty(self):
        from shapely.geometry import LineString
        from d3d_optimization.case import ShapeEvaluator, angle_scope_rows
        wall = np.array([[1.1, 1.], [2.9, 1.], [2.9, -1.], [1.1, -1.], [1.1, 1.]])
        # Both paths hit the right wall first, then cross the target later.
        paths = [LineString([[2., 0.], [3., 0.], [3., -1.8], [2., -1.8], [2., -.8]]),
                 LineString([[2.1, 0.], [3., 0.], [3., -1.8], [2.1, -1.8], [2.1, -.8]])]
        lcfs = np.array([[1.2, .3], [1.3, .4], [1.2, .5]])
        e = ShapeEvaluator(self.field, lcfs, wall, paths, np.ones(2, dtype=bool),
                           np.array([0., .001, .002]), 1.6, 0., 2)
        result = e.baseline
        self.assertIsNone(result["angle"])
        self.assertEqual(result["angle_application"]["status"], "target_unilluminated")
        self.assertEqual(result["constraint_penalties"]["angle_lower_bound"], 0.)
        self.assertFalse(result["feasible"])
        self.assertEqual(result["target_power"], [0.]*5)
        self.assertAlmostEqual(result["maps"][0].mapped_power, 1.)
        self.assertTrue(all(not r["angle_floor_applied"] for r in angle_scope_rows(e, result, "baseline")))

    def test_nonfinite_angle_on_excluded_surface_does_not_pollute_minimum(self):
        from d3d_optimization.case import target_angle_status
        status = target_angle_status([np.nan, 2.], [False, True], [1., 1.])
        self.assertTrue(status["angle_constraint_passed"])
        invalid = target_angle_status([np.nan, 2.], [True, True], [1., 1.])
        self.assertEqual(invalid["status"], "invalid_target_angle")
        self.assertFalse(invalid["angle_constraint_passed"])

    def test_invalid_target_angle_is_not_reported_as_unilluminated(self):
        from d3d_optimization.case import angle_text, case_status
        self.assertIn("invalid", angle_text(None, "invalid_target_angle"))
        result = {"feasible": False, "angle_application": {"status": "invalid_target_angle"}}
        self.assertEqual(case_status(result), "INVALID_TARGET_ANGLE")


if __name__ == '__main__':
    unittest.main()
