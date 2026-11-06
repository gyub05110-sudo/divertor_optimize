"""Selection and finite-setting decision tests; never execute thermal calculations."""
import copy
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


SETTINGS = ("cubic640", "cubic1280", "quintic640", "quintic1280")


def design(value):
    candidate = dict(source="OMP", radius_mm=20, vector=[value] + [0.] * 15)
    identity = hashlib.sha256(json.dumps(candidate, sort_keys=True).encode()).hexdigest()
    return identity, candidate


def state(peak=10., A2=2., eligible=True):
    return dict(eligible=eligible, peak=peak, A2=A2, errors=[] if eligible else ["failed"], ledger={})


def event(value, peak, *, end=100., stage="fine", status="completed", eligible=True):
    identity, candidate = design(value)
    return dict(id=identity, candidate=candidate, stage=stage, start_s=end-1., end_s=end,
                status=status, measurement=state(peak, eligible=eligible))


def checkpoint(e=None):
    if e is None:
        return dict(available=True, kind="original_wall", id=None, candidate=None, peak=10., A2=2.)
    return dict(available=True, kind="new_design", id=e["id"], candidate=e["candidate"],
                peak=e["measurement"]["peak"], A2=e["measurement"]["A2"])


def arm(events, chosen=None, policy="fine_adaptive", seed=20611):
    return dict(case_id="example", seed=seed, policy=policy, budget_s=600,
                baseline=state(), events=events,
                checkpoints={str(t): checkpoint(chosen) for t in (200, 400, 600)})


class UncertaintyTests(unittest.TestCase):
    def setUp(self):
        name = "verification.nf_anytime_uncertainty"
        self.assertIsNotNone(importlib.util.find_spec(name), "uncertainty implementation is absent")
        self.mod = importlib.import_module(name)

    def test_selection_includes_chosen_last_coarse_overrun_and_closest_distinct_native(self):
        chosen, closest, far = event(.001, 8.), event(.002, 8.01), event(.003, 6.)
        late = event(.004, 7., end=601., stage="coarse", status="late")
        invalid = event(.005, 8.001, eligible=False)
        locked = self.mod.lock_candidate_set([arm([chosen, far, invalid, late], chosen),
                                             arm([closest], policy="random_fine", seed=20627)])
        rows = {r["id"]: r for r in locked["candidates"]}
        self.assertEqual(set(rows), {self.mod.BASELINE_ID, chosen["id"], closest["id"], late["id"]})
        self.assertEqual(locked["nearest_observed"][0]["neighbor_id"], closest["id"])
        self.assertTrue(any(r["kind"] == "last_overrun" for r in rows[late["id"]]["reasons"]))
        self.assertEqual(rows[late["id"]]["native_observations"], [])

    def test_nearest_peak_tie_uses_candidate_id_and_excludes_late_fine(self):
        chosen = event(.001, 8.)
        a, b = event(.002, 8.5), event(.003, 7.5)
        late = event(.004, 8., end=601., status="late")
        locked = self.mod.lock_candidate_set([arm([chosen, a, b, late], chosen)])
        self.assertEqual(locked["nearest_observed"][0]["neighbor_id"], min(a["id"], b["id"]))

    def test_late_fine_native_measurement_is_replayed_but_never_enters_nearest_pool(self):
        chosen, neighbor = event(.001, 8.), event(.002, 8.1)
        late = event(.003, 8.001, end=601., status="late")
        # A GPU measurement can exist even when the subsequent CPU guard never
        # finishes, leaving no objective in the terminal event.
        late["error"] = "Deadline: CPU parity guard"
        locked = self.mod.lock_candidate_set([arm([chosen, neighbor, late], chosen)])
        row = next(row for row in locked["candidates"] if row["id"] == late["id"])
        self.assertEqual(len(row["native_observations"]), 1)
        self.assertTrue(row["native_observations"][0]["late"])
        self.assertEqual(row["native_observations"][0]["status"], "late")
        self.assertNotIn(late["id"], locked["observed_native_pool_ids"])
        self.assertEqual(locked["nearest_observed"][0]["neighbor_id"], neighbor["id"])
        self.assertFalse(self.mod.compare_native(state(8.001 * (1 + 2e-10)),
                                                row["native_observations"])["passed"])

    def test_late_fine_geometry_ineligibility_is_retained_for_eligibility_replay(self):
        chosen = event(.001, 8.)
        geometry = event(.002, None, eligible=False)
        geometry["measurement"] = dict(eligible=False, peak=None, A2=None, ledger={},
                                       errors=["continuous_shape_constraint_failed"])
        geometry.update(objective=11., physics_evaluated=False)
        late = copy.deepcopy(geometry)
        late.update(start_s=599., end_s=601., status="late")
        locked = self.mod.lock_candidate_set([arm([chosen, geometry, late], chosen)])
        row = next(row for row in locked["candidates"] if row["id"] == late["id"])
        self.assertEqual(len(row["native_observations"]), 2)
        self.assertFalse(row["native_observations"][0]["late"])
        self.assertTrue(row["native_observations"][1]["late"])
        self.assertNotIn(late["id"], locked["observed_native_pool_ids"])
        self.assertFalse(self.mod.compare_native(state(8.), row["native_observations"])["passed"])
        self.assertTrue(self.mod.compare_native(geometry["measurement"], row["native_observations"])["passed"])

    def test_gpu_measurement_before_parity_deadline_does_not_need_physics_flag(self):
        chosen = event(.001, 8.)
        late = event(.002, 4.582103994277724, end=613.765876, status="late")
        late["start_s"] = 599.287758
        late["measurement"]["A2"] = .39837692279420306
        late["measurement"]["ledger"] = {"balance_error": 0.}
        late["error"] = "Deadline: CPU parity guard"
        self.assertNotIn("physics_evaluated", late)
        self.assertNotIn("objective", late)
        locked = self.mod.lock_candidate_set([arm([chosen, late], chosen)])
        row = next(row for row in locked["candidates"] if row["id"] == late["id"])
        self.assertEqual(len(row["native_observations"]), 1)
        self.assertEqual(row["native_observations"][0]["measurement"], late["measurement"])
        self.assertTrue(row["native_observations"][0]["late"])
        self.assertFalse(self.mod.compare_native(state(4.582103994277724, A2=.39837692279420306, eligible=False),
                                                row["native_observations"])["passed"])

    def test_deadline_placeholder_locks_design_without_inventing_native_measurement(self):
        chosen = event(.001, 8.)
        interrupted = event(.002, None, end=601., status="late", eligible=False)
        interrupted["measurement"] = dict(eligible=False, peak=None, A2=None, ledger={}, errors=["Deadline: "])
        locked = self.mod.lock_candidate_set([arm([chosen, interrupted], chosen)])
        row = next(row for row in locked["candidates"] if row["id"] == interrupted["id"])
        self.assertTrue(any(reason["kind"] == "last_overrun" for reason in row["reasons"]))
        self.assertEqual(row["native_observations"], [])
        self.assertNotIn(interrupted["id"], locked["observed_native_pool_ids"])
        replay = self.mod.compare_native(state(8.), row["native_observations"])
        self.assertTrue(replay["passed"])
        self.assertFalse(replay["previously_scored"])
        self.assertEqual(replay["compared_count"], 0)

    def test_duplicate_candidate_is_locked_once_but_all_reasons_are_preserved(self):
        chosen = event(.001, 8.)
        locked = self.mod.lock_candidate_set([arm([chosen], chosen)])
        rows = [r for r in locked["candidates"] if r["id"] == chosen["id"]]
        self.assertEqual(len(rows), 1)
        self.assertEqual(sum(r["kind"] == "checkpoint" for r in rows[0]["reasons"]), 3)
        self.assertIsNone(locked["nearest_observed"][0]["neighbor_id"])

    def test_malformed_identity_or_ineligible_selected_candidate_is_rejected(self):
        chosen = event(.001, 8.)
        bad = arm([chosen], chosen)
        bad["checkpoints"]["200"]["id"] = "changed"
        with self.assertRaises(ValueError):
            self.mod.lock_candidate_set([bad])
        bad = arm([event(.001, 8., eligible=False)], chosen)
        with self.assertRaises(ValueError):
            self.mod.lock_candidate_set([bad])

    def test_chosen_native_score_requires_on_time_evidence_before_its_checkpoint(self):
        chosen = event(.001, 8., end=201.)
        with self.assertRaises(ValueError):
            self.mod.lock_candidate_set([arm([chosen], chosen)])

    def test_all_paired_improvements_can_coexist_with_overlapping_envelopes(self):
        baseline = {key: state(p) for key, p in zip(SETTINGS, [10., 8., 10., 8.])}
        candidate = {key: state(p, A2=1.9) for key, p in zip(SETTINGS, [9., 7., 9., 7.])}
        result = self.mod.assess_candidate(baseline, candidate)
        self.assertTrue(result["all_settings_paired_improvement"])
        self.assertFalse(result["candidate_max_below_baseline_min"])
        self.assertEqual(result["candidate_peak_range"], [7., 9.])
        self.assertTrue(all(r["A2_decreased"] for r in result["settings"].values()))

    def test_separated_envelopes_and_A2_changes_are_reported(self):
        base = {key: state(10., 2.) for key in SETTINGS}
        cand = {key: state(8., 3.) for key in SETTINGS}
        result = self.mod.assess_candidate(base, cand)
        self.assertTrue(result["candidate_max_below_baseline_min"])
        self.assertAlmostEqual(result["settings"]["cubic1280"]["A2_change_pct"], 50.)
        self.assertFalse(result["local_q_1pct_certified"])

    def test_ineligible_or_missing_setting_cannot_certify_any_all_setting_claim(self):
        base = {key: state() for key in SETTINGS}
        for cand in [{key: state(8.) for key in SETTINGS if key != "quintic1280"},
                     {key: state(8., eligible=key != "quintic1280") for key in SETTINGS}]:
            result = self.mod.assess_candidate(base, cand)
            self.assertFalse(result["all_settings_eligible"])
            self.assertFalse(result["all_settings_paired_improvement"])
            self.assertFalse(result["candidate_max_below_baseline_min"])

    def test_malformed_eligible_peak_fails_closed(self):
        base = {key: state() for key in SETTINGS}
        candidate = {key: state(float("nan")) for key in SETTINGS}
        with self.assertRaises(ValueError):
            self.mod.assess_candidate(base, candidate)

    def test_candidate_range_overlap_is_inclusive_and_unknown_on_failure(self):
        a = {key: state(p) for key, p in zip(SETTINGS, [7., 8., 7., 8.])}
        b = {key: state(p) for key, p in zip(SETTINGS, [8., 9., 8., 9.])}
        self.assertTrue(self.mod.compare_candidate_ranges(a, b)["possible_overlap"])
        b["quintic1280"] = state(None, A2=None, eligible=False)
        self.assertIsNone(self.mod.compare_candidate_ranges(a, b)["possible_overlap"])

    def test_native_replay_detects_lost_eligibility_and_numerical_drift(self):
        observed = [dict(measurement=state(8., 2.))]
        self.assertTrue(self.mod.compare_native(state(8. * (1 + 5e-11), 2.), observed)["passed"])
        self.assertFalse(self.mod.compare_native(state(8. * (1 + 2e-10), 2.), observed)["passed"])
        self.assertFalse(self.mod.compare_native(state(8., eligible=False), observed)["passed"])

    def test_late_ineligible_thermal_replay_detects_peak_and_A2_drift(self):
        chosen = event(.001, 7.)
        late = event(.002, 8., end=601., status="late", eligible=False)
        late["measurement"]["ledger"] = {"balance_error": .01}
        late["measurement"]["errors"] = ["synthetic_physical_constraint_failure"]
        selected = self.mod.lock_candidate_set([arm([chosen, late], chosen)])
        row = next(row for row in selected["candidates"] if row["id"] == late["id"])
        self.assertEqual(len(row["native_observations"]), 1)
        for changed in ("peak", "A2"):
            with self.subTest(observable=changed):
                measured = copy.deepcopy(late["measurement"])
                measured[changed] *= 1 + 2e-10
                replay = self.mod.compare_native(measured, row["native_observations"])
                self.assertFalse(replay["passed"])
                self.assertFalse(replay["all_available_observations_matched"])
                comparison = replay["comparisons"][0]
                self.assertTrue(comparison["eligibility_matches"])
                self.assertTrue(comparison["observable_presence_matches"])
                self.assertGreater(comparison["relative_differences"][changed], self.mod.NATIVE_RTOL)

    def test_ineligible_native_observable_presence_change_fails_in_both_directions(self):
        full = state(8., A2=2., eligible=False)
        for missing_key in ("peak", "A2"):
            missing = dict(full, **{missing_key: None})
            for measured, previous in ((missing, full), (full, missing)):
                with self.subTest(missing_key=missing_key, replayed=measured):
                    replay = self.mod.compare_native(measured, [dict(measurement=previous)])
                    self.assertFalse(replay["passed"])
                    comparison = replay["comparisons"][0]
                    self.assertTrue(comparison["eligibility_matches"])
                    self.assertFalse(comparison["observable_presence_matches"])
                    self.assertEqual(comparison["observable_presence_mismatches"], [missing_key])
                    self.assertEqual(comparison["observable_presence"][missing_key],
                                     dict(stored=previous[missing_key] is not None,
                                          replayed=measured[missing_key] is not None))
                    self.assertEqual(comparison["observable_compared_count"], 1)

    def test_geometry_only_replay_records_eligibility_evidence_without_heat_agreement(self):
        geometry = dict(eligible=False, peak=None, A2=None, ledger={},
                        errors=["continuous_shape_constraint_failed"])
        replay = self.mod.compare_native(geometry, [dict(measurement=geometry)])
        self.assertTrue(replay["passed"])
        self.assertTrue(replay["previously_scored"])
        self.assertEqual(replay["compared_count"], 1)
        self.assertEqual(replay["observable_compared_count"], 0)
        self.assertEqual(replay["eligibility_only_count"], 1)
        self.assertEqual(replay["comparisons"][0]["relative_differences"], {})
        self.assertTrue(replay["comparisons"][0]["eligibility_only"])

    def test_case_replay_coverage_distinguishes_coarse_only_overrun_without_native_history(self):
        chosen = event(.001, 8.)
        late = event(.002, 9., stage="coarse", end=601., status="late")
        selected = self.mod.lock_candidate_set([arm([chosen, late], chosen)])
        values = {row["id"]: {name: state(10. if row["id"] == self.mod.BASELINE_ID else
                                                   8. if row["id"] == chosen["id"] else 9.)
                              for name in SETTINGS} for row in selected["candidates"]}
        quality = self.mod.assess_case(selected, values)
        self.assertTrue(quality["all_available_native_observations_matched"])
        self.assertEqual(quality["native_replay_coverage"],
                         dict(locked_candidate_count=3, previously_scored_candidate_count=2,
                              unobserved_candidate_count=1, compared_observation_count=2,
                              compared_observable_count=4, eligibility_only_observation_count=0))
        no_history = next(row for row in quality["candidates"] if row["id"] == late["id"])
        self.assertFalse(no_history["native_reproducibility"]["previously_scored"])
        self.assertEqual(no_history["native_reproducibility"]["compared_count"], 0)

    def test_quality_compares_chosen_candidate_with_its_locked_nearest_observed_peer(self):
        chosen, neighbor = event(.001, 8.), event(.002, 8.01)
        selection = self.mod.lock_candidate_set([arm([chosen, neighbor], chosen)])
        values = {row["id"]: {name: state(10. if row["id"] == self.mod.BASELINE_ID else
                                                    8. if row["id"] == chosen["id"] else 8.01)
                              for name in SETTINGS} for row in selection["candidates"]}
        quality = self.mod.assess_case(selection, values)
        self.assertIn("nearest_observed_pairs", quality)
        self.assertEqual(len(quality["nearest_observed_pairs"]), 1)
        pair = quality["nearest_observed_pairs"][0]
        self.assertEqual(pair["selected_id"], chosen["id"])
        self.assertEqual(pair["neighbor_id"], neighbor["id"])
        self.assertFalse(pair["possible_overlap"])
        self.assertTrue(pair["first_range_strictly_below_second"])

    def test_all_four_setting_requests_use_locked_ids_and_resume_without_physics(self):
        from verification.nf_anytime_reference import EvidenceStore
        chosen = event(.001, 8.)
        selection = self.mod.lock_candidate_set([arm([chosen], chosen)])
        with tempfile.TemporaryDirectory() as folder:
            store = EvidenceStore(folder)
            store.write("protocol.json", {"test": True})
            store.write("selection.json", {"cases": [selection]})
            prepared = []

            def prepare(*args, **kwargs):
                self.assertEqual(store.read("selection.json"), {"cases": [selection]})
                prepared.append(kwargs)
                return SimpleNamespace(baseline={})

            with patch("d3d_optimization.boundary_evaluator.prepare_boundary", side_effect=prepare), \
                 patch("d3d_optimization.campaign_geometry.evaluator", return_value=object()), \
                 patch("verification.nf_base_budget.measurement", return_value=state()), \
                 patch.object(self.mod, "evaluate_score", return_value=dict(measurement=state(8.), status="completed", error=None, seconds=.01)):
                for setting in SETTINGS:
                    measured, hashes = self.mod.run_setting(store, dict(case_id="example", root=folder), selection,
                                                           setting, object(), "runs/0000")
                    self.assertEqual(set(measured), {self.mod.BASELINE_ID, chosen["id"]})
                    self.assertEqual(len(hashes), 2)
            self.assertEqual([(p["method"], p["count"], p["step"], p["nodes"], p["wall_step"])
                              for p in prepared], [self.mod.SETTINGS[key] for key in SETTINGS])
            with patch("d3d_optimization.boundary_evaluator.prepare_boundary", side_effect=AssertionError("resume prepared")), \
                 patch.object(self.mod, "evaluate_score", side_effect=AssertionError("resume scored")):
                for setting in SETTINGS:
                    self.mod.run_setting(store, dict(case_id="example", root=folder), selection,
                                         setting, object(), "runs/0001")

    def test_failed_preparation_is_retained_and_not_replaced_on_resume(self):
        from verification.nf_anytime_reference import EvidenceStore
        chosen = event(.001, 8.)
        selection = self.mod.lock_candidate_set([arm([chosen], chosen)])
        with tempfile.TemporaryDirectory() as folder:
            store = EvidenceStore(folder)
            store.write("protocol.json", {"test": True})
            store.write("selection.json", {"cases": [selection]})
            with patch("d3d_optimization.boundary_evaluator.prepare_boundary", side_effect=RuntimeError("synthetic preparation failure")), \
                 patch.object(self.mod, "evaluate_score", side_effect=AssertionError("failed preparation scored")):
                measured, _ = self.mod.run_setting(store, dict(case_id="example", root=folder), selection,
                                                  "quintic1280", object(), "runs/0000")
            self.assertTrue(all(not value["eligible"] for value in measured.values()))
            with patch("d3d_optimization.boundary_evaluator.prepare_boundary", side_effect=AssertionError("replaced failure")):
                repeated, _ = self.mod.run_setting(store, dict(case_id="example", root=folder), selection,
                                                   "quintic1280", object(), "runs/0001")
            self.assertEqual(measured, repeated)


if __name__ == "__main__":
    unittest.main()
