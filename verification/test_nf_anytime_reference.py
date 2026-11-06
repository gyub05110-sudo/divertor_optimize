"""Synthetic finite-panel policy and immutable-evidence tests; no transport runs."""
import importlib
import importlib.util
from itertools import combinations
import copy
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


def state(peak=None, A2=2., eligible=True):
    return dict(eligible=eligible, peak=peak, A2=A2 if peak is not None else None,
                errors=[] if eligible else ["not_eligible"], ledger={})


def row(i, coarse, fine, eligible=True):
    return dict(i=i, id=str(i), coarse=state(coarse, eligible=eligible), fine=fine)


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        name = "verification.nf_anytime_reference"
        self.assertIsNotNone(importlib.util.find_spec(name), "finite-panel implementation is absent")
        self.mod = importlib.import_module(name)

    def test_exact_expectation_matches_every_subset_with_ties_and_failures(self):
        values = [2., None, 2., 6., 12.]
        for k in range(len(values) + 1):
            outcomes = [min([10.] + [v for v in sample if v is not None])
                        for sample in combinations(values, k)]
            self.assertAlmostEqual(self.mod.expected_uniform_best(values, k, 10.),
                                   sum(outcomes) / len(outcomes), places=14)

    def test_empty_universe_and_zero_requests_keep_baseline(self):
        self.assertEqual(self.mod.expected_uniform_best([], 0, 10.), 10.)
        self.assertEqual(self.mod.expected_uniform_best([1., None], 0, 10.), 10.)

    def test_malformed_expectation_is_rejected(self):
        for values, k, baseline in [([1.], 2, 10.), ([1.], -1, 10.), ([1.], True, 10.),
                                     ([float("nan")], 1, 10.), ([0.], 1, 10.), ([1.], 1, 0.)]:
            with self.subTest(values=values, k=k, baseline=baseline), self.assertRaises(ValueError):
                self.mod.expected_uniform_best(values, k, baseline)

    def test_groups_follow_coarse_eligible_sequence_and_keep_tail_in_universe(self):
        rows = [row(0, 2., state(8.)), row(1, 1., state(1.), eligible=False),
                row(2, 4., state(7.)), row(3, 2., state(None, eligible=False)),
                row(4, 3., state(6.)), row(5, 1., state(2.))]
        result = self.mod.summarize_panel(rows, state(10.))
        self.assertEqual(result["coarse_eligible_ids"], ["0", "2", "3", "4", "5"])
        self.assertEqual(result["promoted_ids"], ["0"])
        self.assertEqual(result["unused_eligible_tail_ids"], ["5"])
        self.assertEqual(result["fine_requests_k"], 1)
        self.assertEqual(result["selected"]["id"], "0")
        self.assertEqual(result["reference_best"]["id"], "1")
        self.assertEqual(result["regret_relative"], 7.)
        self.assertAlmostEqual(result["uniform_expected_peak"], (8. + 7. + 10. + 6. + 2.) / 5)

    def test_fine_ineligible_promotion_still_uses_one_request(self):
        rows = [row(i, i + 1., state(None, eligible=False) if i == 0 else state(9.))
                for i in range(4)]
        result = self.mod.summarize_panel(rows, state(10.))
        self.assertEqual(result["fine_requests_k"], 1)
        self.assertEqual(result["selected"]["kind"], "original_wall")
        self.assertAlmostEqual(result["uniform_expected_peak"], 9.25)

    def test_no_complete_eligible_group_means_no_promotions(self):
        result = self.mod.summarize_panel([row(i, i + 1., state(1.)) for i in range(3)], state(10.))
        self.assertEqual(result["promoted_ids"], [])
        self.assertEqual(result["uniform_expected_peak"], 10.)

    def test_duplicate_order_or_invalid_eligible_measurement_is_rejected(self):
        for rows in [[row(0, 1., state(1.)), row(0, 2., state(2.))],
                     [row(0, 1., state(float("nan")))]]:
            with self.assertRaises(ValueError):
                self.mod.summarize_panel(rows, state(10.))

    def test_artifact_reads_require_matching_byte_hash_and_writes_are_exclusive(self):
        with tempfile.TemporaryDirectory() as folder:
            store = self.mod.EvidenceStore(Path(folder))
            store.write("panel.json", {"a": 1})
            self.assertEqual(store.read("panel.json"), {"a": 1})
            with self.assertRaises(FileExistsError):
                store.write("panel.json", {"a": 1})
            Path(folder, "panel.json").write_text('{"a": 1}\n', encoding="utf8")
            with self.assertRaisesRegex(ValueError, "hash"):
                store.read("panel.json")

    def test_unsealed_or_unsafe_reuse_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            store = self.mod.EvidenceStore(Path(folder))
            Path(folder, "unsealed.json").write_text("{}", encoding="utf8")
            with self.assertRaises((ValueError, FileNotFoundError)):
                store.read("unsealed.json")
            with self.assertRaises(ValueError):
                store.write("../outside.json", {})

    def test_panel_attempt_order_and_rejection_counts_are_verified(self):
        from d3d_optimization.campaign_geometry import Candidate
        candidate = Candidate("OMP", 20, (.001,) + (0.,) * 15)
        draws = iter([(None, "projection"), (candidate, None)])
        with patch.object(self.mod, "PANEL_SIZE", 1):
            panel = self.mod.generate_panel(None, 20611, draw_fn=lambda *_: next(draws))
            self.mod.validate_panel(panel, 20611)
            panel["rejected"] = {"projection": 2}
            with self.assertRaisesRegex(ValueError, "rejection"):
                self.mod.validate_panel(panel, 20611)
            panel["rejected"] = {"projection": 1}
            panel["attempt_records"][0]["attempt"] = 2
            with self.assertRaisesRegex(ValueError, "attempt"):
                self.mod.validate_panel(panel, 20611)

    def test_scoring_requires_sealed_panel_preserves_failures_and_resumes_without_rescoring(self):
        from d3d_optimization.campaign_geometry import Candidate
        candidates = [Candidate("OMP", 20, (.002*(i+1),) + (0.,)*15) for i in range(4)]
        panel = dict(complete=True, rows=[dict(i=i, id=c.identity, candidate=c.record())
                                         for i, c in enumerate(candidates)])
        with tempfile.TemporaryDirectory() as folder:
            store = self.mod.EvidenceStore(folder)
            store.write("protocol.json", {"test": True})
            calls = []

            def score(level, candidate):
                self.assertEqual(store.read("case/20611/panel.json"), panel)
                calls.append((level, candidate.identity))
                index = candidates.index(candidate)
                if level == "fine" and index == 0:
                    raise RuntimeError("synthetic transport failure")
                return state(index+1. if level == "coarse" else 9.)

            with self.assertRaises(FileNotFoundError):
                self.mod.score_panel(store, "case/20611", panel, {"coarse": "coarse", "fine": "fine"}, state(10.), score_fn=score)
            self.assertEqual(calls, [])
            store.write("case/20611/panel.json", panel)
            first = self.mod.score_panel(store, "case/20611", panel, {"coarse": "coarse", "fine": "fine"}, state(10.), score_fn=score)
            self.assertEqual(len(calls), 8)
            self.assertEqual(first["assessment"]["fine_requests_k"], 1)
            self.assertEqual(first["assessment"]["selected"]["kind"], "original_wall")
            self.assertIn("synthetic transport failure", first["rows"][0]["fine"]["errors"][0])
            second = self.mod.score_panel(store, "case/20611", panel, {}, state(10.), score_fn=lambda *_: self.fail("resume rescored"))
            self.assertEqual(second, first)
            score_file = Path(folder, next(iter(first["score_hashes"])))
            score_file.write_bytes(score_file.read_bytes() + b" ")
            with self.assertRaisesRegex(ValueError, "hash"):
                self.mod.score_panel(store, "case/20611", panel, {}, state(10.))

    def test_posttrial_caches_and_temp_stay_inside_study_runtime_directory(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {}, clear=False):
            self.assertTrue(hasattr(self.mod, "configure_environment"), "post-trial cache placement is absent")
            environment = self.mod.configure_environment(Path(folder), "reference")
            runtime = Path(folder, "posttrial-runtime", "reference").resolve()
            for key in ("CUPY_CACHE_DIR", "CUDA_CACHE_PATH", "TEMP", "TMP", "MPLCONFIGDIR"):
                self.assertTrue(Path(environment[key]).resolve().is_relative_to(runtime))
                self.assertTrue(Path(environment[key]).is_dir())
            for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE"):
                self.assertEqual(environment[key], "1")

    def test_live_panel_overlap_counts_unique_proposals_and_on_time_fine_observations(self):
        def candidate(value):
            record = dict(source="OMP", radius_mm=20, vector=[value] + [0.] * 15)
            return dict(id=self.mod.object_digest(record), candidate=record)

        a, b, c, other = [candidate(value) for value in (.001, .002, .003, .004)]
        panel = dict(case_id="example", seed=20611, rows=[dict(i=i, **row) for i, row in enumerate((a, b, c))])

        def scored(row, *, stage="fine", end=100., status="completed", cached=False, eligible=True):
            return dict(row, stage=stage, start_s=end-1., end_s=end, status=status, cached=cached,
                        measurement=state(8. if eligible else None, eligible=eligible))

        arms = [dict(case_id="example", seed=20611, policy="random_fine", audit_status="completed",
                     budget_s=600, candidates=[a, b, other], candidate_registry_complete=True,
                     events=[scored(a), scored(a), scored(b, stage="coarse"),
                             scored(b, cached=True), scored(b, end=601., status="late"),
                             scored(other, eligible=False)]),
                dict(case_id="example", seed=20627, policy="fine_adaptive", audit_status="completed",
                     budget_s=600, candidates=[c], candidate_registry_complete=True,
                     events=[scored(c, eligible=False)]),
                dict(case_id="different_equilibrium", seed=20611, policy="random_fine",
                     audit_status="completed", budget_s=600, candidates=[b], candidate_registry_complete=True,
                     events=[scored(b)])]
        overlap = self.mod.panel_live_overlap(panel, arms)
        self.assertEqual(set(overlap["live_proposed_ids"]), {r["id"] for r in (a, b, c, other)})
        self.assertEqual(overlap["live_proposed_count"], 4)
        self.assertEqual(set(overlap["live_on_time_fine_observed_ids"]), {r["id"] for r in (a, c, other)})
        self.assertEqual(overlap["live_on_time_fine_observed_count"], 3)
        self.assertEqual(overlap["proposed_intersection_ids"], [r["id"] for r in (a, b, c)])
        self.assertEqual(overlap["proposed_intersection_count"], 3)
        self.assertEqual(overlap["on_time_fine_intersection_ids"], [r["id"] for r in (a, c)])
        self.assertEqual(overlap["on_time_fine_intersection_count"], 2)
        self.assertEqual(len(overlap["by_arm"]), 2)
        self.assertTrue(overlap["complete_live_proposal_coverage"])
        self.assertFalse(overlap["historical_archive_overlap_checked"])
        self.assertFalse(overlap["held_out_validation_claim"])

        # Observed heat values and thermal eligibility do not affect membership
        # or this identity/timing-only provenance calculation.
        original_panel, changed_arms = copy.deepcopy(panel), copy.deepcopy(arms)
        for live_arm in changed_arms:
            for event in live_arm["events"]:
                event["measurement"] = state(1234., A2=9., eligible=not event["measurement"]["eligible"])
        self.assertEqual(self.mod.panel_live_overlap(panel, changed_arms), overlap)
        self.assertEqual(panel, original_panel)

    def test_failed_process_overlap_uses_recorded_ids_and_labels_incomplete_coverage(self):
        candidate = dict(source="OMP", radius_mm=20, vector=[.001] + [0.] * 15)
        identity = self.mod.object_digest(candidate)
        panel = dict(case_id="example", seed=20611, rows=[dict(i=0, id=identity, candidate=candidate)])
        recorded = dict(id=identity, candidate=candidate, stage="draw", status="completed", end_s=100.)
        arms = [dict(case_id="example", seed=20611, policy="random_fine", audit_status="process_failed",
                     budget_s=600, candidates=[], candidate_registry_complete=False, events=[recorded])]
        overlap = self.mod.panel_live_overlap(panel, arms)
        self.assertEqual(overlap["proposed_intersection_ids"], [identity])
        self.assertEqual(overlap["on_time_fine_intersection_ids"], [])
        self.assertFalse(overlap["complete_live_proposal_coverage"])
        self.assertFalse(overlap["complete_live_fine_observation_coverage"])


if __name__ == "__main__":
    unittest.main()
