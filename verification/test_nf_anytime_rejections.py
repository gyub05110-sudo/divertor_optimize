"""Synthetic tests for the separate post-hoc coarse-rejection diagnostic."""
import copy
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


def measurement(peak=8., A2=2., eligible=False):
    return dict(eligible=eligible, peak=peak, A2=A2, ledger={"balance_error": 0.},
                errors=[] if eligible else ["receiving_angle_constraint_failed"])


def event(i, **changes):
    candidate = dict(source="OMP", radius_mm=20, vector=[i/1000.] + [0.] * 15)
    identity = hashlib.sha256(json.dumps(candidate, sort_keys=True).encode()).hexdigest()
    row = dict(id=identity, candidate=candidate, stage="coarse", start_s=99., end_s=100.,
               status="completed", cached=False, physics_evaluated=True, measurement=measurement())
    row.update(changes)
    return row


def arm(events, *, case="case_a", seed=20611, incumbent=9.):
    return dict(case_id=case, seed=seed, policy="coarse_adaptive", budget_s=600,
                audit_status="completed", native_pool_eligible=True, events=events,
                checkpoints={"600": dict(available=True, kind="new_design", id="frozen_live_winner",
                                          peak=incumbent, A2=2.5)}, baseline=measurement(10., eligible=True))


def score(state, status="completed"):
    return dict(measurement=state, status=status, error=None if status == "completed" else "synthetic failure", seconds=.01)


class RejectionTests(unittest.TestCase):
    def setUp(self):
        name = "verification.nf_anytime_rejections"
        self.assertIsNotNone(importlib.util.find_spec(name), "post-hoc rejection audit implementation is absent")
        self.mod = importlib.import_module(name)

    def test_selection_is_all_fresh_on_time_coarse_thermal_rejections_only(self):
        included = event(1)
        other_error = event(2, measurement=measurement(100.))
        other_error["measurement"]["errors"] = ["another_physical_constraint"]
        rows = [included, other_error, event(3, measurement=measurement(1., eligible=True)),
                event(4, status="late", end_s=601.), event(5, cached=True),
                event(6, physics_evaluated=False), event(7, stage="fine"),
                event(8, status="failed"), event(9, end_s=601.), event(10, physics_evaluated=None)]
        selected = self.mod.select_rejections([arm(rows)])
        self.assertEqual(selected["selected_candidate_count"], 2)
        self.assertEqual(selected["qualifying_observation_count"], 2)
        self.assertEqual({r["id"] for r in selected["cases"][0]["candidates"]},
                         {self.mod.BASELINE_ID, included["id"], other_error["id"]})
        self.assertTrue(selected["post_hoc"])
        self.assertFalse(selected["primary_winners_changed"])

    def test_dedup_is_within_case_and_preserves_all_coarse_observations(self):
        rejected = event(1)
        cached = dict(rejected, cached=True)
        selected = self.mod.select_rejections([arm([rejected, cached]),
                                              arm([copy.deepcopy(rejected)], seed=20627, incumbent=7.),
                                              arm([copy.deepcopy(rejected)], case="case_b")])
        self.assertEqual(selected["selected_candidate_count"], 2)
        self.assertEqual(selected["qualifying_observation_count"], 3)
        row = next(r for r in selected["cases"][0]["candidates"] if r["id"] == rejected["id"])
        self.assertEqual(len(row["coarse_observations"]), 3)
        self.assertEqual(len(row["qualifying_observations"]), 2)
        self.assertEqual({o["incumbent_600"]["peak"] for o in row["qualifying_observations"]}, {7., 9.})

    def test_cutoff_is_inclusive_and_selection_never_filters_by_peak(self):
        a, b = event(1, end_s=600.), event(2, end_s=600.0001)
        first = self.mod.select_rejections([arm([a, b])])
        self.assertEqual(first["selected_candidate_count"], 1)
        changed = copy.deepcopy(a)
        changed["measurement"]["peak"] = 1e9
        second = self.mod.select_rejections([arm([changed, b])])
        self.assertEqual([r["id"] for r in first["cases"][0]["candidates"]],
                         [r["id"] for r in second["cases"][0]["candidates"]])

    def test_invalid_identity_and_nonfinite_completion_fail_closed(self):
        for invalid in (dict(event(1), id="changed"), dict(event(1), end_s=float("nan"))):
            with self.assertRaises(ValueError):
                self.mod.select_rejections([arm([invalid])])

    def test_fine_eligible_is_false_negative_even_if_it_does_not_improve_baseline(self):
        rejected = event(1)
        selected = self.mod.select_rejections([arm([rejected])])["cases"][0]
        for peak, improves in ((7., True), (11., False)):
            measured = {self.mod.BASELINE_ID: score(measurement(10., eligible=True)),
                        rejected["id"]: score(measurement(peak, eligible=True))}
            result = self.mod.assess_rejections(selected, measured)
            row = result["rows"][0]
            self.assertTrue(row["fine_eligible"])
            self.assertTrue(row["coarse_eligibility_false_negative"])
            self.assertEqual(row["vs_baseline"]["eligible_peak_improvement"], improves)
            self.assertEqual(row["vs_live_coarse_incumbents"][0]["reference"]["peak"], 9.)
            self.assertFalse(result["primary_winners_changed"])

    def test_fine_ineligibility_and_failed_evaluation_are_distinct(self):
        rejected = event(1)
        selected = self.mod.select_rejections([arm([rejected])])["cases"][0]
        base = {self.mod.BASELINE_ID: score(measurement(10., eligible=True))}
        ineligible = self.mod.assess_rejections(selected, dict(base, **{rejected["id"]: score(measurement(7.))}))["rows"][0]
        self.assertFalse(ineligible["coarse_eligibility_false_negative"])
        self.assertFalse(ineligible["vs_baseline"]["eligible_peak_improvement"])
        self.assertEqual(ineligible["fine"]["measurement"]["errors"], ["receiving_angle_constraint_failed"])
        failure = score(dict(eligible=False, peak=None, A2=None, errors=["synthetic failure"], ledger={}), "failed")
        failed = self.mod.assess_rejections(selected, dict(base, **{rejected["id"]: failure}))["rows"][0]
        self.assertIsNone(failed["fine_eligible"])
        self.assertIsNone(failed["coarse_eligibility_false_negative"])
        self.assertEqual(failed["diagnostic_status"], "fine_evaluation_failed")

    def test_missing_locked_score_is_not_replaced_or_silently_omitted(self):
        selected = self.mod.select_rejections([arm([event(1)])])["cases"][0]
        with self.assertRaises(ValueError):
            self.mod.assess_rejections(selected, {self.mod.BASELINE_ID: score(measurement(10., eligible=True))})

    def test_empty_selection_has_no_baseline_only_physical_case(self):
        selected = self.mod.select_rejections([arm([event(1, measurement=measurement(8., eligible=True))])])
        self.assertEqual(selected["cases"], [])
        self.assertEqual(selected["selected_candidate_count"], 0)

    def test_sealed_selection_precedes_native_only_runner_and_resume_preserves_failure(self):
        from verification.nf_anytime_reference import EvidenceStore
        selected = self.mod.select_rejections([arm([event(1)])])["cases"][0]
        with tempfile.TemporaryDirectory() as folder:
            store = EvidenceStore(folder)
            store.write("protocol.json", {"test": True})
            calls = []

            def runner(evidence, case, chosen, setting, reference, invocation):
                self.assertEqual(evidence.read("selection.json")["cases"], [selected])
                self.assertEqual(setting, "cubic1280")
                calls.append(setting)
                values, hashes = {}, {}
                for row in chosen["candidates"]:
                    failed = row["id"] != self.mod.BASELINE_ID
                    value = (dict(eligible=False, peak=None, A2=None, ledger={}, errors=["synthetic failure"])
                             if failed else measurement(10., eligible=True))
                    name = case["case_id"] + "/scores/" + row["id"] + ".cubic1280.json"
                    record = dict(id=row["id"], candidate=row["candidate"], case_id=case["case_id"],
                                  setting=setting, protocol_sha256=evidence.sha256("protocol.json"),
                                  selection_sha256=evidence.sha256("selection.json"),
                                  **score(value, "failed" if failed else "completed"))
                    evidence.write(name, record)
                    values[row["id"]] = value
                    hashes[name] = evidence.sha256(name)
                return values, hashes

            with self.assertRaises(FileNotFoundError):
                self.mod.run_selected_case(store, {"case_id": "case_a"}, selected, None, "runs/0000", runner=runner)
            self.assertEqual(calls, [])
            store.write("selection.json", {"cases": [selected]})
            first = self.mod.run_selected_case(store, {"case_id": "case_a"}, selected, None, "runs/0000", runner=runner)
            second = self.mod.run_selected_case(store, {"case_id": "case_a"}, selected, None, "runs/0001",
                                                runner=lambda *_: self.fail("failed request was replaced"))
            self.assertEqual(first, second)
            self.assertEqual(calls, ["cubic1280"])
            self.assertIsNone(first["rows"][0]["coarse_eligibility_false_negative"])
            path = store.path(next(iter(first["score_hashes"])))
            path.write_bytes(path.read_bytes() + b" ")
            with self.assertRaisesRegex(ValueError, "hash"):
                self.mod.run_selected_case(store, {"case_id": "case_a"}, selected, None, "runs/0002")

    def test_protocol_locks_the_imported_runner_and_both_shared_source_dependencies(self):
        from verification import nf_anytime_uncertainty
        self.assertIs(self.mod.run_setting, nf_anytime_uncertainty.run_setting)
        self.assertEqual({p.name for p in self.mod.source_paths()},
                         {"nf_anytime_rejections.py", "nf_anytime_reference.py", "nf_anytime_uncertainty.py"})
        details = self.mod.protocol_details()
        self.assertEqual(details["setting"], ["cubic", 1280, .000125, 641, .00025])
        self.assertEqual(details["imported_runner"]["qualified_name"], "verification.nf_anytime_uncertainty.run_setting")
        self.assertTrue(details["post_hoc"])
        self.assertFalse(details["primary_winners_changed"])


if __name__ == "__main__":
    unittest.main()
