"""Post-trial numerical verification uses immutable choices and genuine observations."""
from importlib import import_module, util
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from verification.test_nf_prospective_report import arm, baseline, score, state


class NumericsTests(unittest.TestCase):
    def setUp(self):
        name = "verification.nf_prospective_numerics"
        self.assertIsNotNone(util.find_spec(name), "prospective numerical implementation is absent")
        self.m = import_module(name)

    def test_lock_only_original_primary_and_long_choices_before_any_new_heat(self):
        early = score(.01, 8, 100)
        short = score(.02, 7, 500)
        long = score(.03, 6, 1600)
        late = score(.04, 4, 1802, status="late")
        direct = arm([baseline(), early, short, long, late])
        repeat = arm([baseline(), short], policy="coarse_adaptive")
        row = self.m.lock_case("holdout", [direct, repeat])
        self.assertEqual({c["id"] for c in row["candidates"]}, {"original_wall",short["id"],long["id"]})
        self.assertEqual({c["cutoff_s"] for c in row["choices"]}, {600, 1800})
        locked = next(c for c in row["candidates"] if c["id"] == short["id"])
        self.assertEqual(len(locked["reasons"]), 3)
        self.assertEqual(len(locked["native_observations"]), 2)
        self.assertTrue(all(a["event_index"] >= 0 for a in locked["native_observations"]))

    def test_no_prior_measurement_on_missing_or_killed_arm(self):
        elapsed = arm([baseline(), score(.01, 8, 100)], status="process_failed")
        row = self.m.lock_case("holdout", [elapsed])
        self.assertEqual([c["id"] for c in row["candidates"]], ["original_wall"])
        self.assertFalse(row["candidates"][0]["native_observations"])
        self.assertTrue(all(not c["available"] for c in row["choices"]))

    def test_four_paired_margins_with_overlapping_envelopes_and_ineligible_observable(self):
        a = {s:state(q) for s,q in zip(self.m.SETTINGS,[8,6,8,6])}
        b = {s:state(q) for s,q in zip(self.m.SETTINGS,[9,7,9,7])}
        both = self.m.compare_designs(a,b)
        self.assertEqual(both["all_settings_practical"]["1"], "fast")
        self.assertEqual(both["range_relation"], "overlap_or_touch")
        self.assertEqual(both["eligible_setting_count"], 4)
        b["quintic1280"]["eligible"] = False
        partial = self.m.compare_designs(a,b)
        self.assertEqual(partial["finite_setting_count"], 4)
        self.assertEqual(partial["eligible_setting_count"], 3)
        self.assertEqual(partial["all_settings_practical"]["1"], "incomplete")
        self.assertEqual(partial["range_relation"], "incomplete")
        self.assertIsNotNone(partial["settings"]["quintic1280"]["finite_margin"])
        self.assertIsNone(partial["settings"]["quintic1280"]["eligible_margin"])

    def test_identical_design_stays_tied_in_all_grids(self):
        vals = {s:state(q) for s,q in zip(self.m.SETTINGS,[8,8.2,8.1,8.4])}
        row = self.m.compare_designs(vals,vals,same_design=True)
        self.assertTrue(row["same_design"])
        self.assertEqual(row["eligible_gap_range_pct"], [0.,0.])
        self.assertEqual(row["range_relation"], "overlap_or_touch")
        self.assertEqual(row["all_settings_practical"]["0.5"], "within")

    def test_positive_native_replay_evidence_requires_measurement_and_observations(self):
        observation = [dict(seed=20711, measurement=state())]
        actual = dict(status="completed",measurement=state())
        self.assertTrue(self.m.native_replay(actual, observation)["observables_reproduced"])
        self.assertIsNone(self.m.native_replay(actual, [])["observables_reproduced"])
        self.assertIsNone(self.m.native_replay(actual, [])["passed"])
        failed = dict(status="failed", measurement=state(None, None, False))
        self.assertFalse(self.m.native_replay(failed, observation)["observables_reproduced"])
        ineligible = dict(status="completed",measurement=state(10.,2.,False))
        self.assertFalse(self.m.native_replay(ineligible, observation)["observables_reproduced"])
        changed = dict(status="completed",measurement=state(10.,2.1,True))
        self.assertFalse(self.m.native_replay(changed, observation)["observables_reproduced"])

    def test_fresh_baseline_only_scenario_keeps_original_method_choices(self):
        arms = [arm([baseline()], policy=p, seed=20729, budget=600) for p in ("coarse_adaptive", "fine_adaptive")]
        select = self.m.lock_case("holdout", arms)
        vals = {"original_wall": {k:dict(status="completed",measurement=state(),seconds=0.,error=None)
                                  for k in self.m.SETTINGS}}
        assessed = self.m.assess_case(select, vals)
        self.assertEqual(len(assessed["pairs"]), 1)
        self.assertEqual(assessed["pairs"][0]["eligible_gap_range_pct"], [0.,0.])
        self.assertEqual(assessed["candidates"][0]["native_replay"]["compared_count"], 2)
        self.assertFalse(assessed["candidates"][0]["baseline_comparison"]["all_settings_paired_improvement"])

    def test_scoring_requires_a_completed_separate_lock(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "lock"):
                self.m.require_locked(Path(folder))

    def test_engine_failure_preserves_actual_baseline_and_cannot_overwrite_failed_request(self):
        from verification.nf_anytime_reference import EvidenceStore
        s = score(.006, 7., 100)
        a = arm([baseline(),s], policy="fine_adaptive", seed=20729, budget=600)
        selected = self.m.lock_case("holdout",[a])
        with tempfile.TemporaryDirectory() as folder:
            store = EvidenceStore(folder)
            store.write("protocol.json", {"test": True})
            store.write("selection.json", {"test": True})
            with patch("d3d_optimization.boundary_evaluator.prepare_boundary", return_value=SimpleNamespace(baseline={})), \
                 patch("verification.nf_base_budget.measurement", return_value=state()), \
                 patch("d3d_optimization.campaign_geometry.evaluator", side_effect=RuntimeError("GPU unavailable")):
                first = self.m.score_setting(store,dict(case_id="holdout",root=str(folder)),selected,"cubic1280",None,"runs/0000")
            with patch("d3d_optimization.boundary_evaluator.prepare_boundary", side_effect=AssertionError("should reuse immutable evidence")):
                second = self.m.score_setting(store,dict(case_id="holdout",root=str(folder)),selected,"cubic1280",None,"runs/0001")
            self.assertEqual(first,second)
            self.assertEqual(first["original_wall"]["measurement"],state())
            self.assertEqual(first["original_wall"]["status"],"completed")
            self.assertEqual(first[s["id"]]["status"],"not_evaluated_preparation_failed")
            self.assertIn("GPU unavailable",first[s["id"]]["error"])
            store.verify_all()

    def test_power_accounting_reports_finite_ineligible_measured_ledger_separately(self):
        non = state(8.,1.9,False)
        non["ledger"].update(deposited=.98,unhit=.01,balance_error=.01)
        rows = [dict(status="completed",measurement=state()),dict(status="completed",measurement=non),
                dict(status="failed",measurement=state(None,None,False))]
        p = self.m.power_ledger_summary(rows)
        self.assertEqual((p["registered"],p["measured_ledgers"],p["eligible_ledgers"]),(3,2,1))
        self.assertAlmostEqual(p["max_relative_closure_all"],.01)
        self.assertEqual(p["max_relative_closure_eligible"],0.)

    def test_a_sealed_failure_cannot_be_used_as_an_eligible_measurement(self):
        with self.assertRaisesRegex(ValueError,"failed score"):
            self.m.validate_score(dict(status="failed",seconds=.1,measurement=state(),error="failed"))
        with self.assertRaisesRegex(ValueError,"clock"):
            self.m.validate_score(dict(status="completed",seconds=-.1,measurement=state(),error=None))

    def test_report_from_sealed_scores_rechecks_bytes_and_renders_both_registered_horizons(self):
        from verification.nf_anytime_reference import EvidenceStore
        r = self.m.report
        arms = []
        for seed in (20711,20729,20747):
            for policy in r.POLICIES:
                v,q = (.01,8.) if policy == r.FAST else (.02,9.)
                arms.append(arm([baseline(),score(v,q,100)],seed=seed,policy=policy,budget=1800 if seed == 20711 else 600))
        locked = self.m.lock_case("holdout",arms)
        selection = dict(cases=[locked],primary_evidence_hashes={})
        with tempfile.TemporaryDirectory() as temp:
            store = EvidenceStore(temp)
            store.write("protocol.json",{"test":True})
            store.write("selection.json",selection)
            for design in locked["candidates"]:
                for setting,params in self.m.SETTINGS.items():
                    outcome = design["native_observations"][0]["measurement"]
                    row = dict(case_id="holdout",id=design["id"],candidate=design["candidate"],setting=setting,
                               setting_tuple=list(params),protocol_sha256=store.sha256("protocol.json"),
                               selection_sha256=store.sha256("selection.json"),status="completed",seconds=.05,
                               error=None,measurement=outcome)
                    store.write(self.m.score_name("holdout",design["id"],setting),row)
            n = self.m.summarize(store,selection)
            self.assertEqual((n["registered_design_case_count"],n["registered_setting_requests"]),(3,12))
            live = r.summarize_live(arms)
            live.update(method_baselines={str(t):{key:r.baseline_totals([p[key] for p in live["pairs"] if p["cutoff_s"] == t])
                               for key in ("fast","direct")} for t in (200,400,600,1800)},
                        parent_seconds=0.,parent_costs=[r.summarize_cost(a) for a in arms],registered_arms=6,
                        with_process_receipt=6,closed_arms=6,status_counts={"completed":6})
            dev = dict(registered=6,completed_receipts=6,parent_seconds=3600.,selected={"population":24,"projected_seeds":8},
                       qualified=True,reasons=[],configurations=[])
            rendered = r.korean_report(dict(live=live,development=dev,numerical=n,errors=[],
                cases=[dict(case_id="holdout",time_ms=1400)],input_screening=dict(examined=1)))
            for label in ("600초 주 평가","1800초 첫 시드","1800초 선택의 수치 민감도","원벽 A₂","실제 적격 DE","전력수지"):
                self.assertIn(label,rendered)
            score_path = store.path(self.m.score_name("holdout","original_wall","cubic1280"))
            score_path.write_bytes(score_path.read_bytes()+b" ")
            with self.assertRaisesRegex(ValueError,"hash"):
                self.m.summarize(store,selection)


if __name__ == "__main__":
    unittest.main()
