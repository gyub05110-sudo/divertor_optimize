"""No transport calculations. Adversarial cases for the prospective endpoints."""
from importlib import import_module, util
import json
from pathlib import Path
import tempfile
import unittest

from verification.nf_anytime_reference import object_digest


def state(q=10., a=2., eligible=True):
    return dict(eligible=eligible, peak=q, A2=a, errors=[] if eligible else ["geometry"],
                ledger={} if q is None else dict(prescribed=1., deposited=1., unhit=0., balance_error=0.))


def baseline(end=3.):
    return dict(stage="baseline", start_s=end-1., end_s=end, status="completed", measurement=state())


def score(v, q, end, *, a=2., stage="fine", kind="de_trial", status="completed",
          cached=False, eligible=True):
    candidate = dict(source="OMP", radius_mm=20, vector=[v]+[0.]*15)
    return dict(stage=stage, start_s=end-1., end_s=end, status=status, id=object_digest(candidate),
                candidate=candidate, proposal_kind=kind, cached=cached, physics_evaluated=q is not None,
                measurement=state(q, a, eligible))


def arm(events, *, policy="fine_adaptive", seed=20711, budget=1800, status="completed", elapsed=None):
    return dict(case_id="holdout", seed=seed, policy=policy, budget_s=budget, events=events,
                status=status, elapsed_s=budget if elapsed is None else elapsed, error=None,
                raw_certified=status in ("completed", "arm_failed", "early_stopped"))


def module(case):
    name = "verification.nf_prospective_report"
    case.assertIsNotNone(util.find_spec(name), "prospective reporting implementation is absent")
    return import_module(name)


class ProspectiveReportTests(unittest.TestCase):
    def setUp(self):
        self.m = module(self)

    def test_one_fixed_denominator_in_both_directions_inclusive_boundaries(self):
        for fast, key in [(99., "fast"), (101., "direct"), (99.5, "within"), (100., "within")]:
            row = self.m.margin(fast, 100.)
            self.assertEqual(row["practical"]["1"], key)
            self.assertEqual(row["gap_pct"], 100.-fast)
        self.assertEqual(self.m.margin(99.5, 100.)["practical"]["0.5"], "fast")
        self.assertEqual(self.m.margin(100., 100.)["sign"], "tie")
        self.assertEqual(self.m.margin(101., 100.)["practical"]["2"], "within")

    def test_primary_uses_raw_pre_cutoff_fine_only_and_observes_smaller_A2(self):
        events = [baseline(), score(.01, 8, 40, a=1.8, kind="population_initialization"),
                  score(.02, 1, 50, stage="coarse"), score(.07, 2, 100, cached=True),
                  score(.03, 6, 605, a=3), score(.04, 5, 1801, status="late")]
        a = arm(events)
        first, long = self.m.endpoint(a, 600), self.m.endpoint(a, 1800)
        self.assertEqual(first["selected_id"], events[1]["id"])
        self.assertEqual(first["q_ratio"], .8)
        self.assertTrue(first["A2_decreased"])
        self.assertFalse(first["joint_improvement"])
        self.assertEqual(first["first_5pct"]["time_s"], 40)
        self.assertEqual(first["de"]["eligible_deposition"], 0)
        self.assertEqual(long["selected_id"], events[4]["id"])
        self.assertEqual(long["de"]["eligible_deposition"], 1)
        self.assertEqual(long["power_relative_error"], 0.)

    def test_de_fresh_geometry_and_returned_ineligible_not_eligible_participation(self):
        events = [baseline(), score(.01, None, 30, a=None, eligible=False),
                  score(.02, 7, 40, eligible=False), score(.03, 6, 50, cached=True),
                  score(.04, 8, 60, kind="coordinate_polish"), score(.05, 4, 601)]
        row = self.m.endpoint(arm(events), 600)
        self.assertEqual([row["de"][key] for key in ("requests", "fresh", "deposition", "eligible_deposition")],
                         [3, 2, 1, 0])
        self.assertEqual(row["phases"]["coordinate_polish"]["fresh"], 1)
        self.assertEqual(row["phases"]["coordinate_polish"]["completed_seconds"], 1.)

    def test_uncertified_process_clock_is_exposure_and_never_point_censor_evidence(self):
        a = arm([baseline(),score(.02,7,44)],status="process_failed",elapsed=480.)
        first = self.m.endpoint(a,600)["first_5pct"]
        self.assertEqual(first["reason"],"process_failed")
        self.assertEqual(first["exposure_end_s"],480.)
        self.assertIsNone(first["censor_time_s"])
        self.assertFalse(first["right_censored"])

    def test_geometry_rejection_does_not_subsume_a_measured_ineligible_de_trial(self):
        geom = score(.01,None,40,a=None,eligible=False)
        geom["measurement"]["errors"] = ["continuous_shape_constraint_failed"]
        heat = score(.02,6,70,eligible=False)
        heat["measurement"]["errors"] = ["receiving_angle_constraint_failed"]
        de = self.m.endpoint(arm([baseline(),geom,heat]),600)["de"]
        self.assertEqual((de["geometry_rejections"],de["returned_ineligible"],de["eligible_deposition"]),(1,1,0))

    def test_failed_long_arm_retains_earlier_certified_cutoff_and_censors_later(self):
        fout = dict(stage="progress_io", start_s=700., end_s=701., status="failed", error="disk")
        a = arm([baseline(), score(.01, 9.9, 100), fout], status="arm_failed", elapsed=703.)
        before, after = self.m.endpoint(a, 600), self.m.endpoint(a, 1800)
        self.assertTrue(before["available"])
        self.assertFalse(after["available"])
        self.assertEqual(after["observed_before_stop"]["peak"], 9.9)
        self.assertFalse(after["first_5pct"]["event_observed"])
        self.assertEqual(after["first_5pct"]["censor_time_s"], 701.)
        self.assertEqual(after["first_5pct"]["reason"], "arm_failed")

    def test_failed_process_and_missing_baseline_never_impute_design(self):
        unc = self.m.endpoint(arm([baseline(), score(.01, 2, 20)], status="process_failed"), 600)
        self.assertFalse(unc["available"])
        self.assertIsNone(unc["selected_id"])
        self.assertIsNone(unc["first_5pct"]["time_s"])
        no_base = self.m.endpoint(arm([]), 600)
        self.assertFalse(no_base["available"])
        self.assertEqual(no_base["baseline_status"], "unavailable")
        self.assertEqual(no_base["first_5pct"]["reason"], "baseline_unavailable")

    def test_time_matched_pair_and_eligibility_denominator(self):
        fast = self.m.endpoint(arm([baseline(), score(.01, 8, 40)], policy="coarse_adaptive"), 600)
        direct = self.m.endpoint(arm([baseline(), score(.02, 9, 50)]), 600)
        p = self.m.pair_endpoints(fast, direct)
        self.assertAlmostEqual(p["margin"]["gap_pct"], 100/9)
        self.assertTrue(p["fast_practical_and_direct_de"])
        self.assertFalse(p["same_design"])
        missing = self.m.endpoint(arm([], status="missing"), 600)
        miss = self.m.pair_endpoints(fast, missing)
        self.assertFalse(miss["available"])
        self.assertIsNone(miss["margin"])
        counts = self.m.pair_counts([p, miss])
        self.assertEqual((counts["registered"], counts["comparable"], counts["unavailable"]), (2, 1, 1))
        self.assertEqual(counts["thresholds"]["1"]["fast"], 1)
        with self.assertRaises(ValueError):
            self.m.pair_endpoints(fast, self.m.endpoint(arm([]), 1800))

    def test_geometry_duplicates_and_selection_schedule_are_case_clustered(self):
        single = score(.011, 8, 100)
        arms = [arm([baseline(), single], seed=s, budget=(1800 if s == 20711 else 600), policy=p)
                for s in (20711, 20729, 20747) for p in ("coarse_adaptive", "fine_adaptive")]
        result = self.m.summarize_live(arms)
        groups = {g["cutoff_s"]: g for g in result["case_groups"]}
        self.assertEqual(len(result["pairs"]), 10)  # 3x(200,400,600) and one long pair
        self.assertEqual(groups[600]["registered"], 3)
        self.assertEqual(groups[1800]["registered"], 1)
        self.assertEqual(groups[600]["geometry"]["unique_identities"], 1)
        self.assertEqual(groups[600]["geometry"]["selection_count"], 6)
        self.assertEqual(groups[600]["thresholds"]["1"]["within"], 3)
        self.assertEqual(result["by_cutoff"]["600"]["distinct_equilibria"], 1)

    def test_saved_checkpoint_mutation_invalidates_result(self):
        from verification.nf_anytime_final_audit import checkpoint
        raw = [baseline(), score(.01, 8., 301)]
        rec = dict(events=raw, budget_s=1800, baseline=raw[0]["measurement"], baseline_ready_s=3.,
                   checkpoints={str(t): checkpoint(raw, t, 1800) for t in (200,400,600,1800)},
                   final=checkpoint(raw, 1800, 1800))
        self.m.check_endpoints(rec, raw)
        rec["checkpoints"]["600"]["incumbent"]["peak"] = 7.
        with self.assertRaisesRegex(ValueError, "checkpoint"):
            self.m.check_endpoints(rec, raw)

    def test_corrupt_dead_process_hash_and_unattempted_job_are_distinct(self):
        with tempfile.TemporaryDirectory() as temp:
            study = Path(temp)
            job = dict(phase="holdout", case_id="holdout", seed=20711, policy="fine_adaptive",
                       budget_s=1800, population=24, projected_seeds=8, input_root="staged",
                       relative_dir="live/holdout/20711/fine_adaptive")
            d = study/job["relative_dir"]
            expected = "a"*64
            missing = self.m.load_arm(study, job, expected, {})
            self.assertEqual(missing["status"], "missing")
            d.mkdir(parents=True)
            (d/"job.json").write_text(json.dumps(dict(job,registration_sha256=expected)), encoding="utf8")
            (d/"events.jsonl").write_text(json.dumps(baseline())+"\n" + '{"torn":', encoding="utf8")
            (d/"stdout.log").write_bytes(b"")
            proc = dict(returncode=9, registration_sha256=expected, job_sha256=self.m.audit.digest(d/"job.json"),
                        launch_ns=5, exit_observed_ns=3000000005, elapsed_s=3., result_sha256=None,
                        events_sha256=self.m.audit.digest(d/"events.jsonl"), stdout_sha256=self.m.audit.digest(d/"stdout.log"))
            (d/"process.json").write_text(json.dumps(proc), encoding="utf8")
            failed = self.m.load_arm(study, job, expected, {})
            self.assertEqual(failed["status"], "process_failed")
            self.assertTrue(failed["partial_event_tail"])
            self.assertFalse(self.m.endpoint(failed, 600)["available"])
            (d/"events.jsonl").write_bytes(b"")
            bad = self.m.load_arm(study, job, expected, {})
            self.assertEqual(bad["status"], "audit_failed")
            self.assertIn("hash", bad["error"])

    def test_costs_share_no_budget_credit_across_checkpoints(self):
        a = arm([baseline(), score(.01, 8, 8), score(.02, 7, 1801, status="late")])
        a.update(process=dict(elapsed_s=1805.), elapsed_s=1802., startup_s=.5)
        result = self.m.summarize_cost(a)
        self.assertEqual(result["parent_elapsed_s"], 1805.)
        self.assertEqual(result["budget_s"], 1800)
        self.assertEqual(result["late_scores"], 1)
        self.assertEqual(result["baseline_ready_s"], 3.)

    def test_case_median_requires_complete_seed_set_and_tuned_de_count_is_separate(self):
        arms = []
        for seed in (20711,20729,20747):
            budget = 1800 if seed == 20711 else 600
            arms.append(arm([baseline(),score(.01,8,40)], seed=seed,budget=budget,policy="coarse_adaptive"))
            status = "missing" if seed == 20747 else "completed"
            arms.append(arm([baseline(),score(.02,9,44,kind="population_initialization")],seed=seed,budget=budget,status=status))
        live = self.m.summarize_live(arms)
        reg = live["by_cutoff"]["600"]
        self.assertEqual((reg["registered"],reg["comparable"],reg["complete_equilibria"]),(3,2,0))
        self.assertEqual(reg["equilibria_fast_median_1pct"],0)
        self.assertEqual(reg["fast_practical_and_direct_de"],0)
        result = self.m.baseline_totals([p["direct"] for p in live["pairs"] if p["cutoff_s"] == 600])
        self.assertEqual((result["registered"],result["available"],result["improved_q"],result["q_down_A2_up"]),(3,2,2,0))
        self.assertEqual(result["first_5pct_observed"],2)

    def test_report_exposes_unavailable_setting_and_preserves_published_file(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"REPORT.ko.md"
            self.m.write_once(path,b"first report\n")
            self.m.write_once(path,b"first report\n")
            with self.assertRaisesRegex(ValueError,"different"):
                self.m.write_once(path,b"different\n")
            self.assertEqual(path.read_bytes(),b"first report\n")
        self.assertEqual(self.m.fmt(None),"—")

    def test_unqualified_development_report_does_not_invent_holdout_results(self):
        result = dict(status="development_unqualified", live=None, numerical=dict(status="not_locked"), errors=[],
            development=dict(registered=6,completed_receipts=6,parent_seconds=3660.,selected=None,
                qualified=False,reasons=["no eligible direct comparator"],configurations=[],rows=[]))
        text = self.m.korean_report(result)
        self.assertIn("미등록",text)
        self.assertIn("no eligible direct comparator",text)
        self.assertNotIn("9/9",text)

    def test_primary_original_wall_matches_with_native_replay_precision(self):
        a = self.m.endpoint(arm([baseline()],policy="coarse_adaptive"),600)
        b = self.m.endpoint(arm([baseline()]),600)
        a["baseline"]["peak"] *= (1.+1e-12)
        p = self.m.pair_endpoints(a,b)
        self.assertTrue(p["available"])
        self.assertFalse(p["baseline_exact"])
        a["baseline"]["A2"] *= 1.01
        self.assertFalse(self.m.pair_endpoints(a,b)["available"])

    def test_end_to_end_valid_receipt_reconstructed_from_flushed_worker_events(self):
        from verification.test_nf_prospective_worker import arm as sample_worker
        from verification.nf_anytime_reference import EvidenceStore
        job = dict(phase="holdout",case_id="holdout",seed=20729,policy="fine_adaptive",budget_s=600,
                   population=24,projected_seeds=8,input_root="staged",relative_dir="live/holdout/20729/fine_adaptive")
        registration = "a"*64
        with tempfile.TemporaryDirectory() as temp:
            study = Path(temp)
            d = study/job["relative_dir"]
            e = EvidenceStore(d)
            e.write("job.json",dict(job,registration_sha256=registration))
            rec = sample_worker()
            rec.update(seed=job["seed"],phase="holdout",case_id=job["case_id"],input_root="staged",
                       launch_ns=8000,job_sha256=e.sha256("job.json"),registration_sha256=registration)
            (d/"result.json").write_text(json.dumps(rec),encoding="utf8")
            (d/"events.jsonl").write_text("".join(json.dumps(r)+"\n" for r in rec["events"]),encoding="utf8")
            (d/"stdout.log").write_bytes(b"")
            env = {k:str(d/"cache"/v) for k,v in dict(CUPY_CACHE_DIR="cupy",CUDA_CACHE_PATH="cuda",TEMP="tmp",TMP="tmp",
                                                      MPLCONFIGDIR="matplotlib",HF_HOME="huggingface").items()}
            env.update(NF_PROSPECTIVE_JOB_SHA256=e.sha256("job.json"),OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
            e.write("process.json",dict(environment=env,launch_ns=8000,exit_observed_ns=604000008000,elapsed_s=604.,
                    returncode=0,error=None,registration_sha256=registration,job_sha256=e.sha256("job.json"),
                    events_sha256=self.m.audit.digest(d/"events.jsonl"),result_sha256=self.m.audit.digest(d/"result.json"),
                    stdout_sha256=self.m.audit.digest(d/"stdout.log")))
            hashes = {}
            actual = self.m.load_arm(study,job,registration,hashes,sealed=True)
            self.assertEqual(actual["status"],"completed",actual["error"])
            self.assertEqual(self.m.endpoint(actual,600)["de"]["eligible_deposition"],1)
            self.assertEqual(self.m.endpoint(actual,600)["peak"],8.)
            self.assertEqual(len(hashes),5)


if __name__ == "__main__":
    unittest.main()
