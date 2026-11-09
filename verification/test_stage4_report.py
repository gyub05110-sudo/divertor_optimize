import numpy as np
import pytest


def test_every_regional_peak_is_bound_to_raw_samples_and_segments():
    from verification.stage4_report import verify_regional_samples
    wall = np.array([[0.,0.], [1.,0.], [1.,1.]])
    points = np.array([[.25,0.],[.75,0.],[1.,.4]])
    parent = np.array([0,0,1])
    peaks = np.array([3.,4.,8.])
    names = np.array(["floor","shelf"])
    stored = {"floor":{"peak_MW_m2":4.},"shelf":{"peak_MW_m2":8.}}
    verify_regional_samples(wall,names,points,parent,peaks,stored)
    with pytest.raises(ValueError,match="regional peak"):
        verify_regional_samples(wall,names,points,parent,peaks,dict(stored,floor={"peak_MW_m2":2.}))
    with pytest.raises(ValueError,match="sample coordinate"):
        verify_regional_samples(wall,names,points,np.array([0,1,1]),peaks,stored)


def row(case, arm, peak, *, predicted="floor-only", seed=30713, hit=True, status=None):
    reduction = None if peak is None else 100*(10-peak)/10
    snap=dict(q_native=peak,peak_MW_m2=peak,id=arm if peak is not None else None,reduction_pct=reduction)
    return dict(case_id=case,seed=seed,arm=arm,status=status or ("qualified" if peak else "no_design"),
                checkpoints={str(t):snap for t in (200,400,600)},
                predicted_family=predicted,abstained=(predicted=="full"), predicted_floor_sufficient=(
                    None if predicted=="full" else predicted=="floor-only"),
                baseline_measurable=True,baseline_peak_MW_m2=10.,first_5pct=dict(observed_s=185. if hit else None))


def test_prediction_scored_from_separate_forced_run_and_ties_within_one_percent():
    from verification.stage4_campaign import FOCUSED
    from verification.stage4_report import prediction_case
    arms=["hotspot","full",*FOCUSED]
    rows=[row("x",a,q,predicted="floor-shelf") for a,q in zip(arms,[1.,2.,2.02,2.5,2.3,2.])]
    p=prediction_case(rows)
    assert p["best_families"]==["floor-connector"]
    assert set(p["competitive_families"])=={"floor-connector","floor-only"}
    assert p["selected_forced_regret_pct"]==pytest.approx(15)
    assert p["selected_competitive"] is False and p["actual_floor_sufficient"] is True
    assert p["confusion"]=="false_negative"
    rows[-2]["status"]="no_design"
    assert prediction_case(rows)["scorable"] is False


def test_abstention_not_scored_as_a_prediction_and_misses_no_fake_time():
    from verification.stage4_campaign import FOCUSED
    from verification.stage4_report import prediction_case, product_pair
    rows=[row("x",a,3.,predicted="full") for a in ("hotspot","full",*FOCUSED)]
    assert prediction_case(rows)["scorable"] is False
    assert prediction_case(rows)["selected_competitive"] is None
    auto=row("x","hotspot",None,hit=False)
    full=row("x","full",8.)
    result=product_pair(auto,full,600)
    assert result["gap_pct"] is None and result["availability"]=="full_only"
    assert result["missed_full_5pct"] is True


def test_exact_margin_and_reduction_use_physical_thresholds():
    from verification.stage4_report import product_pair
    auto, full = row("x","hotspot",2.02), row("x","full",2.)
    assert product_pair(auto,full,600)["practical"] == "within_1pct"
    auto = row("x","hotspot",.95*19)
    auto["baseline_peak_MW_m2"] = 19.
    auto["checkpoints"]["600"]["reduction_pct"] = 4.9999999999999964
    assert product_pair(auto,full,600)["auto_5pct"] is True
    auto["status"] = "partial_result"
    assert product_pair(auto,full,600)["auto_q"] == .95*19


def test_control_label_independent_of_automatic_run_jitter():
    from verification.stage4_campaign import FOCUSED
    from verification.stage4_report import prediction_case
    rows=[row("x",a,q,predicted="floor-entry") for a,q in
          zip(("hotspot","full",*FOCUSED),[3.,4.,5.,4.,5.,5.])]
    prediction=prediction_case(rows)
    assert prediction["selected_forced_regret_pct"]==0
    assert prediction["always_floor_competitive"] is False
    assert prediction["confusion"]=="true_negative"
