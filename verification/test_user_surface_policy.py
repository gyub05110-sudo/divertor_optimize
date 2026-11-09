"""Presearch family decisions from the original fine wall's regional measurements."""
import importlib.util

import pytest


def decide(rows, requested="hotspot"):
    assert importlib.util.find_spec("d3d_optimization.user_surface_policy"), "surface policy missing"
    from d3d_optimization.user_surface_policy import choose_surface_policy

    return choose_surface_policy(rows, requested)


def measurements(peaks, watts=None):
    watts = watts or {}
    names = ("fixed", "omp_entry", "shelf", "omp_connector", "floor", "imp_sloped", "imp_vertical")
    return {name: dict(peak_MW_m2=peaks.get(name, 0.), deposited_W=watts.get(name, 1.),
                       first_hit_W=1000. if name == "floor" else 0.) for name in names}


@pytest.mark.parametrize("region,family", [
    ("floor", "floor-only"), ("omp_entry", "floor-entry"), ("shelf", "floor-shelf"),
    ("omp_connector", "floor-connector"),
])
def test_each_adjustable_hotspot_selects_its_family(region, family):
    rows = measurements({region: 2.5})
    policy = decide(rows)
    assert policy["source"] == "original_native_fine_baseline.surface_metrics"
    assert policy["global_peak_MW_m2"] == 2.5
    assert policy["tie_tolerance_MW_m2"] == pytest.approx(2.5e-10)
    assert policy["tied_regions"] == [region]
    assert policy["predicted_region"] == region
    assert policy["predicted_family"] == policy["selected_family"] == family
    assert policy["predicted_floor_sufficient"] is (region == "floor")
    assert policy["fallback_reason"] is None
    assert policy["regional_metrics"] == {
        k: {field: v[field] for field in ("peak_MW_m2", "deposited_W")} for k, v in rows.items()}


@pytest.mark.parametrize("region", ["fixed", "imp_sloped", "imp_vertical"])
def test_even_a_near_tie_with_fixed_or_imp_abstains(region):
    rows = measurements({"floor": 10., region: 10.-5e-10}, {"floor": 900., region: 0.})
    policy = decide(rows)
    assert set(policy["tied_regions"]) == {"floor", region}
    assert policy["predicted_floor_sufficient"] is policy["predicted_region"] is None
    assert policy["selected_family"] == policy["predicted_family"] == "full"
    assert policy["fallback_reason"] == "unadjustable_region_in_peak_tie"


def test_deposited_power_breaks_peak_ties_before_order_and_excludes_values_outside_tie():
    rows = measurements({"floor": 10., "shelf": 10.-5e-10, "omp_entry": 10.-2e-9},
                        {"floor": 5., "shelf": 9., "omp_entry": 1000.})
    policy = decide(rows)
    assert set(policy["tied_regions"]) == {"floor", "shelf"}
    assert policy["predicted_region"] == "shelf"


@pytest.mark.parametrize("ties,winner", [
    (["shelf", "omp_entry", "omp_connector", "floor"], "floor"),
    (["shelf", "omp_entry", "omp_connector"], "omp_connector"),
    (["shelf", "omp_entry"], "shelf"),
])
def test_equal_watts_use_stable_priority_independent_of_input_order(ties, winner):
    rows = measurements(dict.fromkeys(ties, 5.))
    assert decide(dict(reversed(rows.items())))["predicted_region"] == winner


def test_all_forced_controls_preserve_identical_presearch_evidence():
    rows = measurements({"shelf": 3.})
    expected = decide(rows)
    for family in ("full", "floor-only", "floor-entry", "floor-shelf", "floor-connector"):
        assert decide(rows, family) == dict(expected, requested_family=family, selected_family=family)


@pytest.mark.parametrize("bad", [{}, {"floor": {"deposited_W": 10.}},
                                  measurements({"floor": float("nan")}),
                                  measurements({"floor": 1.}, {"floor": -1.})])
def test_incomplete_or_nonfinite_evidence_fails_explicitly(bad):
    with pytest.raises(ValueError, match="regional"):
        decide(bad)
