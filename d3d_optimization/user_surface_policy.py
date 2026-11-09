"""Choose a screened proposal family from the original fine wall's normal heat load."""
from math import isfinite


# This is the final tie-break order, independent of coefficient block order.
REGION_FAMILIES = {
    "floor": "floor-only",
    "omp_connector": "floor-connector",
    "shelf": "floor-shelf",
    "omp_entry": "floor-entry",
}
SEARCH_FAMILIES = ("full", *REGION_FAMILIES.values())
REQUESTED_FAMILIES = (*SEARCH_FAMILIES, "hotspot")
RELATIVE_PEAK_TIE = 1e-10


def choose_surface_policy(rows, requested_family):
    """Use regional peaks in MW/m² and deposited watts as returned by surface_metrics."""
    if requested_family not in REQUESTED_FAMILIES:
        raise ValueError(f"unknown screened family: {requested_family}")
    if not {"fixed", *REGION_FAMILIES}.issubset(rows):
        raise ValueError("missing regional fine baseline measurements")
    metrics = {}
    for name, row in rows.items():
        try:
            peak, power = float(row["peak_MW_m2"]), float(row["deposited_W"])
        except (KeyError, ValueError, TypeError) as exc:
            raise ValueError(f"incomplete regional fine baseline measurements: {name}") from exc
        if not isfinite(peak) or not isfinite(power) or min(peak, power) < 0.:
            raise ValueError(f"invalid regional fine baseline measurements: {name}")
        metrics[name] = {"peak_MW_m2": peak, "deposited_W": power}
    maximum = max(row["peak_MW_m2"] for row in metrics.values())
    if maximum <= 0.:
        raise ValueError("regional fine baseline has no positive peak")
    tolerance = RELATIVE_PEAK_TIE * maximum
    tied = sorted(name for name, row in metrics.items() if maximum - row["peak_MW_m2"] <= tolerance)
    if any(name not in REGION_FAMILIES for name in tied):
        region = None
        predicted = "full"
        reason = "unadjustable_region_in_peak_tie"
    else:
        # max retains the first of equal deposited powers in the frozen priority order.
        region = max((name for name in REGION_FAMILIES if name in tied),
                     key=lambda name: metrics[name]["deposited_W"])
        predicted = REGION_FAMILIES[region]
        reason = None
    return dict(policy_version="regional-hotspot-v1", source="original_native_fine_baseline.surface_metrics",
                global_peak_MW_m2=maximum, tie_tolerance_MW_m2=tolerance, regional_metrics=metrics,
                tied_regions=tied, predicted_region=region, predicted_family=predicted,
                predicted_floor_sufficient=None if region is None else region == "floor",
                fallback_reason=reason, requested_family=requested_family,
                selected_family=predicted if requested_family == "hotspot" else requested_family)
