"""Exercise the Stage 4 figure with in-memory fixtures only."""
import importlib.util
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from verification.stage4_campaign import FOCUSED, SEEDS
from verification.stage4_report import prediction_case, product_pair
from verification.test_stage4_report import row


SCRIPT = Path(__file__).resolve().parents[1]/"results/nf-stage4-20261109/plot_policy_validation.py"


def plotter():
    assert SCRIPT.is_file(), "Stage 4 publication figure script missing"
    spec = importlib.util.spec_from_file_location("policy_validation_figure", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sample(name, paired, forced, forecast):
    group = [row(name, family, value, predicted=forecast) for family, value in zip(FOCUSED, forced)]
    pairs = []
    for seed, (auto, full) in zip(SEEDS, paired):
        a, b = row(name, "hotspot", auto, seed=seed, predicted=forecast), row(name, "full", full, seed=seed)
        group.extend([a, b])
        pairs.extend(product_pair(a, b, t) for t in (200, 400, 600))
    return pairs, prediction_case(group)


def report(*samples):
    return dict(complete=True, schema="nf-stage4-forecast-report.v1",
                pairs=[pair for pairs, _ in samples for pair in pairs],
                predictions=[p for _, p in samples])


def test_exact_sign_ranges_and_independent_forced_forecast():
    data = report(sample("fixture-a", [(9.5, 10), (11, 10), (10.2, 10)], [2.02, 2.5, 2.3, 2.], "floor-shelf"),
                  sample("fixture-b", [(None, None)]*3, [2., 2., 2., 2.], "full"),
                  sample("fixture-c", [(8., 10)]*2+[(None, 10)], [2., 3., 2.5, None], "floor-entry"))
    values = plotter().figure_data(data)
    assert values["total"] == 3 and len(values["rows"]) == 2
    first, second = values["rows"]
    assert first["seed_range"] == pytest.approx((-5., 2., 10.))
    assert first["forced_gaps"] == pytest.approx((1., 25., 15., 0.))
    assert first["forecast"] == "floor-shelf" and first["competitive"] == {"floor-only", "floor-connector"}
    assert second["seed_range"] is None and second["forecast"] is None
    assert second["minima"] == set(FOCUSED)
    assert (values["scorable"], values["forecast_hits"], values["floor_hits"]) == (1, 0, 1)


def test_empty_and_one_class_are_plottable_without_writing_results():
    backend = plotter()
    missing = sample("fixture-a", [(None, None)]*3, [None]*4, None)
    floor = sample("fixture-b", [(10., 10.)]*3, [5., 6., 6., 5.], "floor-only")
    for data in (report(missing), report(floor)):
        figure = backend.make_figure(backend.figure_data(data))
        try:
            if data["predictions"][0]["scorable"]:
                assert any("Floor sufficient for all scored inputs" in item.get_text() for item in figure.texts)
            for extension in ("png", "pdf"):
                output = BytesIO()
                figure.savefig(output, format=extension)
                assert output.tell() > 2000
        finally:
            plt.close(figure)
    with pytest.raises(ValueError, match="complete"):
        backend.figure_data(dict(report(missing), complete=False))
