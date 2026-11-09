"""Whole-grid saddle enclosure must retain pairs hidden from a cell-corner sign test."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

from d3d_optimization.frozen_field import FrozenField


def test_two_stationary_points_hidden_in_one_cubic_span():
    assert importlib.util.find_spec("verification.stage4_xpoint_diagnostics"), "whole-grid diagnostic is missing"
    from verification.stage4_xpoint_diagnostics import stationary_enclosures
    r, z = np.linspace(1., 3., 9), np.linspace(-1., 1., 9)
    rr, zz = np.meshgrid(r, z)
    # Both roots of dpsi/dR lie within one cell bordered by equal signs.
    x = rr - 2.1
    field = FrozenField(r, z, x**3/3 - .03**2*x + (zz-.3)**2, 1.)
    boxes, roots, stats = stationary_enclosures(field)
    assert len(roots) == 2
    np.testing.assert_allclose(sorted(roots, key=lambda p: p[0]), [[2.07, .3], [2.13, .3]], atol=1e-8)
    assert stats["unmatched_terminal_boxes"] == 0
    for root in roots:
        assert np.any(np.all((boxes[:, :2] <= root) & (root <= boxes[:, 2:]), axis=1))


def test_development_census_localizes_all_polished_roots_and_clears_lower_search_ball():
    from verification.stage4_xpoint_diagnostics import run
    study = Path("D:/DivertorOptimization/studies/nf-stage3-20261108")
    if not (study/"inputs"/"raw"/"d3d_shot_8d1575d8c4.parquet").is_file():
        pytest.skip("requires archived Stage 3 development discharges")
    for case in run(study)["cases"]:
        census = case["census"]
        assert census["unmatched_terminal_boxes"] == 0
        assert census["possible_stationary_boxes_within_lower_search_radius"] == 0
        assert case["upper_lcfs_evidence"]["classification"] == "upper_lcfs_saddle"
