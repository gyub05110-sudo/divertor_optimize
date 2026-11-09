"""The companion figure maps each native sample through its own material parent."""
import importlib.util
from pathlib import Path

import numpy as np


def test_native_samples_keep_parent_identity_and_project_to_original_length():
    path = Path(__file__).resolve().parents[1]/"results/nf-stage4-20261109/plot_hotspot_profiles.py"
    assert path.is_file(), "native-profile publication figure is missing"
    spec = importlib.util.spec_from_file_location("native_hotspot_plot", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    base = np.array([[0., 0.], [1., 0.], [2., 0.], [3., 0.]])
    actual = np.array([[0., 0.], [1., 0.], [2., .5], [3., 0.]])
    # Same Cartesian boundary with different limiting samples on the adjoining
    # face. An XY-based nearest-segment guess or a cell-average lookup loses it.
    p = np.array([0, 1, 2, 1, 1])
    points = actual[p] + np.array([1., 0., .3, .5, 1.])[:, None]*(actual[p+1]-actual[p])
    q = np.array([9., 7., 5., 2., 6.])
    x, samples, length = module.native_face_samples(base, actual, np.array(["fixed", "shelf", "fixed"]),
                                                    p, points, q, "shelf")
    np.testing.assert_allclose(x, [0., 500., 1000.])
    np.testing.assert_array_equal(samples, [7., 2., 6.])
    assert length == 1000.
