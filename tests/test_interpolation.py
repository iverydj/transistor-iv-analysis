import numpy as np
import pytest

from transistor_iv import CurveError
from transistor_iv.transfer import fill_gaps, interpolate_curve


def test_fill_gaps_as_in_the_original_scripts():
    np.testing.assert_allclose(fill_gaps(np.array([0.0, 0.25, 0.5]), 0.2), [0, 0.2, 0.25, 0.45, 0.5])
    np.testing.assert_allclose(fill_gaps(np.array([0.0, 1.0]), 0.2), [0, 0.2, 0.4, 0.6, 0.8, 1.0])
    np.testing.assert_allclose(fill_gaps(np.array([0.0, 0.1, 0.2]), 0.2), [0, 0.1, 0.2])


def test_fill_gaps_ignores_single_precision_excess():
    # V_g stored in single precision: steps of 0.2000000x must not create a point 1e-9 V before the next one
    v = np.arange(-30, 40.001, 0.2).astype(np.float32).astype(float)
    g = fill_gaps(v, 0.2)
    np.testing.assert_array_equal(g, v)
    assert np.diff(g).min() > 0.19


def test_spline_passes_through_the_data(device):
    v = np.linspace(-5, 5, 41)
    i = device.sign * 1e-12 * np.exp(device.sign * v)
    c = interpolate_curve(v, i, device.polarity)
    on_grid = np.isin(np.round(c.u, 9), np.round(c.u_data, 9))
    np.testing.assert_allclose(c.log10_i[on_grid], np.log10(np.abs(i))[::1 if device.sign == 1 else -1], rtol=1e-12)
    np.testing.assert_allclose(c.dlog, 1 / np.log(10), rtol=1e-6)    # d log10 e^u / du


def test_too_few_points():
    with pytest.raises(CurveError, match="at least 4"):
        interpolate_curve([0, 1, 2], [1, 2, 3], "n")
