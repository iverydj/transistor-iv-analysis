"""Shared numerical helpers (private)."""

from __future__ import annotations

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .errors import ExtractionError


def sliding_linear_fits(x: np.ndarray, y: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Least-squares lines through every run of ``n`` consecutive points.

    Returns ``(slope, intercept, r)`` for the windows starting at index
    0, 1, ..., len(x) - n. ``r`` is the Pearson correlation coefficient; it
    is 0 for a window in which x or y is constant (as in scipy's linregress).
    """
    if n < 2 or n > x.size:
        raise ValueError(f"window of {n} points does not fit {x.size} points")
    xw = sliding_window_view(x, n)
    yw = sliding_window_view(y, n)
    dx = xw - xw.mean(axis=1, keepdims=True)
    dy = yw - yw.mean(axis=1, keepdims=True)
    sxx = np.sum(dx * dx, axis=1)
    syy = np.sum(dy * dy, axis=1)
    sxy = np.sum(dx * dy, axis=1)
    slope = sxy / sxx
    intercept = yw.mean(axis=1) - slope * xw.mean(axis=1)
    den = np.sqrt(sxx * syy)
    r = np.zeros_like(sxy)
    np.divide(sxy, den, out=r, where=den > 0)
    return slope, intercept, np.clip(r, -1.0, 1.0)


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Least-squares line y = slope * x + intercept; returns ``(slope, intercept, r)``."""
    s, b, r = sliding_linear_fits(np.asarray(x, float), np.asarray(y, float), len(x))
    return float(s[0]), float(b[0]), float(r[0])


def window_bounds(x: np.ndarray, width: float) -> tuple[np.ndarray, np.ndarray]:
    """Index range [lo, hi) of the samples inside a ``width``-wide window around each sample.

    ``x`` must be increasing. The window [x_i - width/2, x_i + width/2] is
    shifted inwards where it would leave the data range.
    """
    if width <= 0:
        raise ValueError(f"window width must be positive, got {width}")
    if x[-1] - x[0] < width:
        raise ExtractionError(f"the sweep spans {x[-1] - x[0]:g} V, less than the {width:g} V window")
    left = np.clip(x - width / 2, x[0], x[-1] - width)
    tol = 1e-4 * width  # instruments often store V_g in single precision (17.8999996)
    lo = np.searchsorted(x, left - tol, side="left")
    hi = np.searchsorted(x, left + width + tol, side="right")
    return lo, hi


def local_quadratic(x: np.ndarray, y: np.ndarray, width: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Smoothed value, first and second derivative of y(x) at every sample.

    At each sample a quadratic is least-squares fitted to the samples inside a
    ``width``-wide window (see :func:`window_bounds`) and evaluated there. On
    a uniform grid this is the Savitzky-Golay filter of polynomial order 2
    (window length = number of samples in ``width``, edges handled like
    scipy's ``mode='interp'``); non-uniform grids are handled directly.
    """
    lo, hi = window_bounds(x, width)
    count = hi - lo
    if count.min() < 5:
        raise ExtractionError(f"a {width:g} V window holds fewer than 5 samples somewhere; widen the window")
    k = np.arange(count.max())
    mask = k[None, :] < count[:, None]
    idx = np.minimum(lo[:, None] + k[None, :], x.size - 1)
    dx = np.where(mask, x[idx] - x[:, None], 0.0)
    yy = np.where(mask, y[idx], 0.0)
    p = [np.sum(mask * dx ** m, axis=1) for m in range(5)]
    gram = np.stack([np.stack(p[i:i + 3], axis=-1) for i in range(3)], axis=-2)  # (n, 3, 3)
    rhs = np.stack([np.sum(yy * dx ** m, axis=1) for m in range(3)], axis=-1)    # (n, 3)
    coef = np.linalg.solve(gram, rhs[..., None])[..., 0]
    return coef[:, 0], coef[:, 1], 2.0 * coef[:, 2]
