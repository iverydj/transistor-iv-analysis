"""Subthreshold region and subthreshold swing.

Definition: SS = dV_g / d log10|I_d|, in mV/dec (millivolts of gate voltage
per decade of drain current), i.e. SS = 1000 / (d log10|I_d| / dV_g) with the
derivative taken in the on-direction (positive for both polarities).

Three selection rules are available:

* :func:`subthreshold_swing` (default, the original scripts' rule): the
  point-wise SS of the interpolated curve inside the subthreshold region,
  minimum among values >= ``floor`` (60 mV/dec).
* :func:`ss_decade_span`: the shortest V_g span over which the (envelope of
  the) measured current rises by a given number of decades.
* :func:`ss_current_window`: a straight-line fit of log10|I_d| vs V_g over an
  explicit current window.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import median_filter

from .._numerics import linear_fit
from ..errors import ExtractionError
from ..sweep import on_direction
from .interpolation import InterpolatedCurve
from .on_off import robust_off_level


@dataclass(frozen=True)
class SubthresholdRegion:
    """Index range [start, end] (inclusive) of the interpolated grid and its V_g range."""

    start: int
    end: int
    v_range: tuple[float, float]


@dataclass(frozen=True)
class SwingResult:
    """A subthreshold swing value and where it was taken.

    ``ss`` in mV/dec; ``v_range`` the (min, max) V_g it refers to (a single
    point for the point-wise rule); ``i_range`` the (low, high) current there (A).
    """

    ss: float
    v_range: tuple[float, float]
    i_range: tuple[float, float]
    method: str


def _v_range(curve_s: int, u0: float, u1: float) -> tuple[float, float]:
    a, b = sorted((curve_s * u0, curve_s * u1))
    return float(a), float(b)


def subthreshold_region(curve: InterpolatedCurve, sensitivity: float = 1.4) -> SubthresholdRegion:
    """Longest run of grid points with d log10|I_d|/dV_g >= 10**(-sensitivity) dec/V.

    With the default ``sensitivity = 1.4`` the threshold is 0.040 dec/V (a
    point-wise swing of at most 25 V/dec); a larger value widens the region.
    Of equally long runs the first (lowest on-direction voltage) is taken.
    """
    idx = np.flatnonzero(np.log10(curve.dlog + 1e-12) >= -sensitivity)
    if idx.size == 0:
        raise ExtractionError(f"no subthreshold region: d log10|I_d|/dV_g never reaches "
                              f"10^-{sensitivity:g} dec/V (raise ss_sensitivity)")
    breaks = np.flatnonzero(np.diff(idx) != 1)
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    k = int(np.argmax(ends - starts))
    start, end = int(starts[k]), int(ends[k])
    return SubthresholdRegion(start, end, _v_range(curve.s, curve.u[start], curve.u[end]))


def swing_curve(curve: InterpolatedCurve) -> np.ndarray:
    """Point-wise SS = 1000 / (d log10|I_d|/dV_g) on the interpolated grid, mV/dec."""
    return 1000.0 / curve.dlog


def subthreshold_swing(curve: InterpolatedCurve, region: SubthresholdRegion,
                       floor: float | None = 60.0) -> SwingResult:
    """Minimum point-wise SS inside ``region`` among the values >= ``floor`` (mV/dec).

    ``floor=None`` takes the plain minimum. The default floor of 60 mV/dec
    (close to the room-temperature thermionic limit ln(10) k_B T/q = 59.5 mV/dec)
    discards steeper values as artefacts of noise or of the interpolation.
    """
    sl = slice(region.start, region.end + 1)
    ss = swing_curve(curve)[sl]
    ok = np.ones(ss.size, bool) if floor is None else ss >= floor
    if not ok.any():
        raise ExtractionError(f"no point-wise SS >= {floor} mV/dec in the subthreshold region")
    k = int(np.flatnonzero(ok)[np.argmin(ss[ok])])
    u, i = curve.u[sl][k], curve.i[sl][k]
    return SwingResult(float(ss[k]), _v_range(curve.s, u, u), (float(i), float(i)), "pointwise")


def ss_current_window(v_g, i_d, polarity: str, i_low: float, i_high: float) -> SwingResult:
    """SS from a straight-line fit of log10|I_d| vs V_g over i_low < |I_d| < i_high (A).

    Uses the measured points of one monotonic segment; needs at least 4 of
    them inside the window. The window is stated, not detected, so the
    result is reproducible; report it together with SS.
    """
    u, a, s = on_direction(v_g, i_d, polarity)
    m = (a > i_low) & (a < i_high)
    if m.sum() < 4:
        raise ExtractionError(f"fewer than 4 points with {i_low:g} A < |I_d| < {i_high:g} A")
    slope, _, _ = linear_fit(u[m], np.log10(a[m]))
    if slope <= 0:
        raise ExtractionError("log10|I_d| does not increase in the on-direction inside the current window")
    return SwingResult(1000.0 / slope, _v_range(s, u[m][0], u[m][-1]), (i_low, i_high), "current_window")


def ss_decade_span(v_g, i_d, polarity: str, decades: float = 1.0, i_min: float | None = None,
                   floor: float | None = None) -> SwingResult:
    """Shortest V_g span over which the measured current rises by ``decades`` decades.

    Rule: (1) log10|I_d| of one monotonic segment is median-filtered over 5
    points; (2) only the branch from the current minimum towards the on-state
    is kept and replaced by its running maximum (the envelope), so noise below
    the subthreshold region cannot mimic a steep slope; (3) the envelope is
    interpolated linearly, and among all spans starting at or above ``i_min``
    (default: ten times :func:`robust_off_level`) over which it rises by
    ``decades``, the shortest is taken; spans steeper than ``floor`` (mV/dec)
    are skipped if a floor is given. SS = span / decades.
    """
    u, a, s = on_direction(v_g, i_d, polarity)
    if np.any(a <= 0):
        raise ExtractionError("|I_d| contains exact zeros; log-scale quantities are undefined")
    y = median_filter(np.log10(a), size=5, mode="mirror")
    if i_min is None:
        i_min = 10.0 * robust_off_level(a)
    k0 = int(np.argmin(y))
    u, env = u[k0:], np.maximum.accumulate(y[k0:])
    lo = np.log10(i_min)

    def crossing(level: np.ndarray) -> np.ndarray:
        """u where the envelope first reaches each level (NaN outside the data)."""
        j = np.searchsorted(env, level, side="left")
        out = np.full(level.shape, np.nan)
        ok = (j > 0) & (j < env.size)
        jj = j[ok]
        f = (level[ok] - env[jj - 1]) / (env[jj] - env[jj - 1])
        out[ok] = u[jj - 1] + f * (u[jj] - u[jj - 1])
        return out

    # candidate spans: starting on a sample, ending on a sample, or starting exactly at i_min;
    # the span length is linear in the start position between these, so they cover all positions
    edge = np.array([lo])
    low = np.concatenate([env, env - decades, edge])
    u0 = np.concatenate([u, crossing(env - decades), crossing(edge)])
    u1 = np.concatenate([crossing(env + decades), u, crossing(edge + decades)])
    ok = np.isfinite(u0) & np.isfinite(u1) & (low >= lo)
    ss = 1000.0 * (u1 - u0) / decades
    if floor is not None:
        ok &= ss >= floor
    if not ok.any():
        raise ExtractionError(f"|I_d| never rises by {decades:g} decade(s) above {i_min:.3g} A"
                              + ("" if floor is None else f" with SS >= {floor} mV/dec"))
    k = int(np.flatnonzero(ok)[np.argmin(ss[ok])])
    return SwingResult(float(ss[k]), _v_range(s, u0[k], u1[k]),
                       (float(10 ** low[k]), float(10 ** (low[k] + decades))), "decade_span")
