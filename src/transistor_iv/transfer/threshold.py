"""Threshold voltage V_th by seven methods.

The six methods of the original scripts (result names in brackets) and one
textbook addition:

==========================================  ===================  =========================================
method                                      result field         V_th is ...
==========================================  ===================  =========================================
constant current                            ``vth_current``      V_g where |I_d| = I_crit
linear fit next to the subthreshold region  ``vth_interpol``     x-intercept of the most linear window
                                                                 of I_d(V_g) starting in the SS region
best linear window                          ``vth_global``       x-intercept of the window maximising
                                                                 |r| * |slope| anywhere in the sweep
maximum log slope                           ``vth_logderivative`` V_g of the largest d log10|I_d|/dV_g
square-root extrapolation (saturation)      ``vth_sqrt``         x-intercept of the most linear window
                                                                 of sqrt|I_d|(V_g)
second derivative                           ``vth_d2``           V_g of the largest d^2|I_d|/dV_g^2
tangent at maximum g_m (addition)           ``vth_gm_tangent``   V_gm - I_d(V_gm)/g_m,max
==========================================  ===================  =========================================

All window-based methods work on the interpolated curve
(:mod:`transistor_iv.transfer.interpolation`); a window of ``window`` volts
is floor(window / mean grid spacing) consecutive grid points (at least 3).
Linear-regime extrapolations report the raw intercept; see
:func:`elr_correction` for the V_d/2 correction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .._numerics import linear_fit, local_quadratic, sliding_linear_fits
from ..errors import ExtractionError
from ..sweep import on_direction
from .interpolation import InterpolatedCurve
from .swing import SubthresholdRegion


@dataclass(frozen=True)
class LineFit:
    """Straight line y = slope * V_g + intercept and its V_g-axis intercept.

    ``y`` is |I_d| (A) or sqrt|I_d| (A^0.5). ``slope`` is in V_g coordinates
    (negative for a p-type curve); ``r`` is the correlation coefficient in the
    on-direction frame (positive for a current rising towards the on-state);
    ``v_range`` is the (min, max) V_g of the fitted grid points.
    """

    vth: float
    slope: float
    intercept: float
    v_range: tuple[float, float]
    r: float
    n_points: int

    def __call__(self, v_g):
        return self.slope * np.asarray(v_g, dtype=float) + self.intercept


def _fit(curve: InterpolatedCurve, y: np.ndarray, start: int, n: int, what: str) -> LineFit:
    """Line through grid points start .. start+n-1 (truncated at the end of the grid)."""
    sl = slice(start, min(start + n, curve.u.size))
    u = curve.u[sl]
    if u.size < 2:
        raise ExtractionError(f"{what}: fewer than 2 grid points in the fit window")
    m, b, r = linear_fit(u, y[sl])
    if m == 0:
        raise ExtractionError(f"{what}: the fitted line is flat")
    v = curve.s * u
    return LineFit(vth=float(curve.s * (-b / m)), slope=float(curve.s * m), intercept=float(b),
                   v_range=(float(v.min()), float(v.max())), r=float(r), n_points=int(u.size))


def _best_window(curve: InterpolatedCurve, y: np.ndarray, n: int, starts: range, score: str) -> int:
    """Start index of the window with the highest score (first of equal ones).

    ``score='r'``: |r|; ``score='r_slope'``: |r| * |slope|.
    """
    if n > curve.u.size or len(starts) == 0:
        return starts.start if len(starts) else 0
    slope, _, r = sliding_linear_fits(curve.u, y, n)
    idx = np.arange(starts.start, min(starts.stop, slope.size))
    if idx.size == 0:
        return starts.start
    val = np.abs(r[idx]) if score == "r" else np.abs(r[idx]) * np.abs(slope[idx])
    return int(idx[np.argmax(val)])


def vth_constant_current(v_g, i_d, polarity: str, i_crit: float = 1e-8,
                         interpolation: str = "linear") -> float:
    """V_g at which |I_d| reaches ``i_crit`` (A), from the measured points of one segment.

    The crossing used is the last one from below: the last point with
    |I_d| < i_crit and the next point (in the on-direction). Between the two,
    V_g is interpolated linearly in |I_d| (``interpolation="linear"``, as the
    original scripts did) or in log10|I_d| (``"log"``).
    """
    u, a, s = on_direction(v_g, i_d, polarity)
    below = np.flatnonzero(a < i_crit)
    if below.size == 0:
        raise ExtractionError(f"constant current: |I_d| never drops below I_crit = {i_crit:g} A")
    k = int(below[-1])
    if k == u.size - 1:
        raise ExtractionError(f"constant current: |I_d| does not reach I_crit = {i_crit:g} A at the on end")
    if interpolation == "linear":
        f = (i_crit - a[k]) / (a[k + 1] - a[k])
    elif interpolation == "log":
        if a[k] <= 0:
            raise ExtractionError("constant current: zero current next to the crossing")
        f = (np.log10(i_crit) - np.log10(a[k])) / (np.log10(a[k + 1]) - np.log10(a[k]))
    else:
        raise ValueError(f"interpolation must be 'linear' or 'log', got {interpolation!r}")
    return float(s * (u[k] + f * (u[k + 1] - u[k])))


def vth_linear_ss_window(curve: InterpolatedCurve, region: SubthresholdRegion, window: float = 5.0,
                         extension: float = 3.0) -> LineFit:
    """``vth_interpol``: linear extrapolation of the most linear I_d window near turn-on.

    Windows of ``window`` volts of the interpolated I_d(V_g) whose first
    point lies between the start of the subthreshold region and ``extension``
    volts past its end (towards the on-state) are fitted with straight lines;
    the one with the largest |r| gives V_th = -intercept / slope. If no full
    window fits before the end of the sweep, the remaining points from the
    region start are used.
    """
    n = curve.n_points(window, minimum=3)
    last = min(region.end + curve.n_points(extension), curve.u.size - 1)
    best = _best_window(curve, curve.i, n, range(region.start, last + 1), "r")
    return _fit(curve, curve.i, best, n, "vth_interpol")


def vth_best_linear_window(curve: InterpolatedCurve, window: float = 5.0) -> LineFit:
    """``vth_global``: linear extrapolation of the window maximising |r| * |slope|.

    All windows of ``window`` volts of the interpolated I_d(V_g) are fitted;
    the score |r| * |slope| prefers straight and steep parts, i.e. the
    linear region near maximum g_m. Independent of the subthreshold region.
    """
    n = curve.n_points(window, minimum=3)
    best = _best_window(curve, curve.i, n, range(curve.u.size - n + 1), "r_slope")
    return _fit(curve, curve.i, best, n, "vth_global")


def vth_sqrt_extrapolation(curve: InterpolatedCurve, window: float = 5.0) -> LineFit:
    """``vth_sqrt``: linear extrapolation of the most linear window of sqrt|I_d|(V_g).

    For a saturated long-channel device, I_d = (W C_i mu / 2L)(V_g - V_th)^2,
    so sqrt|I_d| is linear in V_g with x-intercept V_th and slope
    sqrt(W C_i mu / 2L) (see :func:`transistor_iv.transfer.mobility.saturation_mobility_from_slope`).
    The window with the largest |r| among all windows of ``window`` volts is used.
    """
    n = curve.n_points(window, minimum=3)
    y = np.sqrt(np.maximum(curve.i, 0.0))
    best = _best_window(curve, y, n, range(curve.u.size - n + 1), "r")
    return _fit(curve, y, best, n, "vth_sqrt")


def vth_max_log_slope(curve: InterpolatedCurve, region: SubthresholdRegion) -> float:
    """``vth_logderivative``: V_g of the largest d log10|I_d|/dV_g inside the subthreshold region.

    This is where the point-wise swing is smallest, i.e. the steepest point of
    the subthreshold region; it lies below the strong-inversion threshold.
    """
    sl = slice(region.start, region.end + 1)
    k = int(np.argmax(curve.dlog[sl]))
    return float(curve.s * curve.u[sl][k])


def second_derivative(curve: InterpolatedCurve, window: float | None = None) -> np.ndarray:
    """d^2|I_d|/dV_g^2 on the grid (A/V^2; equal in the V_g and on-direction frames).

    ``window=None``: two successive ``np.gradient`` of the interpolated I_d
    (original scripts); otherwise from local quadratic least-squares fits over
    ``window`` volts (Savitzky-Golay on a uniform grid).
    """
    if window is None:
        return np.gradient(np.gradient(curve.i, curve.u), curve.u)
    return local_quadratic(curve.u, curve.i, window)[2]


def vth_second_derivative(curve: InterpolatedCurve, window: float | None = None) -> float:
    """``vth_d2``: V_g of the maximum of d^2|I_d|/dV_g^2 over the whole sweep.

    The second derivative (= d g_m / dV_g, the "transconductance change")
    peaks where the channel turns on. ``window`` as in :func:`second_derivative`.
    """
    if curve.u.size < 5:
        raise ExtractionError("second derivative: fewer than 5 grid points")
    d2 = second_derivative(curve, window)
    return float(curve.s * curve.u[int(np.argmax(d2))])


def vth_max_gm_tangent(curve: InterpolatedCurve, gm: np.ndarray) -> LineFit:
    """Tangent of I_d(V_g) at the maximum of ``gm`` (textbook linear extrapolation).

    V_th = V_gm - I_d(V_gm) / g_m,max, with V_gm the gate voltage of the
    g_m maximum and I_d the interpolated current there. ``gm`` is the
    on-direction transconductance on the grid
    (:func:`transistor_iv.transfer.transconductance.transconductance`).
    """
    k = int(np.argmax(gm))
    g = float(gm[k])
    if g <= 0:
        raise ExtractionError("gm tangent: g_m is not positive anywhere")
    u_th = curve.u[k] - curve.i[k] / g
    v = float(curve.s * curve.u[k])
    return LineFit(vth=float(curve.s * u_th), slope=curve.s * g, intercept=float(curve.i[k] - g * curve.u[k]),
                   v_range=(v, v), r=1.0, n_points=1)


def elr_correction(vth_intercept: float, v_d: float) -> float:
    """Linear-regime threshold from an extrapolated intercept: V_th = V_gi - V_d/2.

    In the linear regime of the square law, I_d = k (V_g - V_th - V_d/2) V_d,
    so the straight line through I_d(V_g) crosses zero at V_gi = V_th + V_d/2
    (signed V_d, grounded source; valid for both polarities).
    """
    return float(vth_intercept - v_d / 2.0)
