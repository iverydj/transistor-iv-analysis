"""Transconductance g_m = d|I_d|/dV_g in the on-direction (S), on the interpolated grid.

* ``window=None``: ``np.gradient`` of the interpolated current (original
  scripts' numerics). On noisy data this derivative follows every spike of
  the measurement, and its maximum is biased upwards.
* ``window=w`` (V): slope of a local quadratic least-squares fit over w volts
  around each grid point (equal to a Savitzky-Golay derivative of order 2 on
  a uniform grid; w = 1.6 V on a 0.2 V grid is the 9-point filter).

Negative values (current decreasing towards the on-state) are set to zero,
as in the original scripts.
"""

from __future__ import annotations

import numpy as np

from .._numerics import local_quadratic
from .interpolation import InterpolatedCurve


def transconductance(curve: InterpolatedCurve, window: float | None = 1.6) -> np.ndarray:
    """g_m (S) at every grid point of ``curve``; see the module docstring."""
    if window is None:
        gm = np.gradient(curve.i, curve.u)
    else:
        gm = local_quadratic(curve.u, curve.i, window)[1]
    return np.maximum(gm, 0.0)


def gm_at(curve: InterpolatedCurve, gm: np.ndarray, v_g: float) -> float:
    """g_m (S) at gate voltage ``v_g``, linearly interpolated on the grid."""
    u = curve.s * v_g
    if not curve.u[0] <= u <= curve.u[-1]:
        raise ValueError(f"V_g = {v_g} V is outside the sweep")
    return float(np.interp(u, curve.u, gm))
