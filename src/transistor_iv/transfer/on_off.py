"""On-current, off-current and on/off ratio.

Default (original scripts): computed from the raw data of the whole curve
(every sweep segment), I_on = max |I_d| and I_off = the smallest positive
|I_d|, so that the interpolation cannot overshoot. Option ``off="robust"``:
I_off = :func:`robust_off_level`, which a single noise dip cannot lower.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import median_filter

from ..errors import ExtractionError


@dataclass(frozen=True)
class OnOff:
    i_on: float     # A
    i_off: float    # A
    ratio: float


def robust_off_level(i_d) -> float:
    """Off-state current level (A): 10**median of the lowest 10 % of the 5-point
    median-filtered log10|I_d| values."""
    a = np.abs(np.asarray(i_d, dtype=float))
    if np.any(a <= 0):
        raise ExtractionError("|I_d| contains exact zeros; the log-scale off level is undefined")
    y = median_filter(np.log10(a), size=5, mode="mirror")
    return float(10 ** np.median(np.sort(y)[:max(y.size // 10, 1)]))


def on_off_ratio(i_d, off: str = "min") -> OnOff:
    """I_on = max |I_d|; I_off = min positive |I_d| (``off="min"``) or
    :func:`robust_off_level` (``off="robust"``); ratio = I_on / I_off."""
    a = np.abs(np.asarray(i_d, dtype=float))
    if off == "min":
        pos = a[a > 0]
        if pos.size == 0:
            raise ExtractionError("no positive |I_d|")
        i_off = float(pos.min())
    elif off == "robust":
        i_off = robust_off_level(a)
    else:
        raise ValueError(f"off must be 'min' or 'robust', got {off!r}")
    i_on = float(a.max())
    return OnOff(i_on, i_off, i_on / i_off)
