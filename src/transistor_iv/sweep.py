"""Sweep segments and the polarity-independent "on-direction" frame.

A measured V_g sweep can consist of several monotonic segments (a dual sweep
-20 -> 20 -> -20 V has two; a ramp 0 -> -30 V followed by -30 -> 40 -> -30 V
has three). Every transfer-curve method analyses one monotonic segment.

On-direction frame
------------------
With s = +1 for n-type and s = -1 for p-type, the methods work on

    u = s * V_g   (ascending),   a = |I_d|

so that the current increases with u for both polarities. A p-type curve is
therefore analysed exactly like the mirror image of an n-type curve, and
voltages are mapped back with V = s * u.
"""

from __future__ import annotations

import numpy as np

from .errors import CurveError


def monotonic_segments(v_g) -> list[np.ndarray]:
    """Index arrays of the monotonic segments of a V_g sweep, in measurement order.

    A turning point belongs to both adjacent segments; a repeated V_g value
    (zero step, e.g. a doubled end point) separates segments and belongs to
    neither side's interior.

    Example: V_g = 0, -2, -4, -2, 0, 2, 0 gives [0, 1, 2], [2, 3, 4, 5], [5, 6].
    """
    v = np.asarray(v_g, dtype=float)
    if v.ndim != 1 or v.size < 2:
        raise CurveError("V_g must be a 1-D array with at least 2 points")
    trend = np.sign(np.diff(v))
    segments: list[np.ndarray] = []
    start = 0
    i = 0
    while i < trend.size:
        if trend[i] == 0:
            i += 1
            start = i
            continue
        j = i + 1
        while j < trend.size and trend[j] == trend[i]:
            j += 1
        segments.append(np.arange(start, j + 1))
        i = j
        start = j
    if not segments:
        raise CurveError("V_g is constant; this is not a sweep")
    return segments


def select_segment(v_g, index: int = -1) -> np.ndarray:
    """Indices of one monotonic segment: ``index = -1`` the longest (the first of
    equally long ones), ``index >= 0`` the index-th segment in measurement order."""
    segments = monotonic_segments(v_g)
    if index == -1:
        return max(segments, key=len)
    if not 0 <= index < len(segments):
        spans = [(float(np.asarray(v_g)[s[0]]), float(np.asarray(v_g)[s[-1]])) for s in segments]
        raise CurveError(f"segment {index} requested but the sweep has {len(segments)} segments {spans}")
    return segments[index]


def polarity_sign(polarity: str) -> int:
    """+1 for 'n', -1 for 'p'."""
    if polarity not in ("n", "p"):
        raise ValueError(f"polarity must be 'n' or 'p', got {polarity!r}")
    return 1 if polarity == "n" else -1


def on_direction(v_g, i_d, polarity: str) -> tuple[np.ndarray, np.ndarray, int]:
    """``(u, a, s)``: u = s*V_g in ascending order, a = |I_d| in the same order, s = +-1.

    The input must be one monotonic segment (see :func:`select_segment`).
    """
    v = np.asarray(v_g, dtype=float)
    a = np.abs(np.asarray(i_d, dtype=float))
    if v.ndim != 1 or v.shape != a.shape:
        raise CurveError(f"V_g and I_d must be 1-D arrays of equal length, got {v.shape} and {a.shape}")
    if not (np.isfinite(v).all() and np.isfinite(a).all()):
        raise CurveError("V_g and I_d must be finite")
    s = polarity_sign(polarity)
    u = s * v
    du = np.diff(u)
    if np.all(du <= 0):
        u, a = u[::-1], a[::-1]
    elif not np.all(du >= 0):
        raise CurveError("V_g is not monotonic; select one sweep segment first (transistor_iv.select_segment)")
    return u, a, s
