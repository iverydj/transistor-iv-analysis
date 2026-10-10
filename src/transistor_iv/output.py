"""Output-curve (I_d vs V_d at fixed V_g) analysis.

* :func:`symmetry_score`: how similar log|I_d| is at +V_d and -V_d (original
  scripts' contact-symmetry indicator for sweeps through V_d = 0).
* :func:`on_resistance`: low-bias channel resistance (addition).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._numerics import linear_fit
from .io.curves import OutputFamily


def symmetry_score(v_d, i_d) -> float:
    """Symmetry of ln|I_d| between +V_d and -V_d, in percent (100 = perfectly symmetric).

        score = 100 [1 - RMS(y(+V_d) - y(-V_d)) / std(y)],   y = ln|I_d|

    std(y) is taken over all points of the curve. The pairs are formed at the
    positive-bias points whose mirror voltage -V_d lies inside the measured
    negative range, with y(-V_d) interpolated linearly between negative-bias
    points (so the V_d grid need not be symmetric). The score is independent
    of the logarithm base. Zero or negative scores mean the asymmetry is as
    large as the spread of the curve itself.
    """
    v = np.asarray(v_d, dtype=float)
    a = np.abs(np.asarray(i_d, dtype=float))
    if np.any(a == 0):
        raise ValueError("|I_d| has exact zeros; ln|I_d| is undefined")
    y = np.log(a)
    neg, pos = v < 0, v > 0
    if not neg.any() or not pos.any():
        raise ValueError("symmetry_score needs drain voltages of both signs")
    vn, yn = -v[neg], y[neg]
    order = np.argsort(vn)
    vn, yn = vn[order], yn[order]
    vp, yp = v[pos], y[pos]
    inside = (vp >= vn[0]) & (vp <= vn[-1])
    if inside.sum() < 2:
        raise ValueError("fewer than 2 positive-bias points have a mirrored negative-bias counterpart")
    diff = yp[inside] - np.interp(vp[inside], vn, yn)
    return float(100.0 * (1.0 - np.sqrt(np.mean(diff ** 2)) / np.std(y)))


def family_symmetry(family: OutputFamily) -> float:
    """Mean :func:`symmetry_score` over the curves of a family (original scripts' figure title value)."""
    return float(np.mean([symmetry_score(v, i) for v, i in zip(family.v_d, family.i_d)]))


def on_resistance(v_d, i_d, v_d_max: float) -> float:
    """Low-bias resistance R_on = (dI_d/dV_d)^-1 in ohm, from a line fit over |V_d| <= ``v_d_max``.

    Multiply by the channel width W to get the width-normalised value (ohm cm
    with W in cm). Includes the contact resistance.
    """
    v = np.asarray(v_d, dtype=float)
    i = np.asarray(i_d, dtype=float)
    sel = np.abs(v) <= v_d_max
    if sel.sum() < 3:
        raise ValueError(f"fewer than 3 points with |V_d| <= {v_d_max} V")
    slope, _, _ = linear_fit(v[sel], i[sel])
    if slope <= 0:
        raise ValueError(f"I_d does not increase with V_d for |V_d| <= {v_d_max} V (slope {slope:g} S)")
    return 1.0 / slope


def summarize_output(family: OutputFamily, v_d_max: float | None = None) -> pd.DataFrame:
    """One row per curve: gate voltage, max |I_d|, symmetry score (NaN if V_d has one sign)
    and, when ``v_d_max`` is given, R_on."""
    rows = []
    for vg, v, i in zip(family.v_g, family.v_d, family.i_d):
        both = bool((v < 0).any() and (v > 0).any())
        row = {"name": family.name, "v_g": float(vg), "i_d_max": float(np.max(np.abs(i))),
               "symmetry": symmetry_score(v, i) if both else np.nan}
        if v_d_max is not None:
            row["r_on"] = on_resistance(v, i, v_d_max)
        rows.append(row)
    return pd.DataFrame(rows)
