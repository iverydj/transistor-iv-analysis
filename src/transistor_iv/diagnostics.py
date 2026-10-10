"""Optional checks of a transfer analysis (``TransferSettings(diagnostics=True)``).

Regime flags
    ``flag_linear_regime``: |V_d| exceeds the overdrive V_g - V_th at the
    low end of the vth_global fit window (V_th = intercept - V_d/2), so the linear-regime formulas
    (vth_interpol, vth_global, mu_fe_max) are applied outside the linear regime.
    ``flag_saturation``: |V_d| is below the overdrive at the high end of the
    vth_sqrt fit window, so the device is not saturated where the
    square-root fit was made (vth_sqrt, mu_sat_sqrt).
    ``flag_window_at_end``: the vth_global or vth_sqrt window reaches the
    on-end of the sweep (the curve may still be bending there).
Window sensitivities
    ``*_sens``: the largest change of a result when its window is halved
    and doubled (V for threshold voltages, relative for mobilities). A large
    value means the result depends on the window choice (noise, curvature).
    It does not reveal systematic errors such as contact resistance.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .device import Device
from .io.curves import TransferCurve


def _flags(curve: TransferCurve, res) -> dict[str, str]:
    s = 1 if res.polarity == "n" else -1
    vd = abs(curve.v_d) if np.isfinite(curve.v_d) else np.nan
    out = {"flag_linear_regime": "", "flag_saturation": "", "flag_window_at_end": ""}
    on_end = max(s * res.segment_v[0], s * res.segment_v[1])
    ends = []
    g = res.fits.get("vth_global")
    if g is not None and np.isfinite(vd):
        ov_low = min(s * g.v_range[0], s * g.v_range[1]) - s * (g.vth - curve.v_d / 2)
        if vd > ov_low:
            out["flag_linear_regime"] = f"|V_d| = {vd:g} V > V_g - V_th = {ov_low:.3g} V at the fit window start"
    q = res.fits.get("vth_sqrt")
    if q is not None and np.isfinite(vd):
        ov_high = max(s * q.v_range[0], s * q.v_range[1]) - s * q.vth
        if vd < ov_high:
            out["flag_saturation"] = f"|V_d| = {vd:g} V < V_g - V_th = {ov_high:.3g} V at the fit window end"
    for key in ("vth_global", "vth_sqrt"):
        fit = res.fits.get(key)
        if fit is not None and max(s * fit.v_range[0], s * fit.v_range[1]) >= on_end - 1e-9:
            ends.append(key)
    out["flag_window_at_end"] = ", ".join(ends)
    return out


def transfer_diagnostics(curve: TransferCurve, device: Device, settings, res) -> dict[str, object]:
    """Regime flags and window sensitivities of ``res`` (see the module docstring)."""
    from .transfer.analysis import analyze_transfer

    out: dict[str, object] = _flags(curve, res)
    base = replace(settings, diagnostics=False)
    alt = []
    for f in (0.5, 2.0):
        st = replace(base, linear_window=settings.linear_window * f,
                     gm_window=None if settings.gm_window is None else settings.gm_window * f)
        alt.append(analyze_transfer(curve, device, st))
    for key in ("vth_interpol", "vth_global", "vth_sqrt"):
        d = [abs(getattr(a, key) - getattr(res, key)) for a in alt]
        out[f"{key}_sens"] = float(np.nanmax(d)) if d and np.isfinite(d).any() else np.nan
    for key in ("mu_sat_sqrt", "mu_fe_max"):
        d = [abs(getattr(a, key) / getattr(res, key) - 1.0) for a in alt]
        out[f"{key}_sens"] = float(np.nanmax(d)) if d and np.isfinite(d).any() else np.nan
    if settings.gm_window is None:
        out["mu_fe_max_sens"] = np.nan          # no window to vary
    return out
