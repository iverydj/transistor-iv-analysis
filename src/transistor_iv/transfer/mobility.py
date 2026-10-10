"""Field-effect mobilities in cm^2/(V s).

Units: W and L in um (only their ratio enters), C_i in F/cm^2, g_m in S,
I_d in A, voltages in V. |V_d| is used, so the formulas hold for both
polarities. All formulas assume the long-channel gradual-channel model with a
gate-voltage-independent mobility and no contact resistance.

=====================================  ===========================================  ==================
quantity                               formula                                      regime
=====================================  ===========================================  ==================
linear (field-effect) mobility mu_lin  L g_m / (W C_i |V_d|)                        |V_d| << V_g - V_th
effective mobility mu_eff              L |I_d| / (W C_i |V_d| |V_g - V_th|)         |V_d| << V_g - V_th
saturation mobility mu_sat (point)     2 L |I_d| / (W C_i (V_g - V_th)^2)           |V_d| >= V_g - V_th
saturation mobility from sqrt slope    2 L b^2 / (W C_i), b = d sqrt|I_d| / dV_g    |V_d| >= V_g - V_th
mu_FE at V_g from a sqrt fit           2 L b^2 |V_g - V_th| / (W C_i |V_d|)         see function
=====================================  ===========================================  ==================
"""

from __future__ import annotations

import numpy as np

from ..device import Device


def _abs_vd(v_d: float) -> float:
    if not np.isfinite(v_d) or v_d == 0:
        raise ValueError(f"the mobility needs a finite, non-zero drain bias, got V_d = {v_d}")
    return abs(float(v_d))


def linear_mobility(gm, v_d: float, device: Device):
    """mu_lin = L g_m / (W C_i |V_d|)."""
    return np.asarray(gm) / (device.w_over_l * device.c_i * _abs_vd(v_d))


def effective_mobility(v_g, i_d, vth: float, v_d: float, device: Device, min_overdrive: float = 1.0):
    """mu_eff = L |I_d| / (W C_i |V_d| |V_g - V_th|); NaN where |V_g - V_th| <= ``min_overdrive`` (V)."""
    ov = np.abs(np.asarray(v_g, dtype=float) - vth)
    mu = np.abs(np.asarray(i_d)) / (device.w_over_l * device.c_i * _abs_vd(v_d) * np.where(ov > 0, ov, np.nan))
    return np.where(ov > min_overdrive, mu, np.nan)


def saturation_mobility(v_g, i_d, vth: float, device: Device, min_overdrive: float = 1.0):
    """mu_sat = 2 L |I_d| / (W C_i (V_g - V_th)^2); NaN where |V_g - V_th| <= ``min_overdrive`` (V)."""
    ov = np.abs(np.asarray(v_g, dtype=float) - vth)
    mu = 2.0 * np.abs(np.asarray(i_d)) / (device.w_over_l * device.c_i * np.where(ov > 0, ov, np.nan) ** 2)
    return np.where(ov > min_overdrive, mu, np.nan)


def saturation_mobility_from_slope(sqrt_slope: float, device: Device) -> float:
    """mu_sat = 2 L b^2 / (W C_i), with b the slope of sqrt|I_d| vs V_g in A^0.5/V."""
    return float(2.0 * sqrt_slope ** 2 / (device.w_over_l * device.c_i))


def fe_mobility_from_sqrt_fit(sqrt_slope: float, vth: float, v_g: float, v_d: float, device: Device) -> float:
    """mu_FE at gate voltage ``v_g`` from a straight-line fit of sqrt|I_d| (slope b, intercept V_th).

    If sqrt|I_d| = b |V_g - V_th| over the fitted range, then
    g_m = 2 b^2 |V_g - V_th| and the linear-regime field-effect mobility is

        mu_FE = L g_m / (W C_i |V_d|) = 2 L b^2 |V_g - V_th| / (W C_i |V_d|).

    Useful when I_d is quadratic in V_g although the device is in the linear
    regime (|V_d| small), as in 2D-material FETs whose mobility grows with the
    carrier density; then g_m keeps rising to the end of the sweep and its
    maximum is not defined, while the fitted g_m at a stated V_g is.
    """
    return float(2.0 * sqrt_slope ** 2 * abs(v_g - vth) / (device.w_over_l * device.c_i * _abs_vd(v_d)))
