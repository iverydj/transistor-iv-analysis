"""Synthetic transfer and output curves from a smooth long-channel charge model.

Model (n-type; p-type mirrors every voltage and the current), source grounded::

    channel charge   Q(x) = C_i * q(x),   q(x) = S * ln(1 + exp(x / S)),   S = SS / ln(10)
    drain current    I_d = (W mu C_i / L) * [G(V_g - V_th) - G(V_g - V_th - V_d)]
    with             G(x) = integral of q = -S^2 * Li2(-exp(x / S))

i.e. the gradual-channel integral of a charge that is exponential below
threshold and equal to C_i (V_g - V_th) above it. This gives

* linear regime   I_d = (W mu C_i / L) (V_ov - V_d/2) V_d   (exact for V_ov - V_d >> S)
* saturation      I_d = (W mu C_i / 2L) (V_ov^2 + pi^2 S^2 / 3)
* subthreshold    I_d proportional to exp(V_ov / S), i.e. exactly the swing SS

with V_ov = V_g - V_th. Mobility is constant and there is no series
resistance, so g_m rises monotonically to (W mu C_i / L) V_d and the
extraction methods have known answers (see the tests).
"""

from __future__ import annotations

import numpy as np
from scipy.special import spence

from .device import Device
from .io.curves import OutputFamily, TransferCurve


def _g(x: np.ndarray, s_ln: float) -> np.ndarray:
    """-S^2 Li2(-exp(x/S)), using Li2(z) = spence(1 - z) and the inversion formula for x > 0."""
    y = np.asarray(x, dtype=float) / s_ln
    out = np.empty_like(y)
    tiny = y < -5                      # series -Li2(-t) = sum (-1)^(k+1) t^k / k^2, t < 7e-3
    t = np.exp(y[tiny])
    out[tiny] = sum((-1) ** (k + 1) * t ** k / k ** 2 for k in range(1, 10))
    neg = (y <= 0) & ~tiny
    out[neg] = -spence(1.0 + np.exp(y[neg]))
    pos = y > 0
    yp = y[pos]
    out[pos] = np.pi ** 2 / 6 + yp ** 2 / 2 + spence(1.0 + np.exp(-yp))
    return s_ln ** 2 * out


def model_current(v_g, v_d, device: Device, vth: float, ss: float, mu: float,
                  theta: float = 0.0, r_s: float = 0.0) -> np.ndarray:
    """Noise-free drain current in A (signed). ``ss`` in mV/dec, ``mu`` in cm^2/(V s).

    Optional departures from the ideal model (both zero by default):
    ``theta`` (1/V) degrades the mobility as mu / (1 + theta * q_s), q_s being the
    source-side channel charge per C_i (= V_g - V_th far above threshold);
    ``r_s`` (ohm) is a series resistance split equally between source and
    drain, solved self-consistently by bisection.
    """
    s = device.sign
    s_ln = ss * 1e-3 / np.log(10.0)
    k = device.w_over_l * mu * device.c_i
    v_g, v_d = np.broadcast_arrays(np.asarray(v_g, dtype=float), np.asarray(v_d, dtype=float))

    def intrinsic(vgs, vds):
        x = s * (vgs - vth)
        q_s = s_ln * np.logaddexp(0.0, x / s_ln)
        return s * k * (_g(x, s_ln) - _g(x - s * vds, s_ln)) / (1.0 + theta * q_s)

    i0 = intrinsic(v_g, v_d)
    if r_s == 0:
        return i0
    lo, hi = np.zeros_like(i0), i0.copy()       # the root lies between 0 and the intrinsic current
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        above = np.sign(mid - intrinsic(v_g - mid * r_s / 2, v_d - mid * r_s)) == np.sign(i0)
        hi = np.where(above, mid, hi)
        lo = np.where(above, lo, mid)
    return 0.5 * (lo + hi)


def _measure(i: np.ndarray, i_leak: float, noise_abs: float, noise_rel: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    i = i + np.sign(i[np.argmax(np.abs(i))]) * i_leak
    return i * (1.0 + noise_rel * rng.standard_normal(i.size)) + noise_abs * rng.standard_normal(i.size)


def transfer_sweep(v_g, v_d: float, device: Device, vth: float, ss: float, mu: float,
                   i_leak: float = 0.0, noise_abs: float = 0.0, noise_rel: float = 0.0,
                   seed: int = 0, name: str = "synthetic", theta: float = 0.0, r_s: float = 0.0) -> TransferCurve:
    """A synthetic transfer curve at drain bias ``v_d``.

    ``i_leak`` (A) is a constant off-state current in the direction of the
    on-current; ``noise_abs`` (A) and ``noise_rel`` (fraction) are standard
    deviations of additive and multiplicative Gaussian noise; ``theta`` and
    ``r_s`` as in :func:`model_current`.
    """
    v_g = np.asarray(v_g, dtype=float)
    i = model_current(v_g, v_d, device, vth, ss, mu, theta, r_s)
    return TransferCurve(name, v_g, _measure(i, i_leak, noise_abs, noise_rel, seed), v_d)


def hysteresis_sweep(v_start: float, v_stop: float, step: float, v_d: float, device: Device,
                     vth_fwd: float, dvth: float, ss: float, mu: float, i_leak: float = 0.0,
                     noise_abs: float = 0.0, noise_rel: float = 0.0, seed: int = 0,
                     name: str = "synthetic loop") -> TransferCurve:
    """Dual sweep v_start -> v_stop -> v_start whose backward branch has V_th = vth_fwd + dvth.

    The turning point is recorded twice, as in Clarius dual sweeps.
    """
    n = int(round(abs(v_stop - v_start) / step)) + 1
    fwd = np.linspace(v_start, v_stop, n)
    bwd = fwd[::-1]
    i = np.concatenate([model_current(fwd, v_d, device, vth_fwd, ss, mu),
                        model_current(bwd, v_d, device, vth_fwd + dvth, ss, mu)])
    return TransferCurve(name, np.concatenate([fwd, bwd]), _measure(i, i_leak, noise_abs, noise_rel, seed), v_d)


def output_family(v_d, v_g_steps, device: Device, vth: float, ss: float, mu: float,
                  noise_rel: float = 0.0, seed: int = 0, name: str = "synthetic output") -> OutputFamily:
    """Synthetic output curves I_d(V_d), one per gate voltage in ``v_g_steps``."""
    v_d = np.asarray(v_d, dtype=float)
    rng = np.random.default_rng(seed)
    i_d = [model_current(vg, v_d, device, vth, ss, mu) * (1.0 + noise_rel * rng.standard_normal(v_d.size))
           for vg in v_g_steps]
    return OutputFamily(name, np.asarray(v_g_steps, dtype=float), [v_d.copy() for _ in i_d], i_d)
