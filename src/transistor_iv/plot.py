"""Matplotlib figures of measured curves and analysis results.

Every function draws into new axes and returns the figure; save it with
``fig.savefig(...)``.
"""

from __future__ import annotations

from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np

from .hysteresis import HysteresisResult
from .io.curves import OutputFamily, TransferCurve
from .transfer.analysis import TransferResult

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
BAND = "#ecebe7"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")


def style_axes(ax, xlabel: str, ylabel: str, title: str | None = None) -> None:
    """Recessive grid and spines, text in ink colours."""
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel(ylabel, color=INK)
    if title:
        ax.set_title(title, color=INK, fontsize=10, loc="left")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_2, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)


def _vline(ax, x: float, color: str, label: str) -> None:
    if np.isfinite(x):
        ax.axvline(x, color=color, linewidth=1.2, linestyle="--", label=label)


def _legend(ax) -> None:
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=7, frameon=False, labelcolor=INK)


def plot_transfer(curve: TransferCurve):
    """Measured transfer curve: log|I_d| (and |I_g| if measured) and linear I_d, side by side."""
    fig, (ax_log, ax_lin) = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    ax_log.semilogy(curve.v_g, np.abs(curve.i_d), color=SERIES[0], linewidth=1.5, label="$|I_d|$")
    if curve.i_g is not None:
        ax_log.semilogy(curve.v_g, np.abs(curve.i_g), color=SERIES[1], linewidth=1.0, label="$|I_g|$")
    style_axes(ax_log, "$V_g$ (V)", "current (A)", f"{curve.name}  ($V_d$ = {curve.v_d:g} V)")
    _legend(ax_log)
    ax_lin.plot(curve.v_g, curve.i_d, color=SERIES[0], linewidth=1.5)
    style_axes(ax_lin, "$V_g$ (V)", "$I_d$ (A)")
    return fig


def plot_analysis(result: TransferResult):
    """Overview of one :func:`transistor_iv.analyze_transfer` result (2 x 4 panels).

    (a) log|I_d|: data, interpolation, subthreshold region, vth_current, vth_logderivative;
    (b) |I_d|: vth_interpol, vth_global and vth_gm_tangent lines; (c) sqrt|I_d|: vth_sqrt;
    (d) point-wise SS; (e) g_m; (f) d^2|I_d|/dV_g^2 and vth_d2; (g) mobilities; (h) the numbers.
    """
    r, ic, c = result, result.curve, result.curves
    fig, axs = plt.subplots(2, 4, figsize=(16, 7), constrained_layout=True)
    a = axs.ravel()
    v = c["v_g"]
    # (a) log
    if r.region is not None:
        a[0].axvspan(*r.region.v_range, color=BAND, zorder=0, label="subthreshold region")
    a[0].semilogy(ic.s * ic.u_data, ic.a_data, "o", ms=2.5, color=INK_2, label="data")
    a[0].semilogy(v, c["i_d"], color=SERIES[0], linewidth=1.2, label="interpolated")
    _vline(a[0], r.vth_current, SERIES[1], f"vth_current = {r.vth_current:.2f} V")
    _vline(a[0], r.vth_logderivative, SERIES[2], f"vth_logderivative = {r.vth_logderivative:.2f} V")
    style_axes(a[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) log scale")
    _legend(a[0])
    # (b) linear fits
    a[1].plot(v, c["i_d"], color=INK, linewidth=1.2)
    for key, col in (("vth_interpol", SERIES[1]), ("vth_global", SERIES[2]), ("vth_gm_tangent", SERIES[6])):
        fit = r.fits.get(key)
        if fit is None:
            continue
        xs = np.array(sorted((fit.vth, v[np.argmax(c["i_d"])])))
        a[1].plot(xs, fit(xs), color=col, linewidth=1.2, label=f"{key}: {getattr(r, key):.2f} V")
    a[1].set_ylim(bottom=0)
    style_axes(a[1], "$V_g$ (V)", "$|I_d|$ (A)", "(b) linear extrapolations")
    _legend(a[1])
    # (c) sqrt
    sq = np.sqrt(c["i_d"])
    a[2].plot(v, sq, color=INK, linewidth=1.2)
    fit = r.fits.get("vth_sqrt")
    if fit is not None:
        xs = np.array(sorted((fit.vth, v[np.argmax(sq)])))
        a[2].plot(xs, fit(xs), color=SERIES[3], linewidth=1.2, label=f"vth_sqrt: {r.vth_sqrt:.2f} V")
        a[2].axvspan(*fit.v_range, color=BAND, zorder=0)
    a[2].set_ylim(bottom=0)
    style_axes(a[2], "$V_g$ (V)", r"$\sqrt{|I_d|}$ (A$^{1/2}$)", "(c) square-root extrapolation")
    _legend(a[2])
    # (d) SS
    if r.region is not None:
        sl = slice(r.region.start, r.region.end + 1)
        a[3].plot(v[sl], c["ss"][sl], color=SERIES[0], linewidth=1.2)
        if r.swing is not None and r.swing.method == "pointwise":
            a[3].plot(r.swing.v_range[0], r.ss, "o", ms=6, color=SERIES[1], label=f"SS = {r.ss:.1f} mV/dec")
        a[3].set_yscale("log")
    style_axes(a[3], "$V_g$ (V)", "SS (mV/dec)", "(d) point-wise swing in the region")
    _legend(a[3])
    # (e) gm
    if "gm" in c:
        a[4].plot(v, c["gm"], color=SERIES[0], linewidth=1.2)
        a[4].plot(r.v_gm_max, r.gm_max, "o", ms=6, color=SERIES[1], label=f"max {r.gm_max:.3g} S")
    style_axes(a[4], "$V_g$ (V)", "$g_m$ (S)", "(e) transconductance")
    _legend(a[4])
    # (f) d2
    a[5].plot(v, c["d2"], color=SERIES[0], linewidth=1.0)
    _vline(a[5], r.vth_d2, SERIES[1], f"vth_d2 = {r.vth_d2:.2f} V")
    style_axes(a[5], "$V_g$ (V)", r"$d^2|I_d|/dV_g^2$ (A/V$^2$)", "(f) second derivative")
    _legend(a[5])
    # (g) mobility
    for key, col in (("mu_lin", SERIES[0]), ("mu_eff", SERIES[1]), ("mu_sat", SERIES[2])):
        if key in c:
            a[6].plot(v, c[key], color=col, linewidth=1.2, label=key)
    style_axes(a[6], "$V_g$ (V)", r"$\mu$ (cm$^2$ V$^{-1}$ s$^{-1}$)", "(g) point-wise mobilities")
    _legend(a[6])
    # (h) numbers
    a[7].axis("off")
    lines = [f"{k} = {getattr(r, k):.4g}" for k in r.SCALARS]
    a[7].text(0.0, 1.0, f"{r.name}\nV_d = {r.v_d:g} V, {r.polarity}-type\n\n" + "\n".join(lines),
              va="top", ha="left", fontsize=8, color=INK, family="monospace", transform=a[7].transAxes)
    return fig


def plot_overlay(curves: Sequence[TransferCurve], log: bool = True):
    """All curves in one panel (|I_d| on a log axis, or I_d linear)."""
    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis", max(len(curves), 2))
    for k, c in enumerate(curves):
        y = np.abs(c.i_d) if log else c.i_d
        ax.plot(c.v_g, y, color=cmap(k), linewidth=1.0, label=c.name)
    if log:
        ax.set_yscale("log")
    style_axes(ax, "$V_g$ (V)", "$|I_d|$ (A)" if log else "$I_d$ (A)")
    if len(curves) <= 12:
        ax.legend(fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    return fig


def plot_output(family: OutputFamily, title: str | None = None):
    """Output curves of one family: |I_d| on a log axis and I_d linear, coloured by V_g."""
    fig, (ax_log, ax_lin) = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    vmin, vmax = np.nanmin(family.v_g), np.nanmax(family.v_g)
    for vg, v, i in zip(family.v_g, family.v_d, family.i_d):
        col = cmap(0.5 if vmax == vmin else (vg - vmin) / (vmax - vmin))
        ax_log.semilogy(v, np.abs(i), color=col, linewidth=1.2)
        ax_lin.plot(v, i, color=col, linewidth=1.2, label=f"$V_g$ = {vg:g} V")
    style_axes(ax_log, "$V_d$ (V)", "$|I_d|$ (A)", title or family.name)
    style_axes(ax_lin, "$V_d$ (V)", "$I_d$ (A)")
    ax_lin.legend(fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    return fig


def plot_hysteresis(result: HysteresisResult, vth_key: str = "vth_current"):
    """Forward and backward branches on a log axis with the ``vth_key`` threshold of each."""
    fig, ax = plt.subplots(figsize=(6, 4.2), constrained_layout=True)
    for res, col, tag in ((result.forward, SERIES[0], "forward"), (result.backward, SERIES[1], "backward")):
        ic = res.curve
        ax.semilogy(ic.s * ic.u_data, ic.a_data, color=col, linewidth=1.4, label=tag)
        _vline(ax, getattr(res, vth_key), col, f"{vth_key} {tag} = {getattr(res, vth_key):.2f} V")
    ax.text(0.02, 0.98, f"d{vth_key} = {result.dvth[vth_key]:.2f} V\n"
                        f"N_trap = {result.n_trap[vth_key]:.3g} cm$^{{-2}}$\n"
                        f"loop area = {result.loop_area:.3g} dec V",
            transform=ax.transAxes, va="top", fontsize=8, color=INK)
    style_axes(ax, "$V_g$ (V)", "$|I_d|$ (A)", result.name)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    return fig
