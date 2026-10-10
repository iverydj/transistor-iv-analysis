"""Threshold-voltage figures: one figure per method plus a comparison of all methods."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import BAND, DEVICE, INK, INK_2, SERIES, VTH, linear_curve, saturation_curve, save, style_axes
from transistor_iv import TransferSettings, analyze_transfer
from transistor_iv._numerics import sliding_linear_fits
from transistor_iv.transfer import interpolate_curve, second_derivative


@dataclass
class AppConfig:
    noise_rel: float = 0.003      # multiplicative noise of the illustrative curves
    seed: int = 3
    i_crit: float = 1e-8          # A
    d2_window: float = 1.0        # V, smoothed second derivative shown next to the raw one
    theta: float = 0.02           # 1/V, mobility degradation of the vth_interpol pitfall panel


def _line(ax, fit, x_end, color, label):
    xs = np.array([fit.vth, x_end])
    ax.plot(xs, fit(xs), color=color, linewidth=1.3, linestyle="--", label=label)
    ax.plot(fit.vth, 0, "o", color=color, ms=6)


def fig_current(cfg, c, r):
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax[0].semilogy(c.v_g, np.abs(c.i_d), color=SERIES[0], linewidth=1.4)
    ax[0].axhline(cfg.i_crit, color=SERIES[1], linewidth=1, linestyle="--", label=f"$I_{{crit}}$ = {cfg.i_crit:g} A")
    ax[0].axvline(r.vth_current, color=SERIES[1], linewidth=1, label=f"vth_current = {r.vth_current:.3f} V")
    style_axes(ax[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) the criterion")
    ax[0].legend(fontsize=8, frameon=False)
    m = (c.v_g > r.vth_current - 0.35) & (c.v_g < r.vth_current + 0.35)
    ax[1].plot(c.v_g[m], np.abs(c.i_d[m]), "o", color=SERIES[0], ms=5, label="measured points")
    k = np.flatnonzero(np.abs(c.i_d) < cfg.i_crit)[-1]
    ax[1].plot(c.v_g[k:k + 2], np.abs(c.i_d[k:k + 2]), color=SERIES[1], linewidth=1.3,
               label="linear interpolation between\nlast point below and next point")
    ax[1].axhline(cfg.i_crit, color=SERIES[1], linewidth=1, linestyle="--")
    ax[1].plot(r.vth_current, cfg.i_crit, "s", color=INK, ms=6)
    style_axes(ax[1], "$V_g$ (V)", "$|I_d|$ (A)", "(b) zoom: the crossing")
    ax[1].legend(fontsize=8, frameon=False)
    save(fig, "vth_current")


def fig_interpol(cfg, c, r):
    ic, reg, fit = r.curve, r.region, r.fits["vth_interpol"]
    n = ic.n_points(5.0, 3)
    last = min(reg.end + ic.n_points(3.0), ic.u.size - 1)
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8), constrained_layout=True)
    a = ax[0]
    a.axvspan(*reg.v_range, color=BAND, zorder=0, label="subthreshold region")
    a.axvspan(ic.u[reg.start], ic.u[last], ymax=0.06, color=SERIES[2], alpha=0.6, label="allowed window starts")
    a.plot(ic.v_g, ic.i, color=INK, linewidth=1.2)
    sl = slice(int(np.searchsorted(ic.u, fit.v_range[0])), int(np.searchsorted(ic.u, fit.v_range[1])) + 1)
    a.plot(ic.v_g[sl], ic.i[sl], "o", color=SERIES[1], ms=3, label="best window (5 V)")
    _line(a, fit, fit.v_range[1], SERIES[1], f"vth_interpol = {r.vth_interpol:.2f} V")
    a.set_xlim(-2, 14)
    a.set_ylim(0, ic.i[ic.u <= 14].max())
    style_axes(a, "$V_g$ (V)", "$|I_d|$ (A)", "(a) linear fit next to the subthreshold region")
    a.legend(fontsize=8, frameon=False, loc="upper left")
    _, _, rr = sliding_linear_fits(ic.u, ic.i, n)
    starts = np.arange(reg.start, min(last, rr.size - 1) + 1)
    ax[1].plot(ic.v_g[starts], np.abs(rr[starts]), color=SERIES[0], linewidth=1.4)
    k = starts[np.argmax(np.abs(rr[starts]))]
    ax[1].plot(ic.v_g[k], abs(rr[k]), "o", color=SERIES[1], ms=6, label="largest |r|")
    style_axes(ax[1], "window start $V_g$ (V)", "|r| of the window's line fit", "(b) the selection")
    ax[1].legend(fontsize=8, frameon=False)
    # pitfall: mobility degradation makes the curve concave; the most linear window moves into the on-state
    cd = linear_curve(theta=cfg.theta)
    rd = analyze_transfer(cd, DEVICE)
    fd, fg = rd.fits["vth_interpol"], rd.fits["vth_global"]
    a = ax[2]
    a.plot(rd.curves["v_g"], rd.curves["i_d"], color=INK, linewidth=1.2)
    _line(a, fd, fd.v_range[1], SERIES[1], f"vth_interpol = {rd.vth_interpol:.2f} V (window {fd.v_range[0]:.1f}-{fd.v_range[1]:.1f} V)")
    _line(a, fg, fg.v_range[1], SERIES[2], f"vth_global = {rd.vth_global:.2f} V")
    a.axvline(rd.v_gm_max, color=INK_2, linestyle=":", linewidth=0.8, label=f"$g_m$ maximum ({rd.v_gm_max:.1f} V)")
    a.set_ylim(bottom=0)
    style_axes(a, "$V_g$ (V)", "$|I_d|$ (A)", f"(c) pitfall: mobility degradation, theta = {cfg.theta:g} 1/V")
    a.legend(fontsize=7, frameon=False, loc="upper left")
    save(fig, "vth_interpol")


def fig_global(cfg, c, r):
    ic, fit = r.curve, r.fits["vth_global"]
    n = ic.n_points(5.0, 3)
    s, _, rr = sliding_linear_fits(ic.u, ic.i, n)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    a = ax[0]
    a.plot(ic.v_g, ic.i, color=INK, linewidth=1.2)
    sl = slice(int(np.searchsorted(ic.u, fit.v_range[0])), int(np.searchsorted(ic.u, fit.v_range[1])) + 1)
    a.plot(ic.v_g[sl], ic.i[sl], "o", color=SERIES[1], ms=3, label="best window (5 V)")
    _line(a, fit, fit.v_range[1], SERIES[1], f"vth_global = {r.vth_global:.2f} V")
    a.set_ylim(bottom=0)
    style_axes(a, "$V_g$ (V)", "$|I_d|$ (A)", "(a) best linear window anywhere in the sweep")
    a.legend(fontsize=8, frameon=False, loc="upper left")
    score = np.abs(rr) * np.abs(s)
    ax[1].plot(ic.v_g[:score.size], score, color=SERIES[0], linewidth=1.4)
    k = int(np.argmax(score))
    ax[1].plot(ic.v_g[k], score[k], "o", color=SERIES[1], ms=6, label="largest |r| x |slope|")
    style_axes(ax[1], "window start $V_g$ (V)", "|r| x |slope| (A/V)", "(b) the selection")
    ax[1].legend(fontsize=8, frameon=False)
    save(fig, "vth_global")


def fig_logderivative(cfg, c, r):
    ic, reg = r.curve, r.region
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax[0].axvspan(*reg.v_range, color=BAND, zorder=0, label="subthreshold region")
    ax[0].semilogy(ic.v_g, ic.i, color=SERIES[0], linewidth=1.4)
    ax[0].axvline(r.vth_logderivative, color=SERIES[1], linewidth=1, label=f"vth_logderivative = {r.vth_logderivative:.2f} V")
    ax[0].axvline(VTH, color=INK_2, linewidth=0.8, linestyle=":", label=f"model $V_{{th}}$ = {VTH:g} V")
    style_axes(ax[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) log scale")
    ax[0].legend(fontsize=8, frameon=False, loc="lower right")
    ax[1].axvspan(*reg.v_range, color=BAND, zorder=0)
    ax[1].plot(ic.v_g, ic.dlog, color=SERIES[0], linewidth=1.4)
    ax[1].axhline(10 ** -1.4, color=INK_2, linewidth=0.8, linestyle="--", label="region threshold $10^{-1.4}$ dec/V")
    k = int(np.argmax(np.where((ic.v_g >= reg.v_range[0]) & (ic.v_g <= reg.v_range[1]), ic.dlog, -1)))
    ax[1].plot(ic.v_g[k], ic.dlog[k], "o", color=SERIES[1], ms=6, label="maximum")
    style_axes(ax[1], "$V_g$ (V)", r"$d\,\log_{10}|I_d| / dV_g$ (dec/V)", "(b) log slope")
    ax[1].legend(fontsize=8, frameon=False)
    save(fig, "vth_logderivative")


def fig_sqrt(cfg, cs, rs):
    ic, fit = rs.curve, rs.fits["vth_sqrt"]
    y = np.sqrt(ic.i)
    n = ic.n_points(5.0, 3)
    _, _, rr = sliding_linear_fits(ic.u, y, n)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    a = ax[0]
    a.plot(ic.v_g, y, color=INK, linewidth=1.2)
    sl = slice(int(np.searchsorted(ic.u, fit.v_range[0])), int(np.searchsorted(ic.u, fit.v_range[1])) + 1)
    a.plot(ic.v_g[sl], y[sl], "o", color=SERIES[3], ms=3, label="most linear window (5 V)")
    _line(a, fit, 20.0, SERIES[3], f"vth_sqrt = {rs.vth_sqrt:.2f} V, mu_sat_sqrt = {rs.mu_sat_sqrt:.2f}")
    a.set_ylim(bottom=0)
    style_axes(a, "$V_g$ (V)", r"$\sqrt{|I_d|}$ (A$^{1/2}$)", "(a) saturation, $V_d$ = 20 V")
    a.legend(fontsize=8, frameon=False, loc="upper left")
    ax[1].plot(ic.v_g[:rr.size], np.abs(rr), color=SERIES[0], linewidth=1.4)
    k = int(np.argmax(np.abs(rr)))
    ax[1].plot(ic.v_g[k], abs(rr[k]), "o", color=SERIES[3], ms=6, label="largest |r|")
    ax[1].set_ylim(0.9, 1.0005)
    style_axes(ax[1], "window start $V_g$ (V)", r"|r| of the $\sqrt{|I_d|}$ fit", "(b) the selection")
    ax[1].legend(fontsize=8, frameon=False)
    save(fig, "vth_sqrt")


def fig_d2(cfg, c, r):
    ic = interpolate_curve(c.v_g, c.i_d, "n")
    raw, smooth = second_derivative(ic, None), second_derivative(ic, cfg.d2_window)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax[0].plot(ic.v_g, ic.i, color=INK, linewidth=1.2)
    ax[0].axvline(VTH, color=INK_2, linewidth=0.8, linestyle=":", label=f"model $V_{{th}}$ = {VTH:g} V")
    ax[0].set_xlim(-2, 8)
    ax[0].set_ylim(0, ic.i[ic.u <= 8].max())
    style_axes(ax[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) the turn-on")
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].plot(ic.v_g, raw, color=SERIES[0], linewidth=0.8, alpha=0.8, label="np.gradient twice (default)")
    ax[1].plot(ic.v_g, smooth, color=SERIES[1], linewidth=1.6, label=f"local quadratic, {cfg.d2_window:g} V window")
    for d2, col in ((raw, SERIES[0]), (smooth, SERIES[1])):
        k = int(np.argmax(d2))
        ax[1].plot(ic.v_g[k], d2[k], "o", color=col, ms=6)
    style_axes(ax[1], "$V_g$ (V)", r"$d^2|I_d|/dV_g^2$ (A/V$^2$)", "(b) second derivative, 0.3 % noise (dots: maxima)")
    ax[1].legend(fontsize=8, frameon=False)
    save(fig, "vth_d2")


def fig_gm_tangent(cfg, c, r):
    fit = r.fits["vth_gm_tangent"]
    v = r.curves["v_g"]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax[0].plot(v, r.curves["i_d"], color=INK, linewidth=1.2)
    _line(ax[0], fit, 20.0, SERIES[6], f"vth_gm_tangent = {r.vth_gm_tangent:.2f} V")
    ax[0].plot(r.v_gm_max, fit(r.v_gm_max), "s", color=SERIES[6], ms=6)
    ax[0].set_ylim(bottom=0)
    style_axes(ax[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) tangent at the $g_m$ maximum")
    ax[0].legend(fontsize=8, frameon=False, loc="upper left")
    ax[1].plot(v, r.curves["gm"], color=SERIES[0], linewidth=1.4)
    ax[1].plot(r.v_gm_max, r.gm_max, "s", color=SERIES[6], ms=6, label="$g_{m,max}$")
    style_axes(ax[1], "$V_g$ (V)", "$g_m$ (S)", "(b) smoothed transconductance")
    ax[1].legend(fontsize=8, frameon=False)
    save(fig, "vth_gm_tangent")


def fig_comparison(rl, rs):
    keys = ["vth_current", "vth_interpol", "vth_global", "vth_gm_tangent", "vth_logderivative", "vth_d2", "vth_sqrt"]
    fig, ax = plt.subplots(figsize=(8, 3.6), constrained_layout=True)
    y = np.arange(len(keys))
    ax.plot([getattr(rl, k) for k in keys], y + 0.12, "o", color=SERIES[0], ms=7, label="linear regime, $V_d$ = 0.1 V")
    ax.plot([getattr(rs, k) for k in keys], y - 0.12, "s", color=SERIES[1], ms=7, label="saturation, $V_d$ = 20 V")
    ax.axvline(VTH, color=INK_2, linewidth=0.8, linestyle=":")
    ax.text(VTH, -0.75, " model $V_{th}$", color=INK, fontsize=8)
    ax.set_yticks(y, keys)
    style_axes(ax, "$V_{th}$ (V)", "")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.set_ylim(-0.9, len(keys) - 0.5)
    ax.set_title("Noise-free model device at two drain biases: each method answers a different question", fontsize=10,
                 color=INK, loc="left")
    save(fig, "vth_comparison")


def main(cfg: AppConfig) -> None:
    st = TransferSettings(i_crit=cfg.i_crit)
    c = linear_curve(cfg.noise_rel, cfg.seed)
    cs = saturation_curve(cfg.noise_rel, cfg.seed)
    r, rs = analyze_transfer(c, DEVICE, st), analyze_transfer(cs, DEVICE, st)
    fig_current(cfg, c, r)
    fig_interpol(cfg, c, r)
    fig_global(cfg, c, r)
    fig_logderivative(cfg, c, r)
    fig_sqrt(cfg, cs, rs)
    fig_d2(cfg, c, r)
    fig_gm_tangent(cfg, c, r)
    fig_comparison(analyze_transfer(linear_curve(), DEVICE, st), analyze_transfer(saturation_curve(), DEVICE, st))
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
