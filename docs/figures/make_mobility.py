"""Transconductance smoothing and the point-wise mobilities."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import DATA, DEVICE, INK, INK_2, MU, SERIES, linear_curve, saturation_curve, save, style_axes
from transistor_iv import TransferSettings, analyze_transfer, read_clarius_transfer, read_conditions


@dataclass
class AppConfig:
    noise_rel: float = 0.003
    seed: int = 3
    gm_window: float = 1.6       # V, the default
    real_file: str = "TC_example1.xls"
    real_sheet: int = 1


def fig_gm(cfg):
    c = linear_curve(cfg.noise_rel, cfg.seed)
    raw = analyze_transfer(c, DEVICE, TransferSettings(gm_window=None))
    smo = analyze_transfer(c, DEVICE, TransferSettings(gm_window=cfg.gm_window))
    real = read_clarius_transfer(DATA / cfg.real_file)[cfg.real_sheet]
    dev = read_conditions(DATA / "Conditions.xls", eps_r=3.9)[cfg.real_sheet]
    rr = analyze_transfer(real, dev, TransferSettings(gm_window=None))
    rs = analyze_transfer(real, dev, TransferSettings(gm_window=cfg.gm_window))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    for a, r0, r1, title in ((ax[0], raw, smo, f"(a) model device, {cfg.noise_rel * 100:g} % noise, true $\\mu$ = {MU:g}"),
                             (ax[1], rr, rs, f"(b) {cfg.real_file}, sheet '{real.name}'")):
        a.plot(r0.curves["v_g"], r0.curves["mu_lin"], color=INK_2, linewidth=0.7, label="gm_window=None (np.gradient)")
        a.plot(r1.curves["v_g"], r1.curves["mu_lin"], color=SERIES[0], linewidth=1.6, label=f"gm_window={cfg.gm_window:g} V (default)")
        a.plot(r0.v_gm_max, r0.mu_fe_max, "o", color=INK_2, ms=6, label=f"max {r0.mu_fe_max:.3g}")
        a.plot(r1.v_gm_max, r1.mu_fe_max, "o", color=SERIES[1], ms=6, label=f"max {r1.mu_fe_max:.3g}")
        style_axes(a, "$V_g$ (V)", r"$\mu_{lin}$ (cm$^2$ V$^{-1}$ s$^{-1}$)", title)
        a.legend(fontsize=7, frameon=False, loc="upper left")
    save(fig, "gm_smoothing")


def fig_mobilities(cfg):
    rl = analyze_transfer(linear_curve(), DEVICE)
    rs = analyze_transfer(saturation_curve(), DEVICE)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    v = rl.curves["v_g"]
    ax[0].plot(v, rl.curves["mu_lin"], color=SERIES[0], linewidth=1.5, label=r"$\mu_{lin}$ (from $g_m$)")
    ax[0].plot(v, rl.curves["mu_eff"], color=SERIES[1], linewidth=1.5, label=r"$\mu_{eff}$ (from $I_d/(V_g-V_{th})$)")
    ax[0].axvline(rl.vth_interpol, color=INK_2, linestyle=":", linewidth=0.8)
    ax[0].set_ylim(0, 16)
    style_axes(ax[0], "$V_g$ (V)", r"$\mu$ (cm$^2$ V$^{-1}$ s$^{-1}$)", "(a) linear regime, $V_d$ = 0.1 V")
    ax[0].legend(fontsize=8, frameon=False, loc="lower right")
    v = rs.curves["v_g"]
    ax[1].plot(v, rs.curves["mu_sat"], color=SERIES[2], linewidth=1.5, label=r"$\mu_{sat}$ point-wise (with vth_interpol)")
    from transistor_iv.transfer import saturation_mobility
    ax[1].plot(v, saturation_mobility(v, rs.curves["i_d"], rs.vth_sqrt, DEVICE), color=SERIES[6], linewidth=1.5,
               label=r"$\mu_{sat}$ point-wise with vth_sqrt (function call)")
    ax[1].axhline(rs.mu_sat_sqrt, color=SERIES[3], linewidth=1.5, linestyle="--",
                  label=f"mu_sat_sqrt = {rs.mu_sat_sqrt:.2f} (slope of the sqrt fit)")
    ax[1].set_ylim(0, 40)
    style_axes(ax[1], "$V_g$ (V)", r"$\mu$ (cm$^2$ V$^{-1}$ s$^{-1}$)", "(b) saturation, $V_d$ = 20 V")
    ax[1].legend(fontsize=8, frameon=False, loc="upper right")
    save(fig, "mobility_curves")


def main(cfg: AppConfig) -> None:
    fig_gm(cfg)
    fig_mobilities(cfg)
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
