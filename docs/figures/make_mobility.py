"""Transconductance smoothing and the point-wise mobilities."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import DATA, DEVICE, INK, INK_2, MU, SERIES, linear_curve, saturation_curve, save, style_axes
from transistor_iv import TransferSettings, analyze_transfer, read_clarius_transfer


@dataclass
class AppConfig:
    noise_rel: float = 0.003
    seed: int = 3
    gm_window: float = 1.6       # V, the default
    real_file: str = "TC_example5.xls"   # MoS2 FET: smooth on-state, visible derivative noise
    real_sheet: int = 38
    real_segment: int = 1                # forward sweep -30 -> 40 V
    # window study printed by study(): random model devices, 0.2 V steps from -20 to 30 V, V_d = 0.1 V
    study_devices: int = 60
    study_windows: tuple = (None, 1.6, 3.0)
    study_seed: int = 1


def fig_gm(cfg):
    c = linear_curve(cfg.noise_rel, cfg.seed)
    raw = analyze_transfer(c, DEVICE, TransferSettings(gm_window=None))
    smo = analyze_transfer(c, DEVICE, TransferSettings(gm_window=cfg.gm_window))
    real = read_clarius_transfer(DATA / cfg.real_file)[cfg.real_sheet]
    rr = analyze_transfer(real, DEVICE, TransferSettings(segment=cfg.real_segment, gm_window=None))
    rs = analyze_transfer(real, DEVICE, TransferSettings(segment=cfg.real_segment, gm_window=cfg.gm_window))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    a = ax[0]
    a.plot(raw.curves["v_g"], raw.curves["mu_lin"], color=INK_2, linewidth=0.7, label="gm_window=None (np.gradient)")
    a.plot(smo.curves["v_g"], smo.curves["mu_lin"], color=SERIES[0], linewidth=1.6, label=f"gm_window={cfg.gm_window:g} V (default)")
    a.plot(raw.v_gm_max, raw.mu_fe_max, "o", color=INK_2, ms=6, label=f"max {raw.mu_fe_max:.3g}")
    a.plot(smo.v_gm_max, smo.mu_fe_max, "o", color=SERIES[1], ms=6, label=f"max {smo.mu_fe_max:.3g}")
    style_axes(a, "$V_g$ (V)", r"$\mu_{lin}$ (cm$^2$ V$^{-1}$ s$^{-1}$)",
               f"(a) model device, {cfg.noise_rel * 100:g} % noise, true $\\mu$ = {MU:g}")
    a.legend(fontsize=7, frameon=False, loc="upper left")
    a = ax[1]
    a.plot(rr.curves["v_g"], rr.curves["gm"] * 1e9, color=INK_2, linewidth=0.7, label="gm_window=None (np.gradient)")
    a.plot(rs.curves["v_g"], rs.curves["gm"] * 1e9, color=SERIES[0], linewidth=1.6, label=f"gm_window={cfg.gm_window:g} V (default)")
    a.plot(rr.v_gm_max, rr.gm_max * 1e9, "o", color=INK_2, ms=6, label=f"max {rr.gm_max * 1e9:.3g} nS")
    a.plot(rs.v_gm_max, rs.gm_max * 1e9, "o", color=SERIES[1], ms=6, label=f"max {rs.gm_max * 1e9:.3g} nS")
    style_axes(a, "$V_g$ (V)", "$g_m$ (nS)", f"(b) {cfg.real_file}, sheet '{real.name}', forward sweep")
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


def study(cfg):
    """Bias of mu_fe_max vs gm_window: flat g_m with 0.3 % noise, peaked g_m (theta = 0.05 1/V) without noise."""
    from transistor_iv import Device
    from transistor_iv.synthetic import transfer_sweep
    dev = Device.from_dielectric(100, 10, 3.9, 100, "n")
    rng = np.random.default_rng(cfg.study_seed)
    v = np.arange(-20, 30.001, 0.2)
    err = {(case, w): [] for case in ("flat, 0.3 % noise", "peaked, no noise") for w in cfg.study_windows}
    for k in range(cfg.study_devices):
        vth, ss, mu = rng.uniform(-3, 3), rng.uniform(100, 400), rng.uniform(1, 50)
        for case, noise, theta in (("flat, 0.3 % noise", 0.003, 0.0), ("peaked, no noise", 0.0, 0.05)):
            c = transfer_sweep(v, 0.1, dev, vth, ss, mu, i_leak=1e-13, noise_rel=noise, noise_abs=1e-13 if noise else 0,
                               seed=k, theta=theta)
            fine = transfer_sweep(np.arange(-20, 30.001, 0.002), 0.1, dev, vth, ss, mu, theta=theta)
            true = np.gradient(fine.i_d, fine.v_g).max() / (dev.w_over_l * dev.c_i * 0.1)
            for w in cfg.study_windows:
                r = analyze_transfer(c, dev, TransferSettings(gm_window=w))
                err[(case, w)].append(r.mu_fe_max / true - 1)
    for (case, w), e in err.items():
        print(f"mu_fe_max error, {case}, gm_window={w}: median {np.median(e) * 100:+.1f} %")


def main(cfg: AppConfig) -> None:
    fig_gm(cfg)
    fig_mobilities(cfg)
    study(cfg)
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
