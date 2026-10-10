"""Subthreshold-swing figures: region and point-wise rule, the three rules, the 60 mV/dec floor."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import BAND, DEVICE, INK, INK_2, SERIES, V_G, linear_curve, save, style_axes
from transistor_iv import TransferSettings, analyze_transfer
from transistor_iv.synthetic import transfer_sweep
from transistor_iv.transfer import ss_current_window, ss_decade_span


@dataclass
class AppConfig:
    noise_rel: float = 0.003
    seed: int = 3
    current_window: tuple = (1e-11, 1e-9)   # A, for the 'current_window' rule
    steep_ss: float = 45.0                   # mV/dec, for the floor figure


def fig_region(cfg):
    c = linear_curve(cfg.noise_rel, cfg.seed)
    r = analyze_transfer(c, DEVICE)
    v, ss = r.curves["v_g"], r.curves["ss"]
    sl = slice(r.region.start, r.region.end + 1)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax[0].axvspan(*r.region.v_range, color=BAND, zorder=0, label="subthreshold region")
    ax[0].semilogy(v, r.curves["i_d"], color=SERIES[0], linewidth=1.4)
    ax[0].plot(r.swing.v_range[0], r.swing.i_range[0], "o", color=SERIES[1], ms=6, label=f"SS = {r.ss:.1f} mV/dec")
    style_axes(ax[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) interpolated curve")
    ax[0].legend(fontsize=8, frameon=False, loc="lower right")
    ax[1].axvspan(*r.region.v_range, color=BAND, zorder=0)
    ax[1].semilogy(v, ss, color=INK_2, linewidth=0.8, label="point-wise SS = 1000 / (d log$_{10}|I_d|$/d$V_g$)")
    ax[1].semilogy(v[sl], ss[sl], color=SERIES[0], linewidth=1.6, label="inside the region")
    ax[1].axhline(25119, color=INK_2, linestyle="--", linewidth=0.8, label="region limit, 25 V/dec ($10^{-1.4}$ dec/V)")
    ax[1].axhline(60, color=SERIES[1], linestyle=":", linewidth=1, label="floor 60 mV/dec")
    ax[1].plot(r.swing.v_range[0], r.ss, "o", color=SERIES[1], ms=6, label="reported: smallest value >= floor")
    ax[1].set_ylim(30, 1e6)
    style_axes(ax[1], "$V_g$ (V)", "SS (mV/dec)", "(b) the point-wise rule")
    ax[1].legend(fontsize=7, frameon=True, framealpha=0.9, edgecolor="none", loc="center right", bbox_to_anchor=(1.0, 0.42))
    save(fig, "ss_region")


def fig_rules(cfg):
    c = linear_curve(cfg.noise_rel * 3, cfg.seed)
    r = analyze_transfer(c, DEVICE)
    ds = ss_decade_span(c.v_g, c.i_d, "n")
    cw = ss_current_window(c.v_g, c.i_d, "n", *cfg.current_window)
    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    m = (c.v_g > -1) & (c.v_g < 5)
    ax.semilogy(c.v_g[m], np.abs(c.i_d[m]), "o", ms=2.5, color=INK_2, label="measured")
    ax.plot(r.swing.v_range[0], r.swing.i_range[0], "o", color=SERIES[0], ms=7,
            label=f"pointwise (default): {r.ss:.1f} mV/dec")
    ax.plot(ds.v_range, ds.i_range, color=SERIES[1], linewidth=2.5, label=f"decade_span: {ds.ss:.1f} mV/dec")
    ax.axhspan(*cfg.current_window, color=BAND, zorder=0)
    ax.plot(cw.v_range, cw.i_range, color=SERIES[2], linewidth=2.5,
            label=f"current_window {cfg.current_window[0]:g}-{cfg.current_window[1]:g} A: {cw.ss:.1f} mV/dec")
    style_axes(ax, "$V_g$ (V)", "$|I_d|$ (A)", "Three selection rules on one curve (model SS = 200 mV/dec, 1 % noise)")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    save(fig, "ss_rules")


def fig_floor(cfg):
    c = transfer_sweep(V_G, 0.1, DEVICE, 2.0, cfg.steep_ss, 10.0, i_leak=2e-13)
    r = analyze_transfer(c, DEVICE)
    r0 = analyze_transfer(c, DEVICE, TransferSettings(ss_floor=None))
    v, ss = r.curves["v_g"], r.curves["ss"]
    sl = slice(r.region.start, r.region.end + 1)
    fig, ax = plt.subplots(figsize=(6.5, 4.0), constrained_layout=True)
    ax.semilogy(v[sl], ss[sl], color=SERIES[0], linewidth=1.6, label="point-wise SS in the region")
    ax.axhline(60, color=SERIES[1], linestyle=":", linewidth=1, label="floor 60 mV/dec")
    ax.plot(r.swing.v_range[0], r.ss, "o", color=SERIES[1], ms=7, label=f"default (floor 60): {r.ss:.1f} mV/dec")
    ax.plot(r0.swing.v_range[0], r0.ss, "s", color=SERIES[2], ms=7, label=f"ss_floor=None: {r0.ss:.1f} mV/dec")
    ax.set_xlim(r.region.v_range[0] - 0.2, r.region.v_range[0] + 3)
    ax.set_ylim(30, 3000)
    style_axes(ax, "$V_g$ (V)", "SS (mV/dec)", f"A device with SS = {cfg.steep_ss:g} mV/dec (steep-slope or low temperature)")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    save(fig, "ss_floor")


def main(cfg: AppConfig) -> None:
    fig_region(cfg)
    fig_rules(cfg)
    fig_floor(cfg)
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
