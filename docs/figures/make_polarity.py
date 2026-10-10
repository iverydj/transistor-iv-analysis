"""Polarity convention: an n-type curve, its p-type mirror, and the on-direction frame."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import DEVICE, MU, SS, V_G, SERIES, INK_2, save, style_axes
from transistor_iv import Device, TransferCurve, analyze_transfer
from transistor_iv.sweep import on_direction


@dataclass
class AppConfig:
    vth_n: float = 2.0      # V
    v_d: float = 0.1        # V (the p-type mirror uses -v_d)


def main(cfg: AppConfig) -> None:
    from transistor_iv.synthetic import transfer_sweep
    n = transfer_sweep(V_G, cfg.v_d, DEVICE, cfg.vth_n, SS, MU, i_leak=2e-13)
    p_dev = Device(DEVICE.width_um, DEVICE.length_um, DEVICE.c_i, "p")
    p = TransferCurve("p mirror", -n.v_g, -n.i_d, -n.v_d)
    rn, rp = analyze_transfer(n, DEVICE), analyze_transfer(p, p_dev)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax[0].semilogy(n.v_g, np.abs(n.i_d), color=SERIES[0], linewidth=1.5, label=f"n-type, $V_d$ = {n.v_d:g} V: vth_global = {rn.vth_global:.2f} V")
    ax[0].semilogy(p.v_g, np.abs(p.i_d), color=SERIES[1], linewidth=1.5, label=f"p-type, $V_d$ = {p.v_d:g} V: vth_global = {rp.vth_global:.2f} V")
    style_axes(ax[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) measured: p-type turns on towards negative $V_g$")
    ax[0].legend(fontsize=7, frameon=False, loc="upper center")
    for c, pol, col, ls in ((n, "n", SERIES[0], "-"), (p, "p", SERIES[1], "--")):
        u, a, s = on_direction(c.v_g, c.i_d, pol)
        ax[1].semilogy(u, a, color=col, linewidth=1.5, linestyle=ls, label=f"{pol}-type: u = {s:+d} x $V_g$")
    style_axes(ax[1], "on-direction voltage u = s $V_g$ (V)", "$|I_d|$ (A)", "(b) analysed: identical in the frame u = s $V_g$")
    ax[1].legend(fontsize=8, frameon=False, loc="lower right")
    save(fig, "polarity")
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
