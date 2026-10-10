"""Hysteresis: sweep segments, threshold shift and loop area on a synthetic dual sweep."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import BAND, DEVICE, INK, INK_2, SERIES, SS, MU, save, style_axes
from transistor_iv import TransferCurve, analyze_hysteresis, monotonic_segments
from transistor_iv.synthetic import hysteresis_sweep


@dataclass
class AppConfig:
    vth_forward: float = 1.0      # V
    dvth: float = 1.5             # V, backward minus forward
    ramp_start: float = 0.0       # V, the sweep starts here, ramps to v_low, then v_low -> v_high -> v_low
    v_low: float = -10.0
    v_high: float = 20.0
    step: float = 0.2
    noise_floor: float = 1e-12    # A


def main(cfg: AppConfig) -> None:
    loop = hysteresis_sweep(cfg.v_low, cfg.v_high, cfg.step, 0.1, DEVICE, cfg.vth_forward, cfg.dvth, SS, MU,
                            i_leak=2e-13, noise_abs=5e-14, seed=1)
    ramp = np.arange(cfg.ramp_start, cfg.v_low, -2.0)
    i_ramp = np.full(ramp.size, 3e-13)
    curve = TransferCurve("dual sweep", np.r_[ramp, loop.v_g], np.r_[i_ramp, loop.i_d], 0.1)
    h = analyze_hysteresis(curve, DEVICE, segments=(1, 2), noise_floor=cfg.noise_floor)
    segs = monotonic_segments(curve.v_g)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), constrained_layout=True)
    for k, s in enumerate(segs):
        ax[0].plot(s, curve.v_g[s], color=SERIES[k], linewidth=2, label=f"segment {k}")
    style_axes(ax[0], "point index", "$V_g$ (V)", "(a) monotonic_segments: forward = 1, backward = 2")
    ax[0].legend(fontsize=8, frameon=False)
    f, b = h.forward.curve, h.backward.curve
    ax[1].semilogy(f.u_data, f.a_data, color=SERIES[1], linewidth=1.4, label="forward (segment 1)")
    ax[1].semilogy(b.u_data, b.a_data, color=SERIES[2], linewidth=1.4, label="backward (segment 2)")
    m = (f.a_data > cfg.noise_floor) & (np.interp(f.u_data, b.u_data, b.a_data) > cfg.noise_floor)
    ax[1].fill_between(f.u_data[m], f.a_data[m], np.interp(f.u_data[m], b.u_data, b.a_data), color=BAND,
                       label=f"loop area = {h.loop_area:.2f} dec V")
    ax[1].axhline(cfg.noise_floor, color=INK_2, linestyle=":", linewidth=0.8, label="noise_floor")
    y = 1e-8
    ax[1].annotate("", (h.backward.vth_current, y), (h.forward.vth_current, y),
                   arrowprops=dict(arrowstyle="->", color=INK, lw=1.2))
    ax[1].text(4.0, 1e-10, f"dvth_current = {h.dvth['vth_current']:.2f} V\n"
               f"N_trap = {h.n_trap['vth_current']:.2e} cm$^{{-2}}$", fontsize=8, color=INK)
    ax[1].set_xlim(-6, 12)
    style_axes(ax[1], "$V_g$ (V)", "$|I_d|$ (A)", f"(b) model: $V_{{th}}$ shifted by {cfg.dvth:g} V on the way back")
    ax[1].legend(fontsize=7, frameon=False, loc="lower right")
    save(fig, "hysteresis")
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
