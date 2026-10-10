"""On/off ratio: minimum vs robust off level."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import DEVICE, INK_2, SERIES, linear_curve, save, style_axes
from transistor_iv.transfer import on_off_ratio


@dataclass
class AppConfig:
    noise_rel: float = 0.003
    seed: int = 3
    noise_abs: float = 1e-13       # A, instrument noise of the off state


def main(cfg: AppConfig) -> None:
    c = linear_curve(cfg.noise_rel, cfg.seed)
    rng = np.random.default_rng(cfg.seed)
    i = c.i_d + cfg.noise_abs * rng.standard_normal(c.i_d.size)
    lo, rob = on_off_ratio(i, "min"), on_off_ratio(i, "robust")
    fig, ax = plt.subplots(figsize=(6.5, 4.0), constrained_layout=True)
    ax.semilogy(c.v_g, np.abs(i), color=SERIES[0], linewidth=1.2)
    ax.axhline(lo.i_on, color=INK_2, linestyle=":", linewidth=1, label=f"$I_{{on}}$ = max |I_d| = {lo.i_on:.2g} A")
    ax.axhline(lo.i_off, color=SERIES[1], linestyle="--", linewidth=1.2,
               label=f"off='min' (default): {lo.i_off:.2g} A, ratio {lo.ratio:.2g}")
    ax.axhline(rob.i_off, color=SERIES[2], linestyle="--", linewidth=1.2,
               label=f"off='robust': {rob.i_off:.2g} A, ratio {rob.ratio:.2g}")
    style_axes(ax, "$V_g$ (V)", "$|I_d|$ (A)", "Off level of a noisy off state (0.1 pA noise around 0.2 pA)")
    ax.legend(fontsize=8, frameon=False, loc="center right")
    save(fig, "on_off")
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
