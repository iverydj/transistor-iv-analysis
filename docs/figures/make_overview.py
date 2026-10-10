"""README figure: where each quantity comes from on one transfer curve."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import INK, INK_2, SERIES, BAND, linear_curve, DEVICE, save, style_axes
from transistor_iv import TransferSettings, analyze_transfer


@dataclass
class AppConfig:
    noise_rel: float = 0.001
    seed: int = 3
    name: str = "overview"


def main(cfg: AppConfig) -> None:
    c = linear_curve(cfg.noise_rel, cfg.seed)
    r = analyze_transfer(c, DEVICE, TransferSettings())
    v, i = r.curves["v_g"], r.curves["i_d"]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)

    a = ax[0]
    a.axvspan(*r.region.v_range, color=BAND, zorder=0)
    a.semilogy(c.v_g, np.abs(c.i_d), color=SERIES[0], linewidth=1.5)
    a.axhline(r.i_on, color=INK_2, linewidth=0.8, linestyle=":")
    a.axhline(r.i_off, color=INK_2, linewidth=0.8, linestyle=":")
    a.text(-9.5, r.i_on * 1.6, "$I_{on}$", color=INK, fontsize=9)
    a.text(-9.5, r.i_off * 1.6, "$I_{off}$  (on/off ratio = $I_{on}/I_{off}$)", color=INK, fontsize=9)
    a.axhline(1e-8, color=SERIES[1], linewidth=0.8, linestyle="--")
    a.plot(r.vth_current, 1e-8, "o", color=SERIES[1], ms=6)
    a.annotate("constant-current $V_{th}$\n($|I_d| = I_{crit}$)", (r.vth_current, 1e-8), (5.5, 2e-10), color=INK, fontsize=9,
               arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8))
    a.plot(r.swing.v_range[0], r.swing.i_range[0], "s", color=SERIES[2], ms=6)
    a.annotate("subthreshold swing SS\n(steepest point of the shaded\nsubthreshold region)",
               (r.swing.v_range[0], r.swing.i_range[0]), (-9.5, 3e-12), color=INK, fontsize=9,
               arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8))
    style_axes(a, "$V_g$ (V)", "$|I_d|$ (A)", "(a) log scale: off state, switching")

    a = ax[1]
    a.plot(v, i, color=SERIES[0], linewidth=1.5)
    for key, col, lab in (("vth_global", SERIES[1], "linear extrapolation"),
                          ("vth_gm_tangent", SERIES[6], "tangent at max $g_m$")):
        f = r.fits[key]
        xs = np.array([f.vth, 20.0])
        a.plot(xs, f(xs), color=col, linewidth=1.2, linestyle="--", label=f"{lab}: {f.vth:.2f} V")
    a.set_ylim(0, i.max() * 1.05)
    a.legend(fontsize=8, frameon=False, loc="upper left")
    style_axes(a, "$V_g$ (V)", "$|I_d|$ (A)", "(b) linear scale: threshold voltage")

    a = ax[2]
    a.plot(v, r.curves["mu_lin"], color=SERIES[0], linewidth=1.5, label=r"$\mu_{lin} = L g_m / (W C_i |V_d|)$")
    a.plot(r.v_gm_max, r.mu_fe_max, "o", color=SERIES[1], ms=6, label=rf"$\mu_{{FE,max}}$ = {r.mu_fe_max:.2f} cm$^2$/Vs")
    a.legend(fontsize=8, frameon=False, loc="lower right")
    style_axes(a, "$V_g$ (V)", r"$\mu$ (cm$^2$ V$^{-1}$ s$^{-1}$)", "(c) transconductance and mobility")
    fig.suptitle("Synthetic n-type device, $V_d$ = 0.1 V: true $V_{th}$ = 2 V, SS = 200 mV/dec, "
                 r"$\mu$ = 10 cm$^2$/Vs", fontsize=10, color=INK)
    save(fig, cfg.name)


if __name__ == "__main__":
    main(AppConfig())
