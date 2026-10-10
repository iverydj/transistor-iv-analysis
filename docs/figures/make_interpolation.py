"""Interpolation figure: data points, inserted grid points and the spline; the single-precision pitfall."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import DATA, INK, INK_2, SERIES, save, style_axes
from transistor_iv import read_clarius_transfer, select_segment
from transistor_iv.transfer import interpolate_curve, second_derivative


@dataclass
class AppConfig:
    file: str = "TC_example5.xls"      # MoS2 FET, V_g stored in single precision, 0.2 V steps
    sheet: int = 0
    segment: int = 1                   # forward sweep -30 -> 40 V
    zoom: tuple = (-1.0, 1.0)          # V, zoom of panel (a)


def fill_gaps_unguarded(u, step):
    """Gap filling exactly as in the original scripts (no tolerance), for comparison."""
    out = [u[0]]
    for x in u[1:]:
        gap = x - out[-1]
        if gap > step:
            n = int(np.ceil(gap / step)) - 1
            out.extend([out[-1] + step * j for j in range(1, n + 1)])
        out.append(x)
    return np.array(out)


def main(cfg: AppConfig) -> None:
    c = read_clarius_transfer(DATA / cfg.file)[cfg.sheet]
    idx = select_segment(c.v_g, cfg.segment)
    v, i = c.v_g[idx], c.i_d[idx]
    ic = interpolate_curve(v, i, "n")
    coarse = interpolate_curve(v[::5], i[::5], "n")            # 1 V steps: points are inserted
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    m = (coarse.u >= 5) & (coarse.u <= 12)
    md = (coarse.u_data >= 5) & (coarse.u_data <= 12)
    ax[0].semilogy(coarse.u[m], coarse.i[m], "-", color=SERIES[0], linewidth=1.2, label="spline of log$_{10}|I_d|$")
    ax[0].semilogy(coarse.u[m], coarse.i[m], "|", color=SERIES[0], ms=8, label="grid points (0.2 V)")
    ax[0].semilogy(coarse.u_data[md], coarse.a_data[md], "o", color=SERIES[1], ms=6, label="data (every 5th point, 1 V)")
    style_axes(ax[0], "$V_g$ (V)", "$|I_d|$ (A)", "(a) data, inserted points and spline")
    ax[0].legend(fontsize=8, frameon=False, loc="lower right")
    # pitfall: near-duplicate points from single-precision V_g
    g_bad = fill_gaps_unguarded(ic.u_data, 0.2)
    from scipy.interpolate import UnivariateSpline
    sp = UnivariateSpline(ic.u_data, np.log10(np.clip(ic.a_data, 1e-20, None)), s=0)
    i_bad = 10 ** sp(g_bad)
    d2_bad = np.gradient(np.gradient(i_bad, g_bad), g_bad)
    d2_ok = second_derivative(ic)
    ax[1].plot(g_bad, d2_bad, color=INK_2, linewidth=0.8, label=f"original gap filling: {np.sum(np.diff(g_bad) < 1e-6)} "
                                                                f"points < 1 µV apart")
    ax[1].plot(ic.u, d2_ok, color=SERIES[0], linewidth=1.4, label="with the 1 % tolerance")
    kb, ko = int(np.argmax(d2_bad)), int(np.argmax(d2_ok))
    ax[1].plot(g_bad[kb], d2_bad[kb], "o", color=INK_2, ms=7, label=f"maximum (vth_d2) at {g_bad[kb]:.2f} V")
    ax[1].plot(ic.u[ko], d2_ok[ko], "o", color=SERIES[0], ms=7, label=f"maximum (vth_d2) at {ic.u[ko]:.2f} V")
    ax[1].set_yscale("symlog", linthresh=np.max(np.abs(d2_ok)))
    ax[1].set_yticks([-1e-7, -1e-8, 0, 1e-8, 1e-7])
    style_axes(ax[1], "$V_g$ (V)", r"$d^2|I_d|/dV_g^2$ (A/V$^2$)", f"(b) {cfg.file}: V_g = 0.2000000x V steps")
    ax[1].legend(fontsize=8, frameon=False, loc="upper left")
    save(fig, "interpolation")
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
