"""Output curves and the symmetry score."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from common import DATA, INK, INK_2, SERIES, save, style_axes
from transistor_iv import read_clarius_output, symmetry_score


@dataclass
class AppConfig:
    file: str = "OP_ex1.xls"
    sheet: int = 2
    curve: int = 8               # index of the curve used for the symmetry panel


def main(cfg: AppConfig) -> None:
    fam = read_clarius_output(DATA / cfg.file)[cfg.sheet]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for k, (vg, v, i) in enumerate(zip(fam.v_g, fam.v_d, fam.i_d)):
        ax[0].semilogy(v, np.abs(i), color=cmap(k / max(len(fam.v_g) - 1, 1)), linewidth=1.2, label=f"{vg:g} V")
    style_axes(ax[0], "$V_d$ (V)", "$|I_d|$ (A)", f"(a) {cfg.file}, sheet '{fam.name}'")
    ax[0].legend(title="$V_g$", fontsize=7, title_fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))
    v, i = fam.v_d[cfg.curve], fam.i_d[cfg.curve]
    y = np.log(np.abs(i))
    pos = v > 0
    ax[1].plot(v[pos], y[pos], "o-", color=SERIES[0], ms=3, linewidth=1.2, label="$\\ln|I_d|$ at $+V_d$")
    ax[1].plot(-v[v < 0], y[v < 0], "s-", color=SERIES[1], ms=3, linewidth=1.2, label="$\\ln|I_d|$ at $-V_d$, mirrored")
    ax[1].set_title(f"(b) $V_g$ = {fam.v_g[cfg.curve]:g} V: symmetry score = {symmetry_score(v, i):.1f} %",
                    fontsize=10, color=INK, loc="left")
    style_axes(ax[1], "$|V_d|$ (V)", "$\\ln|I_d|$")
    ax[1].legend(fontsize=8, frameon=False, loc="lower right")
    save(fig, "output_curves")
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
