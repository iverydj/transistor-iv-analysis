"""Figures of docs/examples.md: what the analysis of the bundled example files looks like."""

from dataclasses import dataclass
import warnings

import matplotlib.pyplot as plt

from common import DATA, save
from transistor_iv import (Device, TransferSettings, analyze_hysteresis, analyze_transfer, read_clarius_output,
                           read_clarius_transfer, plot)


@dataclass
class AppConfig:
    transfer_file: str = "TC_example2.xls"
    transfer_sheet: int = 0
    loop_file: str = "TC_example5.xls"
    loop_sheet: int = 0
    loop_segments: tuple = (1, 2)
    output_file: str = "OP_ex1.xls"
    output_sheet: int = 2
    # geometry assumed for the illustration (the files do not store it)
    width_um: float = 5.0
    length_um: float = 5.0
    eps_r: float = 3.9
    thickness_nm: float = 300.0


def main(cfg: AppConfig) -> None:
    dev = Device.from_dielectric(cfg.width_um, cfg.length_um, cfg.eps_r, cfg.thickness_nm, "n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = read_clarius_transfer(DATA / cfg.transfer_file)[cfg.transfer_sheet]
    save(plot.plot_analysis(analyze_transfer(c, dev)), "example_analysis", dpi=110)
    loop = read_clarius_transfer(DATA / cfg.loop_file)[cfg.loop_sheet]
    h = analyze_hysteresis(loop, dev, TransferSettings(), segments=cfg.loop_segments)
    save(plot.plot_hysteresis(h, "vth_current"), "example_hysteresis", dpi=120)
    fam = read_clarius_output(DATA / cfg.output_file)[cfg.output_sheet]
    save(plot.plot_output(fam), "example_output", dpi=120)
    plt.close("all")


if __name__ == "__main__":
    main(AppConfig())
