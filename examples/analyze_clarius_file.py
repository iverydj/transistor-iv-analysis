"""Analyse every transfer curve of one Keithley Clarius .xls export.

Edit AppConfig and run:  python examples/analyze_clarius_file.py
Prints a summary table and writes the tables and one analysis figure per curve.
"""

from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

import transistor_iv as tiv  # noqa: E402
from transistor_iv import plot  # noqa: E402

HERE = Path(__file__).resolve().parent


@dataclass
class AppConfig:
    data_file: Path = HERE / "data" / "TC_example2.xls"
    out_dir: Path = HERE / "output" / "analyze_clarius_file"
    # device: geometry is not stored in the file; these values are assumptions for the example
    width_um: float = 5.0
    length_um: float = 5.0
    eps_r: float = 3.9              # SiO2
    thickness_nm: float = 300.0
    polarity: str = "n"
    # analysis settings (see TransferSettings for every option)
    settings: tiv.TransferSettings = field(default_factory=lambda: tiv.TransferSettings(
        segment=-1,                 # longest monotonic segment
        i_crit=1e-8,                # A, constant-current V_th
        linear_window=5.0,          # V
        gm_window=1.6,              # V
    ))
    figures: bool = True


def main(cfg: AppConfig) -> None:
    device = tiv.Device.from_dielectric(cfg.width_um, cfg.length_um, cfg.eps_r, cfg.thickness_nm, cfg.polarity)
    curves = tiv.read_clarius_transfer(cfg.data_file)
    results = tiv.analyze_curves(curves, device, cfg.settings)
    table = tiv.summary_table(curves, results)
    cols = ["name", "v_d", "vth_current", "vth_interpol", "vth_global", "vth_gm_tangent", "vth_sqrt", "ss",
            "on_off", "mu_fe_max", "mu_sat_sqrt"]
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.float_format", "{:.4g}".format):
        print(table[cols].to_string(index=False))
    for p in tiv.write_tables(curves, results, cfg.out_dir):
        print("wrote", p)
    if cfg.figures:
        for r in results:
            if isinstance(r, tiv.CurveError):
                continue
            fig = plot.plot_analysis(r)
            fig.savefig(cfg.out_dir / f"analysis_{r.name}.png", dpi=110)
            plt.close(fig)
        print("figures in", cfg.out_dir)


if __name__ == "__main__":
    main(AppConfig())
