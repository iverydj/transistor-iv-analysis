"""One summary table for several Clarius files, with per-sheet geometry from Conditions.xls.

Edit AppConfig and run:  python examples/batch_summary.py
"""

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

import transistor_iv as tiv

HERE = Path(__file__).resolve().parent


@dataclass
class AppConfig:
    out_dir: Path = HERE / "output" / "batch_summary"
    eps_r: float = 3.9
    # (file, Conditions file or None, fallback geometry W, L in um and thickness in nm)
    files: list = field(default_factory=lambda: [
        (HERE / "data" / "TC_example1.xls", HERE / "data" / "Conditions.xls", None),
        (HERE / "data" / "TC_example2.xls", None, (5.0, 5.0, 300.0)),
    ])
    polarity: str = "n"
    settings: tiv.TransferSettings = field(default_factory=tiv.TransferSettings)


def main(cfg: AppConfig) -> None:
    tables = []
    for path, cond, geom in cfg.files:
        curves = tiv.read_clarius_transfer(path)
        if cond is not None:
            devices = tiv.read_conditions(cond, cfg.eps_r)
        else:
            devices = tiv.Device.from_dielectric(geom[0], geom[1], cfg.eps_r, geom[2], cfg.polarity)
        results = tiv.analyze_curves(curves, devices, cfg.settings)
        tiv.write_tables(curves, results, cfg.out_dir / path.stem)
        tables.append(tiv.summary_table(curves, results).assign(file=path.name))
    table = pd.concat(tables, ignore_index=True)
    table.to_csv(cfg.out_dir / "all_results.csv", index=False)
    summary = table.groupby("file")[["vth_current", "vth_global", "ss", "on_off", "mu_fe_max"]].median()
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print("median per file:")
        print(summary.to_string())
    print("wrote", cfg.out_dir / "all_results.csv")


if __name__ == "__main__":
    main(AppConfig())
