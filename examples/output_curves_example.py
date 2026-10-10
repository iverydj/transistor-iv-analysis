"""Output curves: figures, symmetry score and low-bias resistance.

Edit AppConfig and run:  python examples/output_curves_example.py
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
    data_files: list = field(default_factory=lambda: [HERE / "data" / "OP_ex1.xls", HERE / "data" / "OP_ex2.xls"])
    out_dir: Path = HERE / "output" / "output_curves_example"
    v_d_max: float = 0.5           # V, fit range of the low-bias resistance


def main(cfg: AppConfig) -> None:
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    tables = []
    for path in cfg.data_files:
        for fam in tiv.read_clarius_output(path):
            tab = tiv.summarize_output(fam, cfg.v_d_max).assign(file=path.name)
            tables.append(tab)
            print(f"{path.name} {fam.name}: {len(fam.v_g)} curves, mean symmetry {tiv.family_symmetry(fam):.1f} %")
            fig = plot.plot_output(fam, title=f"{path.name}: {fam.name}")
            fig.savefig(cfg.out_dir / f"output_{path.stem}_{fam.name}.png", dpi=120)
            plt.close(fig)
    table = pd.concat(tables, ignore_index=True)
    table.to_csv(cfg.out_dir / "output_summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(table.head(9).to_string(index=False))
    print("wrote", cfg.out_dir)


if __name__ == "__main__":
    main(AppConfig())
