"""Threshold shift, trapped-charge density and loop area of dual transfer sweeps.

Edit AppConfig and run:  python examples/hysteresis_example.py
"""

from dataclasses import dataclass
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
    data_file: Path = HERE / "data" / "TC_example5.xls"   # MoS2 FETs, 0 -> -30 -> 40 -> -30 V
    out_dir: Path = HERE / "output" / "hysteresis_example"
    segments: tuple = (1, 2)       # forward and backward segment numbers (segment 0 is the ramp to -30 V)
    noise_floor: float = 1e-12     # A, excluded from the loop area
    width_um: float = 5.0          # assumed geometry (N_trap needs only C_i)
    length_um: float = 5.0
    eps_r: float = 3.9
    thickness_nm: float = 300.0
    figure_sheet: int = 0
    vth_method: str = "vth_current"


def main(cfg: AppConfig) -> None:
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    device = tiv.Device.from_dielectric(cfg.width_um, cfg.length_um, cfg.eps_r, cfg.thickness_nm, "n")
    curves = tiv.read_clarius_transfer(cfg.data_file)
    print("segments of the first curve:",
          [(c.v_g[s[0]], c.v_g[s[-1]]) for c in curves[:1] for s in tiv.monotonic_segments(c.v_g)])
    results = [tiv.analyze_hysteresis(c, device, segments=cfg.segments, noise_floor=cfg.noise_floor) for c in curves]
    table = pd.DataFrame([h.summary() for h in results])
    table.to_csv(cfg.out_dir / "hysteresis.csv", index=False)
    m = cfg.vth_method
    cols = ["name", f"{m}_fwd", f"{m}_bwd", f"d{m}", f"n_trap_{m[4:]}", "ss_fwd", "ss_bwd", "loop_area"]
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(table[cols].head(10).to_string(index=False))
        print(table[cols[1:]].describe().loc[["mean", "std", "min", "max"]].to_string())
    fig = plot.plot_hysteresis(results[cfg.figure_sheet], m)
    fig.savefig(cfg.out_dir / "hysteresis.png", dpi=120)
    plt.close(fig)
    print("wrote", cfg.out_dir / "hysteresis.csv", "and hysteresis.png")


if __name__ == "__main__":
    main(AppConfig())
