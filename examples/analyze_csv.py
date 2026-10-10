"""Analyse a transfer curve stored as a CSV file (columns V_g, I_d and optionally V_d, I_g).

Edit AppConfig and run:  python examples/analyze_csv.py
"""

from dataclasses import dataclass
from pathlib import Path

import transistor_iv as tiv

HERE = Path(__file__).resolve().parent


@dataclass
class AppConfig:
    csv_file: Path = HERE / "data" / "transfer_example.csv"
    v_d: float | None = None        # V; give it only if the file has no V_d column
    width_um: float = 5.0           # assumed geometry
    length_um: float = 5.0
    eps_r: float = 3.9
    thickness_nm: float = 300.0
    polarity: str = "n"
    vd_correction: bool = False     # report linear-extrapolation V_th as intercept - V_d/2


def main(cfg: AppConfig) -> None:
    curve = tiv.read_transfer_csv(cfg.csv_file, v_d=cfg.v_d)
    device = tiv.Device.from_dielectric(cfg.width_um, cfg.length_um, cfg.eps_r, cfg.thickness_nm, cfg.polarity)
    r = tiv.analyze_transfer(curve, device, tiv.TransferSettings(vd_correction=cfg.vd_correction))
    print(f"{curve.name}: {curve.v_g.size} points, V_d = {curve.v_d:g} V, "
          f"segment {r.segment_v[0]:g} -> {r.segment_v[1]:g} V")
    for key in tiv.TransferResult.SCALARS:
        print(f"  {key:18s} {getattr(r, key): .5g}")
    for key, fit in r.fits.items():
        print(f"  fit {key:15s} V_g {fit.v_range[0]:.2f} .. {fit.v_range[1]:.2f} V, r = {fit.r:.6f}")
    for key, why in r.notes.items():
        print(f"  note {key}: {why}")


if __name__ == "__main__":
    main(AppConfig())
