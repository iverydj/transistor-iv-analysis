"""What each method returns for a device with known parameters (synthetic compact model).

Edit AppConfig and run:  python examples/synthetic_demo.py
The model (transistor_iv.synthetic) has an exponential subthreshold current with swing SS,
the square law above threshold and a constant mobility, for n- or p-type devices.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

import transistor_iv as tiv
from transistor_iv.synthetic import transfer_sweep


@dataclass
class AppConfig:
    width_um: float = 100.0
    length_um: float = 10.0
    eps_r: float = 3.9
    thickness_nm: float = 100.0
    vth: float = 2.0               # V (n-type value; the p-type device uses -vth)
    ss: float = 200.0              # mV/dec
    mu: float = 10.0               # cm^2/(V s)
    v_d_linear: float = 0.1        # V
    v_d_saturation: float = 20.0   # V
    v_start: float = -10.0
    v_stop: float = 20.0
    step: float = 0.1
    noise_rel: float = 0.003       # multiplicative noise of the second set of rows (the first is exact)
    leak: float = 2e-13            # A
    seed: int = 0


def main(cfg: AppConfig) -> None:
    rows = []
    for noise in (0.0, cfg.noise_rel):
        for pol in ("n", "p"):
            s = 1 if pol == "n" else -1
            dev = tiv.Device.from_dielectric(cfg.width_um, cfg.length_um, cfg.eps_r, cfg.thickness_nm, pol)
            v = s * np.arange(cfg.v_start, cfg.v_stop + cfg.step / 2, cfg.step)
            for regime, v_d in (("linear", cfg.v_d_linear), ("saturation", cfg.v_d_saturation)):
                c = transfer_sweep(v, s * v_d, dev, s * cfg.vth, cfg.ss, cfg.mu, i_leak=cfg.leak,
                                   noise_rel=noise, noise_abs=cfg.leak / 4 if noise else 0.0, seed=cfg.seed)
                r = tiv.analyze_transfer(c, dev)
                row = {"noise": noise, "polarity": pol, "regime": regime}
                row.update({k: getattr(r, k) for k in tiv.transfer.VTH_FIELDS + ("ss", "mu_fe_max", "mu_sat_sqrt")})
                rows.append(row)
    table = pd.DataFrame(rows)
    print(f"model: V_th = +-{cfg.vth} V, SS = {cfg.ss} mV/dec, mu = {cfg.mu} cm2/Vs; "
          f"linear-regime intercepts lie at V_th + V_d/2 = {cfg.vth + cfg.v_d_linear / 2:g} V")
    with pd.option_context("display.width", 250, "display.float_format", "{:.3f}".format):
        print(table.to_string(index=False))


if __name__ == "__main__":
    main(AppConfig())
