"""Shared pieces of the documentation figure scripts: the illustrative device and curves."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402

from transistor_iv import Device  # noqa: E402
from transistor_iv.plot import BAND, INK, INK_2, SERIES, style_axes  # noqa: E402,F401
from transistor_iv.synthetic import hysteresis_sweep, transfer_sweep  # noqa: E402,F401

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "examples" / "data"

# illustrative n-type device: W/L = 100/10 um on 100 nm SiO2 (C_i = 34.5 nF/cm^2)
DEVICE = Device.from_dielectric(100, 10, 3.9, 100, "n")
VTH, SS, MU = 2.0, 200.0, 10.0        # V, mV/dec, cm^2/(V s)
V_G = np.arange(-10, 20.0001, 0.1)


def linear_curve(noise_rel=0.0, seed=0, v_d=0.1, **kw):
    """Linear-regime transfer curve of the illustrative device (0.2 pA leakage, 50 fA noise)."""
    return transfer_sweep(V_G, v_d, DEVICE, VTH, SS, MU, i_leak=2e-13, noise_rel=noise_rel,
                          noise_abs=5e-14, seed=seed, **kw)


def saturation_curve(noise_rel=0.0, seed=0):
    return linear_curve(noise_rel, seed, v_d=20.0)


def save(fig, name, dpi=150):
    fig.savefig(HERE / f"{name}.png", dpi=dpi, metadata={"Software": None})
    print("wrote", HERE / f"{name}.png")
