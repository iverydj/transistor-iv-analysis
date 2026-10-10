"""Many curves at once: one results table plus the original scripts' curve tables.

Tables written by :func:`write_tables` (comma-separated, quoted where
needed, one column pair per curve, shorter columns padded with empty cells):

* ``results.csv``: one row per curve, :meth:`TransferResult.summary`;
* ``transfer_curves.csv``: measured V_g and I_d of every curve;
* ``subthreshold_swing.csv``: V_g and point-wise SS inside the subthreshold region;
* ``mobility_fe.csv``: V_g and mu_lin on the interpolated grid.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from .device import Device
from .errors import CurveError
from .io.curves import TransferCurve
from .transfer.analysis import TransferResult, TransferSettings, analyze_transfer


def analyze_curves(curves: Sequence[TransferCurve], devices: Device | Sequence[Device],
                   settings: TransferSettings = TransferSettings()) -> list[TransferResult | CurveError]:
    """Analyse every curve; ``devices`` is one Device for all or one per curve.

    A curve that cannot be analysed (:class:`CurveError`) is reported with a
    warning and returned as the exception in its place; any other error stops
    the batch.
    """
    if isinstance(devices, Device):
        devices = [devices] * len(curves)
    if len(devices) != len(curves):
        raise ValueError(f"{len(devices)} devices given for {len(curves)} curves")
    out: list[TransferResult | CurveError] = []
    for curve, dev in zip(curves, devices):
        try:
            out.append(analyze_transfer(curve, dev, settings))
        except CurveError as exc:
            warnings.warn(f"{curve.name}: {exc}", stacklevel=2)
            out.append(exc)
    return out


def summary_table(curves: Sequence[TransferCurve], results: Sequence[TransferResult | CurveError]) -> pd.DataFrame:
    """One row per curve; rejected curves have their message in the ``error`` column."""
    rows = []
    for c, r in zip(curves, results):
        rows.append({"name": c.name, "error": str(r)} if isinstance(r, CurveError) else r.summary() | {"error": ""})
    return pd.DataFrame(rows)


def _pairs(columns: list[tuple[str, np.ndarray]]) -> pd.DataFrame:
    n = max((len(v) for _, v in columns), default=0)
    return pd.DataFrame({k: np.pad(np.asarray(v, float), (0, n - len(v)), constant_values=np.nan)
                         for k, v in columns})


def write_tables(curves: Sequence[TransferCurve], results: Sequence[TransferResult | CurveError],
                 directory: str | Path) -> list[Path]:
    """Write the four tables listed in the module docstring into ``directory``; returns their paths."""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    tc, ss, mu = [], [], []
    for c, r in zip(curves, results):
        tc += [(f"{c.name} | V_g (V)", c.v_g), (f"{c.name} | I_d (A)", c.i_d)]
        if isinstance(r, CurveError):
            continue
        if r.region is not None:
            sl = slice(r.region.start, r.region.end + 1)
            ss += [(f"{c.name} | V_g (V)", r.curves["v_g"][sl]), (f"{c.name} | SS (mV/dec)", r.curves["ss"][sl])]
        if "mu_lin" in r.curves:
            mu += [(f"{c.name} | V_g (V)", r.curves["v_g"]), (f"{c.name} | mu_lin (cm2/Vs)", r.curves["mu_lin"])]
    paths = []
    for fname, table in (("results.csv", summary_table(curves, results)), ("transfer_curves.csv", _pairs(tc)),
                         ("subthreshold_swing.csv", _pairs(ss)), ("mobility_fe.csv", _pairs(mu))):
        p = d / fname
        table.to_csv(p, index=False)
        paths.append(p)
    return paths
