"""Plain-text (CSV) input and output.

Transfer curves: one file per curve, one row per point, columns for V_g and
I_d (required), V_d (optional, constant) and I_g (optional). Headers are
matched like Clarius headers (``V_g``, ``Vg``, ``GateV`` ... see
:mod:`transistor_iv.io.columns`)::

    V_g,I_d,V_d,I_g
    -20.0,1.2e-12,0.1,3e-13
    ...

Output curves: one file per family, columns V_g, V_d, I_d, one row per point;
rows with the same V_g form one curve::

    V_g,V_d,I_d
    10.0,0.0,1e-13
    10.0,0.25,2.1e-7
    ...
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .clarius import BIAS_TOLERANCE
from .columns import find_column
from .curves import OutputFamily, TransferCurve


def _columns(df: pd.DataFrame, quantities, headers: dict[str, str] | None) -> dict[str, int | None]:
    headers = headers or {}
    hdr = [str(c) for c in df.columns]
    return {q: find_column(hdr, q, headers.get(q)) for q in quantities}


def transfer_from_dataframe(df: pd.DataFrame, name: str, headers: dict[str, str] | None = None) -> TransferCurve:
    """A :class:`TransferCurve` from a DataFrame (column matching as for CSV files)."""
    idx = _columns(df, ("v_g", "i_d", "v_d", "i_g"), headers)
    if idx["v_g"] is None or idx["i_d"] is None:
        raise ValueError(f"{name}: no gate-voltage or drain-current column among {list(df.columns)}")
    def col(q):
        return df.iloc[:, idx[q]].to_numpy(dtype=float)

    v_d = float("nan")
    if idx["v_d"] is not None:
        vd = col("v_d")
        if np.ptp(vd) > max(BIAS_TOLERANCE, 0.01 * np.max(np.abs(vd))):
            raise ValueError(f"{name}: V_d varies by {np.ptp(vd):g} V; a transfer curve needs a constant drain bias")
        v_d = float(np.median(vd))
    i_g = col("i_g") if idx["i_g"] is not None else None
    return TransferCurve(name, col("v_g"), col("i_d"), v_d, i_g)


def read_transfer_csv(path: str | Path, headers: dict[str, str] | None = None, v_d: float | None = None) -> TransferCurve:
    """Read one transfer curve; its name is the file stem.

    ``v_d`` gives the drain bias when the file has no V_d column (passing it
    for a file that has one is an error).
    """
    path = Path(path)
    curve = transfer_from_dataframe(pd.read_csv(path, float_precision="round_trip"), path.stem, headers)
    if v_d is not None:
        if np.isfinite(curve.v_d):
            raise ValueError(f"{path}: the file has a V_d column; do not pass v_d")
        curve.v_d = float(v_d)
    return curve


def write_transfer_csv(curve: TransferCurve, path: str | Path) -> None:
    """Write a curve in the format read by :func:`read_transfer_csv` (full float precision)."""
    cols = {"V_g": curve.v_g, "I_d": curve.i_d}
    if np.isfinite(curve.v_d):
        cols["V_d"] = np.full(curve.v_g.size, curve.v_d)
    if curve.i_g is not None:
        cols["I_g"] = curve.i_g
    pd.DataFrame(cols).to_csv(path, index=False, float_format="%.17g")


def read_output_csv(path: str | Path, headers: dict[str, str] | None = None) -> OutputFamily:
    """Read an output-curve family; one curve per distinct V_g, in order of first appearance."""
    path = Path(path)
    df = pd.read_csv(path, float_precision="round_trip")
    idx = _columns(df, ("v_g", "v_d", "i_d"), headers)
    if None in idx.values():
        raise ValueError(f"{path}: need gate-voltage, drain-voltage and drain-current columns, got {list(df.columns)}")
    vg_col, vd_col, id_col = (df.iloc[:, idx[q]].to_numpy(dtype=float) for q in ("v_g", "v_d", "i_d"))
    levels = pd.unique(vg_col)
    v_d = [vd_col[vg_col == g] for g in levels]
    i_d = [id_col[vg_col == g] for g in levels]
    return OutputFamily(path.stem, np.asarray(levels, dtype=float), v_d, i_d)
