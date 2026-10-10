"""Reader for Keithley 4200A-SCS Clarius ``.xls`` exports.

Layout of a Clarius export (one workbook per test or per appended run):

* sheet 1: the first data sheet; sheets 2 and 3: ``Calc`` and ``Settings``;
  sheets 4, 5, ...: the remaining data sheets;
* every data sheet has a header row and one column per measured quantity
  (``DrainI``, ``DrainV``, ``GateI``, ``GateV``; numbered ``DrainI(1)``,
  ``DrainI(2)``, ... when one test stores several curves, as output-curve
  tests do);
* the ``Settings`` sheet has one block per data sheet, introduced by a
  line of 34 '=' characters, the sheet name and another '=' line, and
  listing the instrument settings per terminal (``Device Terminal``,
  ``Operation Mode``, ``Start/Bias``, ``Step``, ...).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xlrd

from .columns import find_column, find_numbered_columns
from .curves import OutputFamily, TransferCurve

SETTINGS_DIVIDER = "=" * 34
BIAS_TOLERANCE = 1e-3  # V; a bias column varying less than this (or 1 % of its value) is constant


@dataclass
class ClariusSheet:
    """One data sheet: its header, its cells (empty cells as NaN) and its Settings block."""

    name: str
    header: list[str]
    cells: np.ndarray            # (rows, columns) float
    settings: list[list[str]]    # rows of the sheet's block in the Settings sheet

    def column(self, idx: int) -> np.ndarray:
        """One column without the empty cells Clarius pads shorter columns with."""
        v = self.cells[:, idx]
        finite = np.flatnonzero(np.isfinite(v))
        if finite.size == 0:
            return v[:0]
        v = v[:finite[-1] + 1]
        if not np.isfinite(v).all():
            raise ValueError(f"{self.name}: column {self.header[idx]!r} has empty cells inside the data")
        return v

    def settings_drain_bias(self) -> float:
        """Drain bias from the Settings block (NaN if the drain is swept or not listed).

        Rule of the original scripts: if the drain terminal's ``Step`` entry is
        ``N/A``, the bias is its ``Start/Bias`` entry; NaN if that entry is not a
        number (``N/A`` for a list sweep).
        """
        rows = {str(r[0]).strip(): r for r in self.settings if r}
        if "Device Terminal" not in rows or "Start/Bias" not in rows or "Step" not in rows:
            return float("nan")
        terminals = [str(c).strip() for c in rows["Device Terminal"]]
        if "Drain" not in terminals:
            return float("nan")
        col = terminals.index("Drain")
        if str(rows["Step"][col]).strip() != "N/A":
            return float("nan")
        try:
            return float(rows["Start/Bias"][col])
        except ValueError:
            return float("nan")


def _cell(v) -> float:
    if v == "" or v is None:
        return np.nan
    return float(v)


def _settings_blocks(sheet) -> dict[str, list[list[str]]]:
    rows = [[str(v) for v in sheet.row_values(r)] for r in range(sheet.nrows)]
    starts = [i for i in range(len(rows) - 2)
              if rows[i] and rows[i][0] == SETTINGS_DIVIDER and rows[i + 2][0] == SETTINGS_DIVIDER]
    blocks: dict[str, list[list[str]]] = {}
    for k, i in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(rows)
        name = rows[i + 1][0]
        if name in blocks:
            raise ValueError(f"the Settings sheet lists {name!r} twice")
        blocks[name] = rows[i + 3:end]
    return blocks


def read_clarius_sheets(path: str | Path) -> list[ClariusSheet]:
    """Every data sheet of a Clarius ``.xls`` export, in the order Clarius wrote them."""
    wb = xlrd.open_workbook(str(path))
    names = wb.sheet_names()
    if names[1:3] != ["Calc", "Settings"]:
        raise ValueError(f"{path}: not a Clarius export (sheets 2-3 are {names[1:3]}, expected ['Calc', 'Settings'])")
    blocks = _settings_blocks(wb.sheet_by_name("Settings"))
    sheets = []
    for nm in names[:1] + names[3:]:
        if nm not in blocks:
            raise ValueError(f"{path}: data sheet {nm!r} has no block in the Settings sheet")
        sh = wb.sheet_by_name(nm)
        header = [str(h).strip() for h in sh.row_values(0)]
        cells = np.array([[_cell(v) for v in sh.row_values(r)] for r in range(1, sh.nrows)], dtype=float)
        sheets.append(ClariusSheet(nm, header, cells.reshape(-1, len(header)), blocks[nm]))
    return sheets


def _drain_bias(sheet: ClariusSheet, vd_col: int | None) -> float:
    """Drain bias of a transfer sheet: the DrainV column if present, else the Settings block.

    The data column is what the instrument applied; if the Settings block
    disagrees, a warning is issued and the column value is used.
    """
    bias = sheet.settings_drain_bias()
    if vd_col is None:
        return bias
    vd = sheet.column(vd_col)
    if vd.size == 0:
        return bias
    if np.ptp(vd) > max(BIAS_TOLERANCE, 0.01 * np.max(np.abs(vd))):
        raise ValueError(f"{sheet.name}: drain voltage {sheet.header[vd_col]!r} varies by {np.ptp(vd):g} V; "
                         f"this is not a transfer curve")
    value = float(np.median(vd))
    if np.isfinite(bias) and abs(bias - value) > max(BIAS_TOLERANCE, 0.01 * abs(value)):
        warnings.warn(f"{sheet.name}: the DrainV column says V_d = {value:g} V but the Settings block says "
                      f"{bias:g} V; using the column value", stacklevel=3)
    return value


def read_clarius_transfer(path: str | Path, headers: dict[str, str] | None = None) -> list[TransferCurve]:
    """One :class:`TransferCurve` per data sheet of a Clarius transfer export.

    ``headers`` optionally names the columns, e.g. ``{"i_d": "Id_A", "v_g": "Vgate"}``
    (keys 'i_d', 'v_g', 'v_d', 'i_g'); otherwise the default candidates of
    :mod:`transistor_iv.io.columns` are used. I_d and V_g are required; V_d
    (constant) and I_g are optional. Numbered columns (``DrainI(1)``, ...)
    match by prefix, i.e. the first measurement of the sheet is read.
    """
    headers = headers or {}
    curves = []
    for sh in read_clarius_sheets(path):
        idx = {q: find_column(sh.header, q, headers.get(q)) for q in ("i_d", "v_g", "v_d", "i_g")}
        if idx["i_d"] is None or idx["v_g"] is None:
            raise ValueError(f"{path}: sheet {sh.name!r} has no drain-current or gate-voltage column "
                             f"(header {sh.header}); pass `headers`")
        i_d, v_g = sh.column(idx["i_d"]), sh.column(idx["v_g"])
        if i_d.size != v_g.size:
            raise ValueError(f"{sh.name}: I_d and V_g columns differ in length")
        i_g = sh.column(idx["i_g"]) if idx["i_g"] is not None else None
        if i_g is not None and i_g.size != v_g.size:
            raise ValueError(f"{sh.name}: I_g and V_g columns differ in length")
        curves.append(TransferCurve(sh.name, v_g, i_d, _drain_bias(sh, idx["v_d"]), i_g))
    return curves


def read_clarius_output(path: str | Path, headers: dict[str, str] | None = None) -> list[OutputFamily]:
    """One :class:`OutputFamily` per data sheet of a Clarius output-curve export.

    Each sheet holds numbered columns ``DrainI(n)``, ``DrainV(n)``,
    ``GateV(n)``; curve n has gate voltage mean(GateV(n)). Curves are
    ordered by measurement number.
    """
    headers = headers or {}
    families = []
    for sh in read_clarius_sheets(path):
        cols = {q: find_numbered_columns(sh.header, q, headers.get(q)) for q in ("i_d", "v_d", "v_g")}
        if not cols["i_d"] or not cols["v_d"]:
            raise ValueError(f"{path}: sheet {sh.name!r} has no numbered drain-current/voltage columns "
                             f"(header {sh.header}); pass `headers`")
        v_g, v_d, i_d = [], [], []
        for n in sorted(cols["i_d"]):
            if n not in cols["v_d"]:
                raise ValueError(f"{sh.name}: measurement {n} has no drain-voltage column")
            v_d.append(sh.column(cols["v_d"][n]))
            i_d.append(sh.column(cols["i_d"][n]))
            v_g.append(float(np.mean(sh.column(cols["v_g"][n]))) if n in cols["v_g"] else np.nan)
        families.append(OutputFamily(sh.name, np.array(v_g), v_d, i_d))
    return families
