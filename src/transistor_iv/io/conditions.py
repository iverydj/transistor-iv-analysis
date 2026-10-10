"""Reader for a per-device geometry table (``Conditions.xls``).

One row per data sheet of the measurement file, in the same order, on a
sheet named ``Conditions``. Columns are matched by header (case, spaces and
underscores ignored):

=========================  ==========================================  ========
quantity                   accepted headers                            required
=========================  ==========================================  ========
channel width (um)         ChannelWidth, Ch_Width, Width, W            yes
channel length (um)        ChannelLength, Ch_Length, Length, L         yes
insulator thickness (um)   GiThickness, Gi_Thick, GiThick, Thickness   yes
polarity ('n' / 'p')       DeviceType, Device_Type, Type               no
=========================  ==========================================  ========

The thickness column is in micrometres, as in the original scripts.
"""

from __future__ import annotations

from pathlib import Path

import xlrd

from ..device import Device, gate_capacitance

ALIASES = {
    "width": ("channelwidth", "chwidth", "width", "w"),
    "length": ("channellength", "chlength", "length", "l"),
    "thick": ("githickness", "githick", "thickness"),
    "type": ("devicetype", "type"),
}


def _key(h) -> str:
    return str(h).strip().lower().replace(" ", "").replace("_", "")


def read_conditions(path: str | Path, eps_r: float, polarity: str | None = None) -> list[Device]:
    """Devices from a ``Conditions`` sheet, one per row.

    ``eps_r`` is the insulator permittivity (the table has no column for it).
    ``polarity`` is used for every row when the table has no type column;
    if the table has one, ``polarity`` must be None.
    """
    sh = xlrd.open_workbook(str(path)).sheet_by_name("Conditions")
    header = [_key(h) for h in sh.row_values(0)]
    col = {}
    for key, names in ALIASES.items():
        hits = [i for i, h in enumerate(header) if h in names]
        if len(hits) > 1:
            raise ValueError(f"{path}: several columns match {key}: {[sh.row_values(0)[i] for i in hits]}")
        if hits:
            col[key] = hits[0]
    for req in ("width", "length", "thick"):
        if req not in col:
            raise ValueError(f"{path}: no {req} column in the Conditions header {sh.row_values(0)}")
    if ("type" in col) == (polarity is not None):
        raise ValueError(f"{path}: give `polarity` exactly when the table has no DeviceType column")
    devices = []
    for r in range(1, sh.nrows):
        row = sh.row_values(r)
        pol = str(row[col["type"]]).strip().lower() if "type" in col else polarity
        c_i = gate_capacitance(eps_r, float(row[col["thick"]]) * 1e3)
        devices.append(Device(float(row[col["width"]]), float(row[col["length"]]), c_i, pol))
    return devices
