"""transistor_iv: parameter extraction from transistor transfer and output curves.

Typical use::

    import transistor_iv as tiv
    dev = tiv.Device.from_dielectric(width_um=50, length_um=10, eps_r=3.9, thickness_nm=100, polarity="n")
    curve = tiv.read_clarius_transfer("data.xls")[0]
    res = tiv.analyze_transfer(curve, dev)
    print(res.vth_interpol, res.ss, res.mu_fe_max)
"""

from .batch import analyze_curves, summary_table, write_tables
from .device import Device, gate_capacitance
from .errors import CurveError, ExtractionError
from .hysteresis import HysteresisResult, analyze_hysteresis, loop_area, split_dual_sweep
from .io import (OutputFamily, TransferCurve, read_clarius_output, read_clarius_transfer, read_conditions,
                 read_output_csv, read_transfer_csv, transfer_from_dataframe, write_transfer_csv)
from .output import family_symmetry, on_resistance, summarize_output, symmetry_score
from .sweep import monotonic_segments, select_segment
from .transfer import TransferResult, TransferSettings, analyze_transfer
from .traps import d_it_from_ss, n_trap_from_dvth, thermal_ss_limit

__version__ = "1.0.0"
