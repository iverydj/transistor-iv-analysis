"""Hysteresis of a dual (forward + backward) transfer sweep.

* The sweep is split into monotonic segments (:func:`transistor_iv.monotonic_segments`);
  the forward and backward branches are two of them.
* Each branch is analysed with :func:`transistor_iv.analyze_transfer`.
* Threshold shift, for every V_th method: dV_th = V_th(backward) - V_th(forward).
* Trapped-charge density: N_trap = C_i dV_th / q (cm^-2, signed).
* Loop area in the (V_g, log10|I_d|) plane (dec V, signed), see :func:`loop_area`.

For an n-type device swept up first, dV_th > 0 and a positive loop area mean
that the backward branch carries less current: electrons were trapped during
the forward sweep (a clockwise loop in the I_d-V_g plane).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import trapezoid

from .device import Device
from .errors import CurveError, ExtractionError
from .io.curves import TransferCurve
from .sweep import monotonic_segments
from .transfer.analysis import VTH_FIELDS, TransferResult, TransferSettings, analyze_transfer
from .traps import n_trap_from_dvth


def split_dual_sweep(v_g, segments: tuple[int, int] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Index arrays of the forward and backward branches.

    ``segments=None`` requires exactly two monotonic segments (forward, then
    backward). Otherwise ``segments`` gives their numbers in
    :func:`transistor_iv.monotonic_segments`, e.g. ``(1, 2)`` for a sweep
    0 -> -30 -> 40 -> -30 V whose first segment is the ramp to the start.
    """
    segs = monotonic_segments(v_g)
    if segments is None:
        if len(segs) != 2:
            v = np.asarray(v_g)
            spans = [(float(v[s[0]]), float(v[s[-1]])) for s in segs]
            raise CurveError(f"expected 2 sweep segments, found {len(segs)} {spans}; pass `segments`")
        segments = (0, 1)
    f, b = segments
    if not (0 <= f < len(segs) and 0 <= b < len(segs)) or f == b:
        raise CurveError(f"segments {segments} invalid for a sweep with {len(segs)} segments")
    return segs[f], segs[b]


def loop_area(v_fwd, i_fwd, v_bwd, i_bwd, noise_floor: float = 1e-12) -> float:
    """Signed hysteresis area A = integral [log10|I_fwd| - log10|I_bwd|] dV_g, in dec V.

    Evaluated on the forward branch's V_g points inside the range common to
    both branches (ascending V_g), with log10|I_bwd| interpolated linearly
    there; points where either current is at or below ``noise_floor`` (A) are
    left out so that the noise floor does not contribute, and the trapezoidal
    rule is applied to the remaining points. A > 0: the forward branch carries
    more current.
    """
    vf, af = (np.asarray(x, dtype=float) for x in (v_fwd, i_fwd))
    vb, ab = (np.asarray(x, dtype=float) for x in (v_bwd, i_bwd))
    of, ob = np.argsort(vf, kind="stable"), np.argsort(vb, kind="stable")
    vf, af, vb, ab = vf[of], np.abs(af[of]), vb[ob], np.abs(ab[ob])
    lo, hi = max(vf[0], vb[0]), min(vf[-1], vb[-1])
    if hi <= lo:
        raise CurveError("forward and backward branches do not overlap in V_g")
    keep = (vf >= lo) & (vf <= hi)
    grid, yf = vf[keep], af[keep]
    yb = 10.0 ** np.interp(grid, vb, np.log10(np.maximum(ab, 1e-300)))
    m = (yf > noise_floor) & (yb > noise_floor)
    if m.sum() < 2:
        raise ExtractionError(f"fewer than 2 points above the noise floor {noise_floor:g} A in both branches")
    return float(trapezoid(np.log10(yf[m]) - np.log10(yb[m]), grid[m]))


@dataclass
class HysteresisResult:
    """Forward/backward analysis of one dual sweep (units: V, cm^-2, dec V)."""

    name: str
    forward: TransferResult
    backward: TransferResult
    dvth: dict[str, float]
    n_trap: dict[str, float]
    loop_area: float
    noise_floor: float
    notes: dict[str, str]

    def summary(self) -> dict:
        row: dict[str, object] = {"name": self.name}
        for key in VTH_FIELDS:
            row[f"{key}_fwd"] = getattr(self.forward, key)
            row[f"{key}_bwd"] = getattr(self.backward, key)
            row[f"d{key}"] = self.dvth[key]
            row[f"n_trap_{key[4:]}"] = self.n_trap[key]
        row.update(ss_fwd=self.forward.ss, ss_bwd=self.backward.ss, loop_area=self.loop_area,
                   noise_floor=self.noise_floor, notes="; ".join(f"{k}: {v}" for k, v in self.notes.items()))
        return row


def analyze_hysteresis(curve: TransferCurve, device: Device, settings: TransferSettings = TransferSettings(),
                       segments: tuple[int, int] | None = None, noise_floor: float = 1e-12) -> HysteresisResult:
    """Threshold shift, trapped-charge density and loop area of a dual sweep (see module docstring).

    ``settings.segment`` is ignored: each branch is analysed as a whole.
    """
    f, b = split_dual_sweep(curve.v_g, segments)
    st = replace(settings, segment=-1)

    def branch(idx, tag):
        i_g = None if curve.i_g is None else curve.i_g[idx]
        return TransferCurve(f"{curve.name} {tag}", curve.v_g[idx], curve.i_d[idx], curve.v_d, i_g)

    fwd, bwd = branch(f, "fwd"), branch(b, "bwd")
    rf, rb = analyze_transfer(fwd, device, st), analyze_transfer(bwd, device, st)
    dvth = {k: getattr(rb, k) - getattr(rf, k) for k in VTH_FIELDS}
    n_trap = {k: n_trap_from_dvth(dv, device.c_i) if np.isfinite(dv) else np.nan for k, dv in dvth.items()}
    notes: dict[str, str] = {}
    try:
        area = loop_area(fwd.v_g, fwd.i_d, bwd.v_g, bwd.i_d, noise_floor)
    except ExtractionError as exc:
        area, notes["loop_area"] = np.nan, str(exc)
    return HysteresisResult(curve.name, rf, rb, dvth, n_trap, area, noise_floor, notes)
