"""Full transfer-curve analysis: every quantity of one curve in one call."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from ..device import Device
from ..errors import ExtractionError
from ..io.curves import TransferCurve
from ..sweep import select_segment
from ..traps import d_it_from_ss, thermal_ss_limit
from . import mobility as mob
from .interpolation import InterpolatedCurve, interpolate_curve
from .on_off import on_off_ratio
from .swing import (SubthresholdRegion, SwingResult, ss_current_window, ss_decade_span, subthreshold_region,
                    subthreshold_swing, swing_curve)
from .threshold import (LineFit, elr_correction, second_derivative, vth_best_linear_window, vth_constant_current,
                        vth_linear_ss_window, vth_max_gm_tangent, vth_max_log_slope, vth_second_derivative,
                        vth_sqrt_extrapolation)
from .transconductance import transconductance

VTH_FIELDS = ("vth_current", "vth_interpol", "vth_global", "vth_logderivative", "vth_sqrt", "vth_d2",
              "vth_gm_tangent")


@dataclass(frozen=True)
class TransferSettings:
    """Choices of :func:`analyze_transfer`.

    The defaults follow the original scripts (their documented defaults),
    with the bug fixes listed in ``docs/changes_from_original.md``.

    Sweep and interpolation
        ``segment``: monotonic segment analysed (-1 = longest).
        ``grid_step`` (V): largest spacing of the interpolation grid (None = data points only).
        ``spline_smoothing``: smoothing factor of the log10|I_d| spline (0 = exact interpolation).
    Subthreshold swing
        ``ss_sensitivity``: the subthreshold region is where d log10|I_d|/dV_g >= 10**-ss_sensitivity dec/V.
        ``ss_method``: 'pointwise' (default), 'decade_span' or 'current_window' (see :mod:`.swing`).
        ``ss_floor`` (mV/dec): smallest accepted SS ('pointwise', 'decade_span'); None = no floor.
        ``ss_decades``, ``ss_i_min`` (A): span and lowest start current of 'decade_span'.
        ``ss_current_window`` (A, A): current window of 'current_window'.
    Threshold voltage
        ``i_crit`` (A) and ``cc_interpolation`` ('linear' or 'log'): constant-current method.
        ``linear_window`` (V): window of the linear and sqrt fits.
        ``linear_extension`` (V): how far past the subthreshold region ``vth_interpol`` windows may start.
        ``vd_correction``: report vth_interpol, vth_global, vth_gm_tangent and vth_d2 as intercept - V_d/2.
    Derivatives and mobility
        ``gm_window`` (V): smoothing window of g_m (None = unsmoothed np.gradient).
        ``d2_window`` (V): smoothing window of d^2 I_d/dV_g^2 (None = unsmoothed, original).
        ``min_overdrive`` (V): mu_eff and mu_sat are NaN where |V_g - V_th| <= this.
    Other
        ``off_level``: 'min' (smallest |I_d|, original) or 'robust' (see :mod:`.on_off`).
        ``temperature`` (K): for D_it.
        ``diagnostics``: also compute regime flags and window sensitivities (see :mod:`transistor_iv.diagnostics`).
    """

    segment: int = -1
    grid_step: float | None = 0.2
    spline_smoothing: float = 0.0
    ss_sensitivity: float = 1.4
    ss_method: str = "pointwise"
    ss_floor: float | None = 60.0
    ss_decades: float = 1.0
    ss_i_min: float | None = None
    ss_current_window: tuple[float, float] | None = None
    i_crit: float = 1e-8
    cc_interpolation: str = "linear"
    linear_window: float = 5.0
    linear_extension: float = 3.0
    vd_correction: bool = False
    gm_window: float | None = 1.6
    d2_window: float | None = None
    min_overdrive: float = 1.0
    off_level: str = "min"
    temperature: float = 300.0
    diagnostics: bool = False

    def __post_init__(self) -> None:
        if self.ss_method not in ("pointwise", "decade_span", "current_window"):
            raise ValueError(f"ss_method must be 'pointwise', 'decade_span' or 'current_window', got {self.ss_method!r}")
        if self.ss_method == "current_window" and self.ss_current_window is None:
            raise ValueError("ss_method='current_window' needs ss_current_window=(i_low, i_high)")


@dataclass
class TransferResult:
    """Everything extracted from one transfer curve.

    Scalars are NaN when a quantity could not be extracted; the reason is in
    ``notes``. Units: V, A, S, mV/dec, cm^2/(V s), cm^-2 eV^-1.
    ``fits`` holds the straight lines behind vth_interpol, vth_global,
    vth_sqrt and vth_gm_tangent; ``curve`` the interpolated curve and
    ``curves`` the point-wise quantities on its grid (keys ``v_g``, ``i_d``,
    ``dlog``, ``ss``, ``gm``, ``d2``, ``mu_lin``, ``mu_eff``, ``mu_sat``).
    """

    name: str
    polarity: str
    v_d: float
    segment_v: tuple[float, float]
    vth_current: float = np.nan
    vth_interpol: float = np.nan
    vth_global: float = np.nan
    vth_logderivative: float = np.nan
    vth_sqrt: float = np.nan
    vth_d2: float = np.nan
    vth_gm_tangent: float = np.nan
    ss: float = np.nan
    i_on: float = np.nan
    i_off: float = np.nan
    on_off: float = np.nan
    gm_max: float = np.nan
    v_gm_max: float = np.nan
    mu_fe_max: float = np.nan
    mu_sat_sqrt: float = np.nan
    d_it: float = np.nan
    swing: SwingResult | None = None
    region: SubthresholdRegion | None = None
    fits: dict[str, LineFit] = field(default_factory=dict)
    curve: InterpolatedCurve | None = None
    curves: dict[str, np.ndarray] = field(default_factory=dict)
    notes: dict[str, str] = field(default_factory=dict)
    diagnostics: dict[str, object] = field(default_factory=dict)

    SCALARS = ("vth_current", "vth_interpol", "vth_global", "vth_logderivative", "vth_sqrt", "vth_d2",
               "vth_gm_tangent", "ss", "i_on", "i_off", "on_off", "gm_max", "v_gm_max", "mu_fe_max",
               "mu_sat_sqrt", "d_it")

    def summary(self) -> dict:
        """One flat table row: identification, every scalar, fit windows, notes and diagnostics."""
        row = {"name": self.name, "polarity": self.polarity, "v_d": self.v_d,
               "segment_start": self.segment_v[0], "segment_end": self.segment_v[1]}
        row.update({k: getattr(self, k) for k in self.SCALARS})
        row["ss_vmin"], row["ss_vmax"] = self.swing.v_range if self.swing else (np.nan, np.nan)
        row["ss_region_vmin"], row["ss_region_vmax"] = self.region.v_range if self.region else (np.nan, np.nan)
        for key in ("vth_interpol", "vth_global", "vth_sqrt"):
            fit = self.fits.get(key)
            row[f"{key}_fit_vmin"], row[f"{key}_fit_vmax"] = fit.v_range if fit else (np.nan, np.nan)
            row[f"{key}_fit_r"] = fit.r if fit else np.nan
        row.update(self.diagnostics)
        row["notes"] = "; ".join(f"{k}: {v}" for k, v in self.notes.items())
        return row


def _attempt(notes: dict[str, str], key: str, fn):
    try:
        return fn()
    except ExtractionError as exc:
        notes[key] = str(exc)
        return None


def _swing(ic: InterpolatedCurve, region, v, i, pol: str, st: TransferSettings) -> SwingResult:
    if st.ss_method == "pointwise":
        if region is None:
            raise ExtractionError("no subthreshold region")
        return subthreshold_swing(ic, region, st.ss_floor)
    if st.ss_method == "decade_span":
        return ss_decade_span(v, i, pol, st.ss_decades, st.ss_i_min, st.ss_floor)
    return ss_current_window(v, i, pol, *st.ss_current_window)


def analyze_transfer(curve: TransferCurve, device: Device, settings: TransferSettings = TransferSettings()) -> TransferResult:
    """Extract every transfer-curve quantity from ``curve``.

    One monotonic segment (``settings.segment``) is analysed; the on/off
    ratio uses all points of the curve. A curve that cannot be analysed at
    all raises :class:`transistor_iv.CurveError`; a quantity whose criterion
    is not met is NaN with the reason in ``result.notes``.
    """
    st, pol = settings, device.polarity
    seg = select_segment(curve.v_g, st.segment)
    v, i = curve.v_g[seg], curve.i_d[seg]
    ic = interpolate_curve(v, i, pol, st.grid_step, st.spline_smoothing)
    res = TransferResult(name=curve.name, polarity=pol, v_d=curve.v_d, segment_v=(float(v[0]), float(v[-1])),
                         curve=ic)
    notes = res.notes
    vd_known = np.isfinite(curve.v_d) and curve.v_d != 0

    oo = _attempt(notes, "on_off", lambda: on_off_ratio(curve.i_d, st.off_level))
    if oo is not None:
        res.i_on, res.i_off, res.on_off = oo.i_on, oo.i_off, oo.ratio

    cc = _attempt(notes, "vth_current", lambda: vth_constant_current(v, i, pol, st.i_crit, st.cc_interpolation))
    res.vth_current = np.nan if cc is None else cc

    region = _attempt(notes, "ss_region", lambda: subthreshold_region(ic, st.ss_sensitivity))
    res.region = region
    res.swing = _attempt(notes, "ss", lambda: _swing(ic, region, v, i, pol, st))
    if res.swing is not None:
        res.ss = res.swing.ss
        if res.ss >= thermal_ss_limit(st.temperature):
            res.d_it = d_it_from_ss(res.ss, device.c_i, st.temperature)
        else:
            notes["d_it"] = f"SS below the thermionic limit at {st.temperature:g} K"

    fits = res.fits
    if region is not None:
        fits["vth_interpol"] = _attempt(notes, "vth_interpol",
                                        lambda: vth_linear_ss_window(ic, region, st.linear_window, st.linear_extension))
        res.vth_logderivative = vth_max_log_slope(ic, region)
    else:
        notes["vth_interpol"] = notes["vth_logderivative"] = "no subthreshold region"
    fits["vth_global"] = _attempt(notes, "vth_global", lambda: vth_best_linear_window(ic, st.linear_window))
    fits["vth_sqrt"] = _attempt(notes, "vth_sqrt", lambda: vth_sqrt_extrapolation(ic, st.linear_window))
    d2_peak = _attempt(notes, "vth_d2", lambda: vth_second_derivative(ic, st.d2_window))

    gm = _attempt(notes, "gm", lambda: transconductance(ic, st.gm_window))
    if gm is not None:
        k = int(np.argmax(gm))
        res.gm_max, res.v_gm_max = float(gm[k]), float(ic.s * ic.u[k])
        fits["vth_gm_tangent"] = _attempt(notes, "vth_gm_tangent", lambda: vth_max_gm_tangent(ic, gm))
    fits = {k: f for k, f in fits.items() if f is not None}
    res.fits = fits

    for key in ("vth_interpol", "vth_global", "vth_sqrt", "vth_gm_tangent"):
        if key in fits:
            setattr(res, key, fits[key].vth)
    res.vth_d2 = np.nan if d2_peak is None else d2_peak
    if st.vd_correction:
        if vd_known:
            for key in ("vth_interpol", "vth_global", "vth_gm_tangent", "vth_d2"):
                setattr(res, key, elr_correction(getattr(res, key), curve.v_d))
        else:
            notes["vd_correction"] = "drain bias unknown; raw intercepts reported"

    if "vth_sqrt" in fits:
        res.mu_sat_sqrt = mob.saturation_mobility_from_slope(fits["vth_sqrt"].slope, device)

    res.curves = {"v_g": ic.v_g, "i_d": ic.i, "dlog": ic.dlog, "ss": swing_curve(ic),
                  "d2": second_derivative(ic, st.d2_window) if ic.u.size >= 5 else np.full(ic.u.size, np.nan)}
    if gm is not None:
        res.curves["gm"] = gm
    if not vd_known:
        notes["mobility"] = "drain bias unknown or zero"
    elif gm is not None:
        mu_lin = mob.linear_mobility(gm, curve.v_d, device)
        res.curves["mu_lin"] = mu_lin
        res.mu_fe_max = float(np.max(mu_lin))
        if np.isfinite(res.vth_interpol):
            res.curves["mu_eff"] = mob.effective_mobility(ic.v_g, ic.i, res.vth_interpol, curve.v_d, device,
                                                          st.min_overdrive)
            res.curves["mu_sat"] = mob.saturation_mobility(ic.v_g, ic.i, res.vth_interpol, device, st.min_overdrive)

    if st.diagnostics:
        from ..diagnostics import transfer_diagnostics
        res.diagnostics = transfer_diagnostics(curve, device, st, res)
    return res


def settings_dict(settings: TransferSettings) -> dict:
    """The settings as a plain dict (for writing next to a results table)."""
    return asdict(settings)
