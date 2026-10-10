"""Shared preprocessing of a transfer sweep: the interpolated log-current curve.

Every method except the constant-current V_th and the on/off ratio works on
this curve (the original scripts' "spline" step):

1. one monotonic sweep segment is put in ascending V_g order and duplicate
   V_g values are dropped (first occurrence kept), for both polarities;
2. log10|I_d| (|I_d| clipped at 1e-20 A) is interpolated with a cubic
   smoothing spline, ``scipy.interpolate.UnivariateSpline(V_g, log10|I_d|, s=smoothing)``;
   ``smoothing = 0`` (default) passes exactly through every point;
3. the spline is evaluated on the data points plus extra points inserted
   into every gap wider than ``grid_step`` (default 0.2 V) at multiples of
   ``grid_step`` from the gap's lower-V_g end; a point that would land within
   1 % of a step of the next data point is not inserted (single-precision V_g
   values such as 0.2000000x V would otherwise create points 1e-9 V apart);
4. I = 10**spline on that grid, and the log slope in the on-direction
   s d log10|I_d| / dV_g (``np.gradient``), set to 1e-12 dec/V where it is
   not positive;
5. all arrays are then stored in the on-direction frame u = s V_g
   (s = +1 n-type, -1 p-type; see :mod:`transistor_iv.sweep`), ascending in u.

Steps 1-4 work in ascending V_g for both polarities, as in the original
scripts; a p-type curve and the mirror image of an n-type curve therefore
give mirror-image grids only when no points are inserted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import UnivariateSpline

from ..errors import CurveError
from ..sweep import on_direction

GAP_TOLERANCE = 0.01  # fraction of grid_step below which a gap is not split


def fill_gaps(u: np.ndarray, step: float) -> np.ndarray:
    """``u`` (ascending) with points inserted into gaps wider than ``step``.

    A gap g gets ceil(g/step - 0.01) - 1 points at step, 2*step, ... from its
    left end, e.g. a 0.25 V gap one point (+0.2 V), a 1 V gap four, a
    0.2000003 V gap none.
    """
    if step <= 0:
        raise ValueError(f"grid_step must be positive, got {step}")
    out = [float(u[0])]
    for x in u[1:]:
        base = out[-1]
        n = int(np.ceil((x - base) / step - GAP_TOLERANCE)) - 1
        out.extend(base + step * j for j in range(1, n + 1))
        out.append(float(x))
    return np.array(out)


@dataclass(frozen=True)
class InterpolatedCurve:
    """The interpolated curve of one sweep segment, in the on-direction frame.

    Attributes
    ----------
    s : +1 (n-type) or -1 (p-type); V_g = s * u.
    u_data, a_data : ascending unique on-direction gate voltage (V) and |I_d| (A) of the data.
    u : the evaluation grid (V), ascending.
    log10_i : spline value of log10|I_d| on the grid.
    i : 10**log10_i (A).
    dlog : d log10|I_d| / du (dec/V), values <= 0 replaced by 1e-12.
    """

    s: int
    u_data: np.ndarray
    a_data: np.ndarray
    u: np.ndarray
    log10_i: np.ndarray
    i: np.ndarray
    dlog: np.ndarray

    @property
    def v_g(self) -> np.ndarray:
        """Grid in V_g (ascending in the on-direction, i.e. descending for p-type)."""
        return self.s * self.u

    @property
    def mean_step(self) -> float:
        """Average grid spacing (V)."""
        return float((self.u[-1] - self.u[0]) / (self.u.size - 1))

    def n_points(self, width: float, minimum: int = 0) -> int:
        """Number of grid points that make up a ``width``-volt window: floor(width / mean_step).

        A relative tolerance of 1e-6 absorbs the round-off of single-precision V_g.
        """
        return max(int(np.floor(width / self.mean_step * (1 + 1e-6))), minimum)


def interpolate_curve(v_g, i_d, polarity: str, grid_step: float | None = 0.2,
                      smoothing: float = 0.0) -> InterpolatedCurve:
    """Interpolate one monotonic sweep segment (see the module docstring).

    ``grid_step=None`` evaluates the spline on the data points only.
    Raises :class:`CurveError` for fewer than 4 distinct V_g values.
    """
    u, a, s = on_direction(v_g, i_d, polarity)
    v_asc, a_asc = (u, a) if s == 1 else (-u[::-1], a[::-1])          # ascending V_g
    v_data, first = np.unique(v_asc, return_index=True)
    a_data = a_asc[first]
    if v_data.size < 4:
        raise CurveError(f"only {v_data.size} distinct V_g values; the spline needs at least 4")
    spline = UnivariateSpline(v_data, np.log10(np.clip(a_data, 1e-20, None)), s=smoothing)
    grid = v_data if grid_step is None else fill_gaps(v_data, grid_step)
    log10_i = spline(grid)
    dlog = s * np.gradient(log10_i, grid)
    if s == -1:                                                       # to the on-direction frame
        grid, log10_i, dlog = -grid[::-1], log10_i[::-1], dlog[::-1]
        v_data, a_data = -v_data[::-1], a_data[::-1]
    dlog = np.where(dlog > 0, dlog, 1e-12)
    return InterpolatedCurve(s=s, u_data=v_data, a_data=a_data, u=grid, log10_i=log10_i,
                             i=10.0 ** log10_i, dlog=dlog)
