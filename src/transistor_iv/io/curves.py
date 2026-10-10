"""Containers for measured curves."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class TransferCurve:
    """One transfer measurement I_d(V_g) at a constant drain bias.

    ``v_g`` may hold several sweep segments (e.g. a forward/backward loop).
    ``v_d`` is the signed drain bias with the source grounded, NaN if unknown.
    ``i_g`` is the gate current, ``None`` if not measured.
    """

    name: str
    v_g: np.ndarray
    i_d: np.ndarray
    v_d: float = float("nan")
    i_g: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.v_g = np.asarray(self.v_g, dtype=float)
        self.i_d = np.asarray(self.i_d, dtype=float)
        if self.v_g.ndim != 1 or self.v_g.shape != self.i_d.shape:
            raise ValueError(f"{self.name}: v_g and i_d must be 1-D arrays of equal length, "
                             f"got shapes {self.v_g.shape} and {self.i_d.shape}")
        if self.i_g is not None:
            self.i_g = np.asarray(self.i_g, dtype=float)
            if self.i_g.shape != self.v_g.shape:
                raise ValueError(f"{self.name}: i_g must have the same length as v_g")
        self.v_d = float(self.v_d)


@dataclass
class OutputFamily:
    """A family of output curves I_d(V_d), one per gate voltage.

    ``v_g[k]`` is the gate voltage of curve k, ``v_d[k]`` and ``i_d[k]`` its
    drain-voltage sweep and drain current.
    """

    name: str
    v_g: np.ndarray
    v_d: list[np.ndarray] = field(default_factory=list)
    i_d: list[np.ndarray] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.v_g = np.asarray(self.v_g, dtype=float)
        self.v_d = [np.asarray(v, dtype=float) for v in self.v_d]
        self.i_d = [np.asarray(i, dtype=float) for i in self.i_d]
        if not len(self.v_g) == len(self.v_d) == len(self.i_d):
            raise ValueError(f"{self.name}: v_g, v_d and i_d must have one entry per curve")
        for v, i in zip(self.v_d, self.i_d):
            if v.shape != i.shape:
                raise ValueError(f"{self.name}: a V_d sweep and its I_d differ in length")
