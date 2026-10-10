"""Device description: channel geometry, gate capacitance and carrier polarity.

Units
-----
* channel width W and length L: micrometres (um)
* gate-insulator thickness: nanometres (nm)
* gate capacitance per unit area C_i: F/cm^2
"""

from __future__ import annotations

from dataclasses import dataclass

from scipy.constants import epsilon_0

EPS0_F_PER_CM = epsilon_0 * 1e-2  # vacuum permittivity in F/cm


def gate_capacitance(eps_r: float, thickness_nm: float) -> float:
    """Parallel-plate gate capacitance per unit area, C_i = eps_0 eps_r / t_i, in F/cm^2.

    Example: 300 nm SiO2 (eps_r = 3.9) gives 1.151e-8 F/cm^2 (11.5 nF/cm^2).
    """
    if eps_r <= 0 or thickness_nm <= 0:
        raise ValueError(f"eps_r and thickness_nm must be positive, got {eps_r}, {thickness_nm}")
    return EPS0_F_PER_CM * eps_r / (thickness_nm * 1e-7)


@dataclass(frozen=True)
class Device:
    """A field-effect transistor.

    Parameters
    ----------
    width_um, length_um
        Channel width W and length L in um.
    c_i
        Gate capacitance per unit area in F/cm^2 (see :func:`gate_capacitance`).
    polarity
        ``"n"`` when |I_d| increases with increasing V_g (electron channel),
        ``"p"`` when it increases with decreasing V_g (hole channel).
    """

    width_um: float
    length_um: float
    c_i: float
    polarity: str

    def __post_init__(self) -> None:
        if self.polarity not in ("n", "p"):
            raise ValueError(f"polarity must be 'n' or 'p', got {self.polarity!r}")
        if self.width_um <= 0 or self.length_um <= 0:
            raise ValueError(f"width_um and length_um must be positive, got {self.width_um}, {self.length_um}")
        if self.c_i <= 0:
            raise ValueError(f"c_i must be positive (F/cm^2), got {self.c_i}")

    @classmethod
    def from_dielectric(cls, width_um: float, length_um: float, eps_r: float, thickness_nm: float,
                        polarity: str) -> "Device":
        """Device whose C_i is computed from the insulator permittivity and thickness."""
        return cls(width_um, length_um, gate_capacitance(eps_r, thickness_nm), polarity)

    @property
    def sign(self) -> int:
        """+1 for n-type, -1 for p-type: the direction of V_g that turns the channel on."""
        return 1 if self.polarity == "n" else -1

    @property
    def w_over_l(self) -> float:
        return self.width_um / self.length_um
