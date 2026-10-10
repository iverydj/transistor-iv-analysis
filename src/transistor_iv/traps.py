"""Trap densities from the subthreshold swing and from a threshold-voltage shift.

Constants: q = 1.602176634e-19 C, k_B = 1.380649e-23 J/K (scipy.constants).
"""

from __future__ import annotations

import numpy as np
from scipy.constants import Boltzmann, elementary_charge


def thermal_ss_limit(temperature: float = 300.0) -> float:
    """Thermionic limit of the subthreshold swing, ln(10) k_B T / q, in mV/dec (59.5 at 300 K)."""
    return 1000.0 * np.log(10.0) * Boltzmann * temperature / elementary_charge


def d_it_from_ss(ss: float, c_i: float, temperature: float = 300.0) -> float:
    """Interface-trap density from the subthreshold swing, in cm^-2 eV^-1.

    From SS = ln(10) (k_B T / q) (1 + C_it / C_i) with C_it = q D_it
    (C_it in F/cm^2, D_it in cm^-2 eV^-1):

        D_it = (C_i / q) [SS / (ln(10) k_B T / q) - 1]

    ``ss`` in mV/dec, ``c_i`` in F/cm^2. The depletion capacitance of the
    channel is neglected (adequate for thin, fully depleted channels; an upper
    bound otherwise). Raises ``ValueError`` below the thermionic limit, where
    the formula does not apply.
    """
    limit = thermal_ss_limit(temperature)
    if not ss >= limit:
        raise ValueError(f"SS = {ss:.1f} mV/dec is below the thermionic limit {limit:.1f} mV/dec "
                         f"at T = {temperature:g} K; D_it is undefined")
    return float(c_i / elementary_charge * (ss / limit - 1.0))


def n_trap_from_dvth(dvth: float, c_i: float) -> float:
    """Areal density of trapped charge from a threshold shift, N_trap = C_i dV_th / q, in cm^-2.

    Signed: positive for a positive shift. Use the magnitude when only the
    amount of trapped charge matters.
    """
    return float(c_i * dvth / elementary_charge)
