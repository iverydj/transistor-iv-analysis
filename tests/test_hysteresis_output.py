import numpy as np
import pytest
from scipy.constants import elementary_charge as q
from scipy.integrate import quad

from transistor_iv import (CurveError, Device, OutputFamily, analyze_hysteresis, family_symmetry, loop_area,
                           on_resistance, split_dual_sweep, summarize_output, symmetry_score)
from transistor_iv.synthetic import hysteresis_sweep, model_current, output_family


def test_threshold_shift_and_trap_density(device):
    s = device.sign
    loop = hysteresis_sweep(s * -10, s * 20, 0.1, s * 0.1, device, s * 1.0, s * 0.8, 150, 10, i_leak=1e-13)
    h = analyze_hysteresis(loop, device)
    for key, tol in (("vth_current", 0.01), ("vth_global", 0.01), ("vth_gm_tangent", 0.01), ("vth_d2", 0.11)):
        assert h.dvth[key] == pytest.approx(s * 0.8, abs=tol), key        # d2: limited by the 0.1 V grid
        assert h.n_trap[key] == pytest.approx(device.c_i * h.dvth[key] / q)
    assert h.loop_area > 0          # the forward branch carries more current


def test_loop_area_against_a_fine_integral():
    dev = Device.from_dielectric(100, 10, 3.9, 100, "n")
    v = np.arange(-2.0, 2.0001, 0.01)
    f = model_current(v, 0.1, dev, 0.0, 150, 10)
    b = model_current(v, 0.1, dev, 0.5, 150, 10)

    def integrand(x):
        return (np.log10(model_current(x, 0.1, dev, 0.0, 150, 10))
                - np.log10(model_current(x, 0.1, dev, 0.5, 150, 10)))

    exact = quad(integrand, -2, 2, limit=200)[0]
    assert loop_area(v, f, v[::-1], b[::-1], noise_floor=0.0) == pytest.approx(exact, rel=1e-4)


def test_split_needs_two_segments():
    v = np.r_[0, -1, -2, np.arange(-1, 3), np.arange(1, -3, -1)]
    with pytest.raises(CurveError, match="3"):
        split_dual_sweep(v)
    f, b = split_dual_sweep(v, (1, 2))
    assert v[f][0] == -2 and v[f][-1] == 2 and v[b][-1] == -2


def test_symmetry_score_pairs_by_voltage_not_index():
    v = np.r_[np.linspace(-5, -0.5, 10), np.linspace(0.25, 5, 20)]     # different step on each side
    i = np.sign(v) * np.exp(np.abs(v))
    assert symmetry_score(v, i) == pytest.approx(100.0, abs=1e-9)       # odd I_d(V_d): symmetric
    asym = np.where(v > 0, 2 * i, i)
    assert symmetry_score(v, asym) < 100.0


def test_on_resistance_and_summary():
    v = np.linspace(-1, 1, 21)
    assert on_resistance(v, v / 1e4, 0.5) == pytest.approx(1e4)
    dev = Device.from_dielectric(100, 10, 3.9, 100, "n")
    fam = output_family(np.r_[np.linspace(-5, -0.25, 20), np.linspace(0.25, 5, 20)], [5, 10, 15], dev, 1.0, 150, 10)
    tab = summarize_output(fam, v_d_max=0.5)
    assert list(tab.v_g) == [5, 10, 15] and (tab.r_on > 0).all()
    assert family_symmetry(fam) == pytest.approx(tab.symmetry.mean())
    with pytest.raises(ValueError):
        OutputFamily("x", [1.0], [np.zeros(3)], [np.zeros(4)])
