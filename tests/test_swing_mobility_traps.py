import numpy as np
import pytest
from scipy.constants import elementary_charge as q

from transistor_iv import (Device, ExtractionError, TransferSettings, analyze_transfer, d_it_from_ss,
                           n_trap_from_dvth, thermal_ss_limit)
from transistor_iv.synthetic import transfer_sweep
from transistor_iv.transfer import (effective_mobility, fe_mobility_from_sqrt_fit, interpolate_curve, linear_mobility,
                                    on_off_ratio, saturation_mobility, saturation_mobility_from_slope,
                                    ss_current_window, ss_decade_span, subthreshold_region, subthreshold_swing)


def curve(device, ss=150.0, v_d=0.1, mu=10.0, **kw):
    s = device.sign
    return transfer_sweep(s * np.arange(-10, 20.01, 0.05), s * v_d, device, s * 1.0, ss, mu, **kw)


def test_three_swing_rules_agree_on_an_exponential(device):
    c = curve(device, i_leak=1e-13)
    pol = device.polarity
    ic = interpolate_curve(c.v_g, c.i_d, pol)
    pw = subthreshold_swing(ic, subthreshold_region(ic))
    assert pw.ss == pytest.approx(150, rel=5e-3)
    assert ss_decade_span(c.v_g, c.i_d, pol).ss == pytest.approx(150, rel=1e-2)
    assert ss_current_window(c.v_g, c.i_d, pol, 1e-11, 1e-9).ss == pytest.approx(150, rel=1e-2)


def test_sixty_mv_floor_hides_steeper_swings(device):
    c = curve(device, ss=45.0)                     # e.g. a steep-slope device or low temperature
    default = analyze_transfer(c, device)
    nofloor = analyze_transfer(c, device, TransferSettings(ss_floor=None))
    assert nofloor.ss == pytest.approx(45, rel=1e-2)
    assert default.ss >= 60.0                       # the original rule can only report >= 60 mV/dec
    assert np.isnan(nofloor.d_it) and "d_it" in nofloor.notes


def test_mobility_formulas_by_hand():
    dev = Device(width_um=100, length_um=10, c_i=1e-8, polarity="n")
    assert linear_mobility(1e-6, 0.1, dev) == pytest.approx(10 * 1e-6 / (100 * 1e-8 * 0.1))
    assert linear_mobility(1e-6, -0.1, dev) == pytest.approx(linear_mobility(1e-6, 0.1, dev))   # |V_d|
    with pytest.raises(ValueError):
        linear_mobility(1e-6, 0.0, dev)
    assert effective_mobility([5.0], [1e-6], 1.0, 0.1, dev)[0] == pytest.approx(1e-6 * 10 / (100 * 1e-8 * 0.1 * 4))
    assert saturation_mobility([5.0], [1e-6], 1.0, dev)[0] == pytest.approx(2 * 1e-6 * 10 / (100 * 1e-8 * 16))
    assert np.isnan(saturation_mobility([1.5], [1e-6], 1.0, dev)[0])           # |V_g - V_th| <= 1 V
    assert saturation_mobility_from_slope(1e-3, dev) == pytest.approx(2 * 10 * 1e-6 / (100 * 1e-8))
    assert fe_mobility_from_sqrt_fit(1e-3, 1.0, 11.0, 0.1, dev) == pytest.approx(2 * 10 * 1e-6 * 10 / (100 * 1e-8 * 0.1))


def test_linear_mobility_recovered_with_both_polarities(device):
    r = analyze_transfer(curve(device), device)
    assert r.mu_fe_max == pytest.approx(10.0, rel=2e-3)
    assert r.gm_max > 0


def test_on_off():
    oo = on_off_ratio(np.array([-1e-12, 1e-11, 0.0, 1e-6]))
    assert (oo.i_on, oo.i_off, oo.ratio) == (1e-6, 1e-12, 1e6)
    a = np.r_[np.full(20, 1e-12), 1e-16, np.full(20, 1e-12), np.logspace(-12, -6, 20)]
    assert on_off_ratio(a).i_off == 1e-16                       # the original rule takes the deepest dip
    assert on_off_ratio(a, "robust").i_off == pytest.approx(1e-12)
    with pytest.raises(ExtractionError):
        on_off_ratio(np.zeros(5))


def test_trap_formulas_by_hand():
    assert thermal_ss_limit(300) == pytest.approx(59.53, abs=0.01)
    c_i = 1.151e-8
    expected = c_i / q * (200 / thermal_ss_limit(300) - 1)
    assert d_it_from_ss(200, c_i) == pytest.approx(expected)
    assert d_it_from_ss(200, c_i) == pytest.approx(1.70e11, rel=0.01)
    with pytest.raises(ValueError, match="thermionic"):
        d_it_from_ss(50, c_i)
    assert n_trap_from_dvth(-1.0, c_i) == pytest.approx(-c_i / q)
