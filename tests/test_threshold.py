"""Threshold-voltage methods on the synthetic compact model (known V_th)."""

import numpy as np
import pytest
from scipy.optimize import brentq

from transistor_iv import ExtractionError, TransferCurve, TransferSettings, analyze_transfer
from transistor_iv.synthetic import model_current, transfer_sweep
from transistor_iv.transfer import interpolate_curve, vth_constant_current

VTH, SS, MU = 1.0, 150.0, 10.0


def sweep(device, v_d, step=0.05):
    s = device.sign
    v = s * np.arange(-10, 20 + step / 2, step)
    return transfer_sweep(v, s * v_d, device, s * VTH, SS, MU)


def test_constant_current_matches_model(device):
    c = sweep(device, 0.1)
    s = device.sign
    exact = brentq(lambda v: abs(model_current(v, s * 0.1, device, s * VTH, SS, MU)) - 1e-8, s * -10, s * 20)
    assert vth_constant_current(c.v_g, c.i_d, device.polarity, 1e-8, "log") == pytest.approx(exact, abs=1e-3)
    assert vth_constant_current(c.v_g, c.i_d, device.polarity, 1e-8, "linear") == pytest.approx(exact, abs=0.02)


def test_constant_current_uses_last_crossing_and_fails_loudly(device):
    s = device.sign
    v = s * np.arange(0.0, 6.0)
    i = s * np.array([2e-8, 1e-9, 5e-9, 2e-8, 1e-7, 1e-6])     # leakage above I_crit at the off end
    assert vth_constant_current(v, i, device.polarity, 1e-8) == pytest.approx(s * (2 + 5 / 15))
    with pytest.raises(ExtractionError, match="never drops"):
        vth_constant_current(v, i, device.polarity, 1e-10)


def test_linear_regime_methods(device):
    c = sweep(device, 0.1)
    raw = analyze_transfer(c, device)
    cor = analyze_transfer(c, device, TransferSettings(vd_correction=True))
    s = device.sign
    # straight-line intercepts lie at V_th + V_d/2 (signed); the correction removes V_d/2
    assert raw.vth_global == pytest.approx(s * VTH + s * 0.05, abs=2e-3)
    assert cor.vth_global == pytest.approx(s * VTH, abs=2e-3)
    assert cor.vth_gm_tangent == pytest.approx(s * VTH, abs=2e-3)
    assert cor.vth_interpol == pytest.approx(s * VTH, abs=0.2)
    assert cor.vth_d2 == pytest.approx(s * VTH, abs=0.03)          # d2 peak at V_th + V_d/2
    assert s * raw.vth_logderivative < VTH                         # steepest point lies below V_th
    assert raw.ss == pytest.approx(SS, rel=5e-3)


def test_saturation_sqrt_method(device):
    c = sweep(device, 40.0)
    r = analyze_transfer(c, device)
    assert r.vth_sqrt == pytest.approx(device.sign * VTH, abs=0.02)
    assert r.mu_sat_sqrt == pytest.approx(MU, rel=5e-3)


def test_vd_correction_needs_the_drain_bias(device):
    c = sweep(device, 0.1)
    c = TransferCurve(c.name, c.v_g, c.i_d)            # V_d unknown
    r = analyze_transfer(c, device, TransferSettings(vd_correction=True))
    assert "vd_correction" in r.notes and "mobility" in r.notes
    assert np.isnan(r.mu_fe_max) and np.isfinite(r.vth_global)


def test_p_type_is_the_exact_mirror_of_n_type():
    from transistor_iv import Device
    n = Device.from_dielectric(100, 10, 3.9, 100, "n")
    p = Device(n.width_um, n.length_um, n.c_i, "p")
    c = transfer_sweep(np.arange(-10, 20.01, 0.2), 0.1, n, VTH, SS, MU, i_leak=1e-13, noise_rel=0.01,
                       noise_abs=1e-13, seed=1)
    m = TransferCurve("mirror", -c.v_g, -c.i_d, -c.v_d)
    st = TransferSettings(vd_correction=True, d2_window=1.0)
    rn, rp = analyze_transfer(c, n, st), analyze_transfer(m, p, st)
    for key in ("vth_current", "vth_interpol", "vth_global", "vth_logderivative", "vth_sqrt", "vth_d2",
                "vth_gm_tangent"):
        assert getattr(rp, key) == pytest.approx(-getattr(rn, key), abs=1e-12), key
    for key in ("ss", "on_off", "gm_max", "mu_fe_max", "mu_sat_sqrt", "d_it"):
        assert getattr(rp, key) == pytest.approx(getattr(rn, key), rel=1e-12), key


def test_interpolated_window_points():
    c = interpolate_curve(np.arange(0, 10.01, 0.25), np.exp(np.arange(0, 10.01, 0.25)), "n")
    assert c.mean_step == pytest.approx(0.125)
    assert c.n_points(5.0, minimum=3) == 40       # 5 V, not 25 points of the nominal 0.2 V step


def test_smoothed_gm_is_savitzky_golay_on_a_uniform_grid():
    from scipy.signal import savgol_filter
    from transistor_iv.transfer import transconductance
    v = np.arange(-5, 5.0001, 0.2)
    i = np.exp(v / 2) * (1 + 0.01 * np.sin(7 * v))
    c = interpolate_curve(v, i, "n")
    np.testing.assert_allclose(transconductance(c, 1.6),
                               np.maximum(savgol_filter(c.i, 9, 2, deriv=1, delta=0.2), 0), rtol=1e-9, atol=1e-12)


def test_settings_are_validated():
    with pytest.raises(ValueError, match="ss_method"):
        TransferSettings(ss_method="minimum")
    with pytest.raises(ValueError, match="current_window"):
        TransferSettings(ss_method="current_window")


def test_window_too_small_for_the_data_gives_nan_with_a_note(device):
    s = device.sign
    c = transfer_sweep(s * np.arange(-10, 20.01, 0.5), s * 0.1, device, s * VTH, SS, MU)
    r = analyze_transfer(c, device, TransferSettings(grid_step=None, d2_window=1.0, gm_window=1.0))
    assert np.isnan(r.vth_d2) and "fewer than 5 samples" in r.notes["vth_d2"]
    assert np.isnan(r.mu_fe_max) and "gm" in r.notes
    assert np.isfinite(r.vth_global)
