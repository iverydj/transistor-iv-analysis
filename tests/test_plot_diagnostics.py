import matplotlib.pyplot as plt
import numpy as np

from conftest import read_quiet
from transistor_iv import Device, TransferSettings, analyze_hysteresis, analyze_transfer, read_clarius_output
from transistor_iv import plot
from transistor_iv.synthetic import transfer_sweep
from conftest import DATA


def test_figures_draw():
    dev = Device.from_dielectric(5, 5, 3.9, 300, "n")
    c = read_quiet("TC_example5.xls")[0]
    r = analyze_transfer(c, dev, TransferSettings(segment=1))
    h = analyze_hysteresis(c, dev, segments=(1, 2))
    fam = read_clarius_output(DATA / "OP_ex1.xls")[0]
    for fig in (plot.plot_transfer(c), plot.plot_analysis(r), plot.plot_overlay([c, c]), plot.plot_output(fam),
                plot.plot_hysteresis(h)):
        assert fig.axes
        plt.close(fig)


def test_diagnostics_flags_and_sensitivities():
    dev = Device.from_dielectric(100, 10, 3.9, 100, "n")
    lin = transfer_sweep(np.arange(-10, 20.01, 0.1), 0.1, dev, 1.0, 150, 10)
    sat = transfer_sweep(np.arange(-10, 20.01, 0.1), 40.0, dev, 1.0, 150, 10)
    rl = analyze_transfer(lin, dev, TransferSettings(diagnostics=True))
    rs = analyze_transfer(sat, dev, TransferSettings(diagnostics=True))
    assert rl.diagnostics["flag_linear_regime"] == "" and rl.diagnostics["flag_saturation"] != ""
    assert rs.diagnostics["flag_saturation"] == "" and rs.diagnostics["flag_linear_regime"] != ""
    assert rl.diagnostics["vth_global_sens"] < 0.01
    assert "vth_global_sens" in rl.summary()


def test_options_run_and_change_what_they_should():
    dev = Device.from_dielectric(100, 10, 3.9, 100, "n")
    c = transfer_sweep(np.arange(-10, 20.01, 0.25), 0.1, dev, 1.0, 150, 10, i_leak=1e-13, noise_rel=0.01,
                       noise_abs=3e-14, seed=2)
    base = analyze_transfer(c, dev)
    smooth = analyze_transfer(c, dev, TransferSettings(spline_smoothing=1.0))
    data = analyze_transfer(c, dev, TransferSettings(grid_step=None))
    robust = analyze_transfer(c, dev, TransferSettings(off_level="robust"))
    assert smooth.ss != base.ss and np.isfinite(smooth.ss)
    assert data.curve.u.size == c.v_g.size < base.curve.u.size
    assert robust.i_off > base.i_off and robust.vth_global == base.vth_global
    corrected = analyze_transfer(c, dev, TransferSettings(vd_correction=True))
    np.testing.assert_array_equal(corrected.curves["mu_eff"], base.curves["mu_eff"])   # raw intercept used
    other = analyze_transfer(c, dev, TransferSettings(mobility_vth="vth_global"))
    assert not np.array_equal(other.curves["mu_sat"], base.curves["mu_sat"], equal_nan=True)


def test_plausibility_notes():
    dev = Device.from_dielectric(100, 10, 3.9, 100, "n")
    c = transfer_sweep(np.arange(-10, 20.01, 0.1), 0.1, dev, 2.0, 200, 10, i_leak=2e-13, theta=0.05)
    r = analyze_transfer(c, dev)
    assert "vth_interpol_window" in r.notes          # mobility degradation: window above the g_m maximum
    assert r.vth_interpol < 0 < r.vth_gm_tangent
    wrong = analyze_transfer(c, Device(dev.width_um, dev.length_um, dev.c_i, "p"))
    assert "polarity" in wrong.notes
    convex = transfer_sweep(np.arange(-10, 20.01, 0.1), 0.1, dev, 1.0, 150, 10, theta=-0.02)
    flags = analyze_transfer(convex, dev, TransferSettings(diagnostics=True)).diagnostics
    assert "vth_global" in flags["flag_window_at_end"]      # g_m still rising: steepest window at the end
    assert np.isnan(analyze_transfer(convex, dev, TransferSettings(diagnostics=True, gm_window=None))
                    .diagnostics["mu_fe_max_sens"])
