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
