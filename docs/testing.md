# Tests and validation

Run the tests with

```bash
pip install -e ".[test]"
pytest
```

(about 5 s on one CPU core). Continuous integration runs them on Python 3.10 and 3.12 with
the latest releases of the dependencies, runs every example script and regenerates every
documentation figure. `environment.yml` pins the versions used for development.

## What is verified

**Against the original scripts** (`tests/test_original_scripts.py`). The original
transfer-curve script was run function by function on every bundled example curve (117
segments, as measured and as the exact p-type mirror); its results are stored in
`tests/data/original_scripts_baseline.csv`. The same script with the bug fixes B1-B7 and B9
of [changes_from_original.md](changes_from_original.md) patched in gives
`tests/data/original_scripts_fixed_baseline.csv`. Both files are written by
`tests/tools/make_original_baselines.py`, which needs a copy of the original script (not part
of this repository). The tests check that

* the package with `gm_window=None` (the original's g_m numerics) reproduces the bug-fixed
  original on every curve of both polarities: the six original threshold voltages within
  0.6 mV (the original rounds to 1 mV), SS, on/off, mu_FE,max and mu_sat,sqrt within 1e-5
  relative, and the same drain bias;
* on `TC_example1.xls`, whose 0.25 V steps and consistent drain bias are untouched by the
  grid fixes, the unmodified original's `vth_current`, `vth_logderivative`, `vth_d2`, SS,
  on/off and mu_FE,max are reproduced directly (n-type);
* the original's bugs are visible in its own results (p-type mobility inflated by
  \|V_d\|/1e-6, NaN constant-current V_th on every dual-sweep segment, second-derivative
  thresholds on round-off spikes);
* the default g_m smoothing (B8) changes only the g_m-based results (`gm_max`, `v_gm_max`,
  `mu_fe_max`, `vth_gm_tangent`) on every example curve.

**On a synthetic compact model with known parameters** (`transistor_iv.synthetic`: an
exponential subthreshold current with swing SS joined smoothly to the long-channel square
law, constant mobility, both polarities; `tests/test_threshold.py`,
`tests/test_swing_mobility_traps.py`, `tests/test_hysteresis_output.py`):

* constant-current V_th against a root finder on the model (1 mV with log interpolation,
  20 mV with linear interpolation on 0.05 V steps);
* linear-regime intercepts of `vth_global` and `vth_gm_tangent` at V_th + V_d/2 and the
  V_d/2 correction (2 mV), `vth_interpol` (0.2 V), the second-derivative peak (one grid step),
  `vth_logderivative` below V_th, the square-root method in saturation (20 mV, mu_sat within
  0.5 %);
* SS by all three rules (within 1 %), the effect of the 60 mV/dec floor on a 45 mV/dec
  device, linear mobility (0.2 %);
* hysteresis: threshold shift and N_trap of a model loop, the loop area against an
  adaptive-quadrature integral (1e-4 relative, no noise floor), the effect of the noise
  floor, a missing loop area (NaN with a note);
* exact p/n mirror symmetry of every scalar result on a noisy curve whose grid needs no
  inserted points;
* the smoothed g_m equals the 9-point Savitzky-Golay derivative on a uniform grid;
* the options `spline_smoothing`, `grid_step=None`, `off_level="robust"`, `mobility_vth`,
  `vd_correction` (which leaves mu_eff unchanged), the diagnostics flags (linear regime,
  saturation, window at the end of the sweep) and sensitivities, and the plausibility notes
  (`vth_interpol_window` on a mobility-degradation curve, `polarity`).

**Formulas by hand**: mobilities, D_it (1.70e11 cm⁻² eV⁻¹ for SS = 200 mV/dec on
11.5 nF/cm²), the thermionic limit (59.53 mV/dec at 300 K), N_trap, on/off.

**Readers and tables** (`tests/test_io_batch.py`): sheet and curve counts of the example
files, drain bias from the column and from the Settings block, the mismatch warning, a
list-sweep Settings block, the Conditions table (thickness in um), output families, CSV
round trips with full precision (V_g, I_d, I_g), an extra text column, batch tables with
commas in sheet names, rejected curves.

**Segments, interpolation, figures, settings validation** (`tests/test_sweep.py`,
`tests/test_interpolation.py`, `tests/test_plot_diagnostics.py`, `tests/test_threshold.py`).

## What is not verified

There is no independent reference for real devices: on measured curves the methods differ
by their definitions (see [threshold voltage](threshold_voltage.md)), and the tests only
show that the implementation does what the documentation says. The synthetic model obeys
the assumptions of the extraction formulas (no contact resistance, constant mobility), so it
checks the implementation, not the physics.
