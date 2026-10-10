# Changelog

## 1.0.0 - 2026-10-11

First release.

- Readers: Keithley 4200A-SCS Clarius `.xls` transfer and output exports (drain bias from
  the data column or the Settings block), `Conditions.xls` geometry tables, CSV files.
- Transfer-curve analysis following the author's original Clarius scripts: seven threshold
  voltages, subthreshold region and swing, on/off ratio, transconductance, linear,
  effective and saturation mobility, D_it; explicit n/p polarity convention.
- Options: V_d/2 correction, SS floor and alternative SS rules, robust off level,
  log-interpolated constant-current V_th, smoothed second derivative, regime flags and
  window sensitivities.
- Hysteresis: sweep splitting, threshold shift for every method, N_trap, loop area.
- Output curves: symmetry score, low-bias resistance, summary table.
- Batch tables, matplotlib figures, a synthetic compact model, documentation with
  method figures, example scripts and tests.
- Fixes relative to the original scripts: listed in `docs/changes_from_original.md`.
