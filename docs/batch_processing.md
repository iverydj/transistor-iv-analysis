# Batch processing

## Many curves, one table

```python
import transistor_iv as tiv

curves = tiv.read_clarius_transfer("examples/data/TC_example1.xls")
devices = tiv.read_conditions("examples/data/Conditions.xls", eps_r=3.9)   # one per sheet
settings = tiv.TransferSettings()                                            # defaults
results = tiv.analyze_curves(curves, devices, settings)
table = tiv.summary_table(curves, results)
tiv.write_tables(curves, results, "output/TC_example1")
```

* `analyze_curves(curves, devices, settings)`: `devices` is one `Device` for all curves or
  one per curve. A curve that cannot be analysed (`CurveError`) produces a warning and its
  exception in the result list; any other error stops the batch.
* `summary_table(curves, results)`: a pandas DataFrame with one row per curve,
  `TransferResult.summary()` plus an `error` column.
* `write_tables(curves, results, directory)` writes, with CSV quoting where needed (sheet
  names often contain commas):

| file | content |
|---|---|
| `results.csv` | the summary table |
| `transfer_curves.csv` | measured V_g and I_d of every curve (column pairs `<name> \| V_g (V)`, `<name> \| I_d (A)`) |
| `subthreshold_swing.csv` | V_g and point-wise SS inside each subthreshold region |
| `mobility_fe.csv` | V_g and mu_lin on each interpolated grid |

Shorter columns are padded with empty cells. Keep the settings with the table, e.g.
`pd.Series(tiv.transfer.analysis.settings_dict(settings)).to_csv("settings.csv")`.

## Columns of `results.csv`

| column | meaning |
|---|---|
| `name`, `polarity`, `v_d` | curve, polarity, drain bias used (V) |
| `segment_start`, `segment_end` | V_g at the ends of the analysed segment |
| `vth_current` ... `vth_gm_tangent` | the seven [threshold voltages](threshold_voltage.md) (V) |
| `ss` | [subthreshold swing](subthreshold_swing.md) (mV/dec) |
| `i_on`, `i_off`, `on_off` | [on/off](transconductance_and_on_off.md#onoff-ratio) (A, A, -) |
| `gm_max`, `v_gm_max` | maximum transconductance (S) and its V_g |
| `mu_fe_max`, `mu_sat_sqrt` | [mobilities](mobility.md) (cm²/Vs) |
| `d_it` | [interface-trap density](hysteresis_and_traps.md#interface-trap-density-from-ss) (cm⁻² eV⁻¹) |
| `ss_vmin`, `ss_vmax` | where SS was taken |
| `ss_region_vmin`, `ss_region_vmax` | the subthreshold region |
| `vth_interpol_fit_vmin` ... `vth_sqrt_fit_r` | fit windows (V) and correlation coefficients |
| `flag_*`, `*_sens` | only with `diagnostics=True` |
| `notes` | why a value is NaN |
| `error` | why a curve was rejected |

Correspondence with the original scripts' `results.csv`: `Vth_current` = `vth_current`,
`Vth_interpol` = `vth_interpol`, `Vth_global` = `vth_global`, `Vth_logderivative` =
`vth_logderivative`, `Vth_sqrt` = `vth_sqrt`, `Vth_d2` = `vth_d2`, `onoff ratio` = `on_off`,
`Subthreshold Swing` = `ss`, `u_FE(max)` = `mu_fe_max`, `u_sat_sqrt` = `mu_sat_sqrt`,
`Drain Voltage` = `v_d`.

## Several files

```python
import pandas as pd
rows = []
for path, device in [("a.xls", dev_a), ("b.xls", dev_b)]:
    curves = tiv.read_clarius_transfer(path)
    tab = tiv.summary_table(curves, tiv.analyze_curves(curves, device))
    rows.append(tab.assign(file=path))
pd.concat(rows).to_csv("all_devices.csv", index=False)
```

## Dual sweeps in a batch

```python
loops = tiv.read_clarius_transfer("examples/data/TC_example5.xls")
dev = tiv.Device.from_dielectric(5, 5, 3.9, 300, "n")
tab = pd.DataFrame([tiv.analyze_hysteresis(c, dev, segments=(1, 2)).summary() for c in loops])
```

`HysteresisResult.summary()` gives, per method, the forward and backward V_th, the shift
`dvth_*` and `n_trap_*`, plus `ss_fwd`, `ss_bwd`, `loop_area`.

## Figures

`transistor_iv.plot`: `plot_transfer(curve)` (log and linear), `plot_analysis(result)`
(2 x 4 overview of one analysis), `plot_overlay(curves)`, `plot_hysteresis(result)`,
`plot_output(family)`. See [examples.md](examples.md).
