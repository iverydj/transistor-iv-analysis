# API reference

Everything listed here is importable from `transistor_iv` (top level) or from the module
named in the heading. Docstrings give the formulas and units; the pages linked in each
section explain the methods.

## Devices — `transistor_iv.device`

| name | description |
|---|---|
| `Device(width_um, length_um, c_i, polarity)` | frozen dataclass; `c_i` in F/cm², `polarity` `"n"` or `"p"`; properties `sign` (+1/-1), `w_over_l` |
| `Device.from_dielectric(width_um, length_um, eps_r, thickness_nm, polarity)` | C_i from the insulator |
| `gate_capacitance(eps_r, thickness_nm)` | C_i = eps_0 eps_r / t_i in F/cm² |

## Reading data — `transistor_iv.io` ([input formats](input_formats.md))

| name | description |
|---|---|
| `TransferCurve(name, v_g, i_d, v_d=nan, i_g=None)` | one transfer measurement |
| `OutputFamily(name, v_g, v_d, i_d)` | output curves, one per gate voltage |
| `read_clarius_transfer(path, headers=None)` | list of `TransferCurve`, one per data sheet |
| `read_clarius_output(path, headers=None)` | list of `OutputFamily`, one per data sheet |
| `read_clarius_sheets(path)` | raw sheets: `ClariusSheet(name, header, cells, settings)`, `.column(i)`, `.settings_drain_bias()` |
| `read_conditions(path, eps_r, polarity=None)` | list of `Device` from a `Conditions` sheet |
| `read_transfer_csv(path, headers=None, v_d=None)`, `write_transfer_csv(curve, path)` | CSV transfer curves |
| `transfer_from_dataframe(df, name, headers=None)` | from a pandas DataFrame |
| `read_output_csv(path, headers=None)` | CSV output family |

## Sweeps — `transistor_iv.sweep` ([polarity and units](polarity_and_units.md))

| name | description |
|---|---|
| `monotonic_segments(v_g)` | index arrays of the monotonic segments |
| `select_segment(v_g, index=-1)` | one segment (-1 = longest) |
| `on_direction(v_g, i_d, polarity)` | `(u, abs_i, s)` in the on-direction frame |

## Transfer curves — `transistor_iv.transfer`

**Full analysis**

| name | description |
|---|---|
| `analyze_transfer(curve, device, settings=TransferSettings())` | every quantity of one curve -> `TransferResult` |
| `TransferSettings(...)` | every choice; defaults follow the original scripts except the bug fixes (`gm_window`); fields in the table below |
| `TransferResult` | scalars (`vth_current`, `vth_interpol`, `vth_global`, `vth_logderivative`, `vth_sqrt`, `vth_d2`, `vth_gm_tangent`, `ss`, `i_on`, `i_off`, `on_off`, `gm_max`, `v_gm_max`, `mu_fe_max`, `mu_sat_sqrt`, `d_it`), `swing`, `region`, `fits`, `curve`, `curves`, `notes`, `diagnostics`; `summary()` -> dict |
| `VTH_FIELDS` | the seven threshold-voltage field names |

`TransferSettings` fields: `segment=-1`, `grid_step=0.2`, `spline_smoothing=0.0`,
`ss_sensitivity=1.4`, `ss_method="pointwise"`, `ss_floor=60.0`, `ss_decades=1.0`,
`ss_i_min=None`, `ss_current_window=None`, `i_crit=1e-8`, `cc_interpolation="linear"`,
`linear_window=5.0`, `linear_extension=3.0`, `vd_correction=False`, `gm_window=1.6`,
`d2_window=None`, `min_overdrive=1.0`, `mobility_vth="vth_interpol"`, `off_level="min"`,
`temperature=300.0`, `diagnostics=False`. Invalid choices raise `ValueError`.

`TransferResult.notes` explains every NaN and carries plausibility notes that do not change
any value: `polarity`, `ss_on_state`, `vth_interpol_window`, `<fit>_window` (truncated
window).

**Building blocks** (each works on one monotonic segment)

| name | page |
|---|---|
| `interpolate_curve(v_g, i_d, polarity, grid_step=0.2, smoothing=0.0)` -> `InterpolatedCurve`; `fill_gaps(u, step)` | [interpolation](interpolation.md) |
| `vth_constant_current(v_g, i_d, polarity, i_crit=1e-8, interpolation="linear")` | [threshold voltage](threshold_voltage.md) |
| `vth_linear_ss_window(curve, region, window=5.0, extension=3.0)` -> `LineFit` | " |
| `vth_best_linear_window(curve, window=5.0)` -> `LineFit` | " |
| `vth_sqrt_extrapolation(curve, window=5.0)` -> `LineFit` | " |
| `vth_max_gm_tangent(curve, gm)` -> `LineFit` | " |
| `vth_max_log_slope(curve, region)`, `vth_second_derivative(curve, window=None)`, `second_derivative(curve, window=None)` | " |
| `elr_correction(vth_intercept, v_d)` | " |
| `subthreshold_region(curve, sensitivity=1.4)` -> `SubthresholdRegion` | [subthreshold swing](subthreshold_swing.md) |
| `subthreshold_swing(curve, region, floor=60.0)`, `ss_decade_span(...)`, `ss_current_window(...)` -> `SwingResult`; `swing_curve(curve)` | " |
| `transconductance(curve, window=1.6)`, `gm_at(curve, gm, v_g)` | [transconductance](transconductance_and_on_off.md) |
| `on_off_ratio(i_d, off="min")` -> `OnOff`; `robust_off_level(i_d)` | " |
| `linear_mobility`, `effective_mobility`, `saturation_mobility`, `saturation_mobility_from_slope`, `fe_mobility_from_sqrt_fit` | [mobility](mobility.md) |

`LineFit(vth, slope, intercept, v_range, r, n_points)` is callable: `fit(v_g)` evaluates the line.

## Hysteresis and traps — `transistor_iv.hysteresis`, `transistor_iv.traps` ([page](hysteresis_and_traps.md))

| name | description |
|---|---|
| `analyze_hysteresis(curve, device, settings, segments=None, noise_floor=1e-12)` | -> `HysteresisResult(name, forward, backward, dvth, n_trap, loop_area, noise_floor, notes)`, `.summary()` |
| `split_dual_sweep(v_g, segments=None)` | forward and backward index arrays |
| `loop_area(v_fwd, i_fwd, v_bwd, i_bwd, noise_floor=1e-12)` | dec V |
| `d_it_from_ss(ss, c_i, temperature=300)`, `thermal_ss_limit(temperature)`, `n_trap_from_dvth(dvth, c_i)` | formulas |

## Output curves — `transistor_iv.output` ([page](output_curves.md))

`symmetry_score(v_d, i_d)`, `family_symmetry(family)`, `on_resistance(v_d, i_d, v_d_max)`,
`summarize_output(family, v_d_max=None)`.

## Batches — `transistor_iv.batch` ([page](batch_processing.md))

`analyze_curves(curves, devices, settings)`, `summary_table(curves, results)`,
`write_tables(curves, results, directory)`.

## Diagnostics — `transistor_iv.diagnostics`

Enabled with `TransferSettings(diagnostics=True)`; adds to `TransferResult.diagnostics` (and
the summary row): `flag_linear_regime`, `flag_saturation`, `flag_window_at_end` (empty
string = no problem) and `vth_interpol_sens`, `vth_global_sens`, `vth_sqrt_sens` (V),
`mu_sat_sqrt_sens`, `mu_fe_max_sens` (relative): the largest change when `linear_window` and
`gm_window` are halved and doubled.

## Figures — `transistor_iv.plot`

`plot_transfer(curve)`, `plot_analysis(result)`, `plot_overlay(curves, log=True)`,
`plot_output(family, title=None)`, `plot_hysteresis(result, vth_key="vth_current")`; each
returns a matplotlib figure. Import with `from transistor_iv import plot`.

## Synthetic curves — `transistor_iv.synthetic`

`model_current(v_g, v_d, device, vth, ss, mu, theta=0, r_s=0)`,
`transfer_sweep(v_g, v_d, device, vth, ss, mu, i_leak=0, noise_abs=0, noise_rel=0, seed=0, ...)`,
`hysteresis_sweep(v_start, v_stop, step, v_d, device, vth_fwd, dvth, ss, mu, ...)`,
`output_family(v_d, v_g_steps, device, vth, ss, mu, ...)`. The model: channel charge
C_i S ln(1 + exp((V_g - V_th)/S)) with S = SS/ln 10, integrated along the channel
(gradual-channel approximation); exponential below threshold with swing SS, the square law
above it; optional mobility degradation `theta` (1/V) and series resistance `r_s` (ohm).

## Exceptions — `transistor_iv.errors`

`CurveError` (the curve cannot be analysed), `ExtractionError` (one quantity cannot be
extracted; caught by `analyze_transfer`, which stores NaN and a note).
