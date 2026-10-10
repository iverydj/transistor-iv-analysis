# Changes from the original scripts

This package re-implements the author's two analysis scripts for Keithley Clarius exports:
the transfer-curve script (`TransferCurve_Clarius-xls.py`) and the output-curve script
(`OutputCurve_plot_Clarius-xls.py`). Their definitions, method set and default selection
rules are kept. This page lists **every behavioural difference**, classified as

* **bug fix**: the original did not do what its own definition or documentation says;
* **option**: a different behaviour is available, the original one is the default;
* **design change**: the default behaviour differs for a reason other than a bug;
* **addition**: functionality the original scripts did not have.

"Effect" quotes the change on the bundled example files (n-type as measured, W = L = 5 um on
300 nm SiO2, `TC_example1.xls` with its `Conditions.xls`, `TC_example5.xls` analysed for
segments 1 and 2), obtained by applying each fix alone to the original script. The tests in
`tests/test_original_scripts.py` check that the package reproduces the original script with
all bug fixes applied on every example curve (threshold voltages within the original's
rounding to 1 mV, everything else within 1e-5 relative); see [testing](testing.md).

## Transfer curves: bug fixes

| # | original behaviour | new behaviour | justification | effect on the examples |
|---|---|---|---|---|
| B1 | linear-mobility prefactor L / (W C_i max(V_d, 1e-6)) | L / (W C_i \|V_d\|) | for V_d < 0 (p-type) `max` returned 1e-6 | p-type mu_lin was \|V_d\|/1e-6 too large (1e5 x at 0.1 V, 1e7 x at 10 V); n-type unchanged |
| B2 | constant-current V_th = `np.interp(I_crit, |I_d|, V_g)` over the raw data of the whole sheet | crossing on the analysed segment: last point below I_crit and the next one, linear interpolation in I_d | `np.interp` needs increasing x; on a dual sweep it returned NaN, for p-type wrong values | TC_example1/2 identical (within the original's 1 mV rounding); TC_example5: NaN -> values on all 88 segments |
| B3 | d²I_d/dV_g² negated for p-type before taking the maximum | not negated | d²\|I_d\|/dV_g² keeps its sign under V_g -> -V_g; the negation selected the minimum | p-type `vth_d2` |
| B4 | p-type: the `vth_interpol` window search ran from the subthreshold region towards higher V_g, i.e. into the off state | every method works in the on-direction frame u = s V_g, so p-type curves are analysed as the mirror of n-type ones | the search direction must point to the on-state for both polarities | p-type mirrors of the examples: `vth_interpol` differed from the mirrored n-type value by a median 2.2 V (up to 12 V); now exactly mirrored |
| B5 | the gap filling inserted a point whenever a gap exceeded 0.2 V, also by float round-off | gaps exceeding a multiple of the step by less than 1 % of a step are not split | single-precision V_g (0.2000000x V) created points 3e-9 to 1e-6 V before the next one, whose round-off dominated the derivatives | TC_example5: 182 of 533 grid points affected; `vth_d2` moved by a median 36 V (it sat on round-off spikes), `vth_interpol` by 0.9 to 1.6 V (median), `vth_global` 0.6 to 0.7 V, SS by up to 85 % on single curves, mu_fe_max by a median 3 to 4 % (up to 45 %); TC_example2 (26 points): small changes on several curves, up to 33 V in `vth_sqrt` and 39 V in `vth_d2` on one curve each |
| B6 | fit windows of `LINEAR_WINDOW_FOR_VTH` volts (and the 3 V search extension) converted to points with the nominal 0.2 V step | converted with the mean spacing of the actual grid | the parameter is documented in volts; with 0.25 V data (grid spacing alternating 0.2/0.05 V) the "5 V" window covered 3.1 V, on the 0.1 to 0.5 V list sweep 3 V | TC_example1: `vth_interpol` median 0.38 V (up to 8.3 V), `vth_global` 0.45 V, `vth_sqrt` 0.42 V (up to 68 V on one curve), `mu_sat_sqrt` 13 %; TC_example2: 0.04 to 0.16 V (median) |
| B7 | SS region slice `[start:end]` dropped the last point of the region | `[start:end+1]` | off-by-one (the figure of the region used the inclusive end); a one-point region gave an empty array | none |
| B8 | mu_FE,max = maximum of mu_lin from an unsmoothed `np.gradient` of the exactly interpolated current | g_m from local quadratic fits over `gm_window` = 1.6 V (Savitzky-Golay, 9 points on a 0.2 V grid); `gm_window=None` restores the unsmoothed derivative | the maximum of a noisy derivative selects the largest noise excursion | mu_fe_max: TC_example1 a median 2.7 x lower (up to 5.1 x, spiky on-currents), TC_example2 1.3 %, TC_example5 3 % |
| B9 | V_d from the Settings block whenever the drain `Step` is `N/A`, the DrainV column only as fallback | the DrainV column when present (must be constant), with a warning if the Settings block disagrees | the column is the applied bias | one sheet of TC_example2 (Settings 20 V, column and sheet name 10 V): mu_fe_max x 2 |
| B10 | `results.csv` written with `np.savetxt(delimiter=",")` | written with pandas (CSV quoting) | sheet names with commas (`ND2515(2,6)Vd10 2nd`) shifted the columns | TC_example2 table readable |

## Output curves: bug fix

| # | original behaviour | new behaviour | justification | effect on the examples |
|---|---|---|---|---|
| B11 | symmetry score paired the i-th negative-bias point with the i-th positive-bias point from the end | pairs at mirrored voltage, y(-V_d) interpolated linearly | index pairing is correct only for symmetric, uniform, ordered grids | identical on the bundled files (symmetric 0.25 V grids, max difference 2e-14 %) |

## D_it

| # | | | |
|---|---|---|---|
| B12 | D_it is not in the original scripts; a formula specified earlier for it was not the standard one | D_it = (C_i/q)[SS / (ln(10) k_B T/q) - 1] (cm⁻² eV⁻¹), NaN below the thermionic limit | from SS = ln(10)(k_B T/q)(1 + C_it/C_i), C_it = q D_it | see [hysteresis_and_traps.md](hysteresis_and_traps.md) |

## Options (original behaviour is the default)

| # | option | default (original) | alternative | justification |
|---|---|---|---|---|
| O1 | `vd_correction` | report the extrapolated intercept for `vth_interpol`, `vth_global` (and `vth_gm_tangent`, `vth_d2`) | intercept - V_d/2 | textbook linear-regime threshold, I_d proportional to (V_g - V_th - V_d/2) V_d; see [threshold voltage](threshold_voltage.md#the-v_d2-correction-option) |
| O2 | `ss_floor` | 60 mV/dec | `None` (no floor) or another value | the floor hides genuine sub-60 swings (low temperature, steep-slope devices) |
| O3 | `ss_method` | `"pointwise"`: smallest point-wise SS >= floor in the region | `"decade_span"`: shortest span of one decade of the current envelope; `"current_window"`: line fit in a stated current window | less sensitive to single noisy points; stated windows are reproducible |
| O4 | `off_level` | `"min"`: smallest positive \|I_d\| | `"robust"`: 10^median of the lowest 10 % of the median-filtered log current | a single noise dip sets the minimum |
| O5 | `cc_interpolation` | `"linear"` in I_d | `"log"` in log10\|I_d\| | the current is exponential near a subthreshold criterion |
| O6 | `d2_window` | `None`: `np.gradient` twice | local quadratic fits over the window | noise; see [threshold voltage](threshold_voltage.md#second-derivative-vth_d2) |
| O7 | `diagnostics` | off | regime flags (`flag_linear_regime`, `flag_saturation`, `flag_window_at_end`) and window sensitivities (`*_sens`: change when the windows are halved and doubled) | checks the validity conditions of each formula |
| O8 | `gm_window=None` | (the fix B8 is the default) | restores the original unsmoothed g_m | comparison with old results |

## Design changes

| # | original | now | justification |
|---|---|---|---|
| C1 | `SEGMENT_INDEX_FOR_SPLINE = 2` in the script (the README documents -1) | `segment = -1` (longest segment; the first of equally long ones) | the documented default; 2 was a per-file setting |
| C2 | a segment number that does not exist: warning and longest segment; a segment of fewer than 3 points: the whole sheet | `CurveError` | fail loudly instead of analysing different data than requested |
| C3 | every exception of a dataset caught, the dataset listed as an error | `CurveError` (data problems) recorded per curve in a batch; any other exception stops the run | programming and settings errors must not pass silently |
| C4 | without a subthreshold region, `vth_current`, `vth_d2`, on/off and the mobilities were skipped | only the quantities that need the region (SS, `vth_interpol`, `vth_logderivative`, mu_eff, mu_sat point-wise) are NaN, with a note | the others do not depend on it |
| C5 | V_th and SS rounded to 3 decimals | full precision | rounding belongs to reporting |
| C6 | no positive current at all: I_off = 1e-12 A | `ExtractionError` (NaN with a note) | no invented default values |
| C7 | gate-insulator thickness in um (`GI_THICK`) | `gate_capacitance(eps_r, thickness_nm)` in nm; `Conditions.xls` still in um | nm is the usual unit; the unit is in the argument name |
| C8 | `Conditions.xls` without a type column: `DEVICE_TYPE` used | `read_conditions(..., polarity=...)` must be given exactly when there is no type column | no ambiguity about which polarity applies |
| C9 | the reader assumed the Clarius layout; empty padding cells would fail the float conversion | layout checked (sheets 2-3 `Calc`, `Settings`, a Settings block per sheet), padding removed, an empty cell inside the data is an error | robust reading, explicit failure |
| C10 | header candidates of the scripts | the same plus `I_d`, `V_g`, `V_d`, `I_g` (for CSV files) | CSV input |
| C11 | symmetry score with an exact zero current: NaN or inf propagated | `ValueError` | fail loudly |
| C12 | figures: transfer curve on twin axes (log and linear), a 4 x 3 analysis grid, overlays; result folders with timestamps and a copy of the script | `plot_transfer` (log and linear side by side), `plot_analysis` (2 x 4), `plot_overlay`, `plot_output`, `plot_hysteresis`; the caller chooses the output folder | library use; one y-axis per panel |
| C13 | options `DENOISE_CURRENT` (Gaussian filter and resampling), `REMOVE_OUTLIERS` (g_m outlier removal), `ANALYZE_*` switches, plot settings | not carried over | superseded by `gm_window` and `spline_smoothing`; functions are called individually; matplotlib settings belong to the caller |
| C14 | table names `Vth_current`, ..., `u_FE(max)`; files `TransferCurve.csv`, `SubthresholdSwing.csv`, `MobilityFE.csv` | `vth_current`, ..., `mu_fe_max`; `transfer_curves.csv`, `subthreshold_swing.csv`, `mobility_fe.csv` | Python identifiers; mapping in [batch processing](batch_processing.md#columns-of-resultscsv) |

## Additions

`vth_gm_tangent` (tangent at maximum g_m); D_it; hysteresis analysis (`analyze_hysteresis`:
threshold shift for every method, N_trap, loop area); `fe_mobility_from_sqrt_fit` (mu_FE at a
stated V_g from the square-root fit) and `gm_at`; `ss_decade_span`, `ss_current_window`,
`robust_off_level`; `on_resistance` and `summarize_output`; CSV readers and writer; the
synthetic compact model `transistor_iv.synthetic`.

## For the author's decision

These defaults follow the original scripts but are physically questionable or are choices
made in this rebuild. Each is documented as a pitfall on its page.

| # | item | issue | current default | possible change |
|---|---|---|---|---|
| D1 | `vth_logderivative` | the steepest subthreshold point is not a threshold voltage (it lies below the strong-inversion threshold: 1.6 V for a model V_th of 2 V) | reported as a V_th | rename (e.g. "turn-on voltage") or drop from the V_th set |
| D2 | `ss_floor = 60` mV/dec | hides real sub-60 swings and can report a value far above the true one (45 mV/dec model device: 118.9 reported) | 60 | `None`, or the thermal limit at `temperature` |
| D3 | linear extrapolations report the intercept | in the linear regime the intercept is V_th + V_d/2; with V_d = 5 to 10 V (TC_example1/2) that is 2.5 to 5 V | intercept (`vd_correction=False`) | `vd_correction=True` |
| D4 | regime validity | `vth_sqrt` and `mu_sat_sqrt` are computed in the linear regime too (e.g. -37 V on a TC_example2 curve); `mu_fe_max` applies the linear formula at V_d = 5 to 10 V; the point-wise mu_eff and mu_sat use the linear-regime `vth_interpol` | computed without checks | `diagnostics=True` by default, or NaN outside the regime |
| D5 | I_off = smallest positive \|I_d\| | the deepest noise sample, not the off level (two orders of magnitude in the [on/off example](transconductance_and_on_off.md#onoff-ratio)) | `"min"` | `"robust"` |
| D6 | unsmoothed second derivative | on noisy data the maximum is a noise peak in the on-state (bundled TC_example2 first sheet: a spike at 18.4 V) | `d2_window=None` | `d2_window` of about 1 V |
| D7 | B6 (window in volts via the mean grid spacing) | classified as a bug fix; on list sweeps with changing steps a window is still a fixed number of grid points | mean spacing | please confirm the classification |
| D8 | `gm_window = 1.6` V | the value of the fix B8 was chosen here (the 9-point Savitzky-Golay filter on 0.2 V steps); wider windows remove more noise but flatten peaked g_m (-4 % at 1.6 V, -8 % at 3 V for a mobility-degradation model) | 1.6 V | please confirm |
