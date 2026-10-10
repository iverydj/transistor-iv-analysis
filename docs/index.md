# transistor-iv-analysis documentation

`transistor_iv` extracts the standard figures of merit of field-effect transistors from
measured transfer (I_d vs V_g) and output (I_d vs V_d) curves. Each page below documents one
quantity or one step: its definition with units, the algorithm, the parameters, when it is
valid and what can go wrong. Figures on these pages explain the methods on synthetic or
bundled example curves; results on the example files are collected separately in
[examples](examples.md).

## Data in

* [Input formats](input_formats.md): Keithley 4200A-SCS Clarius `.xls`, `Conditions.xls`
  geometry tables, CSV files, arrays.
* [Polarity, units and sweeps](polarity_and_units.md): n/p convention, units, sweep
  segments, errors and missing values.
* [The interpolated curve](interpolation.md): the spline and grid that most transfer
  methods share.

## Transfer-curve quantities

* [Threshold voltage](threshold_voltage.md): seven methods (constant current, three linear
  extrapolations, maximum log slope, second derivative, square-root extrapolation) and the
  V_d/2 correction.
* [Subthreshold swing](subthreshold_swing.md): the subthreshold region and three selection
  rules.
* [Transconductance and on/off ratio](transconductance_and_on_off.md): g_m and its
  smoothing, I_on, I_off.
* [Mobility](mobility.md): linear, effective and saturation mobility, mu_FE from a
  square-root fit.
* [Hysteresis and trap densities](hysteresis_and_traps.md): forward/backward analysis,
  threshold shift, N_trap, loop area, D_it.

## Output curves and many curves

* [Output curves](output_curves.md): figures, symmetry score, low-bias resistance.
* [Batch processing](batch_processing.md): tables for many curves and files.

## Reference

* [Worked examples](examples.md): the example scripts and what their outputs look like.
* [API reference](api.md).
* [Changes from the original scripts](changes_from_original.md).
* [Tests and validation](testing.md).

## Typical workflow

```python
import transistor_iv as tiv
from transistor_iv import plot

dev = tiv.Device.from_dielectric(width_um=100, length_um=10, eps_r=3.9, thickness_nm=100, polarity="n")
curve = tiv.read_clarius_transfer("my_measurement.xls")[0]          # or read_transfer_csv(...)
res = tiv.analyze_transfer(curve, dev, tiv.TransferSettings())       # defaults = original scripts
print(res.vth_global, res.ss, res.mu_fe_max, res.notes)
plot.plot_analysis(res).savefig("analysis.png")
```

All results come with the information needed to judge them: the fit lines
(`res.fits`), the subthreshold region and the SS position (`res.region`, `res.swing`),
the point-wise curves (`res.curves`) and the reasons for missing values (`res.notes`).
