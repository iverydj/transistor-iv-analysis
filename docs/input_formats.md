# Input formats

Every analysis works on two containers:

* `TransferCurve(name, v_g, i_d, v_d=nan, i_g=None)`: one transfer measurement. `v_g` may
  hold several sweep segments; `v_d` is the signed, constant drain bias with the source
  grounded (NaN if unknown); `i_g` the gate current (optional).
* `OutputFamily(name, v_g, v_d, i_d)`: output curves, one per gate voltage.

Any data source can be used by constructing them from arrays:

```python
import numpy as np
import transistor_iv as tiv
curve = tiv.TransferCurve("my device", v_g=np.array([...]), i_d=np.array([...]), v_d=0.1)
```

Units: V and A.

## Keithley 4200A-SCS Clarius `.xls`

`read_clarius_transfer(path, headers=None)` returns one `TransferCurve` per data sheet;
`read_clarius_output(path, headers=None)` one `OutputFamily` per sheet;
`read_clarius_sheets(path)` the raw sheets with their settings.

**Workbook layout** (checked; anything else raises `ValueError`): sheet 1 is the first data
sheet, sheets 2 and 3 are `Calc` and `Settings`, the remaining sheets are data sheets. Each
data sheet has a header row and one column per quantity. The `Settings` sheet has one block
per data sheet, introduced by a line of 34 `=` characters, the sheet name and another `=`
line; every data sheet must have its block.

**Columns** are found by header, compared after lower-casing and removing everything that
is not a letter or digit (so `Drain I (A)` matches `DrainI`). In this order: the header you
give in `headers` (keys `"i_d"`, `"v_g"`, `"v_d"`, `"i_g"`), the default candidates, then a
prefix match of either:

| quantity | default candidates | required |
|---|---|---|
| drain current | DrainI, DrainI_A, Id, I_d, Drain I, Drain I (A), Id(A) | yes |
| gate voltage | GateV, Vg, V_g, Gate V, Gate V (V), Vg(V) | yes |
| drain voltage | DrainV, Vd, V_d, Drain V, Drain V (V) | no |
| gate current | GateI, Ig, I_g, Gate I, Gate I (A), GateI_A | no |

Clarius numbers the columns of tests that store several curves per sheet: `DrainI(1)`,
`DrainI(2)`, ... For transfer sheets the prefix match reads the first of them; output-curve
sheets are read completely (`read_clarius_output`), one curve per number. Empty cells that
pad shorter columns are removed; an empty cell inside the data raises `ValueError`.

**Drain bias of a transfer sheet.**

* If the sheet has a drain-voltage column, its median is used; it must be constant within
  1 mV (otherwise the sheet is not a transfer curve, `ValueError`).
* Otherwise the `Settings` block is used: if the drain terminal's `Step` entry is `N/A`
  (no drain sweep), its `Start/Bias` entry.
* If both exist and differ by more than 1 mV, the column value is used and a `UserWarning`
  names the sheet (one bundled sheet's Settings block says 20 V while its data column and
  its name say 10 V).
* If neither exists, `v_d` is NaN: the mobilities and the V_d/2 correction are then
  unavailable.

**Per-device geometry: `Conditions.xls`.** `read_conditions(path, eps_r, polarity=None)`
returns one `Device` per row of a sheet named `Conditions`, in data-sheet order. Columns
(matched ignoring case, spaces and underscores): `ChannelWidth` (um), `ChannelLength` (um),
`GiThickness` (gate-insulator thickness in **um**), optional `DeviceType` (`n` / `p`). The
table has no permittivity column: pass `eps_r`. Pass `polarity` exactly when the table has
no `DeviceType` column.

## CSV files

**Transfer curve** (`read_transfer_csv(path, headers=None, v_d=None)`, name = file stem):
one row per point, columns for V_g and I_d, optionally V_d (constant) and I_g, with the same
header matching as above:

```text
V_g,I_d,V_d,I_g
-20.0,1.2e-12,0.1,3e-13
-19.9,1.3e-12,0.1,2e-13
...
```

A file without a V_d column takes the bias from the `v_d` argument. A dual sweep is just the
full point sequence; the segments are found from V_g. `write_transfer_csv(curve, path)`
writes this format with full precision; `transfer_from_dataframe(df, name)` reads from a
pandas DataFrame.

**Output family** (`read_output_csv(path)`): columns V_g, V_d, I_d, one row per point; rows
with the same V_g form one curve (in order of first appearance):

```text
V_g,V_d,I_d
10,0.0,1e-13
10,0.25,2.1e-7
...
```

## Bundled example files

`examples/data/` holds Clarius exports re-saved cell for cell without workbook metadata:

| file | content |
|---|---|
| `TC_example1.xls` + `Conditions.xls` | 19 oxide-TFT transfer curves, -20 to 30 V in 0.25 V steps, V_d = 5 V; geometry per sheet (W/L = 800/200 um, 100 nm insulator) |
| `TC_example2.xls` | 10 oxide-TFT transfer curves, list sweep -40 to 50 V (0.1 to 0.5 V steps), V_d = 10 V, with gate current |
| `TC_example5.xls` | 44 MoS2 dual sweeps 0 -> -30 -> 40 -> -30 V (0.2 V steps), V_d = 0.1 V from the Settings block |
| `OP_ex1.xls`, `OP_ex2.xls` | output families, V_d from -5 to 5 V at nine gate voltages |

The dielectric constant and, except for `TC_example1.xls`, the geometry are not stored; the
examples assume values, so their mobilities are illustrations.
