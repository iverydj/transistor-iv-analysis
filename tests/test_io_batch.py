import warnings

import numpy as np
import pandas as pd
import pytest

from conftest import DATA, read_quiet
from transistor_iv import (Device, TransferCurve, analyze_curves, monotonic_segments, read_clarius_output,
                           read_conditions, read_output_csv, read_transfer_csv, summary_table, write_tables,
                           write_transfer_csv)
from transistor_iv.io import read_clarius_sheets
from transistor_iv.synthetic import transfer_sweep


def test_clarius_transfer_files():
    ex1, ex2, ex5 = (read_quiet(f) for f in ("TC_example1.xls", "TC_example2.xls", "TC_example5.xls"))
    assert (len(ex1), len(ex2), len(ex5)) == (19, 10, 44)
    assert {c.v_d for c in ex1} == {5.0}
    assert ex2[0].name == "ND2515(2,6)Vd10 2nd" and ex2[0].v_g.size == 566 and ex2[0].i_g.size == 566
    assert {c.v_d for c in ex5} == {0.1}                         # no DrainV column: from the Settings block
    v = ex5[0].v_g
    assert [(v[s[0]], v[s[-1]]) for s in monotonic_segments(v)] == [(0, -30), (-30, 40), (40, -30)]


def test_settings_and_data_column_disagree():
    with pytest.warns(UserWarning, match="ND0040.*using the column value"):
        curves = {c.name: c for c in __import__("transistor_iv").read_clarius_transfer(DATA / "TC_example2.xls")}
    assert curves["ND0040(2,5)Vd10 3rd"].v_d == 10.0             # the Settings block says 20 V


def test_settings_blocks_and_layout_check():
    sheets = read_clarius_sheets(DATA / "TC_example5.xls")
    assert sheets[0].settings_drain_bias() == 0.1
    with pytest.raises(ValueError, match="Clarius"):
        read_clarius_sheets(DATA / "Conditions.xls")


def test_conditions_and_output_readers():
    devs = read_conditions(DATA / "Conditions.xls", eps_r=3.9)
    assert len(devs) == 19
    assert devs[0] == Device.from_dielectric(800, 200, 3.9, 100, "n")       # thickness column in um
    with pytest.raises(ValueError, match="polarity"):
        read_conditions(DATA / "Conditions.xls", eps_r=3.9, polarity="n")
    fam = read_clarius_output(DATA / "OP_ex1.xls")
    assert [f.name for f in fam] == ["w4_s102_r1", "w5_s103", "w5_s101"]
    np.testing.assert_allclose(fam[0].v_g, [-40, -30, -20, -10, 0, 10, 20, 30, 40])
    assert fam[0].v_d[0][0] == -5 and fam[0].v_d[0][-1] == 5


def test_csv_round_trip(tmp_path):
    c = read_quiet("TC_example2.xls")[3]
    write_transfer_csv(c, tmp_path / "curve.csv")
    back = read_transfer_csv(tmp_path / "curve.csv")
    assert back.name == "curve" and back.v_d == c.v_d
    np.testing.assert_array_equal(back.v_g, c.v_g)
    np.testing.assert_array_equal(back.i_g, c.i_g)
    pd.DataFrame({"Vg": [0, 1, 2, 3], "Id": [1e-12, 1e-11, 1e-10, 1e-9]}).to_csv(tmp_path / "x.csv", index=False)
    assert read_transfer_csv(tmp_path / "x.csv", v_d=0.1).v_d == 0.1
    with pytest.raises(ValueError, match="V_d column"):
        read_transfer_csv(tmp_path / "curve.csv", v_d=0.1)
    pd.DataFrame({"V_g": [1, 1, 2, 2], "V_d": [0, 1, 0, 1], "I_d": [0, 1, 0, 2]}).to_csv(tmp_path / "o.csv", index=False)
    fam = read_output_csv(tmp_path / "o.csv")
    assert list(fam.v_g) == [1, 2] and list(fam.i_d[1]) == [0, 2]


def test_batch_tables_are_valid_csv(tmp_path):
    dev = Device.from_dielectric(100, 10, 3.9, 100, "n")
    good = transfer_sweep(np.arange(-10, 20.01, 0.2), 0.1, dev, 1.0, 150, 10)
    good.name = "dev(2,5), first"                               # comma in the name
    bad = TransferCurve("flat", np.array([0.0, 0.0, 0.0]), np.ones(3), 0.1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = analyze_curves([good, bad], dev)
    tab = summary_table([good, bad], res)
    assert tab.error.iloc[0] == "" and "constant" in tab.error.iloc[1]
    paths = write_tables([good, bad], res, tmp_path)
    back = pd.read_csv(paths[0], keep_default_na=False)
    assert list(back.name) == ["dev(2,5), first", "flat"]
    assert float(back.vth_global.iloc[0]) == pytest.approx(res[0].vth_global)
    assert pd.read_csv(paths[1]).shape[1] == 4


def test_csv_round_trip_keeps_the_current(tmp_path):
    c = read_quiet("TC_example2.xls")[3]
    write_transfer_csv(c, tmp_path / "c.csv")
    np.testing.assert_array_equal(read_transfer_csv(tmp_path / "c.csv").i_d, c.i_d)
    df = pd.read_csv(tmp_path / "c.csv", float_precision="round_trip")
    df["comment"] = "x"                                          # an extra text column is ignored
    df.to_csv(tmp_path / "d.csv", index=False, float_format="%.17g")
    np.testing.assert_array_equal(read_transfer_csv(tmp_path / "d.csv").v_g, c.v_g)


def test_list_sweep_settings_give_no_drain_bias():
    from transistor_iv.io import ClariusSheet
    rows = [["Device Terminal", "Drain", "Gate"], ["Step", "N/A", "N/A"], ["Start/Bias", "N/A", "N/A"]]
    sheet = ClariusSheet("s", ["DrainI"], np.zeros((0, 1)), rows)
    assert np.isnan(sheet.settings_drain_bias())
