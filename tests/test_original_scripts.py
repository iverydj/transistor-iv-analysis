"""Agreement with the author's original analysis scripts on the bundled example files.

``data/original_scripts_baseline.csv``: results of the original transfer-curve
script, run function by function on every example curve (n-type as measured,
and the exact p-type mirror V_g -> -V_g, I_d -> -I_d, V_d -> -V_d).
``data/original_scripts_fixed_baseline.csv``: the same script with the bug
fixes B1-B7 and B9 of docs/changes_from_original.md patched in (B8, the g_m
smoothing, is compared through gm_window=None). Both are written by
tests/tools/make_original_baselines.py, which needs the original script.
Geometry: W = L = 5 um on 300 nm SiO2, TC_example1 from Conditions.xls.
TC_example5 is analysed for sweep segments 1 (forward) and 2 (backward).
The original rounds threshold voltages and SS to 3 decimals.
"""

import numpy as np
import pandas as pd
import pytest

from conftest import BASELINES, DATA, read_quiet
from transistor_iv import Device, TransferCurve, TransferSettings, analyze_transfer, read_conditions

ORIG = pd.read_csv(BASELINES / "original_scripts_baseline.csv")
FIXED = pd.read_csv(BASELINES / "original_scripts_fixed_baseline.csv")
VTH = ["vth_current", "vth_interpol", "vth_global", "vth_logderivative", "vth_sqrt", "vth_d2"]
REL = ["ss", "on_off", "mu_fe_max", "mu_sat_sqrt"]


def ours(fname, segment, polarity="n", **kw):
    curves = read_quiet(fname)
    if fname == "TC_example1.xls":
        devs = read_conditions(DATA / "Conditions.xls", eps_r=3.9)
    else:
        devs = [Device.from_dielectric(5, 5, 3.9, 300, "n")] * len(curves)
    st = TransferSettings(segment=segment, **kw)
    out = {}
    for c, d in zip(curves, devs):
        if polarity == "p":
            c, d = TransferCurve(c.name, -c.v_g, -c.i_d, -c.v_d), Device(d.width_um, d.length_um, d.c_i, "p")
        out[c.name] = analyze_transfer(c, d, st)
    return out


JOBS = [("TC_example1.xls", -1), ("TC_example2.xls", -1), ("TC_example5.xls", 1), ("TC_example5.xls", 2)]


@pytest.fixture(scope="module")
def unsmoothed():
    return {job: ours(*job, gm_window=None) for job in JOBS}


def rows(table, fname, seg, pol="n"):
    return table[(table.file == fname) & (table.segment == seg) & (table.polarity == pol)].iterrows()


@pytest.fixture(scope="module")
def unsmoothed_p():
    return {job: ours(*job, polarity="p", gm_window=None) for job in JOBS}


@pytest.mark.parametrize("job", JOBS)
@pytest.mark.parametrize("pol", ["n", "p"])
def test_reproduces_the_fixed_original_exactly(job, pol, unsmoothed, unsmoothed_p):
    """Every quantity of every curve equals the bug-fixed original (g_m unsmoothed as there)."""
    res = (unsmoothed if pol == "n" else unsmoothed_p)[job]
    for _, b in rows(FIXED, *job, pol):
        r = res[b["name"]]
        assert r.v_d == b.v_d
        for k in VTH:
            assert getattr(r, k) == pytest.approx(b[k], abs=6e-4), (b["name"], k)
        for k in REL:
            assert getattr(r, k) == pytest.approx(b[k], rel=1e-5), (b["name"], k)


def test_unchanged_quantities_equal_the_unmodified_original(unsmoothed):
    """On TC_example1 (uniform 0.25 V steps, n-type, V_d consistent) only the window-based
    methods (vth_interpol, vth_global, vth_sqrt, mu_sat_sqrt) are affected by a fix."""
    res = unsmoothed[("TC_example1.xls", -1)]
    for _, b in rows(ORIG, "TC_example1.xls", -1):
        r = res[b["name"]]
        for k in ("vth_current", "vth_logderivative", "vth_d2"):
            assert getattr(r, k) == pytest.approx(b[k], abs=6e-4), (b["name"], k)
        for k in ("ss", "on_off", "mu_fe_max"):
            assert getattr(r, k) == pytest.approx(b[k], rel=1e-5), (b["name"], k)


def test_original_bugs_are_visible_in_the_baseline():
    o = ORIG.set_index(["file", "segment", "name", "polarity"])
    n = o.xs("n", level="polarity")
    p = o.xs("p", level="polarity")
    # mobility prefactor max(V_d, 1e-6): p-type mobility inflated by |V_d| / 1e-6
    ratio = (p.mu_fe_max / n.mu_fe_max).to_numpy()
    assert np.all(ratio > 1e4)
    # np.interp on a dual sweep: NaN constant-current V_th on every TC_example5 curve
    assert n.loc["TC_example5.xls"].vth_current.isna().all()
    # near-duplicate grid points from single-precision V_g: round-off spikes decide the d2 peak on TC_example5
    f = FIXED.set_index(["file", "segment", "name", "polarity"]).xs("n", level="polarity")
    shift = (n.loc["TC_example5.xls"].vth_d2 - f.loc["TC_example5.xls"].vth_d2).abs()
    assert shift.median() > 30


@pytest.mark.parametrize("job", JOBS)
def test_default_smoothing_changes_only_the_gm_quantities(job, unsmoothed):
    """gm_window affects gm_max, v_gm_max, mu_fe_max and vth_gm_tangent, nothing else."""
    smooth = ours(*job)
    for name, r in unsmoothed[job].items():
        s = smooth[name]
        for k in VTH + ["ss", "on_off", "mu_sat_sqrt", "d_it"]:
            assert getattr(s, k) == pytest.approx(getattr(r, k), abs=0, nan_ok=True), (name, k)
