"""Regenerate tests/data/original_scripts_*baseline.csv from the original transfer-curve script.

The original script (TransferCurve_Clarius-xls.py) is not part of this repository; set
AppConfig.original_script to a copy of it and run this file. The script is executed with its
final ``main()`` call removed, and its analysis functions are applied to every bundled example
curve, as measured (n-type) and as the exact p-type mirror (V_g, I_d, V_d -> -V_g, -I_d, -V_d).
Geometry: W = L = 5 um, 300 nm, eps_r 3.9; TC_example1 from Conditions.xls; TC_example5 for
sweep segments 1 and 2. Results go to AppConfig.out_dir (a scratch directory: the original
script creates a result/ folder in its working directory when it is loaded).

Variants:
  original : the script as written
  fixed    : with the bug fixes B1-B7 and B9 of docs/changes_from_original.md patched in
             (B8, the g_m smoothing, is not applied: the tests compare with gm_window=None;
             B10 and B11 do not affect these numbers)
"""

import os
import shutil
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


@dataclass
class AppConfig:
    original_script: Path = Path("TransferCurve_Clarius-xls.py")
    data_dir: Path = ROOT / "examples" / "data"
    out_dir: Path = Path("original_baselines_run")
    variants: list = field(default_factory=lambda: ["original", "fixed"])
    jobs: list = field(default_factory=lambda: [("TC_example1.xls", -1, True), ("TC_example2.xls", -1, False),
                                                 ("TC_example5.xls", 1, False), ("TC_example5.xls", 2, False)])


STEP_OLD = "step_base = ctx.Vg_step_interpolate if ctx.Vg_step_interpolate > 0 else max(np.diff(ctx.V_g).min(), 1e-3)"
STEP_OLD2 = "step_base = ctx.Vg_step_interpolate if ctx.Vg_step_interpolate > 0 else max(np.diff(V_g).min(), 1e-3)"
STEP_NEW = "step_base = (ctx.V_g_fine[-1] - ctx.V_g_fine[0]) / (len(ctx.V_g_fine) - 1)"

# (old text, new text, expected number of occurrences)
PATCHES = {
    "B6_window": [(STEP_OLD, STEP_NEW, 1), (STEP_OLD2, STEP_NEW, 2),
                  ("window_size = max(int(LINEAR_WINDOW_FOR_VTH / step_base), 3)",
                   "window_size = max(int(np.floor(LINEAR_WINDOW_FOR_VTH / step_base * (1 + 1e-6))), 3)", 3),
                  ("offset = int(3 / step_base)", "offset = int(np.floor(3 / step_base * (1 + 1e-6)))", 1)],
    "B5_gap": [("        if gap > step:\n            num_elements = int(np.ceil(gap / step)) - 1\n",
                "        num_elements = int(np.ceil(gap / step - 0.01)) - 1\n        if num_elements > 0:\n", 1)],
    "B7_region": [("ctx.V_g_cut = ctx.V_g_fine[ss_start:ss_end]", "ctx.V_g_cut = ctx.V_g_fine[ss_start:ss_end + 1]", 1),
                  ("ctx.I_d_cut = ctx.I_d_fine[ss_start:ss_end]", "ctx.I_d_cut = ctx.I_d_fine[ss_start:ss_end + 1]", 1),
                  ("ctx.dlog_cut = d_logId_dVg[ss_start:ss_end]", "ctx.dlog_cut = d_logId_dVg[ss_start:ss_end + 1]", 1),
                  ("SS_cut = SS_values[ctx.ss_start:ctx.ss_end]", "SS_cut = SS_values[ctx.ss_start:ctx.ss_end + 1]", 1)],
    "B2_cc": [("ctx.vth_current = float(np.round(np.interp(ID_AT_VTH, ctx.I_d_abs, ctx.V_g, left=np.nan, right=np.nan), 3))",
               "ctx.vth_current = _cc_fixed(ctx)", 1)],
    "B3_d2sign": [("        if ctx.device_type == 'p':\n            d2 = -d2\n", "", 1)],
    "B1_absvd": [("prefactor_linear = L / (W * gi_cap * max(ctx.V_d, 1e-6))",
                  "prefactor_linear = L / (W * gi_cap * abs(ctx.V_d))", 1)],
    "B4_interpol_p": [("        ctx.V_lin = ctx.V_g_fine[best_start:best_start + window_size]\n"
                       "        ctx.I_lin = ctx.I_d_fine[best_start:best_start + window_size]\n",
                       "        lo_i, hi_i = (_interpol_window_p(ctx, window_size, offset) if ctx.device_type == 'p'\n"
                       "                      else (best_start, best_start + window_size))\n"
                       "        ctx.V_lin = ctx.V_g_fine[lo_i:hi_i]\n"
                       "        ctx.I_lin = ctx.I_d_fine[lo_i:hi_i]\n", 1)],
}
FIXED = list(PATCHES) + ["B9_vdcol"]

HELPERS = '''
def _cc_fixed(ctx):
    """B2: last crossing on the analysed segment, linear in I_d."""
    seg = get_segment_indices(ctx.V_g, n=SEGMENT_INDEX_FOR_SPLINE)
    v = ctx.V_g[seg]; a = np.abs(ctx.I_d[seg])
    s = 1 if ctx.device_type == 'n' else -1
    if s * (v[-1] - v[0]) < 0:
        v, a = v[::-1], a[::-1]
    below = np.flatnonzero(a < ID_AT_VTH)
    if below.size == 0 or below[-1] == a.size - 1:
        return np.nan
    k = below[-1]
    return float(v[k] + (ID_AT_VTH - a[k]) / (a[k + 1] - a[k]) * (v[k + 1] - v[k]))


def _interpol_window_p(ctx, n, off):
    """B4: p-type windows extend from the subthreshold region towards lower V_g (the on-state)."""
    best_b, best = ctx.ss_end, -1.0
    for b in range(max(ctx.ss_start - off, 0), ctx.ss_end + 1):
        if b - n + 1 < 0:
            continue
        _, _, r_value, *_ = linregress(ctx.V_g_fine[b - n + 1:b + 1], ctx.I_d_fine[b - n + 1:b + 1])
        if abs(r_value) > best:
            best, best_b = abs(r_value), b
    return max(best_b - n + 1, 0), best_b + 1
'''


def patched(text: str, which: list[str]) -> str:
    for name in which:
        for old, new, count in PATCHES.get(name, []):
            assert text.count(old) == count, (name, old)
            text = text.replace(old, new)
    return text + HELPERS


def drain_bias(setting, column, use_column: bool) -> float:
    """The original rule (Settings first, DrainV column as fallback); B9: the column first."""
    v_d = np.nan
    try:
        step_row = np.where(setting[:, 0] == "Step")[0][0]
        terminal_row = np.where(setting[:, 0] == "Device Terminal")[0][0]
        drain_col = np.where(setting[terminal_row, :] == "Drain")[0][0]
        bias_row = np.where(setting[:, 0] == "Start/Bias")[0][0]
        if setting[step_row, drain_col] == "N/A":
            v_d = float(setting[bias_row, drain_col])
    except (IndexError, ValueError):
        pass
    has_column = column is not None and column.size and np.isfinite(column[0])
    if has_column and (use_column or np.isnan(v_d)):
        v_d = float(column[0])
    return v_d


def run_file(ns: dict, fname: str, segment: int, conditions: bool, use_column: bool) -> list[dict]:
    ns.update(DATA_FILE=fname, SEGMENT_INDEX_FOR_SPLINE=segment, CONDITION_DATA=conditions,
              CONDITION_FILE="Conditions.xls")
    names, settings, d_i, g_v, d_v, _, cond = ns["load_workbook_data"]()
    rows, errors = [], []
    for k, name in enumerate(names):
        i_d, v_g = np.array(d_i[k]).flatten(), np.array(g_v[k]).flatten()
        v_d = drain_bias(settings[k], d_v[k], use_column)
        if cond is not None:
            w, l, t = cond[k]["W_um"], cond[k]["L_um"], cond[k]["gi_thick_um"]
        else:
            w, l, t = ns["CH_WIDTH"], ns["CH_LENGTH"], ns["GI_THICK"]
        for pol, sign in (("n", 1), ("p", -1)):
            ctx = ns["AnalysisCtx"](name=name, V_g=sign * v_g, I_d=sign * i_d, V_d=sign * v_d, W_um=w, L_um=l,
                                    gi_thick_um=t, device_type=pol, error_files=errors, I_g_opt=None)
            ns["preprocess_analysis"](ctx)
            if ctx.has_ss_region:
                ns["compute_vth"](ctx)
                ns["compute_ss"](ctx)
                ns["compute_mobility"](ctx)
            if ctx.V_g_fine is not None and ctx.V_g_fine.size >= 4:
                ns["compute_vth_global"](ctx)
                ns["compute_sqrt_sat"](ctx)
            rows.append(dict(file=fname, segment=segment, name=name, polarity=pol, v_d=sign * v_d, W_um=w, L_um=l,
                             t_um=t, vth_current=ctx.vth_current, vth_interpol=ctx.vth_interpol,
                             vth_global=ctx.vth_global, vth_logderivative=ctx.vth_logderivative,
                             vth_sqrt=ctx.vth_sqrt, vth_d2=ctx.vth_d2, on_off=ctx.onoff_ratio,
                             ss=ctx.subthreshold_swing, mu_fe_max=ctx.max_mu_linear, mu_sat_sqrt=ctx.mu_sat_sqrt,
                             ss_v_start=ctx.V_g_fine[ctx.ss_start] if ctx.has_ss_region else np.nan,
                             ss_v_end=ctx.V_g_fine[ctx.ss_end] if ctx.has_ss_region else np.nan,
                             n_fine=ctx.V_g_fine.size))
    return rows


def main(cfg: AppConfig) -> None:
    script = cfg.original_script.resolve()
    src = script.read_text(encoding="utf-8").rstrip()
    assert src.endswith("main()"), "unexpected end of the original script"
    src = src[:-len("main()")]
    out = cfg.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    for f in cfg.data_dir.glob("*.xls"):
        shutil.copy(f, out / f.name)
    os.chdir(out)
    for variant in cfg.variants:
        which = FIXED if variant == "fixed" else []
        ns = {"__file__": str(script), "__name__": "original"}
        exec(compile(patched(src, which), str(script), "exec"), ns)
        rows = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for fname, segment, conditions in cfg.jobs:
                rows += run_file(ns, fname, segment, conditions, use_column="B9_vdcol" in which)
        name = "original_scripts_baseline.csv" if variant == "original" else "original_scripts_fixed_baseline.csv"
        pd.DataFrame(rows).to_csv(out / name, index=False, float_format="%.17g")
        print("wrote", out / name, len(rows), "rows")


if __name__ == "__main__":
    main(AppConfig())
