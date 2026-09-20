#!/usr/bin/env python3
"""Assemble the CO2/MIL-53(Al) lp isotherm from runs/<tag>/p_*/ and plot it.

Every quantity is carried in both conventions (absolute and excess) and both
units (mol/kg and molecules per conventional unit cell, 832.4 g/mol), with
RASPA's 95 % error bars, the accepted insertion/deletion counts and the
relative error per point (see NOTES.md, "What every result carries").

Outputs:
  results/isotherm_MIL53_lp_CO2_304K.csv
  results/isotherm_MIL53_lp_CO2_304K.png

Usage: python scripts/isotherm.py [--tag production] [--no-reference]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import parse_raspa_output as pro  # noqa: E402
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
M_UC = 832.415          # g/mol per conventional unit cell (Al4C32H20O20)
MOLKG_PER_UC = 1000.0 / M_UC   # 1 molecule/uc = 1.2013 mol/kg
THETA_HE = 0.7115       # Widom He, 298 K, eps/k 10.9 K, sigma 2.64 A
REF_CSV = ROOT / "reference" / "bourrelly2005_MIL53Al_CO2_304K.csv"
LP_WINDOW_BAR = 9.0     # only above this is the real material fully lp
STEP_LO_BAR = 5.0       # 5-9 bar: np-lp breathing step


def collect(tag):
    rows = []
    for d in sorted((ROOT / "runs" / tag).glob("p_*"), key=lambda p: float(p.name[2:])):
        files = list((d / "Output" / "System_0").glob("*.data")) if (d / "Output" / "System_0").is_dir() else []
        if not files:
            continue
        r = pro.parse_file(files[0]).iloc[0].to_dict()
        if not r.get("finished"):
            print(f"  (skipping {d.name}: not finished)")
            continue
        rows.append(r)
    if not rows:
        sys.exit("no finished points found")
    df = pd.DataFrame(rows)
    df["p_bar"] = df.p_Pa / 1e5
    df["fugacity_bar"] = df.fugacity_Pa / 1e5
    for kind in ("absolute", "excess"):
        df[f"{kind}_molec_uc_rel_err"] = df[f"{kind}_molec_uc_err95"] / df[f"{kind}_molec_uc"]
    df["theta_He"] = THETA_HE
    df["M_uc_g_per_mol"] = M_UC
    return df.sort_values("p_bar").reset_index(drop=True)


def load_reference():
    if not REF_CSV.exists():
        return None
    r = pd.read_csv(REF_CSV)
    r["in_lp_window"] = r.pressure_bar >= LP_WINDOW_BAR
    return r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="production")
    ap.add_argument("--no-reference", action="store_true")
    a = ap.parse_args(argv)

    df = collect(a.tag)
    ref = None if a.no_reference else load_reference()
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    cols = ["p_bar", "p_Pa", "fugacity_bar", "fugacity_coeff", "T_K", "unit_cells", "n_unit_cells",
            "absolute_mol_kg", "absolute_mol_kg_err95", "absolute_molec_uc", "absolute_molec_uc_err95",
            "absolute_molec_uc_rel_err",
            "excess_mol_kg", "excess_mol_kg_err95", "excess_molec_uc", "excess_molec_uc_err95",
            "excess_molec_uc_rel_err", "bulk_density_kg_m3", "theta_He", "M_uc_g_per_mol",
            "acc_insertion", "accepted_insertion", "acc_deletion", "accepted_deletion",
            "acc_reinsertion", "acc_translation", "acc_rotation", "cycles_init", "cycles_prod",
            "n_warnings", "missing_vdw_pairs", "file"]
    csv = out / "isotherm_MIL53_lp_CO2_304K.csv"
    header = (f"# CO2 in MIL-53(Al) lp, rigid framework, GCMC (RASPA2 2.0.50), T = {df.T_K.iloc[0]:g} K\n"
              f"# absolute and excess loading; 1 molecule/uc = {MOLKG_PER_UC:.4f} mol/kg (M_uc = {M_UC} g/mol)\n"
              f"# excess uses the helium void fraction theta_He = {THETA_HE} (Widom He, 298 K, eps/k 10.9 K, sigma 2.64 A)\n"
              f"# errors are RASPA's 95 % confidence half-widths (5 blocks, t = 2.776)\n")
    with open(csv, "w") as fh:
        fh.write(header)
        df[[c for c in cols if c in df.columns]].to_csv(fh, index=False)
    print(f"wrote {csv.relative_to(ROOT)}")

    show = ["p_bar", "absolute_mol_kg", "absolute_mol_kg_err95", "absolute_molec_uc",
            "excess_mol_kg", "excess_molec_uc", "absolute_molec_uc_rel_err",
            "acc_insertion", "accepted_insertion", "accepted_deletion"]
    with pd.option_context("display.width", 220, "display.float_format", "{:.4g}".format):
        print(df[[c for c in show if c in df.columns]].to_string(index=False))
    bad = df[df.absolute_molec_uc_rel_err > 0.10]
    for _, r in bad.iterrows():
        print(f"  FLAG: {r.p_bar:g} bar has {100 * r.absolute_molec_uc_rel_err:.1f} % relative error on loading (> 10 %)")
    low = df[df.acc_insertion < 0.01]
    for _, r in low.iterrows():
        print(f"  FLAG: {r.p_bar:g} bar insertion acceptance {100 * r.acc_insertion:.2f} % (< 1 %), "
              f"{r.accepted_insertion:.0f} accepted insertions")

    # ---------------- figure ----------------
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    solid = df.absolute_molec_uc_rel_err <= 0.10
    ax.fill_between(df.p_bar, df.excess_mol_kg, df.absolute_mol_kg, color=SIM["color"], alpha=0.16, lw=0,
                    label=f"absolute - excess (theta$_{{He}}$ = {THETA_HE})")
    ax.errorbar(df.p_bar, df.absolute_mol_kg, yerr=df.absolute_mol_kg_err95,
                label="this work, absolute", **SIM)
    ax.errorbar(df.p_bar, df.excess_mol_kg, yerr=df.excess_mol_kg_err95,
                label="this work, excess", **{**SIM, "marker": "v", "ls": "--", "mfc": "white"})
    if (~solid).any():
        ax.plot(df.p_bar[~solid], df.absolute_mol_kg[~solid], "o", mfc="white", mec=SIM["color"], ms=9, lw=0,
                label="rel. error > 10 % (not solid)")
    if ref is not None:
        for flag, style, lab in ((True, REF, f"Bourrelly 2005, excess (assumed), P $\\geq$ {LP_WINDOW_BAR:g} bar"),
                                 (False, {**REF, "color": "#8c8c88", "mfc": "white", "ms": 7},
                                  f"same, {STEP_LO_BAR:g}-{LP_WINDOW_BAR:g} bar: inside the np-lp step, NOT compared")):
            g = ref[ref.in_lp_window == flag]
            if len(g):
                ax.errorbar(g.pressure_bar, g.loading_mmol_per_g, label=lab, **style)
    # Langmuir curves (fitted over P >= 9 bar; see scripts/langmuir.py)
    try:
        from langmuir import fit_langmuir
        w = df[df.p_bar >= LP_WINDOW_BAR]
        if len(w) >= 3:
            fs = fit_langmuir(w.p_bar, w.absolute_mol_kg, w.absolute_mol_kg_err95 / 2.776)
            pp = np.logspace(np.log10(df.p_bar.min()), np.log10(df.p_bar.max()), 300)
            if fs["ok"]:
                ax.plot(pp, fs["N_max"] * fs["b_per_bar"] * pp / (1 + fs["b_per_bar"] * pp),
                        color=SIM["color"], lw=1.1, ls=":", zorder=1,
                        label=f"Langmuir fit, this work (absolute, P $\\geq$ {LP_WINDOW_BAR:g} bar)")
            if ref is not None:
                er = ref[ref.in_lp_window]
                fe = fit_langmuir(er.pressure_bar, er.loading_mmol_per_g, np.full(len(er), 0.05))
                ax.plot(pp, fe["N_max"] * fe["b_per_bar"] * pp / (1 + fe["b_per_bar"] * pp),
                        color=REF["color"], lw=1.1, ls="--", zorder=1,
                        label="Langmuir fit, experiment = Coudert virtual lp curve")
    except Exception as exc:  # never let the figure fail on the fit
        print(f"  (Langmuir overlay skipped: {exc})")
    ax.set_xscale("log")
    ax.set_xlabel("Pressure [bar]")
    ax.set_ylabel("Loading [mol/kg  =  mmol/g]")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x / MOLKG_PER_UC, lambda x: x * MOLKG_PER_UC))
    sec.set_ylabel("Loading [molecules / unit cell]")
    ax.set_title(f"CO$_2$ in MIL-53(Al) lp, rigid framework, {df.T_K.iloc[0]:g} K")
    ax.legend(fontsize=8.5, loc="upper left")
    fig.text(0.01, -0.06,
             "Rigid lp framework: this is the virtual lp branch, so it is comparable with experiment only "
             f"for P $\\geq$ {LP_WINDOW_BAR:g} bar,\nwhere the real material is fully open. Below ~5 bar the real solid is np "
             "(Coudert 2008). Reference points are excess\n(assumed: manometry, convention not stated by the source) and are "
             "never converted; the band shows our absolute-excess gap.\nError bars: 95 % CI.",
             fontsize=7.5, va="top")
    fig.savefig(out / "isotherm_MIL53_lp_CO2_304K.png")
    print(f"wrote {(out / 'isotherm_MIL53_lp_CO2_304K.png').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
