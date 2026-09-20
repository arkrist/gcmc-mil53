#!/usr/bin/env python3
"""Relative deviation simulation / experiment over the lp window (P >= 9 bar).

The experimental reference is the Langmuir fit to the digitised Bourrelly points
(P >= 9 bar), whose K reproduces Coudert's K_lp to a ratio of 1.00, i.e. it IS the
virtual rigid-lp curve. Using the fitted curve rather than the raw points lets the
ratio be evaluated at our own pressures without converting or interpolating the
reference data.

Two series are drawn (our absolute and our excess loading), because the reference
is excess (assumed) while the Langmuir comparison uses our absolute isotherm; the
two bracket the convention mismatch. Points above 29 bar (the last experimental
pressure) rely on extrapolating the experimental fit and are drawn open.

Usage: python scripts/deviation.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from isotherm import LP_WINDOW_BAR, ROOT
from langmuir import fit_langmuir, REF_CSV
from plotstyle import EXTRA, SIM, plt

P_MAX_EXP_BAR = 29.03   # last experimental pressure; beyond this the fit is extrapolated


def main():
    sim = pd.read_csv(ROOT / "results" / "isotherm_MIL53_lp_CO2_304K.csv", comment="#")
    ref = pd.read_csv(REF_CSV)
    e = ref[ref.pressure_bar >= LP_WINDOW_BAR]
    fe = fit_langmuir(e.pressure_bar, e.loading_mmol_per_g, np.full(len(e), 0.05))
    langmuir = lambda p, f: f["N_max"] * f["b_per_bar"] * p / (1 + f["b_per_bar"] * p)

    s = sim[sim.p_bar >= LP_WINDOW_BAR].copy()
    s["exp_fit"] = langmuir(s.p_bar, fe)
    rows = []
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    for col, style, lab in (("absolute", SIM, "absolute (the quantity fitted)"),
                            ("excess", {**SIM, "marker": "v", "ls": "--", "mfc": "white"},
                             f"excess (same convention as the reference, theta$_{{He}}$ = {sim.theta_He.iloc[0]})")):
        ratio = s[f"{col}_mol_kg"] / s.exp_fit
        err = s[f"{col}_mol_kg_err95"] / s.exp_fit
        ax.errorbar(s.p_bar, ratio, yerr=err, label=f"this work, {lab}", **style)
        inside = s.p_bar <= P_MAX_EXP_BAR
        if (~inside).any():
            ax.plot(s.p_bar[~inside], ratio[~inside], style["marker"], ms=11, mfc="none",
                    mec=style["color"], lw=0)
        for p_bar, r in zip(s.p_bar, ratio):
            rows.append({"p_bar": p_bar, "convention": col, "sim_over_exp": r,
                         "extrapolated_reference": p_bar > P_MAX_EXP_BAR})
    ax.axhline(1.0, color="#3d3d3a", lw=1.2)
    ax.axhspan(0.95, 1.05, color="#8c8c88", alpha=0.15, lw=0, label="+/- 5 %")
    ax.annotate("open symbols: experimental fit extrapolated\nbeyond the last data point (29 bar)",
                xy=(0.98, 0.04), xycoords="axes fraction", ha="right", fontsize=8, color="#3d3d3a")
    ax.set_xscale("log")
    ax.set_xticks([10, 20, 30, 50])
    ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(plt.matplotlib.ticker.NullFormatter())
    ax.set_xlabel("Pressure [bar]")
    ax.set_ylabel("simulation / experiment")
    ax.set_title("CO$_2$ in MIL-53(Al) lp, 304 K: deviation over the lp branch\n"
                 "(reference = Langmuir fit to Bourrelly 2005, P $\\geq$ 9 bar)", fontsize=11)
    ax.legend(fontsize=8.5, loc="upper right")
    fig.text(0.01, -0.07,
             "Where the real material is fully open the rigid-lp GCMC converges onto the experiment: the remaining "
             "deviation at 10 bar is\nthe force field over-binding at low coverage (K$_{sim}$/K$_{exp}$ = 3.95), not the "
             "missing framework flexibility. Below 9 bar (not shown)\nthe real solid is np and no comparison is meaningful. "
             "Error bars: 95 % CI on the simulation only.",
             fontsize=7.5, va="top")
    fig.savefig(ROOT / "results" / "deviation_vs_pressure.png")
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "results" / "deviation_vs_pressure.csv", index=False)
    with pd.option_context("display.float_format", "{:.4g}".format):
        print(df.pivot(index="p_bar", columns="convention", values="sim_over_exp").to_string())
    print("wrote results/deviation_vs_pressure.{png,csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
