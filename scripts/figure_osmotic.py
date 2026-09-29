#!/usr/bin/env python3
"""The summary figure: Omega_lp and Omega_np against P, experiment and simulation.

Three panels, which is the argument of 5.8 in one picture:

  (a) EXPERIMENT -- Coudert's Langmuir parameters. Both transitions, 0.27 and 3.97 bar,
      against the measured 0.25-0.3 and 5-6 bar. This panel is the validation.
  (b) SIMULATION, as computed -- the np branch in the experimental CO2-loaded cell and
      the lp branch, both integrated numerically from the measured points. The double
      transition is there (case c), but at 0.013 and 0.27 bar: ~20x too low.
  (c) SIMULATION, corresponding states -- each branch's pressure axis stretched by its
      own K_sim/K_exp. The closing lands on experiment; the reopening does not, and the
      residual factor of ~3 is what 5.8.6 lists three candidate explanations for.

Usage: python scripts/figure_osmotic.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import osmotic_numeric as on  # noqa: E402
from osmotic import (K_LP_EXP, K_NP_EXP, MOLKG_TO_PERCELL, NMAX_NP_EXP,  # noqa: E402
                     V_LP, V_NP)
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KH_LP, KH_NP = 1.86944e-4, 5.6274e-4          # RASPA Widom, both measured
DF = 2.5


def experimental_branches(grid):
    nmax_lp = 12.183 * MOLKG_TO_PERCELL
    b_lp = (K_LP_EXP / 12.183) * 1e5
    b_np = (K_NP_EXP / (NMAX_NP_EXP / MOLKG_TO_PERCELL)) * 1e5
    n_lp = nmax_lp * b_lp * grid / (1 + b_lp * grid)
    n_np = NMAX_NP_EXP * b_np * grid / (1 + b_np * grid)
    return (on.Branch(grid, n_lp, nmax_lp * b_lp / 1e5),
            on.Branch(grid, n_np, NMAX_NP_EXP * b_np / 1e5))


def panel(ax, br_lp, br_np, dF, title, sub, mark_exp=True):
    p, o_lp, o_np, cr = on.transitions_numeric(br_lp, br_np, dF, p_range=(1e2, 2e6), n=8000)
    ax.plot(p / 1e5, o_lp, color=SIM["color"], lw=2.2, label=r"$\Omega_{lp}$")
    ax.plot(p / 1e5, o_np, color=REF["color"], lw=2.2, ls="--", label=r"$\Omega_{np}$")
    if mark_exp:
        for pe, lab in ((0.27, None), (5.5, None)):
            ax.axvspan(pe * 0.9, pe * 1.1, color="#c8c8c2", alpha=0.45, lw=0, zorder=0)
    for pc, w in cr:
        ax.axvline(pc, color=EXTRA[0], lw=1.3, ls=":")
        ax.annotate(f"{pc:.3f} bar\n{w}", xy=(pc, 0.02), xycoords=("data", "axes fraction"),
                    fontsize=7.5, rotation=90, va="bottom", ha="right", color=EXTRA[0],
                    xytext=(-1.5, 0), textcoords="offset points")
    if not cr:
        ax.annotate("no crossing", xy=(0.05, 0.10), xycoords="axes fraction",
                    fontsize=9, color=EXTRA[0])
    ax.set_xscale("log")
    ax.set_xlabel("Pressure [bar]")
    ax.set_title(title, fontsize=10.5)
    # upper right is empty in every panel (Omega falls steeply there), unlike upper
    # left where the near-zero plateau of both curves sits
    ax.annotate(sub, xy=(0.97, 0.955), xycoords="axes fraction", fontsize=8,
                va="top", ha="right", color="#55554f")
    ax.legend(fontsize=9, loc="lower left")
    return cr


def main():
    grid = np.logspace(-2, 1.7, 40)
    ex_lp, ex_np = experimental_branches(grid)
    sim_np, *_ = on.load_branch(ROOT / "results/isotherm_MIL53_npCO2cell_CO2_304K.csv", KH_NP)
    sim_lp, *_ = on.load_branch(ROOT / "results/isotherm_MIL53_lp_CO2_304K.csv", KH_LP)
    f_lp, f_np = KH_LP / K_LP_EXP, KH_NP / K_NP_EXP
    res_lp, res_np = on.rescale_branch(sim_lp, f_lp), on.rescale_branch(sim_np, f_np)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7), gridspec_kw={"wspace": 0.22})
    panel(axes[0], ex_lp, ex_np, DF, "(a) experiment (Coudert 2008 parameters)",
          "validation of the construction")
    panel(axes[1], sim_lp, sim_np, DF, "(b) simulation, as computed",
          "np in the Serre 2007 CO$_2$ cell\nboth branches integrated numerically")
    panel(axes[2], res_lp, res_np, DF, "(c) simulation, corresponding states",
          f"each branch stretched by its own\nK$_{{sim}}$/K$_{{exp}}$ ({f_lp:.1f} lp, {f_np:.1f} np)")
    axes[0].set_ylabel(r"$\Omega_{os}$ [kJ/mol per unit cell]")
    fig.suptitle("Osmotic construction for CO$_2$-induced breathing of MIL-53(Al), 304 K, "
                 rf"$\Delta F_{{host}}$ = {DF:g} kJ/mol per cell", y=1.03, fontsize=12)
    fig.text(0.01, -0.11,
             "Grey bands: the measured transitions, 0.25-0.3 bar (closing) and 5-6 bar (reopening, adsorption branch). "
             "Omega from Coudert eq. 8 integrated\nnumerically on the isotherms themselves (monotone PCHIP in log p, "
             "analytic Henry segment anchored on the measured Widom K$_H$,\nPeng-Robinson V$_m$); the single-site "
             "Langmuir of eq. 11 is not used, because neither simulated branch is Langmuir (chi$^2_{red}$ 318 np, "
             "1368 lp).\nPanel (c) shows that the closing discrepancy is entirely Henry-regime over-binding, while a "
             "factor ~3 in the reopening is not.",
             fontsize=7.6, va="top")
    out = ROOT / "results" / "osmotic_summary.png"
    fig.savefig(out)
    print(f"wrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
