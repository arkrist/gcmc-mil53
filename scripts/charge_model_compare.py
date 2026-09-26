#!/usr/bin/env python3
"""lp-PACMAN vs lp-SABVUN-DDEC: how much does the framework charge set move the
CO2 isotherm of the SAME lp pore?

Why this run exists (NOTES.md 5.2): the np structure only exists as a PACMAN-DDEC6
file, so the osmotic construction needs a PACMAN lp partner to be charge-consistent.
But the PACMAN lp file carries an O-H bond of 0.661 A -- unphysically short, a defect
of that refinement -- and a mu-OH dipole 32 % weaker than its np partner's, at exactly
the site Bourrelly identifies as the first CO2 binding site. So the pair is consistent
in charge METHOD but not in the geometry that matters most.

Four pressures, chosen to straddle the branch (0.1, 1, 10, 30 bar): low coverage where
the mu-OH site dominates, and the saturation plateau where it should not matter. The
question is whether swapping SABVUN-DDEC for PACMAN changes the lp isotherm by more
than the ~5 % that would force a full PACMAN lp isotherm before the construction.

Everything here is ABSOLUTE loading, which no void fraction enters. Both GCMC runs were
given HeliumVoidFraction = 0.7115 (the SABVUN-DDEC value), because the PACMAN helium run
finished at 08:06 and the PACMAN GCMC had started at 06:54. The PACMAN Widom value is
0.709397, 0.3 % lower, so the PACMAN excess column is built on a void fraction 0.3 % too
high -- 0.3 % of a correction that is itself ~3 % of the loading at 30 bar. Both numbers
are carried in results/isotherm_MIL53_lp_pacman_CO2_304K.csv (theta_He, theta_He_widom).

Outputs:
  results/charge_model_comparison.csv
  results/charge_model_comparison.png

Usage: python scripts/charge_model_compare.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from isotherm import MOLKG_PER_UC, ROOT  # noqa: E402
from langmuir import fit_langmuir  # noqa: E402
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

DDEC_CSV = ROOT / "results" / "isotherm_MIL53_lp_CO2_304K.csv"
PACMAN_CSV = ROOT / "results" / "isotherm_MIL53_lp_pacman_CO2_304K.csv"
AGREEMENT_THRESHOLD = 0.05      # 5 %: above this, a full PACMAN lp isotherm is needed


def main():
    ddec = pd.read_csv(DDEC_CSV, comment="#")
    pac = pd.read_csv(PACMAN_CSV, comment="#")

    # keep only the pressures both runs have (the PACMAN run is the 4-point subset)
    common = sorted(set(np.round(pac.p_bar, 10)) & set(np.round(ddec.p_bar, 10)))
    missing = sorted(set(np.round(pac.p_bar, 10)) - set(np.round(ddec.p_bar, 10)))
    if missing:
        print(f"  NOTE: PACMAN pressures with no DDEC counterpart, not compared: {missing}")
    d = ddec[np.round(ddec.p_bar, 10).isin(common)].set_index("p_bar").sort_index()
    p = pac[np.round(pac.p_bar, 10).isin(common)].set_index("p_bar").sort_index()

    rows = []
    for pb in common:
        nd, ed = d.loc[pb, "absolute_mol_kg"], d.loc[pb, "absolute_mol_kg_err95"]
        npc, ep = p.loc[pb, "absolute_mol_kg"], p.loc[pb, "absolute_mol_kg_err95"]
        ratio = npc / nd
        # 95 % band on the ratio, errors added in quadrature (independent runs)
        rel = np.hypot(ed / nd, ep / npc)
        diff = npc - nd
        rows.append({
            "p_bar": pb,
            "ddec_mol_kg": nd, "ddec_err95": ed,
            "pacman_mol_kg": npc, "pacman_err95": ep,
            "diff_mol_kg": diff, "diff_err95": float(np.hypot(ed, ep)),
            "ratio_pacman_over_ddec": ratio,
            "ratio_err95": ratio * rel,
            "pct_change": 100.0 * (ratio - 1.0),
            "significant_95": bool(abs(diff) > np.hypot(ed, ep)),
            "ddec_molec_uc": d.loc[pb, "absolute_molec_uc"],
            "pacman_molec_uc": p.loc[pb, "absolute_molec_uc"],
            "ddec_theta_He": d.loc[pb, "theta_He"],
            "pacman_theta_He": p.loc[pb, "theta_He"],
            "ddec_excess_mol_kg": d.loc[pb, "excess_mol_kg"],
            "pacman_excess_mol_kg": p.loc[pb, "excess_mol_kg"],
            "ddec_acc_insertion": d.loc[pb, "acc_insertion"],
            "pacman_acc_insertion": p.loc[pb, "acc_insertion"],
            "ddec_cycles_prod": d.loc[pb, "cycles_prod"],
            "pacman_cycles_prod": p.loc[pb, "cycles_prod"],
        })
    cmp = pd.DataFrame(rows)

    print("=" * 96)
    print("lp-PACMAN vs lp-SABVUN-DDEC, absolute loading, 304 K, identical force field and box")
    print("=" * 96)
    show = ["p_bar", "ddec_mol_kg", "ddec_err95", "pacman_mol_kg", "pacman_err95",
            "pct_change", "ratio_err95", "significant_95"]
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(cmp[show].to_string(index=False))

    dev = cmp["pct_change"].abs()
    worst = dev.max()
    print(f"\n  largest deviation: {worst:.2f} % (at {cmp.loc[dev.idxmax(), 'p_bar']:g} bar)")
    print(f"  mean |deviation|:  {dev.mean():.2f} %")
    n_sig = int(cmp["significant_95"].sum())
    print(f"  statistically significant at 95 %: {n_sig} of {len(cmp)} points")

    # ---- the decision the construction depends on
    verdict = worst <= 100.0 * AGREEMENT_THRESHOLD
    print("\n" + "-" * 96)
    if verdict:
        print(f"  VERDICT: agreement within {100 * AGREEMENT_THRESHOLD:.0f} % at every pressure.")
        print("  -> Four PACMAN points are enough. The osmotic construction uses the 12-point")
        print("     SABVUN-DDEC Langmuir fit for the lp branch, because the charge set does not")
        print("     move that branch by more than the fit's own uncertainty, and 12 points over")
        print("     0.01-50 bar constrain N_max and b far better than 4 over 0.1-30 bar.")
    else:
        print(f"  VERDICT: deviation exceeds {100 * AGREEMENT_THRESHOLD:.0f} %.")
        print("  -> The remaining PACMAN lp pressures must be run before the construction.")
    print("-" * 96)

    # ---- Langmuir fits, the quantity the construction actually consumes
    print("\nLangmuir fits (absolute, mol/kg; the construction consumes K = N_max b):")
    fits = {}
    for name, df, pmax in (("lp SABVUN-DDEC, 12 pts, 0.01-50 bar", ddec, None),
                           ("lp SABVUN-DDEC, same 4 pressures", d.reset_index(), None),
                           ("lp PACMAN, 4 pts, 0.1-30 bar", pac, None)):
        f = fit_langmuir(df.p_bar, df.absolute_mol_kg, df.absolute_mol_kg_err95 / 2.776)
        fits[name] = f
        print(f"  {name:36s} N_max = {f['N_max']:7.3f} mol/kg, b = {f['b_per_bar']:.4g} /bar, "
              f"K = {f['K_mol_kg_Pa']:.4g} mol/kg/Pa, ok = {f['ok']}")
    k_ddec12 = fits["lp SABVUN-DDEC, 12 pts, 0.01-50 bar"]["K_mol_kg_Pa"]
    k_ddec4 = fits["lp SABVUN-DDEC, same 4 pressures"]["K_mol_kg_Pa"]
    k_pac = fits["lp PACMAN, 4 pts, 0.1-30 bar"]["K_mol_kg_Pa"]
    # The like-for-like number is the one fitted over the SAME four pressures. Comparing
    # the 4-point PACMAN fit with the 12-point DDEC fit mixes the charge effect with a
    # window effect, and the window effect is the larger of the two.
    print(f"  K(PACMAN, 4 pts) / K(DDEC, same 4 pts)  = {k_pac / k_ddec4:.3f}   <- charge effect")
    print(f"  K(PACMAN, 4 pts) / K(DDEC, 12 pts)      = {k_pac / k_ddec12:.3f}   <- charge + fit window")
    print(f"  K(DDEC, 4 pts)   / K(DDEC, 12 pts)      = {k_ddec4 / k_ddec12:.3f}   <- fit window alone")

    out = ROOT / "results"
    csv = out / "charge_model_comparison.csv"
    with open(csv, "w") as fh:
        fh.write("# lp-PACMAN vs lp-SABVUN-DDEC, CO2 in MIL-53(Al) lp, rigid, 304 K, RASPA2 2.0.50\n")
        fh.write("# identical force field, box (4x2x2), cutoff and Ewald precision; only the CIF charges differ\n")
        fh.write("# ABSOLUTE loading is the comparison. Both GCMC runs were given the same\n")
        fh.write(f"# HeliumVoidFraction ({d.theta_He.iloc[0]:g}) in simulation.input, because the PACMAN helium\n")
        fh.write(f"# run (Widom: {p.theta_He_widom.iloc[0]:g}) finished after the GCMC had started. That affects\n")
        fh.write("# the excess column by 0.3 % of the excess correction and the absolute column not at all.\n")
        fh.write("# errors are RASPA's 95 % confidence half-widths; ratio_err95 adds them in quadrature\n")
        fh.write(f"# significant_95: |difference| exceeds the combined 95 % band\n")
        cmp.to_csv(fh, index=False)
    print(f"\nwrote {csv.relative_to(ROOT)}")

    # ---------------- figure ----------------
    # wspace: the left panel carries a secondary axis on its right, which otherwise
    # collides with the right panel's y-label
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.6), gridspec_kw={"wspace": 0.42})
    ax.errorbar(ddec.p_bar, ddec.absolute_mol_kg, yerr=ddec.absolute_mol_kg_err95,
                label="lp SABVUN-DDEC (production, 12 pts)", **SIM)
    ax.errorbar(pac.p_bar, pac.absolute_mol_kg, yerr=pac.absolute_mol_kg_err95,
                label="lp PACMAN-DDEC6 (4 pts)",
                **{**REF, "marker": "s", "ls": "none", "mfc": "white"})
    ax.set_xscale("log")
    ax.set_xlabel("Pressure [bar]")
    ax.set_ylabel("Absolute loading [mol/kg]")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x / MOLKG_PER_UC, lambda x: x * MOLKG_PER_UC))
    sec.set_ylabel("Loading [molecules / unit cell]")
    ax.set_title("Same lp pore, two charge sets", fontsize=10)
    ax.legend(fontsize=8.5, loc="upper left")

    ax2.axhline(0, color="#8c8c88", lw=1.0)
    ax2.axhspan(-100 * AGREEMENT_THRESHOLD, 100 * AGREEMENT_THRESHOLD, color=SIM["color"],
                alpha=0.10, lw=0, label=f"$\\pm${100 * AGREEMENT_THRESHOLD:.0f} % decision band")
    ax2.errorbar(cmp["p_bar"], cmp["pct_change"], yerr=100 * cmp["ratio_err95"],
                 **{**SIM, "marker": "D"}, label="PACMAN - DDEC")
    ax2.set_xscale("log")
    ax2.set_xlabel("Pressure [bar]")
    ax2.set_ylabel("Change in absolute loading [%]")
    ax2.set_ylim(-8, 8)
    ax2.set_title("Deviation, with 95 % bars", fontsize=10)
    ax2.legend(fontsize=8.5, loc="lower right")
    fig.suptitle("Framework charge set: PACMAN-DDEC6 vs SABVUN-DDEC, CO$_2$/MIL-53(Al) lp, 304 K", y=1.02)
    fig.text(0.01, -0.08,
             "Identical UFF + TraPPE force field, 4x2x2 box, 12.0 A cutoff, Ewald 1e-6; only the "
             "_atom_site_charge column differs.\nThe PACMAN lp file has an O-H bond of 0.661 A and a "
             "mu-OH dipole 32 % weaker than its np partner's (NOTES.md 5.2),\nwhich is why this "
             "comparison was run before using the PACMAN pair in the osmotic construction.",
             fontsize=7.5, va="top")
    fig.savefig(out / "charge_model_comparison.png")
    print("wrote results/charge_model_comparison.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
