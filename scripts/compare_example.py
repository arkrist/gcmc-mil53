#!/usr/bin/env python3
"""Sanity check: compare our run of a bundled RASPA example with the reference
output shipped in the RASPA2 source tree.

Pass criterion (agreed in Phase 1): for each component and pressure,
    |L_ours - L_ref| < sqrt(e_ours^2 + e_ref^2)
with L the absolute loading (molecules/unit cell) and e RASPA's 95 % block
error bars. The two runs use different builds, so random streams differ and
only statistical agreement is expected.

Usage: python scripts/compare_example.py <example-name>   (e.g. mfi_ch4)
Writes results/sanity_<name>.csv and results/sanity_<name>.png.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import parse_raspa_output as pro  # noqa: E402
from plotstyle import REF, SIM, plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def load(folder):
    files = sorted(folder.glob("*.data"))
    if not files:
        sys.exit(f"no .data files in {folder}")
    return pd.concat([pro.parse_file(f) for f in files], ignore_index=True)


def main(name):
    ex = ROOT / "examples" / name
    ours = load(ex / "run" / "Output" / "System_0")
    ref = load(ex / "reference")
    keep = ["molecule", "p_Pa", "T_K", "fugacity_coeff", "absolute_molec_uc", "absolute_molec_uc_err95",
            "absolute_mol_kg", "absolute_mol_kg_err95", "excess_mol_kg", "excess_mol_kg_err95",
            "acc_insertion", "acc_deletion", "cycles_init", "cycles_prod", "finished"]
    m = ours[keep].merge(ref[keep], on=["molecule", "p_Pa"], suffixes=("_ours", "_ref"))
    m["delta_molec_uc"] = m.absolute_molec_uc_ours - m.absolute_molec_uc_ref
    m["tol_molec_uc"] = np.hypot(m.absolute_molec_uc_err95_ours, m.absolute_molec_uc_err95_ref)
    m["pass"] = m.delta_molec_uc.abs() < m.tol_molec_uc
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    m.to_csv(out / f"sanity_{name}.csv", index=False)

    show = ["molecule", "p_Pa", "absolute_molec_uc_ours", "absolute_molec_uc_err95_ours",
            "absolute_molec_uc_ref", "absolute_molec_uc_err95_ref", "delta_molec_uc", "tol_molec_uc",
            "acc_insertion_ours", "acc_insertion_ref", "pass"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(m[show].to_string(index=False))

    # Top: loadings (markers only -- 1-2 points per example, a line would suggest an
    # interpolation that was not computed). Bottom: the pass test itself,
    # delta = ours - ref with the combined tolerance as error bar; pass <=> bar crosses 0.
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(5.8, 5.6), sharex=True,
                                  gridspec_kw={"height_ratios": [2.2, 1]})
    for mol, g in m.groupby("molecule"):
        g = g.sort_values("p_Pa")
        p_kpa = g.p_Pa / 1e3
        ax.errorbar(p_kpa, g.absolute_mol_kg_ref, yerr=g.absolute_mol_kg_err95_ref, zorder=2,
                    label=f"{mol}: RASPA reference output (2.0.45)", **{**REF, "ms": 11})
        ax.errorbar(p_kpa, g.absolute_mol_kg_ours, yerr=g.absolute_mol_kg_err95_ours, zorder=3,
                    label=f"{mol}: this work (RASPA 2.0.50, osx-arm64)", **{**SIM, "lw": 0, "elinewidth": 1.5, "ms": 6})
        axd.errorbar(p_kpa, g.delta_molec_uc, yerr=g.tol_molec_uc, zorder=3,
                     **{**SIM, "lw": 0, "elinewidth": 1.5, "ms": 6})
    ax.set_xscale("log")
    if len(m) > 1:
        ax.set_yscale("log")
    ax.set_ylabel("Absolute loading [mol/kg]")
    fw = ours.framework.iloc[0]
    ax.set_title(f"Sanity check: {fw}, T = {ours.T_K.iloc[0]:g} K (error bars: 95 % CI)")
    ax.legend(fontsize=9, loc="upper left")
    axd.axhline(0, color="#3d3d3a", lw=1)
    axd.set_ylabel("ours - ref\n[molec/uc]")
    axd.set_xlabel("Pressure [kPa]")
    axd.margins(x=0.3)
    fig.savefig(out / f"sanity_{name}.png")
    print(f"\n{'PASS' if m['pass'].all() else 'FAIL'}: {int(m['pass'].sum())}/{len(m)} points within tolerance")
    return 0 if m["pass"].all() else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
