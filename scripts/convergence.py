#!/usr/bin/env python3
"""Convergence diagnostics for the isotherm.

1. Loading vs cycle for two pressure points (default: the lowest and the highest),
   initialisation and production shown separately, with the production block
   average and its 95 % band. This is what says whether the run was long enough,
   especially where the insertion acceptance is low.
2. Insertion acceptance rate vs pressure, with the accepted COUNTS annotated
   (the count is what the statistics depend on), and the 1 % flag line.

Usage: python scripts/convergence.py [--tag production] [--points 1000 5000000]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import parse_raspa_output as pro  # noqa: E402
from isotherm import MOLKG_PER_UC, ROOT, collect  # noqa: E402
from plotstyle import EXTRA, SIM, plt  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="production")
    ap.add_argument("--points", nargs="*", type=float, default=None,
                    help="pressures in Pa to trace (default: lowest and highest finished)")
    a = ap.parse_args(argv)
    df = collect(a.tag)
    pts = a.points or [df.p_Pa.min(), df.p_Pa.max()]
    out = ROOT / "results"
    out.mkdir(exist_ok=True)

    fig, axes = plt.subplots(len(pts), 1, figsize=(7.0, 3.1 * len(pts)), squeeze=False)
    for ax, p_pa in zip(axes[:, 0], pts):
        row = df.iloc[(df.p_Pa - p_pa).abs().argmin()]
        tr = pro.loading_trace(row.file)
        n_uc = row.n_unit_cells
        for stage, color, lab in (("init", "#8c8c88", "initialisation (discarded)"), ("prod", SIM["color"], "production")):
            g = tr[tr.stage == stage]
            x = g.cycle if stage == "init" else g.cycle + row.cycles_init
            ax.plot(x, g.N_box / n_uc, lw=1.1, color=color, label=lab)
        ax.axvline(row.cycles_init, color="#3d3d3a", lw=1, ls=":")
        ax.axhline(row.absolute_molec_uc, color=EXTRA[0], lw=1.4,
                   label=f"block average {row.absolute_molec_uc:.3f} molec/uc")
        ax.axhspan(row.absolute_molec_uc - row.absolute_molec_uc_err95,
                   row.absolute_molec_uc + row.absolute_molec_uc_err95, color=EXTRA[0], alpha=0.18, lw=0,
                   label="95 % CI")
        ax.set_title(f"{row.p_bar:g} bar  (insertion acceptance {100 * row.acc_insertion:.2f} %, "
                     f"{row.accepted_insertion:.0f} accepted insertions)", fontsize=10)
        ax.set_ylabel("molecules / uc")
        ax.legend(fontsize=8, loc="lower right", ncol=2)
    axes[-1, 0].set_xlabel("MC cycle (initialisation + production)")
    fig.suptitle(f"Loading vs cycle, CO$_2$ in MIL-53(Al) lp, {df.T_K.iloc[0]:g} K", y=1.0)
    fig.savefig(out / "convergence_loading_vs_cycle.png")

    fig2, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(df.p_bar, 100 * df.acc_insertion, **{**SIM, "ms": 6})
    ax.axhline(1.0, color="#eb6834", lw=1.4, ls="--", label="1 % flag threshold")
    for _, r in df.iterrows():
        ax.annotate(f"{r.accepted_insertion / 1000:.0f}k", (r.p_bar, 100 * r.acc_insertion),
                    textcoords="offset points", xytext=(0, 7), ha="center", fontsize=7.5, color="#3d3d3a")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Pressure [bar]")
    ax.set_ylabel("Insertion acceptance [%]")
    ax.set_title("Insertion acceptance vs pressure\n(labels: accepted insertions, thousands)", fontsize=11)
    ax.legend(fontsize=9)
    fig2.savefig(out / "convergence_acceptance.png")

    cols = ["p_bar", "absolute_molec_uc", "absolute_molec_uc_err95", "absolute_molec_uc_rel_err",
            "acc_insertion", "accepted_insertion", "accepted_deletion", "cycles_prod"]
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(df[cols].to_string(index=False))
    for _, r in df[df.acc_insertion < 0.01].iterrows():
        print(f"  FLAG {r.p_bar:g} bar: acceptance {100 * r.acc_insertion:.2f} % < 1 %, "
              f"{r.accepted_insertion:.0f} accepted insertions / {r.accepted_deletion:.0f} deletions")
    print("wrote results/convergence_loading_vs_cycle.png and results/convergence_acceptance.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
