#!/usr/bin/env python3
"""Assemble the Phase 4 cross-code table from every LAMMPS run family.

Families:
  from_empty   filled from an empty pore (first attempt)
  continued    the same runs, resumed from their final configuration
  seeded       started from RASPA's equilibrated configuration at that pressure
  from_above   started from a configuration with MORE molecules (the 50 bar one)

The point of the last two is that they approach the equilibrium from opposite
sides: if the two ensembles were identical and the sampling ergodic, all families
would meet at the same loading. Where they do not meet, the spread between them
is the honest uncertainty -- not the block error bar of any single run.

Reported per run: mean loading with RASPA's block convention, the drift
(last block - first block), the fluctuation sd(N) against RASPA's (a physical
property of the ensemble), and the pass test |delta| < sqrt(e1^2+e2^2).

Usage: python scripts/crosscheck_summary.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from crosscheck_gcmc import block_stats, lammps_kspace, raspa_fluctuation, read_thermo  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CROSS = ROOT / "runs" / "crosscheck"
FAMILIES = [("from_empty", CROSS / "gcmc", 1500),
            ("continued", CROSS / "continued", 500),
            ("seeded", CROSS / "seeded", 500),
            ("from_above", CROSS / "seeded_from_above", 500)]


def main():
    raspa = pd.read_csv(ROOT / "results" / "isotherm_MIL53_lp_CO2_304K.csv", comment="#")
    rows = []
    for family, base, equil in FAMILIES:
        for d in sorted(base.glob("p_*bar"), key=lambda p: float(p.name[2:-3])):
            log = d / "gcmc.log"
            if not log.exists():
                continue
            p_bar = float(d.name[2:-3])
            th = read_thermo(log)
            if th.empty:
                continue
            prod = th[th.step > th.step.min() + equil]
            mean, err, nsamp, drifting, drift = block_stats(prod.nuc)
            r = raspa[np.isclose(raspa.p_bar, p_bar)].iloc[0]
            sd_raspa, _ = raspa_fluctuation(r.file)
            gvec, kmax1d, nvec = lammps_kspace(log)
            delta = mean - r.absolute_molec_uc
            tol = float(np.hypot(err, r.absolute_molec_uc_err95))
            rows.append({"family": family, "p_bar": p_bar,
                         "start_uc": float(th.nuc.iloc[0]), "end_uc": float(th.nuc.iloc[-1]),
                         "LAMMPS_uc": mean, "LAMMPS_err95": err, "drift": drift,
                         "sd_N_uc_LAMMPS": float(prod.nuc.std()), "sd_N_uc_RASPA": sd_raspa,
                         "fluctuation_ratio": float(prod.nuc.std()) / sd_raspa,
                         "RASPA_uc": r.absolute_molec_uc, "RASPA_err95": r.absolute_molec_uc_err95,
                         "delta": delta, "tolerance": tol,
                         "within_tolerance": bool(abs(delta) < tol),
                         "equilibrated": not drifting,
                         "LAMMPS_G_per_A": gvec, "LAMMPS_kmax1d": kmax1d, "LAMMPS_kvectors": nvec})
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "results" / "crosscheck_gcmc_summary.csv", index=False)
    cols = ["family", "p_bar", "start_uc", "LAMMPS_uc", "LAMMPS_err95", "end_uc", "drift",
            "fluctuation_ratio", "RASPA_uc", "delta", "tolerance", "within_tolerance", "equilibrated"]
    with pd.option_context("display.width", 250, "display.float_format", "{:.3f}".format):
        print(df[cols].to_string(index=False))

    print("\nBracket per pressure (lowest and highest family mean, and where RASPA sits):")
    for p_bar, g in df.groupby("p_bar"):
        lo, hi = g.LAMMPS_uc.min(), g.LAMMPS_uc.max()
        r = g.RASPA_uc.iloc[0]
        inside = lo <= r <= hi
        print(f"  {p_bar:5g} bar: LAMMPS {lo:.3f} ... {hi:.3f} ({len(g)} runs), RASPA {r:.3f} -> "
              + ("RASPA inside the bracket" if inside else "RASPA OUTSIDE the bracket"))
    print("\nwrote results/crosscheck_gcmc_summary.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
