#!/usr/bin/env python3
"""Compare LAMMPS `fix gcmc` loadings with RASPA, point by point.

The LAMMPS loading is block-averaged with RASPA's own convention (5 blocks over
the production part, 95 % confidence half-width = t * s / sqrt(5) with t = 2.776),
so the two error bars mean the same thing.

Pass criterion (same as Phase 1): |L_LAMMPS - L_RASPA| < sqrt(e1^2 + e2^2) on the
absolute loading in molecules per unit cell.

Nothing here tunes anything: if a point fails, the single-point energy comparison
(scripts/crosscheck_energy.py) is where the cause is looked for.

Usage: python scripts/crosscheck_gcmc.py [--equil 1500]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GCMC = ROOT / "runs" / "crosscheck" / "gcmc"
T95 = 2.776


def read_thermo(log):
    """All thermo rows of the run(s) in a LAMMPS log."""
    rows, header = [], None
    for line in Path(log).read_text().splitlines():
        if re.match(r"^\s*Step\s+v_nmol", line):
            header = line.split()
            continue
        if header and re.match(r"^\s*\d+\s", line):
            vals = line.split()
            if len(vals) == len(header):
                rows.append([float(v) for v in vals])
        elif header and line.startswith("Loop time"):
            header = header  # keep parsing further runs
    return pd.DataFrame(rows, columns=["step", "nmol", "nuc", "pe", "evdwl", "ecoul", "elong", "press"])


def block_stats(x, nblocks=5):
    """Mean, 95 % half-width (RASPA's convention) and a DRIFT flag.

    Block statistics describe noise, not a trend. If the series is still drifting
    (last block differs from the first by more than the error bar), the mean is
    biased and the error bar understates the uncertainty -- as happened in the
    first Phase 4 attempt, where the loading was still filling the pore. Such a
    point is flagged rather than reported as a code disagreement.
    """
    x = np.asarray(x, float)
    n = len(x) // nblocks * nblocks
    if n < nblocks:
        return float("nan"), float("nan"), 0, True, float("nan")
    blocks = x[len(x) - n:].reshape(nblocks, -1).mean(axis=1)
    mean = float(blocks.mean())
    err = float(T95 * blocks.std(ddof=1) / np.sqrt(nblocks))
    drift = float(blocks[-1] - blocks[0])
    return mean, err, len(x), bool(abs(drift) > max(err, 1e-12)), drift


def lammps_kspace(log):
    t = Path(log).read_text()
    g = re.search(r"G vector \(1/distance\)\s*=\s*([\d.eE+-]+)", t)
    # LAMMPS prints "KSpace vectors: actual max1d max3d = <nvectors> <kmax1d> <kmax3d>"
    k = re.search(r"KSpace vectors: actual max1d max3d\s*=\s*(\d+)\s+(\d+)", t)
    nvec = int(k.group(1)) if k else None
    kmax1d = int(k.group(2)) if k else None
    return (float(g.group(1)) if g else None, kmax1d, nvec)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--equil", type=int, default=1500)
    a = ap.parse_args(argv)
    raspa = pd.read_csv(ROOT / "results" / "isotherm_MIL53_lp_CO2_304K.csv", comment="#")
    rows = []
    for d in sorted(GCMC.glob("p_*bar"), key=lambda p: float(p.name[2:-3])):
        log = d / "gcmc.log"
        if not log.exists():
            continue
        # LAMMPS buffers its log file; while a run is in flight the captured stdout
        # has more rows. Use whichever currently holds more thermo output.
        screen = d / "screen.txt"
        if screen.exists() and len(read_thermo(screen)) > len(read_thermo(log)):
            log_for_rows = screen
        else:
            log_for_rows = log
        p_bar = float(d.name[2:-3])
        th = read_thermo(log_for_rows)
        prod = th[th.step > a.equil]
        if prod.empty:
            print(f"  ({d.name}: no production rows yet, {len(th)} thermo rows)")
            continue
        mean, err, nsamp, drifting, drift = block_stats(prod.nuc)
        r = raspa[np.isclose(raspa.p_bar, p_bar)].iloc[0]
        delta = mean - r.absolute_molec_uc
        tol = float(np.hypot(err, r.absolute_molec_uc_err95))
        gvec, kmax1d, kmax3d = lammps_kspace(log)
        rows.append({"p_bar": p_bar,
                     "RASPA_molec_uc": r.absolute_molec_uc, "RASPA_err95": r.absolute_molec_uc_err95,
                     "LAMMPS_molec_uc": mean, "LAMMPS_err95": err, "LAMMPS_samples": nsamp,
                     "delta": delta, "tolerance": tol,
                     "pass": bool(abs(delta) < tol and not drifting),
                     "drift_last_minus_first_block": drift, "not_equilibrated": drifting,
                     "RASPA_mol_kg": r.absolute_mol_kg, "LAMMPS_mol_kg": mean * 1000 / 832.415,
                     "LAMMPS_G_per_A": gvec, "LAMMPS_kmax1d": kmax1d, "LAMMPS_kvectors": kmax3d,
                     "RASPA_alpha_per_A": 0.265058, "RASPA_kvec": "7 9 7"})
    if not rows:
        sys.exit("no LAMMPS production data yet")
    df = pd.DataFrame(rows)
    (ROOT / "results").mkdir(exist_ok=True)
    df.to_csv(ROOT / "results" / "crosscheck_gcmc.csv", index=False)
    show = ["p_bar", "RASPA_molec_uc", "RASPA_err95", "LAMMPS_molec_uc", "LAMMPS_err95",
            "delta", "tolerance", "pass", "drift_last_minus_first_block", "not_equilibrated"]
    with pd.option_context("display.width", 220, "display.float_format", "{:.4f}".format):
        print(df[show].to_string(index=False))
    print(f"\nEwald: RASPA alpha = 0.265058 A^-1, kvec 7 9 7; "
          f"LAMMPS G = {df.LAMMPS_G_per_A.iloc[0]} A^-1, kmax1d = {df.LAMMPS_kmax1d.iloc[0]}, "
          f"{df.LAMMPS_kvectors.iloc[0]} vectors")
    for _, r in df[df.not_equilibrated].iterrows():
        print(f"  NOT EQUILIBRATED at {r.p_bar:g} bar: last block - first block = "
              f"{r.drift_last_minus_first_block:+.3f} molec/uc vs error bar {r.LAMMPS_err95:.3f}. "
              "The mean is biased and the error bar understates the uncertainty; this is a sampling "
              "problem, not a code disagreement.")
    print(f"{'PASS' if df['pass'].all() else 'FAIL'}: {int(df['pass'].sum())}/{len(df)} points "
          f"within tolerance and equilibrated")
    print("wrote results/crosscheck_gcmc.csv")
    return 0 if df["pass"].all() else 1


if __name__ == "__main__":
    sys.exit(main())
