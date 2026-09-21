#!/usr/bin/env python3
"""Single-point energy comparison RASPA vs LAMMPS on ONE identical configuration.

Configuration: the rigid MIL-53(Al) lp framework (4x2x2 = 1216 atoms) plus 8 CO2
molecules whose coordinates were taken from a RASPA restart file (12 decimals),
so both codes see exactly the same positions.

Term extraction:
  RASPA  prints the decomposition directly (Host/Adsorbate and Adsorbate/Adsorbate,
         VDW and charge-charge Real/Fourier), plus ONE tail number for the whole
         system; the host-host part of the tail is removed by subtracting a second
         RASPA run with zero molecules.
  LAMMPS gets three single-point runs, A = framework + CO2, B = framework only,
         C = CO2 only, and every cross term is A - B - C. This is exact for the
         pair terms and for the Ewald reciprocal term (a quadratic form in the
         charges: the self terms cancel), and it removes the framework-framework
         energy that RASPA never computes for a rigid host.

Two conventions that must be handled, both verified numerically here:
  * LAMMPS E_vdwl INCLUDES the tail correction (E_tail is printed separately for
    information); RASPA reports its VDW energy without the tail. Verified by
    re-running with `pair_modify tail no`.
  * The real-space / reciprocal-space split of the Coulomb energy depends on the
    Ewald convergence parameter, which the two codes choose by different criteria.
    Only the SUM is physically meaningful, so the sum is what is compared.

Usage: python scripts/crosscheck_energy.py
"""
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CROSS = ROOT / "runs" / "crosscheck"
K_PER_KCAL = 503.2197     # 1 kcal/mol in K


def raspa_terms(data_file):
    t = Path(data_file).read_text()
    t = t[:t.index("Starting simulation")]          # the initial (= only) energy report

    def grab(label):
        m = re.search(rf"^\s*{re.escape(label)}\s+([-\d.]+)", t, re.MULTILINE)
        return float(m.group(1)) if m else None

    return {
        "HA_vdw": grab("Host/Adsorbate VDW energy:"),
        "HA_coul_real": grab("Host/Adsorbate charge-charge Real energy:"),
        "HA_coul_four": grab("Host/Adsorbate charge-charge Fourier energy:"),
        "AA_vdw": grab("Adsorbate/Adsorbate VDW energy:"),
        "AA_coul_real": grab("Adsorbate/Adsorbate charge-charge Real energy:"),
        "AA_coul_four": grab("Adsorbate/Adsorbate charge-charge Fourier energy:"),
        "tail_total": grab("Tail-correction energy:"),
        "alpha": float(re.search(r"Alpha convergence parameter\s*:\s*([\d.]+)", t).group(1)),
        "kvec": re.search(r"kvec \(x,y,z\)\s*:\s*(.+)", t).group(1).strip(),
    }


def lammps_terms(log_file):
    t = Path(log_file).read_text()
    m = re.search(r"Step\s+E_vdwl\s+E_coul\s+E_long\s+E_tail\s+PotEng\s*\n\s*0\s+"
                  r"([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", t)
    keys = ("evdwl", "ecoul", "elong", "etail")
    out = {k: float(m.group(i + 1)) for i, k in enumerate(keys)}
    g = re.search(r"G vector \(1/distance\)\s*=\s*([\d.eE+-]+)\s*\n\s*estimated absolute RMS force accuracy\s*=\s*([\d.eE+-]+)", t)
    if g:
        out["alpha"] = float(g.group(1))
        out["rms_force"] = float(g.group(2))
    k = re.search(r"KSpace vectors: actual max1d max3d\s*=\s*(\d+)\s+(\d+)", t)
    if k:
        out["kmax1d"] = int(k.group(2))
    return out


def main():
    r8 = raspa_terms(next((CROSS / "raspa_singlepoint" / "Output" / "System_0").glob("*.data")))
    r0 = raspa_terms(next((CROSS / "raspa_singlepoint_empty" / "Output" / "System_0").glob("*.data")))
    A = lammps_terms(CROSS / "lammps" / "A_framework_co2.log")
    B = lammps_terms(CROSS / "lammps" / "B_framework.log")
    C = lammps_terms(CROSS / "lammps" / "C_co2.log")

    # LAMMPS cross terms (A - B - C), converted kcal/mol -> K
    x = {k: (A[k] - B[k] - C[k]) * K_PER_KCAL for k in ("evdwl", "ecoul", "elong", "etail")}
    gg = {k: C[k] * K_PER_KCAL for k in ("evdwl", "ecoul", "elong", "etail")}
    rows = [
        ("host-guest LJ (no tail)", r8["HA_vdw"], x["evdwl"] - x["etail"]),
        ("host-guest Coulomb, real space", r8["HA_coul_real"], x["ecoul"]),
        ("host-guest Coulomb, reciprocal", r8["HA_coul_four"], x["elong"]),
        ("host-guest Coulomb, TOTAL", r8["HA_coul_real"] + r8["HA_coul_four"], x["ecoul"] + x["elong"]),
        ("guest-guest LJ (no tail)", r8["AA_vdw"], gg["evdwl"] - gg["etail"]),
        ("guest-guest Coulomb, TOTAL", r8["AA_coul_real"] + r8["AA_coul_four"], gg["ecoul"] + gg["elong"]),
        ("tail correction, host-guest + guest-guest", r8["tail_total"] - r0["tail_total"], x["etail"] + gg["etail"]),
    ]
    df = pd.DataFrame(rows, columns=["term", "RASPA [K]", "LAMMPS [K]"])
    df["diff [K]"] = df["LAMMPS [K]"] - df["RASPA [K]"]
    df["rel. diff"] = df["diff [K]"] / df["RASPA [K]"].abs()
    tot = df[df.term.str.contains("TOTAL|LJ|tail")]
    total_r, total_l = tot["RASPA [K]"].sum(), tot["LAMMPS [K]"].sum()
    df.loc[len(df)] = ["TOTAL interaction energy (LJ + Coulomb + tail)", total_r, total_l,
                       total_l - total_r, (total_l - total_r) / abs(total_r)]
    with pd.option_context("display.width", 200, "display.float_format", "{:.4f}".format):
        print(df.to_string(index=False))
    print(f"\nEwald: RASPA alpha = {r8['alpha']} A^-1, kvec = {r8['kvec']};  "
          f"LAMMPS G = {A.get('alpha')} A^-1, kmax1d = {A.get('kmax1d')}, "
          f"RMS force accuracy = {A.get('rms_force')}")
    print("The real/reciprocal split differs because the two codes pick the Ewald convergence")
    print("parameter by different criteria; only the sums are physically meaningful.")
    (ROOT / "results").mkdir(exist_ok=True)
    df.to_csv(ROOT / "results" / "crosscheck_single_point_energies.csv", index=False)
    print("\nwrote results/crosscheck_single_point_energies.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
