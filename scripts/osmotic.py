#!/usr/bin/env python3
"""Osmotic-ensemble construction for the MIL-53(Al) breathing transition.

Coudert et al., JACS 2008, 130, 14294, eq. 8 and 11. Per unit cell:

    Omega_os(P) = F_host + P V_host - int_0^P N_ads(p) V_m(p) dp

with N_ads from a single-site Langmuir fit of each rigid-phase isotherm. For an
ideal gas, V_m = RT/p, and the integral is analytic (eq. 11):

    int_0^P N_ads V_m dp = N_max R T ln(1 + b P)

Both are computed here: the analytic ideal-gas form and a numerical integral
using the Peng-Robinson molar volume, so the ideal-gas assumption can be judged
rather than assumed.

The stable phase minimises Omega_os; crossings of Omega_lp and Omega_np are the
transition pressures.

VALIDATION FIRST: with Coudert's experimental Langmuir parameters the script must
reproduce his transitions (lp->np ~0.3 bar, np->lp ~5-6 bar). Only then are our
simulated parameters used. Nothing here is tuned to reproduce a target: every
input is stated, and the sensitive ones (Delta F_host, N_max of the np phase) are
scanned.

Usage: python scripts/osmotic.py [--dF 2.5] [--scan 1 4]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from langmuir import fit_langmuir  # noqa: E402
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R = 8.314462618e-3          # kJ/mol/K
T = 304.0
RT = R * T
NA = 6.02214076e23
M_CELL = 832.415            # g/mol, identical for both phases (Z = 4)
MOLKG_TO_PERCELL = M_CELL / 1000.0
V_LP = 1412.0               # A^3 per formula cell, structures/MIL-53_Al_lp.cif
V_NP = 945.8                # A^3 per formula cell, structures/MIL-53_Al_np.cif (1891.6 / 2)
# CO2 critical constants (Span & Wagner), as used by RASPA's Peng-Robinson
TC, PC, OMEGA = 304.1282, 7377300.0, 0.22394
R_SI = 8.314462618          # J/mol/K

# --- Coudert 2008, experimental Langmuir parameters (his Fig. 5b analysis)
K_LP_EXP = 2.6e-5           # mol/kg/Pa
K_NP_EXP = 9.0e-5           # mol/kg/Pa
# N_max of the np branch is NOT in the material we have; Bourrelly 2005 reports
# 0.75 CO2 per structural OH before the step, i.e. 3.0 CO2 per cell (4 OH). This
# is an ASSUMPTION, scanned below, not a fitted quantity.
NMAX_NP_EXP = 3.0           # molecules per cell


def peng_robinson_vm(p_pa, temp=T):
    """Molar volume [m^3/mol] of CO2 from Peng-Robinson (vapour root)."""
    kappa = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA ** 2
    alpha = (1 + kappa * (1 - np.sqrt(temp / TC))) ** 2
    a = 0.45724 * R_SI ** 2 * TC ** 2 / PC * alpha
    b = 0.07780 * R_SI * TC / PC
    out = np.empty_like(np.atleast_1d(p_pa), dtype=float)
    for i, p in enumerate(np.atleast_1d(p_pa)):
        A = a * p / (R_SI * temp) ** 2
        B = b * p / (R_SI * temp)
        roots = np.roots([1.0, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3)])
        real = roots[np.abs(roots.imag) < 1e-9].real
        out[i] = max(real) * R_SI * temp / p        # vapour root
    return out if np.ndim(p_pa) else float(out[0])


def omega_os(p_pa, f_host, v_cell_A3, nmax_percell, b_per_pa, mode="ideal", vm_grid=None):
    """Osmotic potential [kJ/mol per cell] on the pressure array p_pa.

    mode="ideal": analytic Langmuir + ideal gas (Coudert eq. 11).
    mode="pr"   : numerical integral with the Peng-Robinson molar volume. The
                  integral is accumulated once over the whole pressure array
                  (cumulative trapezoid), not recomputed per point.
    """
    p_pa = np.atleast_1d(p_pa)
    pv = p_pa * v_cell_A3 * 1e-30 * NA / 1000.0
    if mode == "ideal":
        ads = nmax_percell * RT * np.log1p(b_per_pa * p_pa)
    else:
        vm = vm_grid if vm_grid is not None else peng_robinson_vm(p_pa)
        n = nmax_percell * b_per_pa * p_pa / (1 + b_per_pa * p_pa)
        integrand = n * vm                                   # m^3/mol per cell
        ads = np.concatenate(([0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(p_pa)))) / 1000.0
    return f_host + pv - ads


def crossings(p_pa, d):
    """Pressures where the sign of d = Omega_np - Omega_lp changes."""
    out = []
    s = np.sign(d)
    for i in np.where(np.diff(s) != 0)[0]:
        lo, hi = p_pa[i], p_pa[i + 1]
        out.append((0.5 * (lo + hi) / 1e5, "lp -> np" if d[i] > 0 else "np -> lp"))
    return out


_VM_CACHE = {}


def transitions(par, dF, mode="ideal", p_range=(1e1, 1e8), n=20000):
    p = np.logspace(np.log10(p_range[0]), np.log10(p_range[1]), n)
    vm = None
    if mode == "pr":
        key = (p_range, n)
        if key not in _VM_CACHE:          # solve the cubic once per grid, then interpolate
            coarse = np.logspace(np.log10(p_range[0]), np.log10(p_range[1]), 400)
            _VM_CACHE[key] = np.interp(p, coarse, peng_robinson_vm(coarse))
        vm = _VM_CACHE[key]
    o_lp = omega_os(p, 0.0, V_LP, par["nmax_lp"], par["b_lp"], mode, vm)
    o_np = omega_os(p, dF, V_NP, par["nmax_np"], par["b_np"], mode, vm)
    return p, o_lp, o_np, crossings(p, o_np - o_lp)


def langmuir_from_isotherm(csv, col="absolute_mol_kg", per_cell_factor=MOLKG_TO_PERCELL, pmax_bar=None):
    """Fit our simulated isotherm; return N_max [per cell] and b [1/Pa]."""
    df = pd.read_csv(csv, comment="#")
    if pmax_bar:
        df = df[df.p_bar <= pmax_bar]
    f = fit_langmuir(df.p_bar, df[col], df[col.replace("mol_kg", "mol_kg_err95")] / 2.776)
    return {"nmax": f["N_max"] * per_cell_factor, "b": f["b_per_bar"] / 1e5, "fit": f}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dF", type=float, default=2.5, help="F_np - F_lp [kJ/mol per cell]")
    ap.add_argument("--scan", type=float, nargs=2, default=[1.0, 4.0])
    a = ap.parse_args(argv)

    # ---------------- 1. validation with Coudert's experimental parameters
    nmax_lp_exp_molkg = 12.183                      # our Langmuir fit to the digitised lp branch
    exp = {"nmax_lp": nmax_lp_exp_molkg * MOLKG_TO_PERCELL,
           "b_lp": K_LP_EXP / nmax_lp_exp_molkg,
           "nmax_np": NMAX_NP_EXP,
           "b_np": K_NP_EXP / (NMAX_NP_EXP / MOLKG_TO_PERCELL)}
    print("=" * 78)
    print("1. VALIDATION -- Coudert's experimental Langmuir parameters")
    print(f"   lp: N_max = {exp['nmax_lp']:.2f} molec/cell (our fit to the digitised branch), "
          f"b = {exp['b_lp']:.3e} /Pa, K = {K_LP_EXP:.1e} mol/kg/Pa")
    print(f"   np: N_max = {exp['nmax_np']:.2f} molec/cell (ASSUMED: Bourrelly 0.75 CO2 per OH), "
          f"b = {exp['b_np']:.3e} /Pa, K = {K_NP_EXP:.1e} mol/kg/Pa")
    for mode in ("ideal", "pr"):
        _, _, _, cr = transitions(exp, a.dF, mode)
        got = ", ".join(f"{p:.3f} bar ({w})" for p, w in cr)
        print(f"   dF = {a.dF} kJ/mol, {mode:5s} V_m -> {got if cr else 'no crossing'}")
    print("   Coudert 2008: lp->np 0.3 bar (microcalorimetry 0.25), np->lp 5-6 bar")

    print("\n   sensitivity (nothing is tuned to the target; both inputs are stated):")
    rows = []
    for dF in np.arange(a.scan[0], a.scan[1] + 0.01, 0.5):
        for nnp in (2.5, 3.0, 3.5, 4.0):
            e = dict(exp, nmax_np=nnp, b_np=K_NP_EXP / (nnp / MOLKG_TO_PERCELL))
            cr = transitions(e, dF)[3]
            rows.append({"dF_kJ_mol": dF, "Nmax_np": nnp,
                         "lp->np_bar": next((p for p, w in cr if w == "lp -> np"), np.nan),
                         "np->lp_bar": next((p for p, w in cr if w == "np -> lp"), np.nan)})
    scan = pd.DataFrame(rows)
    with pd.option_context("display.float_format", "{:.3f}".format):
        print(scan.pivot(index="dF_kJ_mol", columns="Nmax_np", values=["lp->np_bar", "np->lp_bar"]).to_string())
    scan.to_csv(ROOT / "results" / "osmotic_sensitivity.csv", index=False)

    # ---------------- 2. our simulated parameters
    print("\n" + "=" * 78)
    print("2. OUR SIMULATED PARAMETERS")
    sim_lp = langmuir_from_isotherm(ROOT / "results" / "isotherm_MIL53_lp_CO2_304K.csv")
    print(f"   lp (rigid, DDEC): N_max = {sim_lp['nmax']:.2f} molec/cell, b = {sim_lp['b']:.3e} /Pa, "
          f"K = {sim_lp['fit']['K_mol_kg_Pa']:.3e} mol/kg/Pa, fit ok = {sim_lp['fit']['ok']}")
    np_csv = ROOT / "results" / "isotherm_MIL53_np_CO2_304K.csv"
    sets = {"experiment (Coudert)": exp,
            "simulated lp + experimental np": dict(exp, nmax_lp=sim_lp["nmax"], b_lp=sim_lp["b"])}
    if np_csv.exists():
        sim_np = langmuir_from_isotherm(np_csv)
        print(f"   np (rigid, PACMAN): N_max = {sim_np['nmax']:.2f} molec/cell, b = {sim_np['b']:.3e} /Pa, "
              f"K = {sim_np['fit']['K_mol_kg_Pa']:.3e} mol/kg/Pa, fit ok = {sim_np['fit']['ok']}")
        if not sim_np["fit"]["ok"]:
            print("   NOTE: the np Langmuir fit is degenerate (the simulated np isotherm does not "
                  "saturate); K = N_max b is still meaningful, N_max and b separately are not.")
        sets["simulated lp + simulated np"] = {"nmax_lp": sim_lp["nmax"], "b_lp": sim_lp["b"],
                                               "nmax_np": sim_np["nmax"], "b_np": sim_np["b"]}
    for name, par in sets.items():
        cr = transitions(par, a.dF)[3]
        print(f"   {name:32s} -> " + (", ".join(f"{p:.3f} bar ({w})" for p, w in cr) if cr
                                      else "NO TRANSITION: lp stable at every pressure"))

    # ---------------- 3. figure
    fig, axes = plt.subplots(1, len(sets), figsize=(5.2 * len(sets), 4.4), squeeze=False)
    for ax, (name, par) in zip(axes[0], sets.items()):
        p, o_lp, o_np, cr = transitions(par, a.dF, p_range=(1e2, 1e7))   # 0.001 - 100 bar
        ax.plot(p / 1e5, o_lp, color=SIM["color"], lw=2, label="$\\Omega_{lp}$")
        ax.plot(p / 1e5, o_np, color=REF["color"], lw=2, ls="--", label="$\\Omega_{np}$")
        for pc, w in cr:
            ax.axvline(pc, color=EXTRA[0], lw=1.4, ls=":")
            ax.annotate(f"{pc:.2f} bar  {w}", xy=(pc, 0.03), xycoords=("data", "axes fraction"),
                        fontsize=8, rotation=90, va="bottom", ha="right", color=EXTRA[0])
        for pe in (0.25, 6.0):     # Coudert/Bourrelly measured transitions, for reference only
            ax.axvline(pe, color="#8c8c88", lw=1.0, ls="-", alpha=0.6, zorder=0)
        if not cr:
            ax.annotate("no crossing:\nlp stable at every pressure", xy=(0.04, 0.10),
                        xycoords="axes fraction", fontsize=9, color=EXTRA[0])
        ax.set_xscale("log")
        ax.set_xlabel("Pressure [bar]")
        ax.set_title(name, fontsize=10)
        ax.legend(fontsize=9)
    axes[0][0].set_ylabel("$\\Omega_{os}$ [kJ/mol per unit cell]")
    axes[0][0].annotate("grey lines: measured transitions\n(0.25 and ~6 bar)", xy=(0.04, 0.88),
                        xycoords="axes fraction", fontsize=8, color="#6b6b66")
    fig.suptitle(f"Osmotic construction, CO$_2$/MIL-53(Al), {T:g} K, "
                 f"$\\Delta F_{{host}}$ = {a.dF:g} kJ/mol per cell", y=1.02)
    fig.savefig(ROOT / "results" / "osmotic_construction.png")
    print("\nwrote results/osmotic_construction.png and results/osmotic_sensitivity.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
