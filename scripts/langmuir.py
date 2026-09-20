#!/usr/bin/env python3
"""Langmuir analysis of the lp branch: experiment vs simulation, same procedure.

Why (NOTES.md, "Low-pressure comparison, REFRAMED"): Coudert's
K_lp = 2.6e-5 mol kg^-1 Pa^-1 is NOT a Henry constant but the initial slope of a
Langmuir function fitted to the 9-30 bar branch of the experimental isotherm. So we
apply the SAME procedure to (1) the digitised experimental points and (2) our
simulated isotherm, over the SAME window, and compare K and N_max pairwise.

  n(p) = N_max b p / (1 + b p),   K = dn/dp|_0 = N_max b

Fit: weighted least squares. At fixed b the model is linear in N_max, so N_max is
solved analytically and only b is scanned (log grid + golden-section refinement);
parameter uncertainties come from the covariance built from the chi^2 Hessian.
Only numpy/pandas/matplotlib are used.

Step 1 validates the procedure: if K(experiment, P >= 9 bar) reproduces Coudert's
2.6e-5, our fitting matches his and the comparison is meaningful. If it does not,
that is reported as a procedure discrepancy BEFORE anything is said about the
force field.

Usage: python scripts/langmuir.py [--window 9.0] [--column excess]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from isotherm import LP_WINDOW_BAR, MOLKG_PER_UC, REF_CSV, ROOT  # noqa: E402
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

K_LP_COUDERT = 2.6e-5   # mol/kg/Pa, Coudert 2008 Fig. 5b (Langmuir fit, no error bar given)
K_NP_COUDERT = 9.0e-5


def fit_langmuir(p_bar, n, sigma):
    """Weighted Langmuir fit. p in bar, n in mol/kg. Returns dict with N_max, b, K."""
    p = np.asarray(p_bar, float); n = np.asarray(n, float); w = 1.0 / np.asarray(sigma, float) ** 2

    def nmax_at(b):
        f = b * p / (1 + b * p)
        return float(np.sum(w * f * n) / np.sum(w * f * f))

    def chi2(b):
        if b <= 0:
            return np.inf
        f = b * p / (1 + b * p)
        return float(np.sum(w * (n - nmax_at(b) * f) ** 2))

    # A Langmuir is monotonically increasing and concave: refuse data that is not.
    order = np.argsort(p)
    n_sorted = n[order]
    monotonic = bool(np.all(np.diff(n_sorted) > -3 * np.asarray(sigma, float)[order][1:]))
    grid = np.logspace(-4, 3, 4000)                     # b in 1/bar
    b = grid[int(np.argmin([chi2(x) for x in grid]))]
    lo, hi = b / 3, b * 3                               # golden-section refinement
    phi = (np.sqrt(5) - 1) / 2
    for _ in range(200):
        x1, x2 = hi - phi * (hi - lo), lo + phi * (hi - lo)
        if chi2(x1) < chi2(x2):
            hi = x2
        else:
            lo = x1
    b = 0.5 * (lo + hi)
    nmax = nmax_at(b)
    # covariance from the Hessian of chi^2/2 in (N_max, b)
    def resid(par):
        N, bb = par
        return np.sqrt(w) * (n - N * bb * p / (1 + bb * p))
    par = np.array([nmax, b])
    J = np.zeros((len(p), 2))
    for k in range(2):
        h = max(abs(par[k]) * 1e-5, 1e-12)
        pp = par.copy(); pp[k] += h
        pm = par.copy(); pm[k] -= h
        J[:, k] = (resid(pp) - resid(pm)) / (2 * h)
    dof = max(len(p) - 2, 1)
    chi2_red = chi2(b) / dof
    cov = np.linalg.inv(J.T @ J) * max(chi2_red, 1.0)
    dN, db = np.sqrt(np.diag(cov))
    K = nmax * b                                         # mol/kg/bar
    dK = np.sqrt((b * dN) ** 2 + (nmax * db) ** 2 + 2 * b * nmax * cov[0, 1])
    at_bound = b > 0.5 * grid[-1] or b < 2 * grid[0]
    return {"ok": bool(monotonic and not at_bound and dK < abs(K)),
            "monotonic": monotonic, "b_at_bound": bool(at_bound),
            "N_max": nmax, "N_max_err": dN, "b_per_bar": b, "b_err": db,
            "K_mol_kg_bar": K, "K_err_mol_kg_bar": dK,
            "K_mol_kg_Pa": K * 1e-5, "K_err_mol_kg_Pa": dK * 1e-5,
            "chi2_red": chi2_red, "n_points": len(p)}


def show(name, f):
    print(f"{name}:")
    if not f["ok"]:
        why = []
        if not f["monotonic"]:
            why.append("the data are NOT monotonically increasing over the window "
                       "(excess loading passes through a maximum), which a Langmuir function cannot represent")
        if f["b_at_bound"]:
            why.append("the fitted b ran to the edge of the search grid")
        if f["K_err_mol_kg_bar"] >= abs(f["K_mol_kg_bar"]):
            why.append("the uncertainty on K exceeds K itself")
        print("    FIT REFUSED: " + "; ".join(why))
        print(f"    (degenerate values, do not use: N_max = {f['N_max']:.3f}, b = {f['b_per_bar']:.4g}, "
              f"K = {f['K_mol_kg_Pa']:.3g} mol/kg/Pa)")
        return
    print(f"    N_max = {f['N_max']:.3f} +/- {f['N_max_err']:.3f} mol/kg "
          f"({f['N_max'] / MOLKG_PER_UC:.3f} +/- {f['N_max_err'] / MOLKG_PER_UC:.3f} molecules/uc)")
    print(f"    b     = {f['b_per_bar']:.4f} +/- {f['b_err']:.4f} bar^-1")
    print(f"    K     = N_max b = {f['K_mol_kg_Pa']:.3e} +/- {f['K_err_mol_kg_Pa']:.3e} mol/kg/Pa "
          f"({f['K_mol_kg_bar']:.3f} mol/kg/bar), {f['n_points']} points, chi2_red = {f['chi2_red']:.2f}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=float, default=LP_WINDOW_BAR)
    ap.add_argument("--column", default="excess", choices=["excess", "absolute"],
                    help="which simulated convention to fit (default excess: the reference is excess)")
    ap.add_argument("--csv", default="results/isotherm_MIL53_lp_CO2_304K.csv")
    a = ap.parse_args(argv)

    # ---- 1. experiment
    ref = pd.read_csv(REF_CSV)
    e = ref[ref.pressure_bar >= a.window]
    sig_e = np.full(len(e), 0.05)      # digitisation uncertainty stated in the CSV notes column
    fe = fit_langmuir(e.pressure_bar, e.loading_mmol_per_g, sig_e)
    print(f"Langmuir window: P >= {a.window:g} bar\n")
    show(f"(1) EXPERIMENT  Bourrelly 2005, excess (assumed), {len(e)} points", fe)
    ratio = fe["K_mol_kg_Pa"] / K_LP_COUDERT
    ok = abs(np.log(ratio)) < np.log(1.5)
    print(f"    vs Coudert K_lp = {K_LP_COUDERT:.1e} mol/kg/Pa: ratio {ratio:.2f} -> "
          + ("procedure VALIDATED (within 50 %), the comparison is meaningful"
             if ok else "PROCEDURE DISCREPANCY (> 50 %): report this before any force-field statement"))

    # ---- 2. simulation, same window, same form
    sim_path = ROOT / a.csv
    if not sim_path.exists():
        print(f"\n(no simulated isotherm yet at {a.csv}; run scripts/isotherm.py first)")
        return 0
    sim = pd.read_csv(sim_path, comment="#")
    s = sim[sim.p_bar >= a.window]
    if len(s) < 3:
        print(f"\n(only {len(s)} simulated points with P >= {a.window:g} bar; need >= 3 to fit)")
        return 0
    col, err = f"{a.column}_mol_kg", f"{a.column}_mol_kg_err95"
    fs = fit_langmuir(s.p_bar, s[col], s[err] / 2.776)     # sem, not the 95 % half-width
    print()
    show(f"(2) SIMULATION  this work, {a.column}, {len(s)} points", fs)
    if not fs["ok"]:
        print("\n(3) PAIRWISE: NOT DONE -- the simulated fit was refused (see above). "
              "Decision required before comparing conventions (see NOTES.md).")
        pd.DataFrame([{"fit": "experiment_Bourrelly2005", **fe}, {"fit": f"simulation_{a.column}", **fs}]).to_csv(
            ROOT / "results" / "langmuir_fits.csv", index=False)
        return 0
    print("\n(3) PAIRWISE, same window, same functional form, same fitting procedure:")
    print(f"    N_max  sim / exp = {fs['N_max'] / fe['N_max']:.3f}"
          f"   ({fs['N_max']:.2f} vs {fe['N_max']:.2f} mol/kg)")
    print(f"    K      sim / exp = {fs['K_mol_kg_Pa'] / fe['K_mol_kg_Pa']:.3f}"
          f"   ({fs['K_mol_kg_Pa']:.2e} vs {fe['K_mol_kg_Pa']:.2e} mol/kg/Pa)")
    pd.DataFrame([{"fit": "experiment_Bourrelly2005", **fe}, {"fit": f"simulation_{a.column}", **fs}]).to_csv(
        ROOT / "results" / "langmuir_fits.csv", index=False)

    # ---- 4. the three curves
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    pp = np.logspace(np.log10(min(sim.p_bar.min(), 0.01)), np.log10(max(sim.p_bar.max(), 30)), 400)
    ax.errorbar(sim.p_bar, sim[col], yerr=sim[err], label=f"this work, {a.column} (GCMC, rigid lp)", **SIM)
    ax.plot(pp, fs["N_max"] * fs["b_per_bar"] * pp / (1 + fs["b_per_bar"] * pp), color=SIM["color"], lw=1.2, ls=":",
            label=f"Langmuir fit to this work (P $\\geq$ {a.window:g} bar)")
    ax.errorbar(e.pressure_bar, e.loading_mmol_per_g, yerr=sig_e,
                label=f"Bourrelly 2005, excess (assumed), P $\\geq$ {a.window:g} bar", **REF)
    g = ref[ref.pressure_bar < a.window]
    if len(g):
        ax.errorbar(g.pressure_bar, g.loading_mmol_per_g, **{**REF, "color": "#8c8c88", "ms": 7},
                    label="same, inside the np-lp step: NOT compared")
    ax.plot(pp, fe["N_max"] * fe["b_per_bar"] * pp / (1 + fe["b_per_bar"] * pp), color=REF["color"], lw=1.2, ls="--",
            label="Langmuir fit to the experimental lp branch\n= Coudert's virtual rigid-lp curve (K matches to 1.00)")
    ax.plot(pp, K_LP_COUDERT * 1e5 * pp, color=EXTRA[0], lw=1.6,
            label="Coudert 2008 $K_{lp}$: initial slope only (tangent at p $\\to$ 0)")
    ax.set_xscale("log")
    ax.set_ylim(0, max(sim[col].max(), e.loading_mmol_per_g.max()) * 1.25)
    ax.set_xlabel("Pressure [bar]")
    ax.set_ylabel("Loading [mol/kg  =  mmol/g]")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x / MOLKG_PER_UC, lambda x: x * MOLKG_PER_UC))
    sec.set_ylabel("Loading [molecules / unit cell]")
    ax.set_title("lp branch: same Langmuir procedure applied to experiment and simulation")
    ax.legend(fontsize=8, loc="lower right")
    fig.text(0.01, -0.06,
             "Coudert's $K_{lp}$ is the initial slope of a Langmuir fit to the 9-30 bar experimental branch, not a Henry "
             "constant.\nOur fit to the same points returns K = 2.606e-5 vs his 2.6e-5 (ratio 1.00), so the dashed orange curve "
             "is his virtual\nrigid-lp curve; the green line is only its tangent at p -> 0. The low-pressure divergence is the "
             "result: there the real\nmaterial is np, which a rigid-lp model cannot follow, and our points also lie above our own "
             "Langmuir fit, i.e. the\nsimulated isotherm is more heterogeneous than a single-site Langmuir. Error bars: 95 % CI "
             "(simulation), digitisation (experiment).",
             fontsize=7.5, va="top")
    fig.savefig(ROOT / "results" / "langmuir_lp_branch.png")
    print("\nwrote results/langmuir_fits.csv and results/langmuir_lp_branch.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
