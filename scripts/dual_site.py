#!/usr/bin/env python3
"""Dual-site Langmuir fits, and the analytical osmotic construction built on them.

This is the SECONDARY route of the plan: a cross-check on the numerical integration of
eq. 8 (`scripts/osmotic_numeric.py`). A dual-site Langmuir

    n(p) = N1 b1 p / (1 + b1 p) + N2 b2 p / (1 + b2 p)

keeps an analytic integral -- int_0^P n V_m dp = sum_i N_i RT ln(1 + b_i P) for an ideal
gas, exactly Coudert eq. 11 summed over two sites -- while allowing the energetic
heterogeneity that makes a single site fail here (chi2_red 318 np, 1368 lp). If the
dual-site construction agrees with the numerical one, the numerical route is not an
artefact of the interpolation; if it does not, the disagreement measures how much the
functional form still matters.

Fitting: at fixed (b1, b2) the model is LINEAR in (N1, N2), so those are solved exactly
by weighted least squares and only the two b's are searched -- the same trick
scripts/langmuir.py uses for one site. Coarse log grid, then Nelder-Mead refinement
written out here (scipy is not a dependency).

The TERTIARY single-site route is reported alongside, labelled as the poor fit it is.

Usage: python scripts/dual_site.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from langmuir import fit_langmuir  # noqa: E402
from osmotic import MOLKG_TO_PERCELL, NA, RT, V_LP, V_NP  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def solve_N(p, n, w, b1, b2):
    """Exact weighted least squares for (N1, N2) at fixed (b1, b2); N_i >= 0."""
    f1 = b1 * p / (1 + b1 * p)
    f2 = b2 * p / (1 + b2 * p)
    A = np.array([[np.sum(w * f1 * f1), np.sum(w * f1 * f2)],
                  [np.sum(w * f1 * f2), np.sum(w * f2 * f2)]])
    y = np.array([np.sum(w * f1 * n), np.sum(w * f2 * n)])
    try:
        N = np.linalg.solve(A, y)
    except np.linalg.LinAlgError:
        return None
    if np.any(N < 0):                       # clamp to the physical boundary and refit
        N = np.maximum(N, 0.0)
    return N


def chi2(p, n, w, b1, b2):
    if b1 <= 0 or b2 <= 0:
        return np.inf
    N = solve_N(p, n, w, b1, b2)
    if N is None:
        return np.inf
    pred = N[0] * b1 * p / (1 + b1 * p) + N[1] * b2 * p / (1 + b2 * p)
    return float(np.sum(w * (n - pred) ** 2))


def nelder_mead(fn, x0, step=0.4, iters=400):
    """Minimise fn over log-parameters. Small, self-contained, no scipy."""
    x0 = np.asarray(x0, float)
    sim = [x0] + [x0 + step * np.eye(len(x0))[i] for i in range(len(x0))]
    sim = np.array(sim)
    val = np.array([fn(s) for s in sim])
    for _ in range(iters):
        o = np.argsort(val)
        sim, val = sim[o], val[o]
        cen = sim[:-1].mean(axis=0)
        xr = cen + (cen - sim[-1]); fr = fn(xr)
        if fr < val[0]:
            xe = cen + 2 * (cen - sim[-1]); fe = fn(xe)
            sim[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            sim[-1], val[-1] = xr, fr
        else:
            xc = cen + 0.5 * (sim[-1] - cen); fc = fn(xc)
            if fc < val[-1]:
                sim[-1], val[-1] = xc, fc
            else:
                sim[1:] = sim[0] + 0.5 * (sim[1:] - sim[0])
                val[1:] = [fn(s) for s in sim[1:]]
    o = np.argsort(val)
    return sim[o][0], val[o][0]


def fit_dual(p_bar, n_molkg, sigma):
    p = np.asarray(p_bar, float)
    n = np.asarray(n_molkg, float)
    w = 1.0 / np.asarray(sigma, float) ** 2
    best = (np.inf, None)
    grid = np.logspace(-3, 3, 25)
    for i, b1 in enumerate(grid):
        for b2 in grid[i + 1:]:
            c = chi2(p, n, w, b1, b2)
            if c < best[0]:
                best = (c, (b1, b2))
    x, c = nelder_mead(lambda z: chi2(p, n, w, np.exp(z[0]), np.exp(z[1])),
                       np.log(np.array(best[1])))
    b1, b2 = np.exp(x)
    if b1 > b2:                              # site 1 is the weaker-binding one
        b1, b2 = b2, b1
    N = solve_N(p, n, w, b1, b2)
    dof = max(len(p) - 4, 1)
    return {"N1": N[0], "b1": b1, "N2": N[1], "b2": b2,
            "N_total": N[0] + N[1], "chi2_red": c / dof,
            "K_mol_kg_Pa": (N[0] * b1 + N[1] * b2) / 1e5}


def omega_dual(p_pa, par, f_host, v_cell_A3):
    """Analytic Omega for a dual-site Langmuir, ideal-gas V_m (Coudert eq. 11 x 2)."""
    ads = np.zeros_like(p_pa)
    for Nm, b in ((par["N1_cell"], par["b1_pa"]), (par["N2_cell"], par["b2_pa"])):
        ads = ads + Nm * RT * np.log1p(b * p_pa)
    pv = p_pa * v_cell_A3 * 1e-30 * NA / 1000.0
    return f_host + pv - ads


def crossings(p_pa, d):
    out = []
    s = np.sign(d)
    for i in np.where(np.diff(s) != 0)[0]:
        f = d[i] / (d[i] - d[i + 1])
        out.append((p_pa[i] * (p_pa[i + 1] / p_pa[i]) ** f / 1e5,
                    "lp -> np" if d[i] > 0 else "np -> lp"))
    return out


def as_cell(fit):
    return {"N1_cell": fit["N1"] * MOLKG_TO_PERCELL, "b1_pa": fit["b1"] / 1e5,
            "N2_cell": fit["N2"] * MOLKG_TO_PERCELL, "b2_pa": fit["b2"] / 1e5}


def main():
    branches = {
        "np CO2 cell": ROOT / "results/isotherm_MIL53_npCO2cell_CO2_304K.csv",
        "lp SABVUN-DDEC": ROOT / "results/isotherm_MIL53_lp_CO2_304K.csv",
        "lp PACMAN": ROOT / "results/isotherm_MIL53_lp_pacman_CO2_304K.csv",
    }
    print("=" * 94)
    print("SECONDARY -- dual-site Langmuir fits (TERTIARY single-site alongside)")
    print("=" * 94)
    fits = {}
    for name, csv in branches.items():
        d = pd.read_csv(csv, comment="#")
        sig = d.absolute_mol_kg_err95 / 2.776
        du = fit_dual(d.p_bar, d.absolute_mol_kg, sig)
        si = fit_langmuir(d.p_bar, d.absolute_mol_kg, sig)
        fits[name] = (du, si)
        print(f"\n  {name}  ({len(d)} points)")
        print(f"    dual-site  : N1 {du['N1']:7.3f} mol/kg, b1 {du['b1']:9.4g} /bar | "
              f"N2 {du['N2']:7.3f}, b2 {du['b2']:9.4g} /bar")
        print(f"                 N_total {du['N_total']:.3f} mol/kg = "
              f"{du['N_total'] * MOLKG_TO_PERCELL:.3f} CO2/cell, K {du['K_mol_kg_Pa']:.4e} mol/kg/Pa, "
              f"chi2_red {du['chi2_red']:.2f}")
        print(f"    single-site: N_max {si['N_max']:7.3f} mol/kg = "
              f"{si['N_max'] * MOLKG_TO_PERCELL:.3f} CO2/cell, K {si['K_mol_kg_Pa']:.4e}, "
              f"chi2_red {si['chi2_red']:.2f}   <-- POOR FIT")
        print(f"    improvement in chi2_red: {si['chi2_red'] / du['chi2_red']:.1f}x")

    print("\n" + "=" * 94)
    print("Analytical construction on the DUAL-SITE fits (cross-check on the numerical route)")
    print("=" * 94)
    p = np.logspace(2, 7, 20000)
    np_par = as_cell(fits["np CO2 cell"][0])
    for lp_name in ("lp PACMAN", "lp SABVUN-DDEC"):
        lp_par = as_cell(fits[lp_name][0])
        print(f"\n  {lp_name} + np CO2 cell")
        for dF in (1.0, 1.5, 2.0, 2.5, 3.0):
            o_lp = omega_dual(p, lp_par, 0.0, V_LP)
            o_np = omega_dual(p, np_par, dF, V_NP)
            cr = crossings(p, o_np - o_lp)
            print(f"    dF={dF:.1f}: " + (", ".join(f"{q:.3f} bar ({w})" for q, w in cr)
                                          if cr else "no crossing"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
