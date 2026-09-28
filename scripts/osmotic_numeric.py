#!/usr/bin/env python3
"""Osmotic construction by NUMERICAL integration of the simulated isotherms.

WHY THIS REPLACES THE LANGMUIR ROUTE. Coudert et al. (JACS 2008, 130, 14294) eq. 8 is

    Omega_os(P) = F_host + P V_host - int_0^P N_ads(p) V_m(p) dp

and eq. 11 is only its analytical specialisation for a single-site Langmuir isotherm,
used for his taxonomy; the method works with any adequate description of N(p). Our
simulated branches are NOT single-site Langmuir -- weighted chi2_red is 318 for the np
CO2-cell branch and 1368 for the lp branch, with systematic S-shaped residuals, because
the simulated isotherms are more energetically heterogeneous than one site (NOTES.md,
"What a rigid-framework GCMC can and cannot reproduce"). Forcing that form imposes an
N_max the data do not support, and N_max is precisely the parameter the high-pressure
upturn corrupts.

Integrating the data directly removes the problem entirely, and it also removes the
N_max question: the integral runs from 0 to the transition pressure, both transitions
lie below 10 bar, so points above 10 bar never enter. N_max mattered only because a
Langmuir extrapolates through it.

HOW N(p) IS REPRESENTED
  * between data points: PCHIP -- monotone piecewise cubic Hermite (Fritsch-Carlson) in
    log p. Monotone matters: an ordinary spline overshoots between widely spaced points
    and can make N(p) non-monotonic, which is unphysical for adsorption and would put
    spurious structure into Omega.
  * below the lowest data point: a Henry-law segment N = K' p. Its contribution is
    analytic, int_0^p0 K' p (RT/p) dp = N(p0) RT, and it is NOT negligible -- about
    1.0 kJ/mol per cell for the np branch, against a Delta F_host of 1-4 kJ/mol.
    Two anchors are carried: the Widom K_H, and the lowest isotherm point itself.
  * V_m(p): Peng-Robinson, the same routine the analytical script uses.

VALIDATION FIRST (step 2 of the plan, and it runs before anything else): synthetic
points are generated from Coudert's experimental Langmuir parameters on OUR pressure
grid, integrated by this same numerical machinery, and must reproduce the analytical
transitions of 0.273 and 3.97 bar. If they do not, the integration is wrong and nothing
downstream means anything.

Uncertainties come from a bootstrap over the loading error bars.

Usage:
  python scripts/osmotic_numeric.py [--dF 2.5] [--scan 1 4] [--boot 400]
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from osmotic import (MOLKG_TO_PERCELL, NA, RT, T, V_LP, V_NP,  # noqa: E402
                     K_LP_EXP, K_NP_EXP, NMAX_NP_EXP, peng_robinson_vm)
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------- PCHIP
def pchip_slopes(x, y):
    """Fritsch-Carlson monotone derivatives. numpy only; scipy is not a dependency."""
    h = np.diff(x)
    delta = np.diff(y) / h
    n = len(x)
    d = np.zeros(n)
    for i in range(1, n - 1):
        if delta[i - 1] * delta[i] <= 0:
            d[i] = 0.0
        else:
            w1 = 2 * h[i] + h[i - 1]
            w2 = h[i] + 2 * h[i - 1]
            d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])
    # one-sided ends with the standard monotonicity clamps
    def end(hh0, hh1, d0, d1):
        dd = ((2 * hh0 + hh1) * d0 - hh0 * d1) / (hh0 + hh1)
        if np.sign(dd) != np.sign(d0):
            return 0.0
        if np.sign(d0) != np.sign(d1) and abs(dd) > 3 * abs(d0):
            return 3 * d0
        return dd
    d[0] = end(h[0], h[1], delta[0], delta[1]) if n > 2 else delta[0]
    d[-1] = end(h[-1], h[-2], delta[-1], delta[-2]) if n > 2 else delta[-1]
    return d


def pchip_eval(x, y, d, xq):
    xq = np.atleast_1d(xq)
    idx = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    h = x[idx + 1] - x[idx]
    t = (xq - x[idx]) / h
    h00 = 2 * t ** 3 - 3 * t ** 2 + 1
    h10 = t ** 3 - 2 * t ** 2 + t
    h01 = -2 * t ** 3 + 3 * t ** 2
    h11 = t ** 3 - t ** 2
    return h00 * y[idx] + h10 * h * d[idx] + h01 * y[idx + 1] + h11 * h * d[idx + 1]


class Branch:
    """N(p) [molecules per formula cell] from tabulated points plus a Henry tail."""

    def __init__(self, p_bar, n_percell, k_henry_percell_per_pa=None):
        o = np.argsort(p_bar)
        self.p = np.asarray(p_bar, float)[o]
        self.n = np.asarray(n_percell, float)[o]
        self.lp = np.log(self.p)
        self.d = pchip_slopes(self.lp, self.n)
        self.p0_pa = self.p[0] * 1e5
        self.n0 = self.n[0]
        # Henry slope: the Widom value if given, else the lowest point through the origin
        self.k = k_henry_percell_per_pa if k_henry_percell_per_pa else self.n0 / self.p0_pa

    def __call__(self, p_pa):
        p_pa = np.atleast_1d(np.asarray(p_pa, float))
        out = np.empty_like(p_pa)
        lo = p_pa <= self.p0_pa
        out[lo] = self.k * p_pa[lo]
        if (~lo).any():
            q = np.clip(p_pa[~lo] / 1e5, self.p[0], self.p[-1])
            out[~lo] = pchip_eval(self.lp, self.n, self.d, np.log(q))
        return out


# ----------------------------------------------------------------- Omega
def omega_numeric(p_pa, branch, f_host, v_cell_A3, vm):
    """Omega_os [kJ/mol per cell] on p_pa (ascending, p_pa[0] = the Henry crossover)."""
    n = branch(p_pa)
    integrand = n * vm                      # m^3/mol per cell
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(p_pa))))
    # analytic Henry piece from 0 to p_pa[0]: int_0^p0 K p (RT/p) dp = N(p0) RT
    head = branch(p_pa[0])[0] * RT          # kJ/mol per cell
    ads = head + cum / 1000.0
    pv = p_pa * v_cell_A3 * 1e-30 * NA / 1000.0
    return f_host + pv - ads


def crossings(p_pa, d):
    out = []
    s = np.sign(d)
    for i in np.where(np.diff(s) != 0)[0]:
        # linear interpolation on the sign change gives a much better pressure than
        # the midpoint of a log grid
        f = d[i] / (d[i] - d[i + 1])
        pc = p_pa[i] * (p_pa[i + 1] / p_pa[i]) ** f
        out.append((pc / 1e5, "lp -> np" if d[i] > 0 else "np -> lp"))
    return out


def transitions_numeric(br_lp, br_np, dF, p_range=(1e2, 1e6), n=6000):
    p = np.logspace(np.log10(p_range[0]), np.log10(p_range[1]), n)
    vm = peng_robinson_vm(np.logspace(np.log10(p_range[0]), np.log10(p_range[1]), 400))
    vm = np.interp(p, np.logspace(np.log10(p_range[0]), np.log10(p_range[1]), 400), vm)
    o_lp = omega_numeric(p, br_lp, 0.0, V_LP, vm)
    o_np = omega_numeric(p, br_np, dF, V_NP, vm)
    return p, o_lp, o_np, crossings(p, o_np - o_lp)


def fmt(cr):
    return ", ".join(f"{p:.3f} bar ({w})" for p, w in cr) if cr else "no crossing"


# ----------------------------------------------------------------- validation
def validate(grid_bar, dF, verbose=True):
    """Synthetic Coudert-Langmuir points -> numerical route -> must give 0.273 / 3.97."""
    nmax_lp = 12.183 * MOLKG_TO_PERCELL
    b_lp = (K_LP_EXP / 12.183) * 1e5                      # per bar
    nmax_np = NMAX_NP_EXP
    b_np = (K_NP_EXP / (NMAX_NP_EXP / MOLKG_TO_PERCELL)) * 1e5
    n_lp = nmax_lp * b_lp * grid_bar / (1 + b_lp * grid_bar)
    n_np = nmax_np * b_np * grid_bar / (1 + b_np * grid_bar)
    br_lp = Branch(grid_bar, n_lp, k_henry_percell_per_pa=nmax_lp * b_lp / 1e5)
    br_np = Branch(grid_bar, n_np, k_henry_percell_per_pa=nmax_np * b_np / 1e5)
    _, _, _, cr = transitions_numeric(br_lp, br_np, dF)
    got = {w: p for p, w in cr}
    if verbose:
        print("=" * 84)
        print("STEP 2 -- VALIDATION OF THE NUMERICAL ROUTE")
        print("=" * 84)
        print(f"  synthetic Coudert-Langmuir points on our own pressure grid "
              f"({len(grid_bar)} points, {grid_bar[0]:g}-{grid_bar[-1]:g} bar)")
        print(f"  numerical  -> {fmt(cr)}")
        print(f"  analytical -> 0.273 bar (lp -> np), 3.972 bar (np -> lp)")
        for w, target in (("lp -> np", 0.273), ("np -> lp", 3.972)):
            if w in got:
                err = 100 * (got[w] - target) / target
                print(f"    {w}: {got[w]:.3f} vs {target:.3f} bar  -> {err:+.2f} %")
        ok = all(w in got and abs(got[w] - t) / t < 0.05
                 for w, t in (("lp -> np", 0.273), ("np -> lp", 3.972)))
        print(f"  -> {'VALIDATED' if ok else 'NOT VALIDATED'} (both within 5 %)")
        return ok, got
    return None, got


# ----------------------------------------------------------------- data
def load_branch(csv, k_henry_mol_kg_pa=None):
    d = pd.read_csv(csv, comment="#")
    k = k_henry_mol_kg_pa * MOLKG_TO_PERCELL if k_henry_mol_kg_pa else None
    return (Branch(d.p_bar.to_numpy(), d.absolute_molec_uc.to_numpy(), k),
            d.p_bar.to_numpy(), d.absolute_molec_uc.to_numpy(),
            (d.absolute_molec_uc_err95 / 2.776).to_numpy())


def bootstrap(p_lp, n_lp, s_lp, k_lp, p_np, n_np, s_np, k_np, dF, nboot, rng):
    lo, hi = [], []
    for _ in range(nboot):
        bl = Branch(p_lp, np.maximum(n_lp + rng.normal(0, s_lp), 1e-9), k_lp)
        bn = Branch(p_np, np.maximum(n_np + rng.normal(0, s_np), 1e-9), k_np)
        cr = transitions_numeric(bl, bn, dF)[3]
        g = {w: p for p, w in cr}
        lo.append(g.get("lp -> np", np.nan))
        hi.append(g.get("np -> lp", np.nan))
    return np.array(lo), np.array(hi)


def ci(a):
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return "none"
    return f"{np.median(a):.3f} [{np.percentile(a, 2.5):.3f}, {np.percentile(a, 97.5):.3f}]"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dF", type=float, default=2.5)
    ap.add_argument("--scan", type=float, nargs=2, default=[1.0, 4.0])
    ap.add_argument("--boot", type=int, default=400)
    ap.add_argument("--np-csv", default="results/isotherm_MIL53_npCO2cell_CO2_304K.csv")
    ap.add_argument("--kh-np", type=float, default=None,
                    help="Widom K_H of the np branch [mol/kg/Pa]; default anchors on the lowest point")
    ap.add_argument("--kh-lp", type=float, default=1.86944e-4)
    ap.add_argument("--seed", type=int, default=20260928)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)

    np_csv = ROOT / a.np_csv
    d_np_raw = pd.read_csv(np_csv, comment="#")
    grid = d_np_raw.p_bar.to_numpy()
    ok, _ = validate(grid, a.dF)
    if not ok:
        print("\n  numerical route did not reproduce the analytical result; stopping.")
        return 1

    print("\n" + "=" * 84)
    print("STEP 1 -- PRIMARY: numerical integration of the SIMULATED isotherms")
    print("=" * 84)
    lp_sets = {
        "lp PACMAN (consistent pair, 4 pts)": ROOT / "results/isotherm_MIL53_lp_pacman_CO2_304K.csv",
        "lp SABVUN-DDEC (12 pts)": ROOT / "results/isotherm_MIL53_lp_CO2_304K.csv",
    }
    br_np, p_np, n_np, s_np = load_branch(np_csv, a.kh_np)
    print(f"  np branch: {np_csv.name}, {len(p_np)} points, "
          f"{p_np[0]:g}-{p_np[-1]:g} bar, Henry anchor "
          + (f"Widom K_H = {a.kh_np:.3e} mol/kg/Pa" if a.kh_np else
             f"lowest point ({br_np.k / MOLKG_TO_PERCELL:.3e} mol/kg/Pa equivalent)"))

    rows = []
    for name, csv in lp_sets.items():
        br_lp, p_lp, n_lp, s_lp = load_branch(csv, a.kh_lp)
        p, o_lp, o_np, cr = transitions_numeric(br_lp, br_np, a.dF)
        print(f"\n  {name}  ({len(p_lp)} points)")
        print(f"    dF = {a.dF} kJ/mol -> {fmt(cr)}")
        blo, bhi = bootstrap(p_lp, n_lp, s_lp, a.kh_lp * MOLKG_TO_PERCELL,
                             p_np, n_np, s_np, br_np.k, a.dF, a.boot, rng)
        print(f"    bootstrap ({a.boot} replicates, 95 % CI):")
        print(f"      lp -> np : {ci(blo)} bar   (experiment 0.25-0.3)")
        print(f"      np -> lp : {ci(bhi)} bar   (experiment 5-6)")
        for dF in np.arange(a.scan[0], a.scan[1] + 1e-9, 0.25):
            cr2 = transitions_numeric(br_lp, br_np, float(dF))[3]
            g = {w: q for q, w in cr2}
            rows.append({"lp_branch": name, "dF_kJ_mol": float(dF),
                         "lp->np_bar": g.get("lp -> np", np.nan),
                         "np->lp_bar": g.get("np -> lp", np.nan)})
    scan = pd.DataFrame(rows)
    print("\n  Delta F_host scan, numerical route:")
    with pd.option_context("display.float_format", "{:.3f}".format, "display.width", 200):
        print(scan.pivot(index="dF_kJ_mol", columns="lp_branch",
                         values=["lp->np_bar", "np->lp_bar"]).to_string())
    out = ROOT / "results" / "osmotic_numeric_scan.csv"
    scan.to_csv(out, index=False)
    print(f"\n  wrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
