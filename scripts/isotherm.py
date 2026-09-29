#!/usr/bin/env python3
"""Assemble the CO2/MIL-53(Al) lp isotherm from runs/<tag>/p_*/ and plot it.

Every quantity is carried in both conventions (absolute and excess) and both
units (mol/kg and molecules per conventional unit cell, 832.4 g/mol), with
RASPA's 95 % error bars, the accepted insertion/deletion counts and the
relative error per point (see NOTES.md, "What every result carries").

Outputs:
  results/isotherm_MIL53_lp_CO2_304K.csv
  results/isotherm_MIL53_lp_CO2_304K.png

Usage: python scripts/isotherm.py [--tag production] [--no-reference]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import parse_raspa_output as pro  # noqa: E402
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
M_UC = 832.415          # g/mol per conventional unit cell (Al4C32H20O20)
MOLKG_PER_UC = 1000.0 / M_UC   # 1 molecule/uc = 1.2013 mol/kg
THETA_HE = 0.7115       # Widom He, 298 K, eps/k 10.9 K, sigma 2.64 A
REF_CSV = ROOT / "reference" / "bourrelly2005_MIL53Al_CO2_304K.csv"
LP_WINDOW_BAR = 9.0     # only above this is the real material fully lp
STEP_LO_BAR = 5.0       # 5-9 bar: np-lp breathing step


def collect(tag):
    rows = []
    for d in sorted((ROOT / "runs" / tag).glob("p_*"), key=lambda p: float(p.name[2:])):
        files = list((d / "Output" / "System_0").glob("*.data")) if (d / "Output" / "System_0").is_dir() else []
        if not files:
            continue
        r = pro.parse_file(files[0]).iloc[0].to_dict()
        r.update(convergence_guards(files[0]))
        # store the path relative to the repository, not the machine it ran on
        r["file"] = str(Path(r["file"]).resolve().relative_to(ROOT))
        if not r.get("finished"):
            print(f"  (skipping {d.name}: not finished)")
            continue
        rows.append(r)
    if not rows:
        sys.exit("no finished points found")
    df = pd.DataFrame(rows)
    df["p_bar"] = df.p_Pa / 1e5
    df["fugacity_bar"] = df.fugacity_Pa / 1e5
    for kind in ("absolute", "excess"):
        df[f"{kind}_molec_uc_rel_err"] = df[f"{kind}_molec_uc_err95"] / df[f"{kind}_molec_uc"]
    df["theta_He"] = THETA_HE
    df["M_uc_g_per_mol"] = M_UC
    return df.sort_values("p_bar").reset_index(drop=True)


def convergence_guards(data_file, nblocks=5):
    """Drift and fluctuation guards on one RASPA production run.

    Both were made permanent after Phase 4, where a still-filling pore produced
    smooth-looking block error bars that described a TREND, not noise (NOTES.md,
    "Step 2, second attempt"). They are applied here too, because the CO2-cell np
    points run at high loading where insertion acceptance is lowest and sampling is
    hardest.

    drift       last block mean minus first block mean, against the 95 % error bar.
                A point whose loading is still climbing has a biased mean.
    sd(N)       standard deviation of the instantaneous molecule count. <dN^2> is a
                physical property of the grand-canonical ensemble; an sd near zero
                means N is effectively frozen and the samples are correlated, so the
                block error bar understates the true uncertainty.
    """
    tr = pro.loading_trace(data_file)
    # an in-flight or trace-less run gives an empty frame with no columns at all
    prod = tr[tr.stage == "prod"] if ("stage" in tr.columns and not tr.empty) else tr
    if prod.empty or len(prod) < nblocks:
        return {"drift_molec_uc": np.nan, "drifting": True, "sd_N_box": np.nan,
                "sd_N_per_uc": np.nan, "n_trace": len(prod), "frozen_N": True}
    N = prod.N_box.to_numpy()
    n = len(prod) // nblocks * nblocks
    blocks = N[len(N) - n:].reshape(nblocks, -1).mean(axis=1)
    drift = float(blocks[-1] - blocks[0])
    err95 = float(2.776 * blocks.std(ddof=1) / np.sqrt(nblocks))
    # The endpoint test alone is MIS-CALIBRATED: |last - first| has a 1-sigma scale of
    # sqrt(2)*sd(blocks) = 1.41 sd, while err95 = 2.776/sqrt(5) sd = 1.24 sd, so the
    # threshold sits below the natural scale of the statistic and a perfectly converged
    # run trips it about 46 % of the time (measured, 200,000 synthetic runs). Phase 4's
    # real drift was monotonic filling, so a trend test is what distinguishes the two.
    # A point is called drifting only if the endpoint test fires AND the loading trend
    # across 10 blocks is genuinely directional; joint false-positive rate ~6 %.
    n10 = len(N) // 10 * 10
    first_gap = np.nan
    under_equilibrated = False
    if n10 >= 10:
        b10 = N[len(N) - n10:].reshape(10, -1).mean(axis=1)
        trend = float(np.corrcoef(np.arange(10), b10)[0, 1]) if b10.std() > 0 else 0.0
        # Equilibration check, independent of the trend test. A run started from an
        # empty box at high loading may still be filling during the first production
        # block; the trend statistic can miss that, because one low block among ten
        # barely moves a correlation. So the first block is compared directly with
        # the other nine: if it sits more than 2 sd below them, initialisation was too
        # short for THAT point and it needs more init cycles, not more production.
        rest = b10[1:]
        sd_rest = float(rest.std(ddof=1))
        first_gap = float((b10[0] - rest.mean()) / sd_rest) if sd_rest > 0 else 0.0
        under_equilibrated = bool(first_gap < -2.0)
    else:
        trend = 0.0
    sd = float(prod.N_box.std())
    mean_n = float(prod.N_box.mean())
    # a Poisson-like floor: sd(N) should be of order sqrt(N) in the grand canonical
    # ensemble; far below that means the particle number is not really moving
    return {"drift_N_box": drift, "drift_err95_N_box": err95, "trend_r_10blocks": trend,
            "first_block_sd_below_rest": first_gap, "under_equilibrated": under_equilibrated,
            "drifting": bool(abs(drift) > max(err95, 1e-12) and abs(trend) >= 0.5),
            "drift_endpoint_only": bool(abs(drift) > max(err95, 1e-12)),
            "sd_N_box": sd, "mean_N_box": mean_n,
            "sd_over_sqrtN": sd / np.sqrt(mean_n) if mean_n > 0 else np.nan,
            "n_trace": int(len(prod)),
            "frozen_N": bool(mean_n > 0 and sd < 0.2 * np.sqrt(mean_n))}


def load_reference():
    if not REF_CSV.exists():
        return None
    r = pd.read_csv(REF_CSV)
    r["in_lp_window"] = r.pressure_bar >= LP_WINDOW_BAR
    return r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="production")
    ap.add_argument("--no-reference", action="store_true")
    ap.add_argument("--theta-he", type=float, default=THETA_HE,
                    help="helium void fraction used for the excess column (phase-specific)")
    ap.add_argument("--z", type=int, default=1,
                    help="formula cells per RASPA unit cell (np file cell holds 2)")
    ap.add_argument("--out-prefix", default="isotherm_MIL53_lp_CO2_304K")
    ap.add_argument("--phase", choices=("lp", "np"), default="lp",
                    help="which MIL-53 phase this run is, for the title and the caption")
    ap.add_argument("--fit-window", type=float, default=LP_WINDOW_BAR,
                    help="lowest pressure [bar] included in the Langmuir overlay. The default is "
                         "the lp window (P >= 9 bar, where the real material is fully open); the np "
                         "branch has no such window, so pass 0 to fit all of it.")
    a = ap.parse_args(argv)

    df = collect(a.tag)
    if a.z != 1:   # RASPA reports per FILE cell; convert to per formula cell
        for c in [c for c in df.columns if "molec_uc" in c and not c.endswith("rel_err")]:
            df[c] = df[c] / a.z
    # theta_He: RASPA computed the excess column with the value in simulation.input.
    # --theta-he is the Widom value we want on record. They can differ when the helium
    # run finished after the GCMC started, so carry BOTH and say which did what.
    df["theta_He"] = df["theta_He_used"]          # what the excess column is built on
    df["theta_He_widom"] = a.theta_he             # what we measured
    mism = df[(df.theta_He_used - a.theta_he).abs() > 1e-6]
    if len(mism):
        used = sorted(set(mism.theta_He_used.round(8)))
        print(f"  NOTE: RASPA built the excess column with theta_He = {used} (from simulation.input), "
              f"not the Widom value {a.theta_he} passed here.")
        print(f"        Absolute loading is unaffected. The excess correction shifts by "
              f"{100 * abs(used[0] - a.theta_he) / a.theta_he:.2f} % of itself; both values are in the CSV.")
    df["formula_cells_per_raspa_cell"] = a.z
    ref = None if a.no_reference else load_reference()
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    cols = ["p_bar", "p_Pa", "fugacity_bar", "fugacity_coeff", "T_K", "unit_cells", "n_unit_cells",
            "absolute_mol_kg", "absolute_mol_kg_err95", "absolute_molec_uc", "absolute_molec_uc_err95",
            "absolute_molec_uc_rel_err",
            "excess_mol_kg", "excess_mol_kg_err95", "excess_molec_uc", "excess_molec_uc_err95",
            "excess_molec_uc_rel_err", "bulk_density_kg_m3", "theta_He", "theta_He_widom", "M_uc_g_per_mol",
            "acc_insertion", "accepted_insertion", "acc_deletion", "accepted_deletion",
            "acc_reinsertion", "acc_translation", "acc_rotation", "cycles_init", "cycles_prod",
            "drift_N_box", "drift_err95_N_box", "trend_r_10blocks", "drifting",
            "drift_endpoint_only", "first_block_sd_below_rest", "under_equilibrated",
            "sd_N_box", "mean_N_box",
            "sd_over_sqrtN", "frozen_N", "n_trace",
            "n_warnings", "missing_vdw_pairs", "file"]
    csv = out / f"{a.out_prefix}.csv"
    header = (f"# CO2 in MIL-53(Al) {a.phase}, rigid framework, GCMC (RASPA2 2.0.50), T = {df.T_K.iloc[0]:g} K\n"
              f"# absolute and excess loading; 1 molecule/uc = {MOLKG_PER_UC:.4f} mol/kg (M_uc = {M_UC} g/mol)\n"
              f"# excess uses theta_He = {df.theta_He.iloc[0]:g}, the value RASPA was given in simulation.input\n"
              f"# theta_He_widom = {a.theta_he:g} is the measured Widom value (He, 298 K, eps/k 10.9 K, sigma 2.64 A)\n"
              f"# loadings per FORMULA cell (M = {M_UC} g/mol, Z = 4); RASPA unit cell holds {a.z} formula cell(s)\n"
              f"# errors are RASPA's 95 % confidence half-widths (5 blocks, t = 2.776)\n")
    with open(csv, "w") as fh:
        fh.write(header)
        df[[c for c in cols if c in df.columns]].to_csv(fh, index=False)
    print(f"wrote {csv.relative_to(ROOT)}")

    show = ["p_bar", "absolute_mol_kg", "absolute_mol_kg_err95", "absolute_molec_uc",
            "excess_mol_kg", "excess_molec_uc", "absolute_molec_uc_rel_err",
            "acc_insertion", "accepted_insertion", "accepted_deletion"]
    with pd.option_context("display.width", 220, "display.float_format", "{:.4g}".format):
        print(df[[c for c in show if c in df.columns]].to_string(index=False))
    bad = df[df.absolute_molec_uc_rel_err > 0.10]
    for _, r in bad.iterrows():
        print(f"  FLAG: {r.p_bar:g} bar has {100 * r.absolute_molec_uc_rel_err:.1f} % relative error on loading (> 10 %)")
    for _, r in df[df.get("drifting", False) == True].iterrows():
        print(f"  GUARD FAIL: {r.p_bar:g} bar is still DRIFTING -- last block minus first = "
              f"{r.drift_N_box:+.2f} molecules against a {r.drift_err95_N_box:.2f} error bar, "
              f"with a directional 10-block trend r = {r.trend_r_10blocks:+.2f}; "
              f"the mean is biased and the error bar describes a trend, not noise")
    for _, r in df[(df.get("drift_endpoint_only", False) == True) & (df.get("drifting", False) == False)].iterrows():
        print(f"  (info: {r.p_bar:g} bar trips the endpoint drift test "
              f"[{r.drift_N_box:+.2f} vs {r.drift_err95_N_box:.2f}] but its 10-block trend is "
              f"r = {r.trend_r_10blocks:+.2f}, not directional -- endpoint noise, not filling)")
    for _, r in df[df.get("under_equilibrated", False) == True].iterrows():
        print(f"  GUARD FAIL: {r.p_bar:g} bar is UNDER-EQUILIBRATED -- the first production "
              f"block sits {abs(r.first_block_sd_below_rest):.1f} sd below the mean of the other "
              f"nine; extend INITIALISATION for this point, not production")
    for _, r in df[df.get("frozen_N", False) == True].iterrows():
        print(f"  GUARD FAIL: {r.p_bar:g} bar has sd(N) = {r.sd_N_box:.2f} against sqrt(<N>) = "
              f"{np.sqrt(r.mean_N_box):.2f}; the particle number is effectively frozen, so the "
              f"block error bar understates the true uncertainty")
    for _, r in df[df.absolute_molec_uc_rel_err > 0.05].iterrows():
        print(f"  NOTE: {r.p_bar:g} bar has {100 * r.absolute_molec_uc_rel_err:.1f} % relative "
              f"error (> 5 %, the threshold set when the cycle count was cut to 50,000)")
    low = df[df.acc_insertion < 0.01]
    for _, r in low.iterrows():
        print(f"  FLAG: {r.p_bar:g} bar insertion acceptance {100 * r.acc_insertion:.2f} % (< 1 %), "
              f"{r.accepted_insertion:.0f} accepted insertions")

    # ---------------- figure ----------------
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    solid = df.absolute_molec_uc_rel_err <= 0.10
    ax.fill_between(df.p_bar, df.excess_mol_kg, df.absolute_mol_kg, color=SIM["color"], alpha=0.16, lw=0,
                    label=f"absolute - excess (theta$_{{He}}$ = {a.theta_he})")
    ax.errorbar(df.p_bar, df.absolute_mol_kg, yerr=df.absolute_mol_kg_err95,
                label="this work, absolute", **SIM)
    ax.errorbar(df.p_bar, df.excess_mol_kg, yerr=df.excess_mol_kg_err95,
                label="this work, excess", **{**SIM, "marker": "v", "ls": "--", "mfc": "white"})
    if (~solid).any():
        ax.plot(df.p_bar[~solid], df.absolute_mol_kg[~solid], "o", mfc="white", mec=SIM["color"], ms=9, lw=0,
                label="rel. error > 10 % (not solid)")
    if ref is not None:
        for flag, style, lab in ((True, REF, f"Bourrelly 2005, excess (assumed), P $\\geq$ {LP_WINDOW_BAR:g} bar"),
                                 (False, {**REF, "color": "#8c8c88", "mfc": "white", "ms": 7},
                                  f"same, {STEP_LO_BAR:g}-{LP_WINDOW_BAR:g} bar: inside the np-lp step, NOT compared")):
            g = ref[ref.in_lp_window == flag]
            if len(g):
                ax.errorbar(g.pressure_bar, g.loading_mmol_per_g, label=lab, **style)
    # Langmuir curves (fitted over P >= 9 bar; see scripts/langmuir.py)
    try:
        from langmuir import fit_langmuir
        w = df[df.p_bar >= a.fit_window]
        if len(w) >= 3:
            fs = fit_langmuir(w.p_bar, w.absolute_mol_kg, w.absolute_mol_kg_err95 / 2.776)
            pp = np.logspace(np.log10(df.p_bar.min()), np.log10(df.p_bar.max()), 300)
            if fs["ok"]:
                ax.plot(pp, fs["N_max"] * fs["b_per_bar"] * pp / (1 + fs["b_per_bar"] * pp),
                        color=SIM["color"], lw=1.1, ls=":", zorder=1,
                        label=(f"Langmuir fit, this work (absolute, P $\\geq$ {a.fit_window:g} bar):\n"
                               f"N$_{{max}}$ = {fs['N_max']:.3g} mol/kg, K = {fs['K_mol_kg_Pa']:.3g} mol/kg/Pa"))
            if ref is not None:
                er = ref[ref.in_lp_window]
                fe = fit_langmuir(er.pressure_bar, er.loading_mmol_per_g, np.full(len(er), 0.05))
                ax.plot(pp, fe["N_max"] * fe["b_per_bar"] * pp / (1 + fe["b_per_bar"] * pp),
                        color=REF["color"], lw=1.1, ls="--", zorder=1,
                        label="Langmuir fit, experiment = Coudert virtual lp curve")
    except Exception as exc:  # never let the figure fail on the fit
        print(f"  (Langmuir overlay skipped: {exc})")
    ax.set_xscale("log")
    ax.set_xlabel("Pressure [bar]")
    ax.set_ylabel("Loading [mol/kg  =  mmol/g]")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x / MOLKG_PER_UC, lambda x: x * MOLKG_PER_UC))
    sec.set_ylabel("Loading [molecules / unit cell]")
    ax.set_title(f"CO$_2$ in MIL-53(Al) {a.phase}, rigid framework, {df.T_K.iloc[0]:g} K")
    ax.legend(fontsize=8.5, loc="upper left")
    captions = {
        "lp": ("Rigid lp framework: this is the virtual lp branch, so it is comparable with experiment only "
               f"for P $\\geq$ {LP_WINDOW_BAR:g} bar,\nwhere the real material is fully open. Below ~5 bar the real solid is np "
               "(Coudert 2008). Reference points are excess\n(assumed: manometry, convention not stated by the source) and are "
               "never converted; the band shows our absolute-excess gap.\nError bars: 95 % CI."),
        "np": ("Rigid np framework, the EMPTY dehydrated np cell. This is NOT the physical np branch: the real np phase under "
               "CO$_2$ is expanded,\nwhich is what breathing means. The cell is too tight for CO$_2$ (LCD 2.83 A vs 3.3 A kinetic "
               "diameter, theta$_{He}$ = 0.121),\nso the isotherm stays linear to 10 bar and reaches 0.21 molecules/cell against a "
               "literature np capacity of ~3.0.\nLoadings are per FORMULA cell: the file cell is a genuine superstructure holding 2 "
               "(NOTES.md 5.1). Error bars: 95 % CI."),
    }
    fig.text(0.01, -0.06, captions[a.phase], fontsize=7.5, va="top")
    fig.savefig(out / f"{a.out_prefix}.png")
    print(f"wrote results/{a.out_prefix}.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
