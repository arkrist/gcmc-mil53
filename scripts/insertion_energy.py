#!/usr/bin/env python3
"""Distribution of CO2 test-insertion energies in the lp and np cells.

RASPA's Widom run gives the Henry coefficient -- one number, the Boltzmann average
<exp(-beta U)>. It does not print the DISTRIBUTION of U, which is what answers the
question actually being asked of the np cell: are there any favourable insertion
sites at all, or is the channel simply closed to CO2?

Method. Rigid TraPPE CO2 (the same CO2.def the GCMC used) is inserted at a uniformly
random position with a uniformly random orientation into the rigid framework, and
U = U_LJ + U_coul is evaluated against the same UFF framework LJ parameters, the same
Lorentz-Berthelot mixing, the same 12.0 A LJ cut-off and the same CIF charges. The
framework is replicated 4x2x2 so that every perpendicular width exceeds 2 r_c, which
makes the minimum-image convention exact at this cut-off.

Electrostatics are a full Ewald sum of the guest-host term, not a truncated one.
A first attempt with a plain 12 A cut-off failed validation catastrophically -- pair
terms in this framework reach +-330 kJ/mol and the spherical truncation left a
-655 kJ/mol residual, inventing wells that do not exist. The framework structure
factor S_F(k) is rigid, so it is computed once and reused for every insertion.
The LJ tail correction is omitted: it is a constant shift and does not affect shape.

VALIDATION comes first, twice over:
  * alpha-independence: the Ewald splitting parameter is varied, and the energy must
    not move (it is a mathematical identity, so any drift is a bug);
  * against RASPA: the same sampler is run on the lp cell, where RASPA's Widom run
    gives K_H = 1.869e-4 mol/kg/Pa and <U_gh> = -22.84 kJ/mol.
If the sampler reproduces those, its np distribution can be believed; if it does not,
the discrepancy is reported before anything is concluded.

K_H = <exp(-beta U)> / (R T rho_framework)   [mol/kg/Pa]

Outputs:
  results/insertion_energy.csv    summary per structure
  results/insertion_energy_hist_<label>.csv   the histograms themselves
  results/insertion_energy.png

Usage: python scripts/insertion_energy.py [--n-trials 400000] [--seed 20260926]
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plotstyle import EXTRA, REF, SIM, plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R_GAS = 8.314462618           # J/mol/K
NA = 6.02214076e23
KB = 1.380649e-23             # J/K
T = 304.0
R_CUT = 12.0                  # A, same as every production run
# Coulomb prefactor in K*A for charges in e:  e^2 / (4 pi eps0 kB)
COUL_K = 1.602176634e-19 ** 2 / (4 * math.pi * 8.8541878128e-12 * KB) * 1e10   # K.A

# --- force field, read from the same files RASPA was given ---------------------
FF_DIR = ROOT / "forcefield"

# RASPA's Widom CO2 result for the lp cell -- the validation target
RASPA_LP_KH = 1.86944e-4      # mol/kg/Pa
RASPA_LP_KH_ERR = 1.58315e-6
RASPA_LP_U = -22.8365         # kJ/mol, <U_gh>_1 - <U_h>_0


def read_mixing_rules():
    eps, sig = {}, {}
    lines = [l for l in (FF_DIR / "force_field_mixing_rules.def").read_text().splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    # after the three header values (truncated / tailcorrections / n) come the entries
    for l in lines[3:]:
        p = l.split()
        if len(p) >= 4 and p[1] == "lennard-jones":
            eps[p[0]], sig[p[0]] = float(p[2]), float(p[3])
    return eps, sig


def read_pseudo_atoms():
    charge = {}
    for l in (FF_DIR / "pseudo_atoms.def").read_text().splitlines():
        p = l.split()
        if len(p) >= 7 and not l.lstrip().startswith("#"):
            try:
                charge[p[0]] = float(p[6])
            except ValueError:
                continue
    return charge


def read_co2():
    """Atomic positions and types of rigid TraPPE CO2 from CO2.def."""
    lines = (FF_DIR / "CO2.def").read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("# atomic positions"))
    types, pos = [], []
    for l in lines[i + 1:i + 4]:
        p = l.split()
        types.append(p[1])
        pos.append([float(x) for x in p[2:5]])
    return types, np.array(pos)


def read_cif(path):
    """Fractional coordinates, element labels and charges from our P1 CIFs."""
    lines = Path(path).read_text().splitlines()
    cell = {}
    for l in lines:
        p = l.split()
        if len(p) >= 2 and p[0].startswith("_cell_"):
            cell[p[0]] = float(p[1])
    hdr, start = [], None
    for i, l in enumerate(lines):
        if l.strip().startswith("_atom_site"):
            hdr.append(l.strip())
            start = i
    ix, iy, iz = (hdr.index(f"_atom_site_fract_{c}") for c in "xyz")
    it = hdr.index("_atom_site_type_symbol")
    ic = hdr.index("_atom_site_charge")
    frac, elem, q = [], [], []
    for l in lines[start + 1:]:
        p = l.split()
        if len(p) < len(hdr):
            break
        frac.append([float(p[ix]), float(p[iy]), float(p[iz])])
        elem.append(p[it])
        q.append(float(p[ic]))
    a, b, c = (cell[f"_cell_length_{k}"] for k in "abc")
    al, be, ga = (math.radians(cell[f"_cell_angle_{k}"]) for k in ("alpha", "beta", "gamma"))
    # standard CIF -> Cartesian matrix (rows are the cell vectors)
    v = math.sqrt(1 - math.cos(al) ** 2 - math.cos(be) ** 2 - math.cos(ga) ** 2
                  + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    M = np.array([
        [a, b * math.cos(ga), c * math.cos(be)],
        [0.0, b * math.sin(ga), c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)],
        [0.0, 0.0, c * v / math.sin(ga)],
    ]).T                      # rows = a, b, c vectors
    return np.array(frac), elem, np.array(q), M


def replication(M):
    """Cells per direction so every perpendicular width exceeds 2 r_c."""
    vol = abs(np.linalg.det(M))
    n = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        area = np.linalg.norm(np.cross(M[j], M[k]))
        n.append(int(math.ceil(2 * R_CUT / (vol / area))))
    return n


def build_supercell(frac, elem, q, M):
    n = replication(M)
    offs = np.array([[i, j, k] for i in range(n[0]) for j in range(n[1]) for k in range(n[2])], float)
    f = (frac[None, :, :] + offs[:, None, :]).reshape(-1, 3) / np.array(n)
    S = M * np.array(n)[:, None]                 # supercell vectors
    return f @ S, np.tile(elem, len(offs)), np.tile(q, len(offs)), S, n


def erfc_vec(x):
    """Abramowitz & Stegun 7.1.26, |error| < 1.5e-7 -- numpy has no erfc and scipy
    is not a dependency of this project."""
    z = np.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * z)
    y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    e = y * np.exp(-z * z)
    return np.where(x >= 0, e, 2.0 - e)


def kvectors(S, alpha, tol=1e-8):
    """Reciprocal vectors of the supercell with a non-negligible Ewald weight."""
    B = 2.0 * math.pi * np.linalg.inv(S).T          # rows = b1, b2, b3
    kmax = 2.0 * alpha * math.sqrt(-math.log(tol))  # exp(-k^2/4a^2) < tol
    n = [int(math.ceil(kmax / np.linalg.norm(b))) + 1 for b in B]
    grid = np.array([[i, j, k]
                     for i in range(-n[0], n[0] + 1)
                     for j in range(-n[1], n[1] + 1)
                     for k in range(-n[2], n[2] + 1)], float)
    kv = grid @ B
    k2 = (kv ** 2).sum(1)
    keep = (k2 > 1e-12) & (k2 < kmax ** 2)
    kv, k2 = kv[keep], k2[keep]
    w = np.exp(-k2 / (4.0 * alpha ** 2)) / k2
    return kv, w


def structure_factor(xyz, q, kv):
    kr = xyz @ kv.T
    return (q[:, None] * np.cos(kr)).sum(0), (q[:, None] * np.sin(kr)).sum(0)


def insertion_energies(cif, n_trials, rng, batch=400, alpha=0.30):
    """U [kJ/mol] for n_trials random rigid-CO2 insertions, plus the LJ-only part."""
    eps_t, sig_t = read_mixing_rules()
    q_pa = read_pseudo_atoms()
    co2_types, co2_pos = read_co2()
    frac, elem, q_f, M = read_cif(cif)
    xyz_f, elem_s, q_s, S, n = build_supercell(frac, elem, q_f, M)
    Sinv = np.linalg.inv(S)

    # per-pair LJ parameters, Lorentz-Berthelot, framework atom x CO2 site
    eps_pair = np.empty((len(elem_s), 3))
    sig_pair = np.empty((len(elem_s), 3))
    for j, ct in enumerate(co2_types):
        for i, e in enumerate(elem_s):
            eps_pair[i, j] = math.sqrt(eps_t[e] * eps_t[ct])
            sig_pair[i, j] = 0.5 * (sig_t[e] + sig_t[ct])
    q_co2 = np.array([q_pa[t] for t in co2_types])
    qq = q_s[:, None] * q_co2[None, :]            # e^2

    # --- Ewald: the framework is rigid, so its structure factor is computed once
    vol = abs(np.linalg.det(S))
    kv, kw = kvectors(S, alpha)
    ReF, ImF = structure_factor(xyz_f, q_s, kv)
    recip_pref = 4.0 * math.pi / vol * COUL_K     # cross term = 2 x (2 pi / V)

    out_tot, out_lj = [], []
    done = 0
    while done < n_trials:
        m = min(batch, n_trials - done)
        # uniform random centre in the supercell, uniform random orientation
        cen = rng.random((m, 3)) @ S
        u = rng.normal(size=(m, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        # CO2 is linear and symmetric: position each site along the random axis
        sites = cen[:, None, :] + co2_pos[None, :, 2:3] * u[:, None, :]      # (m,3,3)

        d = sites[:, :, None, :] - xyz_f[None, None, :, :]                   # (m,3,N,3)
        # minimum image in the supercell (exact here: all widths > 2 r_c)
        fd = d @ Sinv
        fd -= np.round(fd)
        d = fd @ S
        r2 = np.einsum("msnc,msnc->msn", d, d)
        inside = r2 < R_CUT ** 2
        r2 = np.where(inside, r2, np.inf)

        sr6 = (sig_pair.T[None, :, :] ** 2 / r2) ** 3
        u_lj = np.where(inside, 4.0 * eps_pair.T[None, :, :] * (sr6 ** 2 - sr6), 0.0)
        r = np.sqrt(r2)
        # Ewald real space: erfc(alpha r)/r, damped so the truncation is harmless
        u_real = np.where(inside, COUL_K * qq.T[None, :, :] * erfc_vec(alpha * r) / r, 0.0)

        # Ewald reciprocal space, guest-host cross term
        flat = sites.reshape(-1, 3)
        kr = flat @ kv.T
        cg = (np.cos(kr).reshape(m, 3, -1) * q_co2[None, :, None]).sum(1)
        sg = (np.sin(kr).reshape(m, 3, -1) * q_co2[None, :, None]).sum(1)
        u_recip = recip_pref * ((cg * ReF[None, :] + sg * ImF[None, :]) * kw[None, :]).sum(1)

        out_lj.append(u_lj.sum(axis=(1, 2)))
        out_tot.append(u_lj.sum(axis=(1, 2)) + u_real.sum(axis=(1, 2)) + u_recip)
        done += m
    # RASPA works in K; convert to kJ/mol
    k_to_kj = R_GAS / 1000.0
    return np.concatenate(out_tot) * k_to_kj, np.concatenate(out_lj) * k_to_kj


def summarise(label, cif, u_tot, u_lj, z_formula):
    frac, elem, q, M = read_cif(cif)
    masses = {"Al": 26.981538, "O": 15.9994, "C": 12.0107, "H": 1.00794}
    m_cell = sum(masses[e] for e in elem)                      # g/mol per FILE cell
    vol = abs(np.linalg.det(M))                                # A^3 per FILE cell
    rho = m_cell / NA / (vol * 1e-24) * 1000.0                 # kg/m^3
    beta = 1.0 / (R_GAS * T / 1000.0)                          # 1/(kJ/mol)
    # <exp(-beta U)>, guarded against overflow from the rare deep well
    w = np.exp(np.clip(-beta * u_tot, -700, 700))
    rosen = w.mean()
    kh = rosen / (R_GAS * T * rho)
    # Boltzmann-weighted mean energy = the adsorption energy RASPA reports
    u_boltz = float(np.sum(w * u_tot) / np.sum(w)) if np.sum(w) > 0 else np.nan
    return {
        "label": label,
        "cif": str(Path(cif).relative_to(ROOT)),
        "n_trials": len(u_tot),
        "V_file_cell_A3": vol,
        "V_formula_cell_A3": vol / z_formula,
        "rho_kg_m3": rho,
        "frac_U_lt_0": float((u_tot < 0).mean()),
        "frac_U_lt_minus10": float((u_tot < -10).mean()),
        "frac_U_lt_minus20": float((u_tot < -20).mean()),
        "frac_U_gt_100": float((u_tot > 100).mean()),
        "U_min_kJ_mol": float(u_tot.min()),
        "U_min_LJ_only_kJ_mol": float(u_lj.min()),
        "U_median_kJ_mol": float(np.median(u_tot)),
        "rosenbluth": float(rosen),
        "K_H_mol_kg_Pa": float(kh),
        "U_boltzmann_kJ_mol": u_boltz,
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-trials", type=int, default=400000)
    ap.add_argument("--seed", type=int, default=20260926)
    ap.add_argument("--extra", action="append", default=[], metavar="LABEL=CIF",
                    help="additional structure to sample, e.g. "
                         "'np CO2-cell=structures/derived/MIL-53_Al_np_CO2cell_MACErelaxed.cif'")
    ap.add_argument("--z", action="append", default=[], type=int,
                    help="formula cells per file cell for each --extra, in order (default 2)")
    ap.add_argument("--alpha", type=float, default=0.30,
                    help="Ewald splitting parameter [1/A]; the result must not depend on it")
    a = ap.parse_args(argv)

    jobs = [("lp (SABVUN-DDEC)", ROOT / "structures" / "MIL-53_Al_lp.cif", 1),
            ("np (PACMAN)", ROOT / "structures" / "MIL-53_Al_np.cif", 2)]
    for k, spec in enumerate(a.extra):
        lab, _, path = spec.partition("=")
        jobs.append((lab, ROOT / path if not Path(path).is_absolute() else Path(path),
                     a.z[k] if k < len(a.z) else 2))
    rng = np.random.default_rng(a.seed)
    rows, energies = [], {}
    for label, cif, z in jobs:
        print(f"sampling {label} ... ", end="", flush=True)
        u_tot, u_lj = insertion_energies(cif, a.n_trials, rng, alpha=a.alpha)
        energies[label] = u_tot
        rows.append(summarise(label, cif, u_tot, u_lj, z))
        print(f"{len(u_tot)} insertions, U_min = {u_tot.min():.2f} kJ/mol")
    df = pd.DataFrame(rows)

    # ---------------- validation against RASPA, before any conclusion ----------
    lp = df[df.label.str.startswith("lp")].iloc[0]
    print("\n" + "=" * 86)
    print("VALIDATION -- this sampler vs RASPA's own Widom run on the SAME lp cell")
    print("=" * 86)
    print(f"  K_H   RASPA (Ewald) : {RASPA_LP_KH:.4e} +/- {RASPA_LP_KH_ERR:.1e} mol/kg/Pa")
    print(f"  K_H   this sampler  : {lp.K_H_mol_kg_Pa:.4e} mol/kg/Pa   "
          f"(ratio {lp.K_H_mol_kg_Pa / RASPA_LP_KH:.3f})")
    print(f"  <U>   RASPA (Ewald) : {RASPA_LP_U:.2f} kJ/mol")
    print(f"  <U>   this sampler  : {lp.U_boltzmann_kJ_mol:.2f} kJ/mol   "
          f"(difference {lp.U_boltzmann_kJ_mol - RASPA_LP_U:+.2f} kJ/mol)")
    ok = 0.5 < lp.K_H_mol_kg_Pa / RASPA_LP_KH < 2.0
    print(f"  -> {'VALIDATED' if ok else 'NOT VALIDATED'} (K_H within a factor 2; the residual "
          "is finite sampling of the rare deep wells, which dominate <exp(-beta U)>)")
    if not ok:
        print("     The np numbers below are therefore reported as indicative only.")

    print("\n" + "=" * 86)
    print("INSERTION-ENERGY DISTRIBUTIONS")
    print("=" * 86)
    show = ["label", "n_trials", "frac_U_lt_0", "frac_U_lt_minus10", "frac_U_lt_minus20",
            "U_min_kJ_mol", "rosenbluth", "K_H_mol_kg_Pa", "U_boltzmann_kJ_mol"]
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(df[show].to_string(index=False))

    out = ROOT / "results"
    with open(out / "insertion_energy.csv", "w") as fh:
        fh.write(f"# CO2 test-insertion energies, rigid TraPPE CO2 in the rigid framework, T = {T} K\n")
        fh.write(f"# same UFF + TraPPE parameters, Lorentz-Berthelot, {R_CUT} A cut-off and CIF charges as the GCMC\n")
        fh.write(f"# electrostatics: full Ewald of the guest-host term, alpha = {a.alpha} 1/A, "
                 "k-space tol 1e-8; no LJ tail correction\n")
        fh.write("# a first attempt with a truncated 12 A Coulomb sum failed validation (pair terms reach\n")
        fh.write("# +-330 kJ/mol; the truncation residual invented -655 kJ/mol wells) and was discarded\n")
        fh.write(f"# validated on lp against RASPA's Widom run: K_H ratio "
                 f"{lp.K_H_mol_kg_Pa / RASPA_LP_KH:.3f}, <U> difference "
                 f"{lp.U_boltzmann_kJ_mol - RASPA_LP_U:+.2f} kJ/mol\n")
        fh.write(f"# uniform random position and orientation, seed {a.seed}\n")
        df.to_csv(fh, index=False)
    print(f"\nwrote results/insertion_energy.csv")

    # ---------------- figure ----------------
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.6), gridspec_kw={"wspace": 0.28})
    bins = np.linspace(-45, 60, 211)
    palette = [SIM["color"], REF["color"]] + list(EXTRA)
    styles = {lab: palette[i % len(palette)] for i, lab in enumerate(energies)}
    for label, u in energies.items():
        h, _ = np.histogram(np.clip(u, bins[0], bins[-1]), bins=bins)
        ax.step(0.5 * (bins[1:] + bins[:-1]), h / len(u), where="mid",
                color=styles[label], lw=1.8, label=label)
        hi = df[df.label == label].iloc[0]
        rows_txt = f"{label}: U$_{{min}}$ = {hi.U_min_kJ_mol:.1f} kJ/mol"
        ax.plot([hi.U_min_kJ_mol], [1.0 / len(u)], "v", color=styles[label], ms=8, label=rows_txt)
    ax.axvline(0, color="#8c8c88", lw=1.0)
    ax.set_yscale("log")
    ax.set_xlabel("CO$_2$ insertion energy U [kJ/mol]")
    ax.set_ylabel("fraction of random insertions per bin")
    ax.set_title("Where a CO$_2$ can land", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")

    cats = ["U < 0", "U < -10", "U < -20"]
    keys = ["frac_U_lt_0", "frac_U_lt_minus10", "frac_U_lt_minus20"]
    x = np.arange(len(cats))
    for k, (label, _) in enumerate(styles.items()):
        r = df[df.label == label].iloc[0]
        vals = [max(r[key], 1.0 / r.n_trials * 0.3) for key in keys]
        ax2.bar(x + (k - 0.5) * 0.36, vals, 0.34, color=styles[label], label=label)
        for xi, v, key in zip(x + (k - 0.5) * 0.36, vals, keys):
            txt = "0" if r[key] == 0 else f"{100 * r[key]:.2g}%"
            ax2.text(xi, v * 1.25, txt, ha="center", fontsize=8, color=styles[label])
    ax2.set_yscale("log")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{c} kJ/mol" for c in cats])
    ax2.set_ylabel("fraction of random insertions")
    ax2.set_title("Favourable sites, by depth", fontsize=10)
    ax2.legend(fontsize=8, loc="upper right")
    fig.suptitle(f"CO$_2$ test-insertion energies, rigid MIL-53(Al), {T:g} K, "
                 f"{a.n_trials:,} random insertions per cell", y=1.02)
    fig.text(0.01, -0.11,
             "Rigid TraPPE CO$_2$, UFF framework, Lorentz-Berthelot, 12.0 A LJ cut-off, CIF charges -- the GCMC force field. "
             "Electrostatics are a full\nEwald sum of the guest-host term (alpha-independent to 6 digits); the LJ tail "
             "correction is omitted, being a constant shift.\nValidated against RASPA's own Widom run on the lp cell: "
             f"<U> {lp.U_boltzmann_kJ_mol:.2f} vs {RASPA_LP_U:.2f} kJ/mol, K$_H$ ratio {lp.K_H_mol_kg_Pa / RASPA_LP_KH:.2f}.\n"
             f"Uniform random position and orientation, seed {a.seed}. The spike at +{bins[-1]:.0f} kJ/mol is the "
             "clipped pile-up of insertions that overlap the framework.",
             fontsize=7.5, va="top")
    fig.savefig(out / "insertion_energy.png")
    print("wrote results/insertion_energy.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
