#!/usr/bin/env python3
"""Variable-cell relaxation of MIL-53(Al) with the MACE-MP-0 foundation potential.

WHY. The rigid np isotherm is near-empty because the np file is the dehydrated lt
geometry, whose pore-limiting diameter is 2.52 A against CO2's 3.30 A (NOTES.md 5.6b).
The narrow-pore form CO2 actually induces is wider, and no CO2-loaded MIL-53 structure
is downloadable. A flexible potential can relax the cell WITH CO2 inside and predict
that structure -- if the potential can be trusted for this material.

VALIDATION FIRST, AND IT IS NOT A FORMALITY. MACE-MP-0 is trained on Materials Project
PBE data, which has no dispersion correction, and the lp/np balance of MIL-53 is
famously dispersion-controlled: the np form is held closed by van der Waals attraction
between the terephthalate linkers. A potential without dispersion should over-expand the
np cell badly. So both variants are run and reported:

    plain MACE-MP-0        (PBE-like, no dispersion)
    MACE-MP-0 + D3         (Grimme D3 through torch-dftd, as mace_mp(dispersion=True))

Targets, both measured from our own CIFs:
    empty lp   1412.0 A^3 per formula cell   (SABVUN, the production isotherm cell)
    empty np    945.8 A^3 per formula cell   (CoRE 2024 SI, the hydrated-lt geometry)

Pass criterion, set BEFORE running: both cells within a few percent. The np cell is the
demanding one -- reproducing lp alone proves nothing, because almost any potential keeps
an open framework open.

EXTERNAL BENCHMARK (Stavitski et al., Langmuir 2011, 27, 3970, SI). Their MIL-53(Al)
cells were optimised at PLAIN PBE, with dispersion added only for internal relaxation at
fixed cell -- so they are a published "PBE without dispersion" reference, not truth:
    empty lp 1437.7    empty np  970.2    CO2@lp 1509.5    CO2@np 1232.9 A^3
Plain MACE-MP-0 is trained on Materials Project PBE, also without dispersion, so it
should land near THOSE numbers, not near experiment. MACE-MP-0 + D3 should instead
approach experiment. That contrast is the cleanest test of whether D3 is doing its job,
and it is reported alongside the experimental comparison.

Anything produced here is a MODEL-DERIVED STRUCTURE, never a measurement, and is written
to structures/derived/ with that in its header.

Usage:
  python scripts/mace_relax.py --validate                  # both cells, both variants
  python scripts/mace_relax.py --cif <path> [--co2-per-formula-cell 3] [--dispersion]
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
V_LP_TARGET = 1412.0          # A^3 per formula cell
V_NP_TARGET = 945.8
PASS_PCT = 5.0                # "within a few percent", fixed before the run
# Stavitski 2011 Langmuir 27, 3970 SI -- plain-PBE (dispersion-free) cell volumes
PBE_REF = {"lp": 1437.7, "np": 970.2, "CO2@lp": 1509.5, "CO2@np": 1232.9}


def formula_cells(atoms):
    """How many Al4C32H20O20 units the cell holds (framework atoms only)."""
    n_al = sum(1 for s in atoms.get_chemical_symbols() if s == "Al")
    return n_al // 4


def load(cif):
    from ase.io import read
    a = read(str(cif))
    a.set_pbc(True)
    return a


def make_calc(dispersion, dtype="float64"):
    from mace.calculators import mace_mp
    return mace_mp(model="medium", default_dtype=dtype, device="cpu",
                   dispersion=dispersion)


def relax_guests_only(atoms, calc, n_frame, fmax=0.2, steps=200, logfile=None):
    """Relax the inserted molecules at FIXED cell and FIXED framework.

    Inserting CO2 into a channel whose pore-limiting diameter is 2.5 A necessarily
    starts with hard overlaps. Releasing everything at once from such a start lets the
    optimiser take one enormous step and leaves MACE far outside anything it was
    trained on. So the guests are relieved first, with the framework held, and only
    then is the cell released.
    """
    from ase.constraints import FixAtoms
    from ase.optimize import FIRE
    atoms.calc = calc
    atoms.set_constraint(FixAtoms(indices=list(range(n_frame))))
    f0 = float(np.abs(atoms.get_forces()).max())
    opt = FIRE(atoms, logfile=logfile)
    opt.run(fmax=fmax, steps=steps)
    f1 = float(np.abs(atoms.get_forces()).max())
    atoms.set_constraint()
    return f0, f1, opt.get_number_of_steps()


def relax(atoms, calc, fmax=0.05, steps=400, logfile=None, traj=None,
          constant_volume=False, fixed_cell=False):
    """Relax. fixed_cell=True moves atoms only; otherwise the cell is a degree of freedom.

    fixed_cell is what 5.7 needs: the cell comes from experiment (Serre 2007) and must
    not move, while the internal coordinates -- strained 2-4 % by the fractional-
    coordinate transplant -- have to be repaired.
    """
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    atoms.calc = calc
    v0 = atoms.get_volume()
    e0 = atoms.get_potential_energy()
    f0 = float(np.abs(atoms.get_forces()).max())
    t0 = time.time()
    target = atoms if fixed_cell else FrechetCellFilter(atoms, constant_volume=constant_volume)
    opt = FIRE(target, logfile=logfile, trajectory=traj)
    opt.run(fmax=fmax, steps=steps)
    return {
        "V0_A3": v0, "V_A3": atoms.get_volume(),
        "E0_eV": e0, "E_eV": atoms.get_potential_energy(),
        "fmax0_eV_A": f0,
        "fmax_eV_A": float(np.abs(atoms.get_forces()).max()),
        "steps": opt.get_number_of_steps(),
        "converged": bool(opt.converged()),
        "wall_s": time.time() - t0,
    }


def _min_image_dist(points, others, cell):
    """Minimum-image distance from each point to the nearest of `others`."""
    d = points[:, None, :] - others[None, :, :]
    fd = np.linalg.solve(cell.T, d.reshape(-1, 3).T).T.reshape(len(points), len(others), 3)
    fd -= np.round(fd)
    return np.linalg.norm(fd @ cell, axis=2).min(axis=1)


def add_co2(atoms, per_formula_cell, grid=28):
    """Insert CO2 greedily at the points of greatest clearance.

    Spacing the molecules along one cell axis is wrong: the np cell is only 6.57 A long
    in the channel direction, so N molecules laid along it would sit ~6.6/N A apart --
    closer than CO2 is long. Instead each molecule is placed at the grid point furthest
    from EVERY atom already present, framework and previously-placed CO2 alike, and its
    axis is turned to the direction of greatest local clearance. That finds whatever
    channels the cell actually has without assuming how many there are.

    Note on units: this clearance is measured to the nearest NUCLEUS, so it is larger
    than the Zeo++ free-sphere radius, which is measured to the nearest vdW SURFACE.
    A 2.8 A nucleus clearance against an O (r_vdW 1.52 A) is a ~1.3 A free radius,
    consistent with Zeo++'s 1.41 A for this cell (NOTES.md 5.6b).
    """
    from ase import Atoms
    z = formula_cells(atoms)
    n = per_formula_cell * z
    cell = np.array(atoms.get_cell())
    g = (np.arange(grid) + 0.5) / grid
    fr = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    pts = fr @ cell

    # candidate CO2 axes: the 13 distinct directions of a cubic +-1/0 lattice
    axes = []
    for v in np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 0, 1],
                       [1, 0, -1], [0, 1, 1], [0, 1, -1], [1, 1, 1], [1, 1, -1],
                       [1, -1, 1], [-1, 1, 1]], float):
        axes.append(v / np.linalg.norm(v))
    axes = np.array(axes)

    placed, clearances = [], []
    occupied = atoms.get_positions().copy()
    for _ in range(n):
        d = _min_image_dist(pts, occupied, cell)
        c = pts[int(np.argmax(d))]
        # orient so both O atoms are as far from everything as possible
        best_u, best_score = axes[0], -1e9
        for u in axes:
            ends = np.array([c + 1.16 * u, c - 1.16 * u])
            score = _min_image_dist(ends, occupied, cell).min()
            if score > best_score:
                best_score, best_u = score, u
        co2 = Atoms("OCO", positions=[c + 1.16 * best_u, c, c - 1.16 * best_u])
        atoms += co2
        occupied = np.vstack([occupied, co2.get_positions()])
        placed.append(c)
        clearances.append(float(d.max()))
    return atoms, n, clearances


BOND_MAX = {("Al", "O"): 2.20, ("H", "O"): 1.20, ("C", "O"): 1.60,   # keys MUST be sorted
            ("C", "C"): 1.70, ("C", "H"): 1.25}


def bond_stats(atoms):
    """Bond-length statistics by element pair, for the charge-transfer decision."""
    sym = atoms.get_chemical_symbols()
    n = len(atoms)
    out = {}
    for i in range(n - 1):                     # i = n-1 has no partners to its right
        d = atoms.get_distances(i, list(range(i + 1, n)), mic=True)
        for off, r in enumerate(d):
            j = i + 1 + off
            key = tuple(sorted((sym[i], sym[j])))
            lim = BOND_MAX.get(key)
            if lim and r < lim:
                out.setdefault("-".join(key), []).append(float(r))
    return {k: {"n": len(v), "mean": float(np.mean(v)), "min": float(np.min(v)),
                "max": float(np.max(v))} for k, v in out.items()}


def describe(tag, res, target, z):
    v = res["V_A3"] / z
    v0 = res["V0_A3"] / z
    err = 100.0 * (v - target) / target if target else float("nan")
    print(f"  {tag:34s} V {v0:7.1f} -> {v:7.1f} A^3/formula cell", end="")
    if target:
        print(f"   target {target:6.1f}   error {err:+6.2f} %", end="")
    print(f"   [{res['steps']} steps, fmax {res['fmax_eV_A']:.3f} eV/A, "
          f"{'converged' if res['converged'] else 'NOT converged'}, {res['wall_s'] / 60:.1f} min]")
    return err


def write_cif(atoms, path, header_lines):
    """Write a relaxed structure -- WITHOUT charges, and say so loudly.

    ASE's CIF writer has no _atom_site_charge column, so anything written here is
    unusable for GCMC: RASPA reads framework charges from that column and would
    silently run an uncharged framework. To keep the charges, recover the structure
    from the trajectory with scripts/traj_to_cif.py, which re-attaches labels and
    charges by index after asserting the atom order is unchanged.
    """
    from ase.io import write
    path = path.with_name(path.stem + "_NOCHARGES" + path.suffix)
    path.parent.mkdir(parents=True, exist_ok=True)
    write(str(path), atoms, format="cif")
    body = path.read_text()
    path.write_text("".join(f"# {l}\n" for l in header_lines)
                    + "# NO _atom_site_charge COLUMN -- NOT USABLE FOR GCMC.\n"
                    + "# Use scripts/traj_to_cif.py on the trajectory to keep the charges.\n"
                    + body)
    print(f"  wrote {path} (geometry only)")
    print(f"  *** this file has NO CHARGES; for GCMC use scripts/traj_to_cif.py on "
          f"the trajectory instead ***")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--cif")
    ap.add_argument("--co2-per-formula-cell", type=int, default=0)
    ap.add_argument("--dispersion", action="store_true")
    ap.add_argument("--fmax", type=float, default=0.05)
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--label", default="")
    ap.add_argument("--fixed-cell", action="store_true",
                    help="relax internal coordinates only; the cell is held exactly")
    ap.add_argument("--no-stage", action="store_true",
                    help="skip the fixed-cell guest-only pre-relaxation (not advised when inserting)")
    a = ap.parse_args(argv)

    out = ROOT / "results"
    out.mkdir(exist_ok=True)

    if a.validate:
        jobs = [("lp", ROOT / "structures" / "MIL-53_Al_lp.cif", V_LP_TARGET),
                ("np", ROOT / "structures" / "MIL-53_Al_np.cif", V_NP_TARGET)]
        rows = []
        print("=" * 100)
        print("VALIDATION -- MACE-MP-0 variable-cell relaxation of the two EMPTY cells")
        print("=" * 100)
        for disp in (False, True):
            name = "MACE-MP-0 + D3" if disp else "MACE-MP-0 (no dispersion)"
            print(f"\n{name}")
            try:
                calc = make_calc(disp)
            except Exception as exc:
                print(f"  could not build calculator: {exc}")
                continue
            for tag, cif, target in jobs:
                atoms = load(cif)
                z = formula_cells(atoms)
                res = relax(atoms, calc, fmax=a.fmax, steps=a.steps,
                            logfile=str(ROOT / "logs" / f"mace_{tag}_{'d3' if disp else 'nod3'}.log"))
                err = describe(f"{tag} (empty, Z={z})", res, target, z)
                v = res["V_A3"] / z
                pbe = PBE_REF[tag]
                err_pbe = 100.0 * (v - pbe) / pbe
                print(f"  {'':34s} vs Stavitski plain-PBE {pbe:6.1f} A^3: {err_pbe:+6.2f} %")
                rows.append({"variant": name, "cell": tag, "dispersion": disp,
                             "V_per_formula_cell_A3": v, "target_A3": target,
                             "error_pct": err, "pbe_ref_A3": pbe, "error_vs_pbe_pct": err_pbe,
                             "z": z, **res})
        import pandas as pd
        df = pd.DataFrame(rows)
        csv = out / "mace_validation.csv"
        with open(csv, "w") as fh:
            fh.write("# MACE-MP-0 variable-cell relaxation of the EMPTY MIL-53(Al) cells\n")
            fh.write("# targets are our own CIFs: lp 1412.0, np 945.8 A^3 per formula cell\n")
            fh.write(f"# FIRE + FrechetCellFilter, fmax = {a.fmax} eV/A, float64, CPU\n")
            fh.write(f"# pass criterion fixed before the run: both cells within {PASS_PCT} %\n")
            fh.write("# pbe_ref_A3: Stavitski 2011 Langmuir 27, 3970 SI, plain-PBE cells (no dispersion)\n")
            fh.write("# plain MACE-MP-0 is PBE-trained and should track pbe_ref; +D3 should track target\n")
            df.to_csv(fh, index=False)
        print(f"\nwrote {csv.relative_to(ROOT)}")

        print("\n" + "=" * 100)
        for disp in (False, True):
            g = df[df.dispersion == disp]
            if not len(g):
                continue
            worst = g.error_pct.abs().max()
            name = g.variant.iloc[0]
            ok = worst <= PASS_PCT
            print(f"  {name:28s} worst error {worst:5.2f} %  ->  "
                  f"{'PASS' if ok else 'FAIL'} (criterion {PASS_PCT} %)")
        print("=" * 100)
        return 0

    if not a.cif:
        ap.error("give --validate or --cif")
    atoms = load(ROOT / a.cif if not Path(a.cif).is_absolute() else a.cif)
    z = formula_cells(atoms)
    label = a.label or Path(a.cif).stem
    n_co2, clear = 0, []
    if a.co2_per_formula_cell:
        atoms, n_co2, clear = add_co2(atoms, a.co2_per_formula_cell)
        print(f"inserted {n_co2} CO2 ({a.co2_per_formula_cell} per formula cell, Z = {z})")
        print("  clearance of each CO2 centre to the nearest framework atom [A]: "
              + ", ".join(f"{c:.2f}" for c in clear))
    calc = make_calc(a.dispersion)
    stage = None
    if n_co2 and not a.no_stage:
        n_frame = len(atoms) - 3 * n_co2
        f0, f1, ns = relax_guests_only(atoms, calc, n_frame,
                                       logfile=str(ROOT / "logs" / f"mace_{label}_stage1.log"))
        stage = {"fmax_before_eV_A": f0, "fmax_after_eV_A": f1, "steps": ns}
        print(f"  stage 1 (guests only, cell and framework fixed): "
              f"fmax {f0:.1f} -> {f1:.3f} eV/A in {ns} steps")
    res = relax(atoms, calc, fmax=a.fmax, steps=a.steps, fixed_cell=a.fixed_cell,
                logfile=str(ROOT / "logs" / f"mace_{label}.log"),
                traj=str(ROOT / "logs" / f"mace_{label}.traj"))
    describe(f"{label} + {n_co2} CO2", res, None, z)
    print(f"  V per formula cell: {res['V0_A3'] / z:.1f} -> {res['V_A3'] / z:.1f} A^3")
    print(f"  bond lengths after relaxation (compare with the original np geometry):")
    for kind, r in sorted(bond_stats(atoms).items()):
        print(f"    {kind:6s} n={r['n']:3d}  mean {r['mean']:.3f} A  "
              f"min {r['min']:.3f}  max {r['max']:.3f}")

    disp_tag = "D3" if a.dispersion else "noD3"
    write_cif(atoms, ROOT / "structures" / "derived" / f"{label}_{n_co2}CO2_MACE_{disp_tag}.cif", [
        "MODEL-DERIVED STRUCTURE -- NOT A MEASUREMENT.",
        f"Variable-cell relaxation with MACE-MP-0 ({'with' if a.dispersion else 'without'} D3 dispersion),",
        f"FIRE + FrechetCellFilter, fmax = {a.fmax} eV/A, float64, CPU.",
        f"Start: {Path(a.cif).name}, {n_co2} CO2 inserted ({a.co2_per_formula_cell} per formula cell).",
        f"V per formula cell {res['V0_A3'] / z:.1f} -> {res['V_A3'] / z:.1f} A^3 in {res['steps']} steps,",
        f"final fmax {res['fmax_eV_A']:.4f} eV/A, converged = {res['converged']}.",
        "Do not compare this with experiment as if it were a refined structure.",
    ])
    (out / f"mace_{label}_{n_co2}CO2_{disp_tag}.json").write_text(
        json.dumps({**res, "z": z, "n_co2": n_co2, "stage1": stage,
                    "V_per_formula_cell_A3": res["V_A3"] / z,
                    "clearance_A": clear, "dispersion": a.dispersion}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
