#!/usr/bin/env python3
"""Build LAMMPS data files for the cross-code check against RASPA.

Writes three systems with the SAME box and force field:
  A = framework + CO2, B = framework only, C = CO2 only
so every energy term can be split as cross = A - B - C. That isolates the
host-guest interaction exactly, including the Ewald reciprocal part (the
reciprocal energy is a quadratic form in the charges, so the self terms cancel),
and it removes the framework-framework energy, which RASPA never computes for a
rigid host.

Coordinates come from a RASPA restart file (12 decimals), so both codes see the
identical configuration. The framework is the same CIF, replicated as in RASPA.

Force field (identical to forcefield/*.def):
  UFF framework LJ, TraPPE CO2, Lorentz-Berthelot, 12.0 A cut-off, tail
  corrections, Ewald. Energies are converted K -> kcal/mol for LAMMPS `real`
  units with kB = 1/503.2197 kcal/mol/K.

Usage:
  python scripts/make_lammps_inputs.py --restart runs/crosscheck/raspa_singlepoint/Restart/System_0/restart_* \
      --cif structures/MIL-53_Al_lp.cif --cells 4 2 2 --out runs/crosscheck/lammps
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from check_cif import atom_site_loop  # noqa: E402
from make_inputs import cell_matrix, read_cell  # noqa: E402

K_TO_KCAL = 1.0 / 503.2197
# type order used everywhere below (1-based in LAMMPS)
TYPES = ["Al", "O", "C", "H", "C_co2", "O_co2"]
MASS = {"Al": 26.981538, "O": 15.9994, "C": 12.0107, "H": 1.00794, "C_co2": 12.0107, "O_co2": 15.9994}
EPS_K = {"Al": 254.126, "O": 30.193, "C": 52.838, "H": 22.142, "C_co2": 27.0, "O_co2": 79.0}
SIGMA = {"Al": 4.0082, "O": 3.1181, "C": 3.4309, "H": 2.5711, "C_co2": 2.80, "O_co2": 3.05}
Q_CO2 = {"C_co2": 0.70, "O_co2": -0.35}
# RASPA's CO2.def atom order within a molecule: 0 = O, 1 = C, 2 = O
CO2_ORDER = ["O_co2", "C_co2", "O_co2"]


def read_framework(cif, cells):
    text = Path(cif).read_text()
    cols, rows = atom_site_loop(text)
    li, qi = cols.index("_atom_site_label"), cols.index("_atom_site_charge")
    si = cols.index("_atom_site_type_symbol")
    xi, yi, zi = (cols.index(f"_atom_site_fract_{c}") for c in "xyz")
    h = cell_matrix(*read_cell(cif))
    atoms = []
    for r in rows:
        f = np.array([float(r[xi]), float(r[yi]), float(r[zi])])
        for i in range(cells[0]):
            for j in range(cells[1]):
                for k in range(cells[2]):
                    atoms.append((r[si], (f + [i, j, k]) @ h, float(r[qi])))
    return atoms, h


def read_restart_adsorbates(path):
    mols = {}
    for m in re.finditer(r"Adsorbate-atom-position:\s+(\d+)\s+(\d+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)",
                         Path(path).read_text()):
        mols.setdefault(int(m.group(1)), {})[int(m.group(2))] = np.array([float(m.group(i)) for i in (3, 4, 5)])
    return [[mols[i][j] for j in sorted(mols[i])] for i in sorted(mols)]


def write_data(path, box, framework, co2, note):
    """box = 3x3 cell matrix (rows a, b, c) with a along x, b in the xy plane."""
    lx, xy, ly, xz, yz, lz = box[0, 0], box[1, 0], box[1, 1], box[2, 0], box[2, 1], box[2, 2]
    atoms, bonds = [], []
    for el, xyz, q in framework:
        atoms.append((TYPES.index(el) + 1, 0, q, xyz))
    for mi, mol in enumerate(co2, start=1):
        first = len(atoms) + 1
        for el, xyz in zip(CO2_ORDER, mol):
            atoms.append((TYPES.index(el) + 1, mi, Q_CO2[el], xyz))
        # bonds only define the 1-2 / 1-3 exclusions (bond_style zero)
        bonds += [(first, first + 1), (first + 1, first + 2)]
    lines = [f"# {note}", "",
             f"{len(atoms)} atoms", f"{len(bonds)} bonds", "",
             f"{len(TYPES)} atom types", "1 bond types", "",
             f"0.0 {lx:.12f} xlo xhi", f"0.0 {ly:.12f} ylo yhi", f"0.0 {lz:.12f} zlo zhi",
             f"{xy:.12f} {xz:.12f} {yz:.12f} xy xz yz", "", "Masses", ""]
    lines += [f"{i + 1} {MASS[t]:.6f}   # {t}" for i, t in enumerate(TYPES)]
    lines += ["", "Atoms # full", ""]
    for i, (t, mol, q, xyz) in enumerate(atoms, start=1):
        lines.append(f"{i} {mol} {t} {q: .8f} {xyz[0]: .8f} {xyz[1]: .8f} {xyz[2]: .8f}")
    if bonds:
        lines += ["", "Bonds", ""]
        lines += [f"{i} 1 {a} {b}" for i, (a, b) in enumerate(bonds, start=1)]
    Path(path).write_text("\n".join(lines) + "\n")
    return len(atoms)


def write_in(path, data, note):
    eps = {t: EPS_K[t] * K_TO_KCAL for t in TYPES}
    lines = [f"# {note}", "units real", "atom_style full", "boundary p p p",
             "bond_style zero", "special_bonds lj/coul 0.0 0.0 0.0   # exclude intramolecular 1-2 and 1-3",
             "pair_style lj/cut/coul/long 12.0 12.0",
             "kspace_style ewald 1.0e-6", "", f"read_data {data}", "",
             "pair_modify mix arithmetic tail yes   # Lorentz-Berthelot + analytic tail correction", ""]
    lines += [f"pair_coeff {i + 1} {i + 1} {eps[t]:.10f} {SIGMA[t]:.6f}   # {t}" for i, t in enumerate(TYPES)]
    lines += ["", "bond_coeff 1", "neighbor 2.0 bin", "",
              "thermo_style custom step evdwl ecoul elong etail pe",
              "thermo_modify format float %20.8f", "run 0", ""]
    Path(path).write_text("\n".join(lines) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--restart", required=True)
    ap.add_argument("--cif", required=True)
    ap.add_argument("--cells", nargs=3, type=int, default=[4, 2, 2])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    framework, h = read_framework(a.cif, a.cells)
    box = np.array([h[0] * a.cells[0], h[1] * a.cells[1], h[2] * a.cells[2]])
    co2 = read_restart_adsorbates(a.restart)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"framework atoms: {len(framework)}, CO2 molecules: {len(co2)}")
    print("box (LAMMPS): lx ly lz = {:.6f} {:.6f} {:.6f}, xy xz yz = {:.3e} {:.3e} {:.3e}".format(
        box[0, 0], box[1, 1], box[2, 2], box[1, 0], box[2, 0], box[2, 1]))
    for tag, fw, gas in (("A_framework_co2", framework, co2), ("B_framework", framework, []),
                         ("C_co2", [], co2)):
        n = write_data(out / f"{tag}.data", box, fw, gas, f"cross-check {tag}")
        write_in(out / f"{tag}.in", f"{tag}.data", f"single-point energy, {tag}")
        print(f"  {tag}: {n} atoms")
    return 0


if __name__ == "__main__":
    sys.exit(main())
