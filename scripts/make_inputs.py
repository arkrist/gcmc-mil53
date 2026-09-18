#!/usr/bin/env python3
"""Generate one RASPA run directory per pressure from templates/simulation.input.template.

  runs/<tag>/p_<pressure in Pa>/
      simulation.input, <framework>.cif, pseudo_atoms.def,
      force_field_mixing_rules.def, force_field.def, CO2.def

Force-field files and the CIF are copied, not linked, so every run directory is
a complete record of what was simulated.

Replication: the minimum-image convention requires every perpendicular width of
the simulation box to exceed 2 x cutoff. For a triclinic cell,
    width_i = V / |a_j x a_k|,
and n_i = ceil(2 * r_cut / width_i).

Usage:
  python scripts/make_inputs.py --cif structures/MIL-53_Al_lp.cif --helium-vf 0.xx \
         [--tag production] [--cycles 50000 --init 10000] [--pressures 1e4 5e4 ...]
  python scripts/make_inputs.py --cif ... --replication-only     # just print the replication
  python scripts/make_inputs.py --cif ... --helium-run [--cycles 500000]   # runs/helium_void_fraction/
  python scripts/make_inputs.py --cif ... --widom-run  [--cycles 50000]    # runs/henry_widom_CO2/
"""
import argparse
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FF_FILES = ["pseudo_atoms.def", "force_field_mixing_rules.def", "force_field.def", "CO2.def", "helium.def"]
PRESSURES_PA = [1e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 3e6, 5e6]  # 0.01 ... 5 MPa
R_CUT = 12.0  # A, LJ and real-space Coulomb cutoff (12.8 left only 0.03 A margin along c; see NOTES.md)


def read_cell(cif):
    txt = Path(cif).read_text()

    def get(key):
        m = re.search(rf"^\s*{key}\s+([-+]?\d*\.?\d+)", txt, re.MULTILINE)
        if not m:
            sys.exit(f"{key} not found in {cif}")
        return float(m.group(1))

    return [get(f"_cell_length_{x}") for x in "abc"], [get(f"_cell_angle_{x}") for x in ("alpha", "beta", "gamma")]


def cell_matrix(lengths, angles_deg):
    a, b, c = lengths
    al, be, ga = np.radians(angles_deg)
    va = np.array([a, 0, 0])
    vb = np.array([b * np.cos(ga), b * np.sin(ga), 0])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    vc = np.array([cx, cy, np.sqrt(c**2 - cx**2 - cy**2)])
    return np.array([va, vb, vc])


def replication(cif, r_cut=R_CUT):
    lengths, angles = read_cell(cif)
    h = cell_matrix(lengths, angles)
    vol = abs(np.linalg.det(h))
    widths = [vol / np.linalg.norm(np.cross(h[(i + 1) % 3], h[(i + 2) % 3])) for i in range(3)]
    n = [math.ceil(2 * r_cut / w) for w in widths]
    # ceil of an exact multiple would give width*n == 2 r_cut; require strictly larger
    n = [k + 1 if k * w <= 2 * r_cut else k for k, w in zip(n, widths)]
    return lengths, angles, widths, n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cif", required=True)
    ap.add_argument("--helium-vf", type=float, help="helium void fraction (from runs/helium_void_fraction)")
    ap.add_argument("--tag", default="production")
    ap.add_argument("--cycles", type=int, default=50000)
    ap.add_argument("--init", type=int, default=10000)
    ap.add_argument("--print-every", type=int, default=500)
    ap.add_argument("--temperature", type=float, default=303.0)
    ap.add_argument("--pressures", type=float, nargs="+", default=PRESSURES_PA)
    ap.add_argument("--no-cif-charges", action="store_true",
                    help="do not read _atom_site_charge (only after an explicit decision, see NOTES.md)")
    ap.add_argument("--replication-only", action="store_true")
    ap.add_argument("--helium-run", action="store_true",
                    help="write runs/helium_void_fraction/ (Widom He insertion, 298 K) instead of the isotherm")
    ap.add_argument("--widom-run", action="store_true",
                    help="write runs/henry_widom_CO2/ (Widom CO2 insertion, Henry coefficient) instead of the isotherm")
    args = ap.parse_args(argv)

    cif = Path(args.cif)
    lengths, angles, widths, n = replication(cif)
    print(f"cell a,b,c = {lengths} A, angles = {angles} deg")
    print("perpendicular widths [A]: " + ", ".join(f"{w:.3f}" for w in widths))
    print(f"replication for r_cut = {R_CUT} A: {n[0]} x {n[1]} x {n[2]}  ->  box widths "
          + ", ".join(f"{k * w:.3f}" for k, w in zip(n, widths)) + f"  (need > {2 * R_CUT:.1f})")
    if args.replication_only:
        return 0
    if args.helium_run:
        d = ROOT / "runs" / "helium_void_fraction"
        d.mkdir(parents=True, exist_ok=True)
        tmpl = (ROOT / "templates" / "helium_void_fraction.input.template").read_text()
        cycles = args.cycles if args.cycles != 50000 else 500000  # RASPA example default for He
        (d / "simulation.input").write_text(tmpl.format(cycles=cycles, framework=cif.stem, cutoff=R_CUT,
                                                        unit_cells=" ".join(map(str, n))))
        shutil.copy(cif, d / cif.name)
        for f in FF_FILES:
            shutil.copy(ROOT / "forcefield" / f, d / f)
        print(f"wrote {d.relative_to(ROOT)}")
        return 0
    if args.widom_run:
        d = ROOT / "runs" / "henry_widom_CO2"
        d.mkdir(parents=True, exist_ok=True)
        tmpl = (ROOT / "templates" / "henry_widom.input.template").read_text()
        (d / "simulation.input").write_text(tmpl.format(
            cycles=args.cycles, framework=cif.stem, cutoff=R_CUT, temperature=args.temperature,
            unit_cells=" ".join(map(str, n))))
        shutil.copy(cif, d / cif.name)
        for f in FF_FILES:
            shutil.copy(ROOT / "forcefield" / f, d / f)
        print(f"wrote {d.relative_to(ROOT)}")
        return 0
    if args.helium_vf is None:
        sys.exit("--helium-vf is required (run the helium void-fraction calculation first)")

    template = (ROOT / "templates" / "simulation.input.template").read_text()
    for p in args.pressures:
        d = ROOT / "runs" / args.tag / f"p_{p:.0f}"
        d.mkdir(parents=True, exist_ok=True)
        text = template.format(
            cycles=args.cycles, init_cycles=args.init, print_every=args.print_every,
            use_cif_charges="no" if args.no_cif_charges else "yes",
            framework=cif.stem, unit_cells=" ".join(map(str, n)), helium_vf=args.helium_vf,
            temperature=args.temperature, pressure=f"{p:.6g}", cutoff=R_CUT)
        (d / "simulation.input").write_text(text)
        shutil.copy(cif, d / cif.name)
        for f in FF_FILES:
            shutil.copy(ROOT / "forcefield" / f, d / f)
        print(f"wrote {d.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
