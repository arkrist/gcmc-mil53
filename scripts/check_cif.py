#!/usr/bin/env python3
"""Pre-flight check of a framework CIF before any GCMC run.

Reports: cell parameters and volume, symmetry (P1 or not), atom count and
composition, whether partial charges (_atom_site_charge) are present, their net
charge and range, and the unit-cell molar mass. Then it checks that every atom
label (with digits stripped, as RASPA does with RemoveAtomNumberCodeFromLabel)
has a Lennard-Jones entry in forcefield/force_field_mixing_rules.def.

Why the last check matters: RASPA does NOT stop on a missing LJ pair; it prints
"WARNING: THERE ARE ATOM-PAIRS WITH NO VDW INTERACTION" and sets that
interaction to zero. This script exits with status 1 instead.

Usage: python scripts/check_cif.py structures/MIL-53_Al_lp.cif
"""
import re
import shlex
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from make_inputs import R_CUT, replication  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MASS = {"Al": 26.981538, "O": 15.9994, "C": 12.0107, "H": 1.00794}


def atom_site_loop(text):
    """Return (column names, rows) of the loop that contains _atom_site_label."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            cols, j = [], i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                cols.append(lines[j].strip().split()[0])
                j += 1
            rows = []
            while j < len(lines) and lines[j].strip() and not lines[j].strip().startswith(("_", "loop_", "#")):
                rows.append(shlex.split(lines[j]))
                j += 1
            if "_atom_site_label" in cols:
                return cols, rows
            i = j
        else:
            i += 1
    sys.exit("no _atom_site_label loop found")


def ff_labels():
    txt = (ROOT / "forcefield" / "force_field_mixing_rules.def").read_text().splitlines()
    return {ln.split()[0] for ln in txt if re.search(r"\blennard-jones\b", ln, re.I)}


def main(cif):
    text = Path(cif).read_text()
    ok = True
    lengths, angles, widths, n = replication(cif)
    h_vol = np.prod(lengths) * np.sqrt(1 - sum(np.cos(np.radians(angles)) ** 2)
                                       + 2 * np.prod(np.cos(np.radians(angles))))
    sg = re.search(r"_symmetry_space_group_name_H-M\s+'?([^'\n]+)'?", text) or \
        re.search(r"_space_group_name_H-M_alt\s+'?([^'\n]+)'?", text)
    print(f"file            : {cif}")
    print(f"cell a b c [A]  : {lengths[0]:.4f} {lengths[1]:.4f} {lengths[2]:.4f}")
    print(f"angles [deg]    : {angles[0]:.3f} {angles[1]:.3f} {angles[2]:.3f}")
    print(f"volume [A^3]    : {h_vol:.2f}")
    print(f"space group     : {sg.group(1).strip() if sg else 'not given'}")

    cols, rows = atom_site_loop(text)
    lab_i = cols.index("_atom_site_label")
    sym_i = cols.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in cols else None
    q_i = cols.index("_atom_site_charge") if "_atom_site_charge" in cols else None
    labels = [r[lab_i] for r in rows]
    elements = [r[sym_i] if sym_i is not None else re.sub(r"[^A-Za-z]", "", r[lab_i]) for r in rows]
    comp = {e: elements.count(e) for e in sorted(set(elements))}
    print(f"atoms in loop   : {len(rows)}  composition {comp}")
    p1 = sg and sg.group(1).strip().replace(" ", "") in ("P1", "P1(1)")
    if not p1:
        print("NOTE            : not P1 -- atom count above is the asymmetric unit; RASPA expands "
              "symmetry, so check 'Number of framework atoms' in a RASPA output.")
    mass = sum(MASS.get(e, float("nan")) for e in elements)
    print(f"unit-cell mass  : {mass:.3f} g/mol  (only meaningful for a P1 cell)")
    print(f"formula units   : {mass / 208.1045:.3f} x Al(OH)(O2C-C6H4-CO2), M = 208.1045 g/mol")

    if q_i is None:
        print("CHARGES         : ABSENT (no _atom_site_charge column)")
        ok = False
    else:
        q = np.array([float(r[q_i]) for r in rows])
        print(f"CHARGES         : present, net = {q.sum():+.6f} e, range [{q.min():+.4f}, {q.max():+.4f}]")
        for e in comp:
            qe = q[[k for k, el in enumerate(elements) if el == e]]
            print(f"   {e:3s} mean {qe.mean():+.4f}  min {qe.min():+.4f}  max {qe.max():+.4f}")
        if abs(q.sum()) > 1e-3:
            print("WARNING         : cell is not neutral; Ewald with a net charge adds a uniform background.")

    stripped = sorted({re.sub(r"\d+", "", lab) for lab in labels})
    known = ff_labels()
    missing = [s for s in stripped if s not in known]
    print(f"labels (no digits): {stripped}")
    if missing:
        print(f"FF CHECK        : FAIL -- no LJ parameters for {missing}")
        ok = False
    else:
        print("FF CHECK        : every label has an LJ entry")
    print(f"replication     : {n[0]} x {n[1]} x {n[2]} for r_cut = {R_CUT} A "
          f"(box widths {', '.join(f'{k * w:.2f}' for k, w in zip(n, widths))} A)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
