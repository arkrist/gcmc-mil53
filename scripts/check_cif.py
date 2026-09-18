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

FIRST GATE -- phase identity by cell (not by refcode). Accepted: MIL-53(Al) lp,
Loiseau et al., Chem. Eur. J. 2004, 10, 1373: Imma (No. 74), a = 6.608,
b = 16.675, c = 12.813 A, V ~ 1412 A^3, Al4C32H20O20 (76 atoms), 832.4 g/mol.
Rejected: as-synthesised Pnma (17.129, 6.628, 12.182 A; bdc in the pores) and
the monoclinic narrow-pore hydrated forms (Cc / P21/c). The comparison is done on
the SORTED axis lengths, so a P1 cell with permuted axes is still recognised
(and reported as permuted). Any mismatch -> exit status 1, nothing else is trusted.

Usage: python scripts/check_cif.py structures/MIL-53_Al_lp.cif
"""
import itertools
import re
import shlex
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from make_inputs import R_CUT, replication  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LP_CELL = (6.608, 16.675, 12.813)            # A, Imma, Loiseau 2004 (lp / ht form)
LP_COMPOSITION = {"Al": 4, "C": 32, "H": 20, "O": 20}
LP_MASS = 832.4                               # g/mol per cell
AS_SYNTH_CELL = (17.129, 6.628, 12.182)      # Pnma, as-synthesised -- must be rejected
REL_TOL_LENGTH = 0.01                         # 1 % per axis
TOL_ANGLE = 0.5                               # deg
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


def implied_conventional(lengths, angles, tol_deg=0.3):
    """Smallest orthogonal cell spanned by short lattice vectors (n_i in -2..2).

    Recognises a primitive / reduced setting of an orthorhombic lattice (the
    CoRE 2014 DDEC files store I-centred MIL-53 as a 38-atom primitive cell).
    Returns (sorted lengths, volume, multiplicity vs file cell) or None.
    """
    from make_inputs import cell_matrix
    h = cell_matrix(lengths, angles)
    vfile = abs(np.linalg.det(h))
    vecs = sorted(((np.linalg.norm(np.array(n) @ h), np.array(n) @ h)
                   for n in itertools.product(range(-2, 3), repeat=3) if any(n)), key=lambda x: x[0])[:80]
    best = None
    for (l1, v1), (l2, v2), (l3, v3) in itertools.combinations(vecs, 3):
        cosines = [abs(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b)) for a, b in ((v2, v3), (v1, v3), (v1, v2))]
        if all(c < np.sin(np.radians(tol_deg)) for c in cosines):
            vol = abs(np.linalg.det(np.array([v1, v2, v3])))
            if vol > 1 and (best is None or vol < best[1] - 1e-6):
                best = (sorted([l1, l2, l3]), vol, vol / vfile)
    return best


def ff_labels():
    txt = (ROOT / "forcefield" / "force_field_mixing_rules.def").read_text().splitlines()
    return {ln.split()[0] for ln in txt if re.search(r"\blennard-jones\b", ln, re.I)}


def main(cif):
    text = Path(cif).read_text()
    ok = True
    lengths, angles, widths, n = replication(cif)
    h_vol = np.prod(lengths) * np.sqrt(1 - sum(np.cos(np.radians(angles)) ** 2)
                                       + 2 * np.prod(np.cos(np.radians(angles))))
    sg = re.search(r"_symmetry_space_group_name_H-M\s+['\"]?([^'\"\n]+)['\"]?", text) or \
        re.search(r"_space_group_name_H-M_alt\s+'?([^'\n]+)'?", text)
    print(f"file            : {cif}")
    print(f"cell a b c [A]  : {lengths[0]:.4f} {lengths[1]:.4f} {lengths[2]:.4f}")
    print(f"angles [deg]    : {angles[0]:.3f} {angles[1]:.3f} {angles[2]:.3f}")
    print(f"volume [A^3]    : {h_vol:.2f}")
    print(f"space group     : {sg.group(1).strip() if sg else 'not given'}")

    # ---- gate 1: is this the lp cell?
    rel = [abs(x - y) / y for x, y in zip(sorted(lengths), sorted(LP_CELL))]
    rel_as = [abs(x - y) / y for x, y in zip(sorted(lengths), sorted(AS_SYNTH_CELL))]
    ortho = all(abs(g - 90.0) < TOL_ANGLE for g in angles)
    print("expected lp    : a b c = 6.608 16.675 12.813 A, 90/90/90, V ~ 1412 A^3")
    print("rel. deviation : " + ", ".join(f"{100 * r:.2f} %" for r in rel) + " (sorted axes)")
    cell_ok = ortho and max(rel) < REL_TOL_LENGTH
    if cell_ok:
        order = "same order" if all(abs(x - y) / y < REL_TOL_LENGTH for x, y in zip(lengths, LP_CELL)) \
            else "axes permuted with respect to Imma a,b,c (harmless, reported for the record)"
        print(f"CELL CHECK     : PASS -- MIL-53(Al) lp ({order})")
    else:
        conv = None if ortho else implied_conventional(lengths, angles)
        if conv and max(abs(x - y) / y for x, y in zip(conv[0], sorted(LP_CELL))) < REL_TOL_LENGTH:
            print(f"CELL CHECK     : LATTICE MATCHES lp, but the file is a PRIMITIVE/REDUCED setting "
                  f"(implied conventional cell {' '.join(f'{x:.4f}' for x in conv[0])} A, V = {conv[1]:.1f} A^3, "
                  f"{conv[2]:.0f}x the file cell). Using it as-is or transforming it is a decision. STOP.")
            return 1
        why = "not orthorhombic (narrow-pore monoclinic?)" if not ortho else \
            "matches the AS-SYNTHESISED Pnma cell" if max(rel_as) < REL_TOL_LENGTH else \
            "does not match the lp cell"
        if ortho and abs(h_vol - 1412 / 2) < 15:
            why += "; volume ~ half of lp -> possibly a primitive cell of the body-centred lattice"
        print(f"CELL CHECK     : FAIL -- {why}. STOP.")
        return 1

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
    if comp != LP_COMPOSITION or abs(mass - LP_MASS) > 0.5:
        print(f"COMPOSITION    : FAIL -- expected {LP_COMPOSITION} ({LP_MASS} g/mol), got {comp}. STOP.")
        ok = False
    else:
        print("COMPOSITION    : PASS -- Al4C32H20O20, 76 atoms (includes the mu-OH hydrogens)")

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
