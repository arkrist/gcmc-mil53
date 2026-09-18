#!/usr/bin/env python3
"""Exact change of basis: primitive (reduced) cell of the I-centred MIL-53(Al) lp
lattice -> conventional orthorhombic cell. Same atoms, same charges.

No atom is moved, added or re-derived from cell parameters: the conventional
cell vectors are integer combinations of the primitive ones,
    C = M P,  M integer, det(M) = 2,
and the 76 conventional atoms are the 38 primitive atoms plus their images
under primitive lattice translations that fall inside the conventional cell.

The output file is RE-READ and verified against the input (decision of
2026-09-19; if any check fails we fall back to the primitive cell, we do not
patch):
  1. det(M) = 2
  2. 76 atoms, no duplicates (min distance), no losses (each primitive atom
     appears exactly det(M) times)
  3. composition Al4C32H20O20 and 832.4 g/mol
  4. net charge 0 within 1e-4 e, and every atom keeps its source charge
  5. crystal density identical to the primitive file
  6. nearest-neighbour distance distribution unchanged

Axes are ordered as in Loiseau's Imma setting: a (6.6 A) < c (12.8 A) < b (16.7 A),
written as a = 6.6, b = 16.7, c = 12.8, right-handed.

Usage: python scripts/primitive_to_conventional.py IN_primitive.cif OUT_conventional.cif
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from check_cif import MASS, atom_site_loop  # noqa: E402
from make_inputs import cell_matrix, read_cell  # noqa: E402

TARGET = {"Al": 4, "C": 32, "H": 20, "O": 20}
NA_A3 = 1.66053906660  # g/mol per A^3 -> g/cm^3


def read_atoms(path):
    text = Path(path).read_text()
    cols, rows = atom_site_loop(text)
    li = cols.index("_atom_site_label")
    si = cols.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in cols else None
    xi, yi, zi = (cols.index(f"_atom_site_fract_{c}") for c in "xyz")
    qi = cols.index("_atom_site_charge")
    el = [r[si] if si is not None else "".join(ch for ch in r[li] if ch.isalpha()) for r in rows]
    frac = np.array([[float(r[xi]), float(r[yi]), float(r[zi])] for r in rows])
    q = np.array([float(r[qi]) for r in rows])
    lengths, angles = read_cell(path)
    return el, frac, q, cell_matrix(lengths, angles), lengths, angles


def find_basis(h):
    """Integer M (rows) with C = M h orthogonal, smallest volume, axes ordered a<b>c as Imma."""
    cands = [(np.array(n), np.array(n) @ h) for n in itertools.product(range(-2, 3), repeat=3) if any(n)]
    cands.sort(key=lambda x: np.linalg.norm(x[1]))
    cands = cands[:60]
    best = None
    for (n1, v1), (n2, v2), (n3, v3) in itertools.combinations(cands, 3):
        cos = [abs(a @ b) / np.linalg.norm(a) / np.linalg.norm(b) for a, b in ((v1, v2), (v1, v3), (v2, v3))]
        if max(cos) > 1e-4:
            continue
        vol = abs(np.linalg.det(np.array([v1, v2, v3])))
        if vol < 1e-6 or (best and vol >= best[0] - 1e-6):
            continue
        best = (vol, [(n1, v1), (n2, v2), (n3, v3)])
    trip = sorted(best[1], key=lambda x: np.linalg.norm(x[1]))  # short, medium, long
    a, c, b = trip                                                # Imma: a=6.6, b=16.7, c=12.8
    M = np.array([a[0], b[0], c[0]])
    if np.linalg.det(M) < 0:                                      # right-handed
        M[2] = -M[2]
    return M


def build(el, frac_p, q, h, M):
    C = M @ h
    Cinv = np.linalg.inv(C)
    cart_p = frac_p @ h
    out = []  # (element, frac_conv, charge, source index)
    for i in range(len(el)):
        for n in itertools.product(range(-3, 4), repeat=3):
            f = (cart_p[i] + np.array(n) @ h) @ Cinv
            f = np.where(np.abs(f - np.round(f)) < 1e-9, np.round(f), f)  # snap exact boundary values
            if np.all(f >= 0) and np.all(f < 1):
                out.append((el[i], f, q[i], i))
    return C, out


def write_cif(path, C, atoms, src_note):
    L = [np.linalg.norm(v) for v in C]
    ang = [np.degrees(np.arccos(C[j] @ C[k] / L[j] / L[k])) for j, k in ((1, 2), (0, 2), (0, 1))]
    counter = {}
    lines = [f"# {line}" for line in src_note.splitlines()]
    lines += ["data_MIL-53_Al_lp", "_symmetry_space_group_name_H-M   'P 1'", "_symmetry_Int_Tables_number   1",
              f"_cell_length_a   {L[0]:.8f}", f"_cell_length_b   {L[1]:.8f}", f"_cell_length_c   {L[2]:.8f}",
              f"_cell_angle_alpha   {ang[0]:.6f}", f"_cell_angle_beta   {ang[1]:.6f}", f"_cell_angle_gamma   {ang[2]:.6f}",
              "loop_", " _symmetry_equiv_pos_site_id", " _symmetry_equiv_pos_as_xyz", "  1  'x, y, z'",
              "loop_", " _atom_site_label", " _atom_site_type_symbol",
              " _atom_site_fract_x", " _atom_site_fract_y", " _atom_site_fract_z", " _atom_site_charge"]
    for e, f, qq, _ in sorted(atoms, key=lambda a: ("Al", "O", "C", "H").index(a[0])):
        counter[e] = counter.get(e, 0) + 1
        lines.append(f" {e}{counter[e]:<3d} {e:2s} {f[0]:.8f} {f[1]:.8f} {f[2]:.8f} {qq: .6f}")
    Path(path).write_text("\n".join(lines) + "\n")


def nn_distances(el, frac, h, rmax=5.0):
    """Per-atom sorted list of all neighbour distances < rmax (periodic)."""
    cart = frac @ h
    widths = [abs(np.linalg.det(h)) / np.linalg.norm(np.cross(h[(i + 1) % 3], h[(i + 2) % 3])) for i in range(3)]
    nrep = [int(np.ceil(rmax / w)) + 1 for w in widths]
    shifts = np.array([np.array(n) @ h for n in itertools.product(*[range(-k, k + 1) for k in nrep])])
    res = []
    for i in range(len(el)):
        d = np.linalg.norm(cart[None, :, :] + shifts[:, None, :] - cart[i], axis=2).ravel()
        res.append(np.sort(d[(d > 1e-6) & (d < rmax)]))
    return res


def main(src, dst):
    el_p, fr_p, q_p, h_p, _, _ = read_atoms(src)
    M = find_basis(h_p)
    detM = round(np.linalg.det(M))
    C, atoms = build(el_p, fr_p, q_p, h_p, M)
    note = (f"MIL-53(Al) lp, conventional orthorhombic cell. Exact change of basis of {Path(src).name}\n"
            f"(CoRE MOF 2014 DDEC, SABVUN_clean; DDEC charges unchanged). Basis matrix M (C = M P), rows:\n"
            + "\n".join(" ".join(f"{x:+d}" for x in row) for row in M) + f"\ndet(M) = {detM}. Written by scripts/primitive_to_conventional.py")
    write_cif(dst, C, atoms, note)

    # ---------------- verification on the RE-READ output ----------------
    el_c, fr_c, q_c, h_c, Lc, Ac = read_atoms(dst)
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        ok &= bool(cond)
        print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")

    print(f"M (rows = conventional a, b, c in primitive units):\n{M}")
    check("1 det(M) = 2", abs(np.linalg.det(M) - 2) < 1e-9, f"det = {np.linalg.det(M):.6f}")
    check("  cell", True, f"{Lc[0]:.4f} {Lc[1]:.4f} {Lc[2]:.4f} A, angles {Ac[0]:.4f} {Ac[1]:.4f} {Ac[2]:.4f}")
    check("  orthogonal", max(abs(a - 90) for a in Ac) < 1e-3, f"max |angle-90| = {max(abs(a - 90) for a in Ac):.2e} deg")

    counts = np.bincount([a[3] for a in atoms], minlength=len(el_p))
    dmin = min(np.min(d) for d in nn_distances(el_c, fr_c, h_c, rmax=2.0) if len(d)) if len(el_c) else 0
    check("2 atom count", len(el_c) == 76 and len(el_c) == detM * len(el_p), f"{len(el_c)} atoms = {detM} x {len(el_p)}")
    check("  no losses", np.all(counts == detM), f"every primitive atom appears exactly {detM}x (min {counts.min()}, max {counts.max()})")
    check("  no duplicates", dmin > 0.5, f"shortest interatomic distance {dmin:.4f} A")

    comp = {e: el_c.count(e) for e in sorted(set(el_c))}
    mass_c = sum(MASS[e] for e in el_c)
    check("3 composition", comp == TARGET, f"{comp}")
    check("  molar mass", abs(mass_c - 832.4) < 0.5, f"{mass_c:.3f} g/mol per conventional cell")

    check("4 net charge", abs(q_c.sum()) < 1e-4, f"{q_c.sum():+.2e} e (primitive file: {q_p.sum():+.2e} e)")
    same_q = all(abs(qq - q_p[i]) < 5e-7 for (_, _, qq, i) in atoms)
    check("  charges carried over", same_q, "each atom has its source atom's charge (to the 6 written decimals)")

    mass_p = sum(MASS[e] for e in el_p)
    rho_p = mass_p / abs(np.linalg.det(h_p)) * NA_A3
    rho_c = mass_c / abs(np.linalg.det(h_c)) * NA_A3
    check("5 density", abs(rho_c - rho_p) / rho_p < 1e-6, f"conventional {rho_c:.6f} vs primitive {rho_p:.6f} g/cm^3")

    nn_p = nn_distances(el_p, fr_p, h_p)
    nn_c = nn_distances(el_c, fr_c, h_c)
    # match each conventional atom to its source primitive atom and compare environments
    src = [a[3] for a in sorted(atoms, key=lambda a: ("Al", "O", "C", "H").index(a[0]))]
    maxdev, nmatch = 0.0, 0
    for k, i in enumerate(src):
        if len(nn_c[k]) != len(nn_p[i]):
            maxdev = float("inf")
            break
        maxdev = max(maxdev, float(np.max(np.abs(nn_c[k] - nn_p[i]))) if len(nn_p[i]) else 0.0)
        nmatch += 1
    check("6 NN distances", maxdev < 1e-5, f"all {nmatch} atoms: identical neighbour lists within 5 A, max |delta| = {maxdev:.2e} A")
    first_p = np.array([d[0] for d in nn_p])
    first_c = np.array([nn_c[k][0] for k in range(len(src))])
    print("  nearest-neighbour distance histogram (primitive x2 vs conventional), bins of 0.1 A:")
    bins = np.arange(0.8, 2.01, 0.1)
    hp, _ = np.histogram(np.repeat(first_p, detM), bins)
    hc, _ = np.histogram(first_c, bins)
    for lo, a, b in zip(bins[:-1], hp, hc):
        if a or b:
            print(f"    [{lo:.1f},{lo + 0.1:.1f}) A   primitive x2: {a:3d}   conventional: {b:3d}")
    check("  histogram equal", np.array_equal(hp, hc), "")
    print("\nALL CHECKS PASS" if ok else "\nA CHECK FAILED -> fall back to the primitive cell (option a); do not patch")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
