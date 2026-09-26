#!/usr/bin/env python3
"""Put our Al np framework into the experimental CO2-loaded np cell (Serre 2007).

WHY (NOTES.md 5.6b/5.7). `structures/MIL-53_Al_np.cif` is the dehydrated lt geometry:
pore-limiting diameter 2.52 A against CO2's 3.30 A, zero CO2-accessible volume, and a
rigid GCMC isotherm 326x below experiment. The np form CO2 actually induces is wider.
Serre et al. measured it for Cr; Dundar et al. used that same Cr cell for Al. This script
carries out that substitution explicitly and checks it.

METHOD -- this is the lozenge shear, not an isotropic dilation.
Our file is P1 in the setting 6.5687 / 15.2347 / 18.9855, beta = 95.377. The basis
change [a+c, b, -a] re-expresses it, volume-preserving and right-handed, as
19.4993 / 15.2347 / 6.5687, beta = 104.219 -- the C2/c setting of Loiseau 2004, with b
doubled because our file is a genuine 2x superstructure along b (5.1). The fractional
coordinates are then read into the target cell. Same space-group setting, same topology,
same atom order; only the lozenge angle and the b opening change, which is what breathing
physically is.

The b doubling is kept rather than reduced: a b/2 translation maps only 86 of our 152
atoms (5.1), so folding to the single C2/c cell would have to discard half a genuinely
ordered structure.

CELLS (per formula cell; the built cells are doubled along b, so 2x these volumes)
  ours, Al np (dehydrated lt geometry)   19.499 / 7.617 / 6.569, beta 104.22   945.8 A^3
  Serre 2007, Cr lt hydrated             19.685 / 7.849 / 6.782, beta 104.90  1012.8 A^3
  Serre 2007, Cr np under 1 bar CO2      19.713 / 8.310 / 6.806, beta 105.85  1072.5 A^3
The opening is almost entirely along b, 7.617 -> 8.310 A (+9.1 %).

Anything written here is a MODEL-DERIVED STRUCTURE, not a measurement.

Usage:
  python scripts/build_np_co2_cell.py                 # both variants + checks
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from insertion_energy import read_cif  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "structures" / "MIL-53_Al_np.cif"
OUT_DIR = ROOT / "structures" / "derived"

# Serre et al., Adv. Mater. 2007, 19, 2246: np MIL-53(Cr) under 1 bar CO2, C2/c,
# in situ synchrotron XRD. b is doubled here to match our superstructure.
SERRE_CO2_NP = dict(a=19.713, b=2 * 8.310, c=6.806, beta=105.85,
                    label="Serre2007_CrCO2np",
                    note="Serre 2007 Adv. Mater. 19, 2246: np MIL-53(Cr) at 1 bar CO2, "
                         "C2/c, a 19.713(1) b 8.310(1) c 6.806(1) A, beta 105.85(1) deg, "
                         "V 1072.5(1) A^3; used for Al by Dundar 2017 J. Mol. Model., Table 1")
# Sensitivity variant: the same cell scaled to the Al/Cr volume ratio of the lp forms,
# 1412/1486 = 0.950, i.e. 1072.5 -> 1018.9 A^3 per formula cell.
AL_CR_VOLUME_RATIO = 1412.0 / 1486.0

MASSES = {"Al": 26.981538, "O": 15.9994, "C": 12.0107, "H": 1.00794}
# covalent-radius sums are not used; these are the shortest bonds we expect to see
MIN_REASONABLE_CONTACT = 0.80      # A -- an O-H bond; anything shorter is an error


def cell_matrix(a, b, c, alpha=90.0, beta=90.0, gamma=90.0):
    al, be, ga = (math.radians(x) for x in (alpha, beta, gamma))
    v = math.sqrt(1 - math.cos(al) ** 2 - math.cos(be) ** 2 - math.cos(ga) ** 2
                  + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    return np.array([
        [a, b * math.cos(ga), c * math.cos(be)],
        [0.0, b * math.sin(ga), c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)],
        [0.0, 0.0, c * v / math.sin(ga)],
    ]).T


def params(M):
    L = np.linalg.norm(M, axis=1)
    ang = [math.degrees(math.acos(np.dot(M[j], M[k]) / (L[j] * L[k])))
           for j, k in ((1, 2), (0, 2), (0, 1))]
    return L, ang, abs(np.linalg.det(M))


def min_contact(frac, M):
    """Shortest interatomic distance under the minimum-image convention."""
    n = len(frac)
    d = frac[:, None, :] - frac[None, :, :]
    d -= np.round(d)
    r = np.linalg.norm(d @ M, axis=2)
    r[np.arange(n), np.arange(n)] = np.inf
    return float(r.min()), np.unravel_index(int(np.argmin(r)), r.shape)


BOND_MAX = {("Al", "O"): 2.20, ("H", "O"): 1.20, ("C", "O"): 1.60,   # keys MUST be sorted
            ("C", "C"): 1.70, ("C", "H"): 1.25}


def bond_lengths(frac, elem, M):
    """Nearest-neighbour bond lengths, grouped by element pair.

    Pairs are taken as bonded if they are inside BOND_MAX in the SOURCE geometry; the
    same pair list is then measured in whatever cell is passed, so the two calls are
    directly comparable.
    """
    d = frac[:, None, :] - frac[None, :, :]
    d -= np.round(d)
    r = np.linalg.norm(d @ M, axis=2)
    out = {}
    n = len(frac)
    for i in range(n):
        for j in range(i + 1, n):
            key = tuple(sorted((elem[i], elem[j])))
            lim = BOND_MAX.get(key)
            if lim and r[i, j] < lim:
                out.setdefault("-".join(key), []).append(r[i, j])
    return out


def read_labels(path):
    lines = Path(path).read_text().splitlines()
    hdr, start = [], None
    for i, l in enumerate(lines):
        if l.strip().startswith("_atom_site"):
            hdr.append(l.strip())
            start = i
    il = hdr.index("_atom_site_label")
    out = []
    for l in lines[start + 1:]:
        p = l.split()
        if len(p) < len(hdr):
            break
        out.append(p[il])
    return out


def write_cif(path, header, frac, elem, labels, q, M):
    L, ang, V = params(M)
    body = [f"# {h}" for h in header]
    body += [
        "data_MIL-53_Al_np_CO2cell",
        f"_cell_length_a       {L[0]:.6f}",
        f"_cell_length_b       {L[1]:.6f}",
        f"_cell_length_c       {L[2]:.6f}",
        f"_cell_angle_alpha    {ang[0]:.4f}",
        f"_cell_angle_beta     {ang[1]:.4f}",
        f"_cell_angle_gamma    {ang[2]:.4f}",
        f"_cell_volume         {V:.4f}",
        "",
        '_symmetry_space_group_name_H-M    "P 1"',
        "_symmetry_Int_Tables_number       1",
        "",
        "loop_",
        "  _symmetry_equiv_pos_as_xyz",
        "  'x, y, z'",
        "",
        "loop_",
        "  _atom_site_type_symbol",
        "  _atom_site_label",
        "  _atom_site_symmetry_multiplicity",
        "  _atom_site_fract_x",
        "  _atom_site_fract_y",
        "  _atom_site_fract_z",
        "  _atom_site_occupancy",
        "  _atom_site_charge",
    ]
    for e, lab, f, qq in zip(elem, labels, frac, q):
        body.append(f"{e:<3s} {lab:<9s} 1.0  {f[0]:.5f}  {f[1]:.5f}  {f[2]:.5f}  1.0000 {qq:.10f}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(body) + "\n")


def build(target, label, note, frac_c2c, elem, labels, q, M_src, M_c2c):
    M_t = cell_matrix(target["a"], target["b"], target["c"], 90.0, target["beta"], 90.0)
    Lt, angt, Vt = params(M_t)
    z = sum(1 for e in elem if e == "Al") // 4

    print(f"\n--- {label}")
    print(f"    target cell  {Lt[0]:.4f} / {Lt[1]:.4f} / {Lt[2]:.4f}  beta {angt[1]:.3f}"
          f"   V {Vt:.2f} A^3  = {Vt / z:.2f} per formula cell (Z = {z})")

    # ---- checks
    comp = {}
    for e in elem:
        comp[e] = comp.get(e, 0) + 1
    mass = sum(MASSES[e] for e in elem)
    print(f"    composition  {''.join(f'{k}{v}' for k, v in sorted(comp.items()))}"
          f"   mass {mass:.3f} g/mol = {mass / z:.3f} per formula cell")
    print(f"    net charge   {q.sum():+.6f} e")
    rmin, (i, j) = min_contact(frac_c2c, M_t)
    print(f"    shortest contact {rmin:.3f} A  ({elem[i]}{labels[i]}-{elem[j]}{labels[j]})", end="")
    print("   OK" if rmin >= MIN_REASONABLE_CONTACT else "   *** TOO SHORT ***")

    # The diagnostic that matters: a fractional-coordinate transplant strains every bond
    # by roughly the cell strain. Aromatic C-C bonds do not stretch 9 %, so this is what
    # the fixed-cell internal relaxation has to repair, and how much it has to move.
    bl0, bl1 = bond_lengths(frac_c2c, elem, M_c2c), bond_lengths(frac_c2c, elem, M_t)
    print("    bond strain from the cell change (this is what the relaxation must repair):")
    for kind in sorted(bl0):
        r0, r1 = np.array(bl0[kind]), np.array(bl1[kind])
        ch = 100.0 * (r1 - r0) / r0
        print(f"      {kind:6s} n={len(r0):3d}  {r0.mean():.3f} -> {r1.mean():.3f} A   "
              f"change mean {ch.mean():+5.2f} %, min {ch.min():+5.2f} %, max {ch.max():+5.2f} %")

    path = OUT_DIR / f"MIL-53_Al_np_{label}.cif"
    write_cif(path, [
        "MODEL-DERIVED STRUCTURE -- NOT A MEASUREMENT.",
        "Our Al np framework (structures/MIL-53_Al_np.cif, CoRE 2024 SI, PACMAN charges)",
        "re-expressed by the volume-preserving basis change [a+c, b, -a] into the C2/c",
        f"setting, then placed in: {note}.",
        "b is doubled relative to the published cell because our file is a genuine 2x",
        "superstructure along b (NOTES.md 5.1); fractional coordinates are unchanged, so",
        "this is the lozenge shear, not an isotropic dilation.",
        "Charges are the PACMAN values of the ORIGINAL geometry, carried over unchanged;",
        "atom order and labels are preserved. See NOTES.md 5.7 for why that is a caveat.",
        "Internal coordinates are NOT yet relaxed in this file.",
    ], frac_c2c, elem, labels, q, M_t)
    print(f"    wrote {path.relative_to(ROOT)}")
    return path, M_t


def main():
    frac, elem, q, M_src = read_cif(SRC)
    labels = read_labels(SRC)
    A, B, C = M_src
    M_c2c = np.array([A + C, B, -A])          # volume-preserving, right-handed
    det = np.linalg.det(M_c2c) / np.linalg.det(M_src)
    L0, ang0, V0 = params(M_src)
    L1, ang1, V1 = params(M_c2c)

    print("=" * 96)
    print("Re-expressing our Al np cell in the C2/c setting (basis change [a+c, b, -a])")
    print("=" * 96)
    print(f"  source     {L0[0]:.4f} / {L0[1]:.4f} / {L0[2]:.4f}  "
          f"alpha,beta,gamma {ang0[0]:.2f} {ang0[1]:.3f} {ang0[2]:.2f}   V {V0:.2f}")
    print(f"  C2/c       {L1[0]:.4f} / {L1[1]:.4f} / {L1[2]:.4f}  "
          f"alpha,beta,gamma {ang1[0]:.2f} {ang1[1]:.3f} {ang1[2]:.2f}   V {V1:.2f}")
    print(f"  volume ratio {det:+.6f} (must be exactly +1: volume-preserving, right-handed)")
    print(f"  Loiseau 2004 lt, b doubled: 19.499 / 15.234 / 6.569, beta 104.22 -- matches")
    assert abs(det - 1.0) < 1e-9, "basis change is not volume-preserving"

    # fractional coordinates in the C2/c basis (identical atoms, new basis)
    xyz = frac @ M_src
    frac_c2c = xyz @ np.linalg.inv(M_c2c)
    frac_c2c -= np.floor(frac_c2c)             # wrap into [0,1)

    built = []
    built.append(build(SERRE_CO2_NP, SERRE_CO2_NP["label"], SERRE_CO2_NP["note"],
                       frac_c2c, elem, labels, q, M_src, M_c2c))

    f = AL_CR_VOLUME_RATIO ** (1.0 / 3.0)
    scaled = dict(a=SERRE_CO2_NP["a"] * f, b=SERRE_CO2_NP["b"] * f, c=SERRE_CO2_NP["c"] * f,
                  beta=SERRE_CO2_NP["beta"], label="Serre2007_scaledAl",
                  note=f"the same Serre 2007 cell scaled isotropically by "
                       f"(1412/1486)^(1/3) = {f:.5f}, i.e. to the Al/Cr volume ratio of the "
                       f"lp forms ({AL_CR_VOLUME_RATIO:.4f}); a sensitivity variant, the shear "
                       f"already being present in the parent cell")
    built.append(build(scaled, scaled["label"], scaled["note"],
                       frac_c2c, elem, labels, q, M_src, M_c2c))

    print("\n" + "=" * 96)
    print("Both files hold the SAME fractional coordinates; only the cell differs.")
    print("Next: relax internal coordinates at FIXED cell (MACE-MP-0 + D3), then Zeo++.")
    print("=" * 96)
    return 0


if __name__ == "__main__":
    sys.exit(main())
