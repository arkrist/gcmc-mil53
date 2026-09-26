#!/usr/bin/env python3
"""Write the last frame of an ASE trajectory back to a CIF, keeping labels and charges.

An ASE trajectory carries elements, positions and the cell, but not the CIF's
_atom_site_label or _atom_site_charge. Atom ORDER is preserved through a relaxation,
so those two columns are taken from the source CIF and re-attached by index. The
script asserts that the elements still match position by position, which is the check
that makes the re-attachment safe -- and the same check the charge-transfer decision
in NOTES.md 5.7 depends on.

Usage:
  python scripts/traj_to_cif.py <traj> <source.cif> <out.cif> [--header-file f]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from build_np_co2_cell import bond_lengths, params, write_cif  # noqa: E402
from insertion_energy import read_cif  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("traj")
    ap.add_argument("source_cif")
    ap.add_argument("out_cif")
    ap.add_argument("--header", action="append", default=[])
    a = ap.parse_args(argv)

    from ase.io import read
    frames = read(a.traj, index=":")
    last = frames[-1]
    print(f"{a.traj}: {len(frames)} frames, using the last one")

    frac0, elem0, q, M0 = read_cif(a.source_cif)
    labels = read_labels(a.source_cif)
    sym = last.get_chemical_symbols()
    assert len(sym) == len(elem0), f"atom count changed: {len(sym)} vs {len(elem0)}"
    bad = [i for i, (x, y) in enumerate(zip(sym, elem0)) if x != y]
    assert not bad, f"atom ORDER changed at indices {bad[:10]} -- charges cannot be re-attached"
    print(f"  atom order verified: {len(sym)} atoms, element-by-element identical to the source")

    M = np.array(last.get_cell())
    frac = last.get_scaled_positions(wrap=True)
    L, ang, V = params(M)
    L0, ang0, V0 = params(M0)
    print(f"  source cell   {L0[0]:.4f} / {L0[1]:.4f} / {L0[2]:.4f}  beta {ang0[1]:.3f}  V {V0:.2f}")
    print(f"  relaxed cell  {L[0]:.4f} / {L[1]:.4f} / {L[2]:.4f}  beta {ang[1]:.3f}  V {V:.2f}")

    # bond lengths: transplanted (strained) vs relaxed, both in the relaxed cell
    bl_in = bond_lengths(frac0, elem0, M0)
    bl_out = bond_lengths(frac, elem0, M)
    print("  bond lengths, transplanted -> relaxed:")
    for kind in sorted(bl_in):
        r0, r1 = np.array(bl_in[kind]), np.array(bl_out.get(kind, []))
        if not len(r1):
            print(f"    {kind:6s} n={len(r0):3d}  -> NONE FOUND after relaxation (bond broken?)")
            continue
        print(f"    {kind:6s} n={len(r0):3d}->{len(r1):3d}  {r0.mean():.3f} -> {r1.mean():.3f} A "
              f"(min {r1.min():.3f}, max {r1.max():.3f})")

    # displacement of each atom inside the fixed cell
    d = frac - frac0
    d -= np.round(d)
    disp = np.linalg.norm(d @ M, axis=1)
    print(f"  atom displacement during relaxation: mean {disp.mean():.3f} A, "
          f"max {disp.max():.3f} A ({elem0[int(np.argmax(disp))]}{labels[int(np.argmax(disp))]})")
    print(f"  net charge carried over: {q.sum():+.6f} e")

    write_cif(Path(a.out_cif), a.header or ["relaxed structure"], frac, elem0, labels, q, M)
    print(f"  wrote {a.out_cif}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
