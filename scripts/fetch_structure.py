#!/usr/bin/env python3
"""Download the MIL-53(Al) lp structure, verify it, and build the working CIF.

Source: CoRE MOF 2014 DDEC database (Nazarian, Camp & Sholl, Chem. Mater. 2016,
28, 785), Zenodo record 3986573, file core-mof-1.0-ddec.tar, member
./re-labeled/SABVUN_clean.cif. The download is checked against the md5 published
by Zenodo before anything is extracted.

Then, in order:
  structures/original/SABVUN_clean.cif      the untouched original (read-only)
  structures/SABVUN_clean_primitive.cif     copy with only the non-CIF metadata
                                            line commented out
  structures/MIL-53_Al_lp.cif               exact change of basis to the 76-atom
                                            conventional cell (verified; see
                                            scripts/primitive_to_conventional.py)

Phase identity is checked by cell, never by refcode: see scripts/check_cif.py and
NOTES.md ("Structure provenance"). Nothing here is silent -- every step prints.

Usage: python scripts/fetch_structure.py [--keep-tar]
"""
import argparse
import hashlib
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
URL = "https://zenodo.org/records/3986573/files/core-mof-1.0-ddec.tar?download=1"
MD5 = "a2e3578003968739640b6faeed361299"          # as published by Zenodo
MEMBER = "./re-labeled/SABVUN_clean.cif"
SRC_MD5 = "7b79c613f3c5da8c31a051df553575ff"      # of the extracted CIF


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-tar", action="store_true")
    a = ap.parse_args(argv)
    struct = ROOT / "structures"
    (struct / "original").mkdir(parents=True, exist_ok=True)
    original = struct / "original" / "SABVUN_clean.cif"

    if original.exists() and md5sum(original) == SRC_MD5:
        print(f"{original.relative_to(ROOT)} already present (md5 verified)")
    else:
        tar_path = struct / "core-mof-1.0-ddec.tar"
        print(f"downloading {URL}")
        urllib.request.urlretrieve(URL, tar_path)
        got = md5sum(tar_path)
        if got != MD5:
            sys.exit(f"md5 mismatch: expected {MD5}, got {got}. Refusing to use this download.")
        print(f"md5 {got} matches Zenodo")
        with tarfile.open(tar_path) as tf:
            member = tf.getmember(MEMBER)
            with tf.extractfile(member) as src, open(original, "wb") as dst:
                dst.write(src.read())
        if not a.keep_tar:
            tar_path.unlink()
        original.chmod(0o444)
        print(f"extracted {MEMBER} -> {original.relative_to(ROOT)} (md5 {md5sum(original)})")

    primitive = struct / "SABVUN_clean_primitive.cif"
    text = original.read_text().splitlines()
    primitive.write_text("\n".join([
        "# Copy of structures/original/SABVUN_clean.cif (CoRE MOF 2014 DDEC, Zenodo 3986573). ONLY change: the next line",
        "# (a Python-dict metadata line, not CIF) is commented out. Its _srcid 'CAKYAQ_clean' is a database artefact; see NOTES.md.",
        "# " + text[0], *text[1:]]) + "\n")
    print(f"wrote {primitive.relative_to(ROOT)} (metadata line commented out)")

    print("\nchange of basis to the conventional cell:")
    rc = subprocess.call([sys.executable, str(ROOT / "scripts" / "primitive_to_conventional.py"),
                          str(primitive), str(struct / "MIL-53_Al_lp.cif")])
    if rc != 0:
        sys.exit("transformation checks FAILED -- see the output above; do not patch, use the primitive cell")
    print("\nnow run: python scripts/check_cif.py structures/MIL-53_Al_lp.cif")
    return 0


if __name__ == "__main__":
    sys.exit(main())
