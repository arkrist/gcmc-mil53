#!/usr/bin/env python3
"""Parse a RASPA2 output file (Output/System_0/*.data).

Extracts, per adsorbed component:
  * absolute and excess loading in molecules/unit cell and mol/kg framework,
    each with RASPA's block-average error bar;
  * acceptance rates of the swap insertion, swap deletion, reinsertion,
    translation and rotation moves;
  * the run conditions (T, p, fugacity coefficient, unit cells, framework mass).
It also exposes `loading_trace()` (instantaneous N vs cycle from the periodic
status prints) for convergence plots.

Error bars
----------
RASPA splits the production run into 5 blocks and prints
    avg +/- 2.776 * s / sqrt(5)
where s is the sample standard deviation of the 5 block averages and 2.776 is
Student's t (95 %, 4 d.o.f.) -- see RASPA2 src/statistics.h,
ERROR_CONFIDENCE_INTERVAL_95. So "+/-" is a 95 % confidence half-width, NOT a
standard error. We keep RASPA's number in the `*_err95` columns and also
report `*_sem` = err95 / 2.776.

Usage:
    python scripts/parse_raspa_output.py path/to/output.data [...]  [--csv out.csv]
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

T95 = 2.776  # Student t, 95 %, 4 dof (RASPA uses 5 blocks)
FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def _first(pattern, text, cast=float, default=None):
    m = re.search(pattern, text, re.MULTILINE)
    return cast(m.group(1)) if m else default


def parse_conditions(text):
    ncells = [_first(rf"Number of unitcells \[{ax}\]:\s*(\d+)", text, int, 1) for ax in "abc"]
    n_uc = ncells[0] * ncells[1] * ncells[2]
    box_mass = _first(rf"Framework Mass:\s*({FLOAT})", text)  # g/mol, whole simulation box
    return {
        "framework": _first(r"Framework name:\s*(\S+)", text, str),
        "T_K": _first(rf"External temperature:\s*({FLOAT})", text),
        "unit_cells": "x".join(map(str, ncells)),
        "n_unit_cells": n_uc,
        "framework_mass_box_g_per_mol": box_mass,
        "framework_mass_uc_g_per_mol": box_mass / n_uc if box_mass else None,
        "cycles_init": _first(r"Number of initializing cycles:\s*(\d+)", text, int),
        "cycles_prod": _first(r"^Number of cycles:\s*(\d+)", text, int),
        "finished": "Simulation finished" in text,
    }


def parse_component_meta(text):
    """Fugacity coefficient and partial fugacity per component (from the header)."""
    meta = {}
    for m in re.finditer(r"^Component (\d+) \[(.+?)\] \(Adsorbate molecule\)", text, re.MULTILINE):
        idx, name = int(m.group(1)), m.group(2)
        chunk = text[m.end(): m.end() + 6000]
        meta[name] = {
            "component": idx,
            "partial_pressure_Pa": _first(rf"Partial pressure:\s*({FLOAT})", chunk),
            "fugacity_coeff": _first(rf"Fugacity coefficient:\s*({FLOAT})", chunk),
            "fugacity_Pa": _first(rf"Partial fugacity:\s*({FLOAT})", chunk),
            "bulk_density_kg_m3": _first(rf"Density of the bulk fluid phase:\s*({FLOAT})", chunk),
        }
    return meta


def parse_loadings(text):
    """Loadings from the final 'Number of molecules:' block-average section."""
    start = text.rfind("Number of molecules:\n====")
    if start < 0:
        return {}
    section = text[start:]
    out = {}
    for m in re.finditer(r"^Component \d+ \[(.+?)\]\n-+\n(.*?)(?=^Component \d+ \[|\n\n\n)",
                         section, re.MULTILINE | re.DOTALL):
        name, body = m.group(1), m.group(2)
        rec = {}
        for kind in ("absolute", "excess"):
            for unit_tag, key in (("molecules/unit cell", "molec_uc"), ("mol/kg framework", "mol_kg")):
                pat = rf"Average loading {kind} \[{re.escape(unit_tag)}\]\s+({FLOAT})\s+\+/-\s+({FLOAT})"
                mm = re.search(pat, body)
                if mm:
                    val, err = float(mm.group(1)), float(mm.group(2))
                    rec[f"{kind}_{key}"] = val
                    rec[f"{kind}_{key}_err95"] = err
                    rec[f"{kind}_{key}_sem"] = err / T95
        out[name] = rec
    return out


def parse_acceptance(text):
    """Acceptance ratios (fraction, 0-1) of MC moves per component."""
    out = {}

    def put(name, key, val):
        out.setdefault(name, {})[key] = val

    # Swap / reinsertion: "Component [CO2] total tried: N ... accepted: M (P [%])"
    for header, key in (("swap addition", "acc_insertion"),
                        ("swap deletion", "acc_deletion"),
                        ("Reinsertion", "acc_reinsertion")):
        sec = re.search(rf"Performance of the {header} move:\n=+\n(.*?)\n\n", text, re.DOTALL)
        if not sec:
            continue
        for m in re.finditer(rf"Component \[(.+?)\] total tried:\s*({FLOAT}).*?accepted:\s*({FLOAT})",
                             sec.group(1)):
            tried, acc = float(m.group(2)), float(m.group(3))
            put(m.group(1), key, acc / tried if tried > 0 else float("nan"))
            put(m.group(1), key.replace("acc_", "tried_"), tried)

    # Translation / rotation: per-direction totals and successes -> overall ratio
    for header, key in (("translation", "acc_translation"), ("rotation", "acc_rotation")):
        sec = re.search(rf"Performance of the {header} move:\n=+\n(.*?)\n\n", text, re.DOTALL)
        if not sec:
            continue
        for m in re.finditer(rf"Component \d+ \[(.+?)\]\n\s*total\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\n"
                             rf"\s*succesfull\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})", sec.group(1)):
            tot = sum(float(m.group(i)) for i in (2, 3, 4))
            suc = sum(float(m.group(i)) for i in (5, 6, 7))
            put(m.group(1), key, suc / tot if tot > 0 else float("nan"))
    return out


def parse_file(path):
    """One row per component with conditions, loadings and acceptance rates."""
    text = Path(path).read_text(errors="replace")
    cond = parse_conditions(text)
    meta = parse_component_meta(text)
    loads = parse_loadings(text)
    acc = parse_acceptance(text)
    # Total pressure = sum of partial pressures. (The header "External Pressure"
    # line lists every pressure of a multi-pressure input, so it is ambiguous.)
    cond["p_Pa"] = sum(m.get("partial_pressure_Pa") or 0.0 for m in meta.values()) or None
    rows = []
    for name in meta or loads:
        row = {"file": str(path), "molecule": name, **cond, **meta.get(name, {}),
               **loads.get(name, {}), **acc.get(name, {})}
        rows.append(row)
    return pd.DataFrame(rows)


def loading_trace(path):
    """Instantaneous number of molecules in the box vs cycle, from status prints.

    Returns DataFrame(stage, cycle, molecule, N_box). stage is 'init' or 'prod'.
    Resolution is RASPA's PrintEvery.
    """
    text = Path(path).read_text(errors="replace")
    rows = []
    for m in re.finditer(r"^(\[Init\] )?Current cycle: (\d+) out of \d+\n(.*?)(?=^(?:\[Init\] )?Current cycle:|^Finishing)",
                         text, re.MULTILINE | re.DOTALL):
        stage = "init" if m.group(1) else "prod"
        for c in re.finditer(r"^Component \d+ \((.+?)\), current number of integer/fractional/reaction molecules: (\d+)/",
                             m.group(3), re.MULTILINE):
            rows.append({"stage": stage, "cycle": int(m.group(2)), "molecule": c.group(1),
                         "N_box": int(c.group(2))})
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--csv", help="write the table to this CSV")
    args = ap.parse_args(argv)
    df = pd.concat([parse_file(f) for f in args.files], ignore_index=True)
    if args.csv:
        df.to_csv(args.csv, index=False)
    cols = ["framework", "molecule", "T_K", "p_Pa", "fugacity_coeff",
            "absolute_molec_uc", "absolute_molec_uc_err95", "absolute_mol_kg", "absolute_mol_kg_err95",
            "excess_mol_kg", "excess_mol_kg_err95", "acc_insertion", "acc_deletion", "acc_reinsertion",
            "acc_translation", "acc_rotation", "finished"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(df[[c for c in cols if c in df.columns]].to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
