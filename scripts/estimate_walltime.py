#!/usr/bin/env python3
"""Estimate the wall time of the full isotherm from the short timing runs.

For each pressure point of runs/<tag>/ (default: timing) we know the total wall
time of a short run (driver.log) and N(cycle) from RASPA's status prints. One
RASPA cycle is max(20, N) MC moves, so the cost of a cycle is taken as
proportional to max(20, N). The production half of the test is used, where
the loading is closest to equilibrium:
    s_per_cycle(prod) = t_total * sum_prod max(20,N) / sum_all max(20,N) / n_prod
    t_full            = s_per_cycle(prod) * (init_full + prod_full)
If N is still rising in the last quarter of the test, the equilibrium loading
(and hence the cost) is higher than seen, and the estimate is flagged as a
LOWER bound. Assumption to keep in mind: the cost per move is taken as
constant, while in reality it grows mildly with N (more guest-guest pairs,
Ewald), so the high-pressure estimates lean low.

Usage: python scripts/estimate_walltime.py [tag] [--init 10000 --prod 50000]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import parse_raspa_output as pro  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("tag", nargs="?", default="timing")
    ap.add_argument("--init", type=int, default=10000)
    ap.add_argument("--prod", type=int, default=50000)
    a = ap.parse_args(argv)
    base = ROOT / "runs" / a.tag
    wall = {m.group(1): int(m.group(2)) for m in
            re.finditer(r"(p_\d+): done in (\d+) s", (base / "driver.log").read_text())}
    rows = []
    for name, t in sorted(wall.items(), key=lambda kv: float(kv[0][2:])):
        f = next((base / name / "Output" / "System_0").glob("*.data"))
        info = pro.parse_file(f).iloc[0]
        tr = pro.loading_trace(f)
        n_init, n_prod = int(info.cycles_init), int(info.cycles_prod)
        tr["x"] = np.where(tr.stage == "init", tr.cycle, n_init + tr.cycle)
        tr = tr.sort_values("x")
        w = np.maximum(20, tr.N_box.to_numpy())
        prod = (tr.stage == "prod").to_numpy()
        s_cycle = t * w[prod].sum() / w.sum() / n_prod
        last = tr.N_box.to_numpy()[-max(2, len(tr) // 4):]
        rising = last[-1] > last[0] + 2 * np.sqrt(max(last[0], 1))
        rows.append({"point": name, "p_bar": float(name[2:]) / 1e5, "test_wall_s": t,
                     "N_box_end": int(tr.N_box.iloc[-1]),
                     "molec_uc_end": tr.N_box.iloc[-1] / info.n_unit_cells,
                     "acc_insertion": info.acc_insertion,
                     "s_per_cycle_prod": s_cycle,
                     "full_h": s_cycle * (a.init + a.prod) / 3600,
                     "flag": "N still rising -> lower bound" if rising else ""})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 200, "display.float_format", "{:.3g}".format):
        print(df.to_string(index=False))
    print(f"\nfull isotherm, sequential ({a.init} + {a.prod} cycles per point): "
          f"{df.full_h.sum():.1f} h = {df.full_h.sum() / 24:.2f} days")
    df.to_csv(base / "walltime_estimate.csv", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
