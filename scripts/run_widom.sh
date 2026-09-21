#!/usr/bin/env bash
# Zero-coverage Henry coefficient of CO2 by Widom test insertion (304 K, 20k cycles).
# This is the true zero-coverage limit OF THE FORCE FIELD; see NOTES.md for why it
# must not be compared directly with Coudert's Langmuir K_lp.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RASPA_DIR="${RASPA_DIR:-$HOME/micromamba/envs/gcmc-mil53}"; export RASPA_DIR
python "$ROOT/scripts/make_inputs.py" --cif "$ROOT/structures/MIL-53_Al_lp.cif" --widom-run --cycles 20000
cd "$ROOT/runs/henry_widom_CO2"
"$RASPA_DIR/bin/simulate" simulation.input > simulate.log 2>&1
grep -E "\[CO2\] Average (Henry coefficient|Widom Rosenbluth-weight)|<U_gh>_1-<U_h>_0:" Output/System_0/*.data
