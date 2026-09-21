#!/usr/bin/env bash
# Helium void fraction of the rigid framework (Widom He insertion, 298 K, 500k cycles).
# Result goes into the isotherm inputs as HeliumVoidFraction; see NOTES.md.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RASPA_DIR="${RASPA_DIR:-$HOME/micromamba/envs/gcmc-mil53}"; export RASPA_DIR
python "$ROOT/scripts/make_inputs.py" --cif "$ROOT/structures/MIL-53_Al_lp.cif" --helium-run
cd "$ROOT/runs/helium_void_fraction"
"$RASPA_DIR/bin/simulate" simulation.input > simulate.log 2>&1
grep -A7 "Average Widom Rosenbluth factor" Output/System_0/*.data | tail -3
