#!/usr/bin/env bash
# Run one of the RASPA2 bundled examples that ship with reference output, in
# examples/<name>/run/, and keep a copy of the reference output in
# examples/<name>/reference/ for comparison (scripts/compare_example.py).
#
# The simulation.input is copied UNCHANGED from the RASPA2 v2.0.50 source tree
# (the conda package does not include the examples/ directory), so cycle
# counts, force field and moves are exactly those of the reference run.
#
# Usage: scripts/run_example.sh <mfi_ch4|cubtc_co2|irmof1_mix>
# Env:   RASPA_SRC  path to the RASPA2 source checkout (default ~/src/RASPA2-2.0.50)
#        RASPA_DIR  conda env prefix holding bin/simulate and share/raspa
set -euo pipefail

case "${1:-}" in
  mfi_ch4)    ex="Basic/7_Adsorption_of_Methane_in_MFI" ;;
  cubtc_co2)  ex="Basic/8_Adsorption_of_CO2_in_CU-BTC" ;;
  irmof1_mix) ex="Non-Basic/1_Adsorption_of_Binary_Mixture_CO2_CH4_in_IRMOF-1" ;;
  *) echo "usage: $0 <mfi_ch4|cubtc_co2|irmof1_mix>" >&2; exit 1 ;;
esac

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RASPA_SRC="${RASPA_SRC:-$HOME/src/RASPA2-2.0.50}"
RASPA_DIR="${RASPA_DIR:-$HOME/micromamba/envs/gcmc-mil53}"
export RASPA_DIR   # simulate looks up $RASPA_DIR/share/raspa/{forcefield,molecules,structures}

dest="$ROOT/examples/$1"
mkdir -p "$dest/run" "$dest/reference"
cp "$RASPA_SRC/examples/$ex/simulation.input" "$dest/run/"
cp "$RASPA_SRC/examples/$ex/Output/System_0/"*.data "$dest/reference/"

cd "$dest/run"
echo "[$(date '+%F %T')] running $ex in $dest/run"
/usr/bin/time -p "$RASPA_DIR/bin/simulate" simulation.input > simulate.log 2> time.log
echo "[$(date '+%F %T')] done: $(grep real time.log)"
