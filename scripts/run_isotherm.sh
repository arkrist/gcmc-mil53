#!/usr/bin/env bash
# Run all pressure points of runs/<tag>/ sequentially, lowest pressure first.
#
# Resume: a point is "done" when its Output/System_0/*.data contains
# "Simulation finished". Done points are skipped. An unfinished point (killed,
# crashed, laptop asleep) is rerun FROM SCRATCH: its partial Output/ is moved
# to Output.incomplete.<timestamp>/. Restarting mid-run from a RASPA binary
# restart is possible but not used, to keep every result from one continuous
# run.
#
# Usage: scripts/run_isotherm.sh [tag]            (default tag: production)
# Env:   RASPA_DIR  conda env prefix (default ~/micromamba/envs/gcmc-mil53)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TAG="${1:-production}"
RASPA_DIR="${RASPA_DIR:-$HOME/micromamba/envs/gcmc-mil53}"
export RASPA_DIR
LOG="$ROOT/runs/$TAG/driver.log"

[ -d "$ROOT/runs/$TAG" ] || { echo "no runs/$TAG -- run scripts/make_inputs.py first" >&2; exit 1; }

# p_<Pa> directories, sorted numerically by pressure
dirs=$(ls -d "$ROOT/runs/$TAG"/p_* | awk -F'p_' '{print $NF" "$0}' | sort -g | cut -d' ' -f2-)

for d in $dirs; do
  name=$(basename "$d")
  if grep -qs "Simulation finished" "$d"/Output/System_0/*.data; then
    echo "[$(date '+%F %T')] $name: already finished, skipping" | tee -a "$LOG"
    continue
  fi
  if [ -d "$d/Output" ]; then
    mv "$d/Output" "$d/Output.incomplete.$(date +%Y%m%d%H%M%S)"
    echo "[$(date '+%F %T')] $name: found incomplete output, rerunning from scratch" | tee -a "$LOG"
  fi
  echo "[$(date '+%F %T')] $name: start" | tee -a "$LOG"
  start=$(date +%s)
  (cd "$d" && "$RASPA_DIR/bin/simulate" simulation.input > simulate.log 2>&1)
  echo "[$(date '+%F %T')] $name: done in $(( $(date +%s) - start )) s" | tee -a "$LOG"
  if grep -qs "NO VDW INTERACTION" "$d"/Output/System_0/*.data; then
    echo "[$(date '+%F %T')] $name: ABORT -- RASPA reports atom pairs with no VDW interaction (see NOTES.md)" | tee -a "$LOG"
    exit 2
  fi
done
echo "[$(date '+%F %T')] all points of runs/$TAG finished" | tee -a "$LOG"
