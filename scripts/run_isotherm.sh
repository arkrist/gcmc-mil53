#!/usr/bin/env bash
# Run the pressure points of runs/<tag>/, optionally several at a time.
#
# Resume, per point: a point is "done" when its Output/System_0/*.data contains
# "Simulation finished". Done points are skipped, so re-running the script after
# a crash, a kill or a reboot only costs the points that were in flight. An
# unfinished point is rerun FROM SCRATCH (its partial Output/ is moved aside to
# Output.incomplete.<timestamp>/), so every result comes from one continuous run.
# A .lock directory per point (mkdir, atomic) keeps parallel workers, or a second
# invocation of this script, from starting the same point twice; a lock left by a
# killed run is stale and is removed if its point is not finished.
#
# Usage: scripts/run_isotherm.sh [tag] [jobs]      (defaults: production, 1)
#   jobs > 1 runs that many points concurrently (the M2 has 4 performance cores).
# Env:   RASPA_DIR  conda env prefix (default ~/micromamba/envs/gcmc-mil53)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TAG="${1:-production}"
JOBS="${2:-1}"
RASPA_DIR="${RASPA_DIR:-$HOME/micromamba/envs/gcmc-mil53}"
export RASPA_DIR ROOT TAG
LOG="$ROOT/runs/$TAG/driver.log"

[ -d "$ROOT/runs/$TAG" ] || { echo "no runs/$TAG -- run scripts/make_inputs.py first" >&2; exit 1; }

run_point() {                      # $1 = point directory
  local d="$1" name start rc
  name=$(basename "$d")
  if grep -qs "Simulation finished" "$d"/Output/System_0/*.data; then
    echo "[$(date '+%F %T')] $name: already finished, skipping" | tee -a "$LOG"; return 0
  fi
  if ! mkdir "$d/.lock" 2>/dev/null; then
    echo "[$(date '+%F %T')] $name: locked by another worker, skipping" | tee -a "$LOG"; return 0
  fi
  trap 'rmdir "$d/.lock" 2>/dev/null || true' RETURN
  if [ -d "$d/Output" ]; then
    mv "$d/Output" "$d/Output.incomplete.$(date +%Y%m%d%H%M%S)"
    echo "[$(date '+%F %T')] $name: found incomplete output, rerunning from scratch" | tee -a "$LOG"
  fi
  echo "[$(date '+%F %T')] $name: start" | tee -a "$LOG"
  start=$(date +%s)
  rc=0
  (cd "$d" && "$RASPA_DIR/bin/simulate" simulation.input > simulate.log 2>&1) || rc=$?
  if [ $rc -ne 0 ]; then
    echo "[$(date '+%F %T')] $name: FAILED (exit $rc) after $(( $(date +%s) - start )) s" | tee -a "$LOG"; return 0
  fi
  echo "[$(date '+%F %T')] $name: done in $(( $(date +%s) - start )) s" | tee -a "$LOG"
  if grep -qs "NO VDW INTERACTION" "$d"/Output/System_0/*.data; then
    echo "[$(date '+%F %T')] $name: WARNING -- RASPA reports atom pairs with no VDW interaction (see NOTES.md)" | tee -a "$LOG"
  fi
}
export -f run_point
export LOG

batch_start=$(date +%s)
echo "[$(date '+%F %T')] start: tag=$TAG jobs=$JOBS" | tee -a "$LOG"
# lowest pressure first; with several jobs the long high-pressure points overlap
ls -d "$ROOT/runs/$TAG"/p_* | awk -F'p_' '{print $NF" "$0}' | sort -g | cut -d' ' -f2- \
  | xargs -P "$JOBS" -I{} bash -c 'run_point "$@"' _ {}
echo "[$(date '+%F %T')] all points of runs/$TAG finished; batch wall time $(( $(date +%s) - batch_start )) s" | tee -a "$LOG"
