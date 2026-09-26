#!/usr/bin/env bash
# Print the state of every GCMC run in runs/ without needing Claude Code.
#
# A pressure point is DONE when its Output/System_0/*.data contains
# "Simulation finished" -- the same test scripts/run_isotherm.sh uses to decide
# what to skip on resume. Anything else is either RUNNING (a live simulate
# process has that directory open) or STALE (output exists, no process, not
# finished -- run_isotherm.sh will redo it from scratch; see NOTES.md).
#
# Usage: scripts/status.sh [tag ...]      (default: every runs/* directory)
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# pids of live RASPA processes, and the directories they are working in
live_dirs() {
  local pid
  for pid in $(pgrep -x simulate 2>/dev/null); do
    lsof -a -p "$pid" -d cwd -Fn 2>/dev/null | sed -n 's/^n//p'
  done
}
LIVE="$(live_dirs)"

is_live() { [ -n "$LIVE" ] && printf '%s\n' "$LIVE" | grep -qxF "$1"; }

point_state() {                       # $1 = point dir -> DONE|RUNNING|STALE|EMPTY
  local d="$1" data
  data=$(ls "$d"/Output/System_0/*.data 2>/dev/null | head -1)
  if [ -n "$data" ] && grep -qs "Simulation finished" "$data"; then echo DONE; return; fi
  if is_live "$d"; then echo RUNNING; return; fi
  if [ -n "$data" ]; then echo STALE; else echo EMPTY; fi
}

printf '=== %s  --  %s\n' "GCMC MIL-53 job status" "$(date '+%F %T %Z')"

# ---------------------------------------------------------------- processes
echo
echo "-- processes"
n_sim=$(pgrep -x simulate 2>/dev/null | wc -l | tr -d ' ')
echo "   RASPA (simulate):  $n_sim running"
if [ "$n_sim" -gt 0 ]; then
  ps -o pid=,etime=,%cpu=,command= -p "$(pgrep -x simulate | tr '\n' ',' | sed 's/,$//')" \
    | sed 's/^/     /'
fi
n_caf=$(pgrep -x caffeinate 2>/dev/null | wc -l | tr -d ' ')
echo "   caffeinate:        $n_caf running (sleep is $( [ "$n_caf" -gt 0 ] && echo inhibited || echo NOT inhibited))"
if [ -d "$ROOT/logs" ]; then
  for L in "$ROOT"/logs/*.log; do
    [ -e "$L" ] || continue
    printf '   log %-28s %s  (last line: %s)\n' "$(basename "$L")" \
      "$(stat -f '%Sm' -t '%F %T' "$L")" "$(tail -1 "$L" | cut -c1-70)"
  done
fi

# ---------------------------------------------------------------- runs
TAGS=("$@")
if [ ${#TAGS[@]} -eq 0 ]; then
  TAGS=()
  for d in "$ROOT"/runs/*/; do TAGS+=("$(basename "$d")"); done
fi

echo
echo "-- runs"
for tag in "${TAGS[@]}"; do
  dir="$ROOT/runs/$tag"
  [ -d "$dir" ] || { printf '   %-34s (no such run)\n' "$tag"; continue; }
  pts=("$dir"/p_*)
  if [ ! -e "${pts[0]}" ]; then                       # no p_* points
    if [ ! -d "$dir/Output/System_0" ]; then          # e.g. runs/crosscheck: its own layout
      printf '   %-34s %-12s not a pressure-point run, not tracked here\n' "$tag" "n/a"
    else
      printf '   %-34s %-12s single run\n' "$tag" "$(point_state "$dir")"
    fi
    continue
  fi
  done_n=0 run_n=0 stale_n=0 empty_n=0 bad=""
  for d in "${pts[@]}"; do
    case "$(point_state "$d")" in
      DONE)    done_n=$((done_n+1));;
      RUNNING) run_n=$((run_n+1));    bad="$bad ${d##*/}:running";;
      STALE)   stale_n=$((stale_n+1)); bad="$bad ${d##*/}:STALE";;
      EMPTY)   empty_n=$((empty_n+1)); bad="$bad ${d##*/}:empty";;
    esac
  done
  total=${#pts[@]}
  overall=OK
  [ "$run_n"   -gt 0 ] && overall=RUNNING
  [ "$stale_n" -gt 0 ] && overall=NEEDS-RERUN
  [ "$empty_n" -gt 0 ] && [ "$run_n" -eq 0 ] && overall=NOT-STARTED
  [ "$done_n" -eq "$total" ] && overall=COMPLETE
  printf '   %-34s %-12s %2d/%2d done' "$tag" "$overall" "$done_n" "$total"
  [ -n "$bad" ] && printf '  --%s' "$bad"
  printf '\n'
  # stale locks: a lock whose point is not finished and has no live process
  for l in "$dir"/p_*/.lock; do
    [ -d "$l" ] || continue
    p="${l%/.lock}"
    [ "$(point_state "$p")" = RUNNING ] || printf '        stale lock: %s (remove it before re-running)\n' "${p#$ROOT/}"
  done
  if [ -f "$dir/driver.log" ]; then
    printf '        driver.log: %s\n' "$(tail -1 "$dir/driver.log" | cut -c1-96)"
  fi
done

# ---------------------------------------------------------------- artifacts
echo
echo "-- derived artifacts (results/)"
for f in isotherm_MIL53_lp_CO2_304K.csv \
         isotherm_MIL53_np_CO2_304K.csv \
         isotherm_MIL53_lp_pacman_CO2_304K.csv \
         charge_model_comparison.csv \
         langmuir_fits.csv \
         osmotic_sensitivity.csv \
         osmotic_construction.png; do
  if [ -f "$ROOT/results/$f" ]; then
    printf '   %-42s %s\n' "$f" "$(stat -f '%Sm' -t '%F %T' "$ROOT/results/$f")"
  else
    printf '   %-42s MISSING\n' "$f"
  fi
done
echo
