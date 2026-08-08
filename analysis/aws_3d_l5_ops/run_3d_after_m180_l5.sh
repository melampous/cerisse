#!/usr/bin/env bash

# Run the validated L5/SFC 3D hierarchy after the current m180 job exits.
# This script never signals, edits, or restarts the 2D calculation.

set -uo pipefail

D3=/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_3d_gpu
D2=/shared/cerisse_reflux/exm/underexpanded_jet/2d/m180_L3seed
INPUT=inputs_pilot_l5_prod
GRID=profile_20260715/fixed_grids_v5_L5.dat
SEED=v5_staged/chk04500
EX=./main3d.gnu.MPI.CUDA.ex
BIND=./bind8_v5.sh
MPI=/opt/amazon/openmpi/bin/mpirun
RUNROOT="$D3/l5_prod_cb_20260715"
CHAIN_LOG="$RUNROOT/chain.log"
DEADLINE_UTC="2026-07-16 10:50:00"

mkdir -p "$RUNROOT"
exec 9>"$D3/run_3d_after_m180_l5.lock"
if ! flock -n 9; then
  printf '[%s] another L5 chain already holds the lock\n' "$(date -u +%FT%TZ)" >&2
  exit 2
fi

log() {
  printf '[%s] %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$CHAIN_LOG"
}

fail() {
  log "ERROR: $*"
  exit 1
}

checkpoint_time() {
  sed -n '3p' "$1/Header"
}

latest_checkpoint() {
  find "$RUNROOT" -maxdepth 1 -type d -name 'chk[0-9][0-9][0-9][0-9][0-9]' \
    -print | sort -V | tail -1
}

audit_inputs() {
  cd "$D3" || return 1
  test -x "$EX" || fail "missing executable $D3/$EX"
  test -x "$BIND" || fail "missing GPU binding script $D3/$BIND"
  test -x "$MPI" || fail "missing OpenMPI 4 launcher $MPI"
  test -f "$INPUT" || fail "missing input $D3/$INPUT"
  test -f "$GRID" || fail "missing fixed grid $D3/$GRID"
  test -f "$SEED/Header" || fail "missing seed checkpoint $D3/$SEED"

  grep -Eq '^amr\.max_level[[:space:]]*=[[:space:]]*5$' "$INPUT" ||
    fail "input is not capped at L5"
  grep -Eq '^amr\.regrid_int[[:space:]]*=[[:space:]]*100000$' "$INPUT" ||
    fail "input does not freeze regridding"
  grep -Eq '^DistributionMapping\.strategy[[:space:]]*=[[:space:]]*SFC$' "$INPUT" ||
    fail "input does not use SFC distribution"
  grep -Fq "amr.regrid_file = $GRID" "$INPUT" ||
    fail "input does not select the audited L5 grid"
  grep -Eq '^cns\.stats_start_time[[:space:]]*=[[:space:]]*6\.5e-3$' "$INPUT" ||
    fail "statistics start time is not 6.5 ms"

  local grid_levels
  grid_levels=$(sed -n '1p' "$GRID")
  test "$grid_levels" = 5 || fail "L5 grid file declares $grid_levels levels"

  log "audit OK: seed_t=$(checkpoint_time "$SEED") L5/SFC fixed hierarchy selected"
  sha256sum "$INPUT" "$GRID" "$EX" "$BIND" > "$RUNROOT/sha256_manifest.txt"
}

m180_is_running() {
  pgrep -f '[m]ain2d\.gnu\.MPI\.CUDA\.ex.*inputs_pB_h100' >/dev/null
}

wait_for_m180() {
  local deadline last_report=0 now
  deadline=$(date -ud "$DEADLINE_UTC" +%s)
  while m180_is_running; do
    now=$(date -u +%s)
    if (( now >= deadline )); then
      fail "m180 still running at the 10:50 UTC safety deadline"
    fi
    if (( now - last_report >= 300 )); then
      local progress
      progress=$(grep -E '^STEP = ' "$D2/m180_phaseB_h100.log" 2>/dev/null | tail -1)
      log "waiting for m180 to finish; ${progress:-progress unavailable}"
      last_report=$now
    fi
    sleep 60
  done

  local final_time
  final_time=$(awk '/^STEP = /{t=$6} END{print t}' "$D2/m180_phaseB_h100.log" 2>/dev/null)
  test -n "$final_time" || fail "m180 process exited but no final time was found"
  python3 - "$final_time" <<'PY' || fail "m180 exited before the requested 10.0 ms stop time"
import sys
sys.exit(0 if float(sys.argv[1]) >= 0.009999 else 1)
PY
  grep -q 'AMReX .* finalized' "$D2/m180_phaseB_h100.log" ||
    fail "m180 exited without a clean AMReX finalization marker"

  local occupied
  occupied=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits |
    sed '/^[[:space:]]*$/d' | wc -l)
  test "$occupied" -eq 0 || fail "$occupied GPU compute contexts remain after m180"
  log "m180 completed cleanly at t=$final_time; all eight H100 GPUs are free"
}

GUARD_PID=
start_deadline_guard() {
  local deadline now delay
  deadline=$(date -ud "$DEADLINE_UTC" +%s)
  now=$(date -u +%s)
  (( now < deadline )) || fail "3D start requested after the safety deadline"
  delay=$((deadline - now))
  (
    sleep "$delay"
    touch "$D3/dump_and_stop"
    printf '[%s] deadline reached; requested a graceful 3D checkpoint and stop\n' \
      "$(date -u +%FT%TZ)" >> "$CHAIN_LOG"
  ) &
  GUARD_PID=$!
  log "3D-only deadline guard armed for ${DEADLINE_UTC} UTC"
}

stop_deadline_guard() {
  if test -n "${GUARD_PID:-}" && kill -0 "$GUARD_PID" 2>/dev/null; then
    kill "$GUARD_PID" 2>/dev/null || true
    wait "$GUARD_PID" 2>/dev/null || true
  fi
}
trap stop_deadline_guard EXIT

run_leg() {
  local name=$1 restart=$2 cfl=$3 stop_time=$4 max_step=$5 stats=$6 perturb=$7
  local log_file="$RUNROOT/${name}.log"

  rm -f "$D3/dump_and_stop"
  log "starting $name: restart=$restart cfl=$cfl stop=$stop_time max_step=$max_step stats=$stats"

  "$MPI" -np 8 --bind-to none "$BIND" "$EX" "$INPUT" \
    amr.restart="$restart" \
    amr.max_level=5 \
    amr.regrid_file="$GRID" \
    amr.regrid_on_restart=1 \
    amr.regrid_int=100000 \
    DistributionMapping.strategy=SFC \
    cfl="$cfl" stop_time="$stop_time" max_step="$max_step" \
    cns.do_reflux=1 cns.check_state=1 \
    cns.record_stats="$stats" cns.stats_start_time=6.5e-3 \
    prob.perturb_enabled="$perturb" \
    prob.perturb_t_start=5.0e-3 prob.perturb_duration=2.0e-4 \
    amr.checkpoint_files_output=1 amr.check_file="$RUNROOT/chk" \
    amr.plot_files_output=1 amr.plot_file="$RUNROOT/plt" \
    amr.check_int=500 amr.plot_int=500 \
    amr.data_log="$RUNROOT/probes_${name}.log" > "$log_file" 2>&1
  local rc=$?

  if grep -q 'Stopped by user' "$log_file"; then
    log "$name stopped cleanly by the 10:50 UTC deadline"
    return 9
  fi
  test "$rc" -eq 0 || fail "$name exited with MPI status $rc"
  grep -q 'AMReX .* finalized' "$log_file" || fail "$name lacks a clean finalization marker"
  if grep -qE 'STATE CHECK FAILED|amrex::Abort|SIGABRT|Segfault|NaN encountered' "$log_file"; then
    fail "$name contains a state or runtime failure"
  fi
  log "$name completed; checkpoint=$(latest_checkpoint)"
  return 0
}

check_ramp_gate() {
  local log_file="$RUNROOT/ramp.log"
  grep -Eq 'Level 5[[:space:]]+26 grids[[:space:]]+12582912 cells' "$log_file" ||
    fail "ramp did not construct the audited L5 hierarchy"
  if grep -qE 'Level 6[[:space:]]+[0-9]+ grids' "$log_file"; then
    fail "ramp unexpectedly constructed L6"
  fi

  local median
  median=$(python3 - "$log_file" <<'PY'
import re, statistics, sys
values = [float(x) for x in re.findall(r'Coarse TimeStep time: ([0-9.]+)', open(sys.argv[1]).read())]
if len(values) < 10:
    raise SystemExit(1)
print(f'{statistics.median(values[-10:]):.6f}')
PY
  ) || fail "fewer than ten ramp timings were recorded"
  python3 - "$median" <<'PY' || fail "performance gate failed: median ${median}s > 6.5s"
import sys
sys.exit(0 if float(sys.argv[1]) <= 6.5 else 1)
PY
  log "performance gate PASS: last-10 median=${median}s <= 6.5s"
}

main() {
  audit_inputs
  if test "${1:-}" = --audit; then
    log "audit-only request complete"
    return 0
  fi

  log "waiting for the existing m180 process; no action will be taken on the 2D case"
  wait_for_m180
  test -z "$(pgrep -f '[m]ain3d\.gnu\.MPI\.CUDA\.ex' || true)" ||
    fail "another 3D process is already running"
  start_deadline_guard

  run_leg ramp "$SEED" 0.15 6.5e-3 4520 0 1 || {
    test "$?" -eq 9 && return 0
    return 1
  }
  check_ramp_gate

  local checkpoint
  checkpoint=$(latest_checkpoint)
  test -n "$checkpoint" || fail "ramp produced no checkpoint"
  run_leg develop "$checkpoint" 0.30 6.5e-3 9999999 0 1 || {
    test "$?" -eq 9 && return 0
    return 1
  }

  checkpoint=$(latest_checkpoint)
  test -n "$checkpoint" || fail "development leg produced no checkpoint"
  python3 - "$(checkpoint_time "$checkpoint")" <<'PY' || fail "development did not end at 6.5 ms"
import sys
sys.exit(0 if abs(float(sys.argv[1]) - 0.0065) <= 1.0e-10 else 1)
PY

  run_leg statistics "$checkpoint" 0.30 10.16e-3 9999999 1 0 || {
    test "$?" -eq 9 && return 0
    return 1
  }
  log "3D L5 production chain completed at checkpoint $(latest_checkpoint)"
}

main "$@"
