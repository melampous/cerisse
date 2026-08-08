#!/usr/bin/env bash

# Request a normal checkpoint-and-stop as soon as m180 reaches t = 10 ms.

set -euo pipefail

D2=/shared/cerisse_reflux/exm/underexpanded_jet/2d/m180_L3seed
RUN_LOG="$D2/m180_phaseB_h100.log"
STOP_FILE="$D2/dump_and_stop"
WATCH_LOG="$D2/stop_at_10ms.log"
TARGET=0.010

exec 9>"$D2/stop_at_10ms.lock"
if ! flock -n 9; then
  printf '[%s] another 10 ms stop monitor is already active\n' \
    "$(date -u +%FT%TZ)" >&2
  exit 2
fi

log() {
  printf '[%s] %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$WATCH_LOG"
}

m180_is_running() {
  pgrep -f '[m]ain2d\.gnu\.MPI\.CUDA\.ex.*inputs_pB_h100' >/dev/null
}

latest_state() {
  awk '/^STEP = /{step=$3; time=$6} END{if (time != "") print step, time}' \
    "$RUN_LOG"
}

target_reached() {
  awk -v time="$1" -v target="$TARGET" 'BEGIN{exit !(time >= target)}'
}

log "monitoring m180; target physical time is ${TARGET} s"

while m180_is_running; do
  read -r step time < <(latest_state)
  if test -n "${time:-}" && target_reached "$time"; then
    touch "$STOP_FILE"
    log "target reached at step=$step time=$time; requested checkpoint and normal exit"

    for _ in $(seq 1 900); do
      m180_is_running || break
      sleep 1
    done
    m180_is_running && {
      log "ERROR: m180 did not exit within 15 minutes of the stop request"
      exit 1
    }

    grep -q 'AMReX .* finalized' "$RUN_LOG" || {
      log "ERROR: m180 exited without a clean AMReX finalization marker"
      exit 1
    }

    read -r final_step final_time < <(latest_state)
    target_reached "$final_time" || {
      log "ERROR: final time $final_time is below the requested target"
      exit 1
    }
    log "m180 completed cleanly at step=$final_step time=$final_time"
    exit 0
  fi
  sleep 1
done

read -r final_step final_time < <(latest_state)
log "ERROR: m180 exited before the monitor requested a stop (step=${final_step:-unknown} time=${final_time:-unknown})"
exit 1
