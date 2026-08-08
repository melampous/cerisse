#!/bin/bash
# =============================================================================
# queue_after_shift100.sh — RUNS ON THE A100 BOX (10.0.135.240).
# Waits for the m142_L5_shift100 chain to finish, then builds the NPR-sweep
# executable with the box's native CUDA (/usr/local/cuda) and runs the sweep:
# smoke (0.5 ms pace test) -> 16-leg chain (NPR 1.893 -> 5.746).
# Queued per user authorization 2026-07-20 ~01:55Z ("设计好直接提交,排队").
#
# Gate logic:
#   - primary : m142_L5_shift100_TRUESTAT_DONE appears in shift100/bcvar.status
#   - fallback: no main3d process for 30 consecutive minutes (chain presumed
#               dead) -> proceed, loudly logged
#   - timeout : 12 h -> abort (nothing launched)
# All progress lines -> NPRSWEEP.status (on /shared, visible from head node).
# =============================================================================
set -u
CASE=/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_npr_sweep
GATE=/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_L5_shift100/bcvar.status
S=$CASE/NPRSWEEP.status
log(){ echo "[$(date -u +%FT%TZ)] $*" >> "$S"; }

log "QUEUE_WATCHER_START pid=$$ host=$(hostname) gate=shift100_TRUESTAT_DONE"

idle=0; waited=0
while true; do
  if grep -q "m142_L5_shift100_TRUESTAT_DONE" "$GATE" 2>/dev/null; then
    log "GATE_OPEN shift100 TRUESTAT_DONE"
    break
  fi
  if pgrep -f "main3d.*\.ex" > /dev/null 2>&1; then
    idle=0
  else
    idle=$((idle+1))
    if [ "$idle" -ge 30 ]; then
      log "GATE_OPEN_FALLBACK no main3d process for 30 min, TRUESTAT_DONE absent"
      break
    fi
  fi
  waited=$((waited+1))
  if [ "$waited" -ge 720 ]; then log "ABORT_TIMEOUT gate not open after 12 h"; exit 1; fi
  sleep 60
done

# small settle margin so the finishing leg's final checkpoint write completes
sleep 120

cd "$CASE" || { log "ABORT cd $CASE failed"; exit 1; }

log "BUILD_START native /usr/local/cuda, make -j24"
rm -rf tmp_build_dir
export PATH=/usr/local/cuda/bin:$PATH
make -j24 > build_a100.log 2>&1
if [ ! -f main2d.gnu.MPI.CUDA.ex ]; then
  log "BUILD_FAIL see build_a100.log — nothing launched"
  exit 1
fi
cp -f main2d.gnu.MPI.CUDA.ex main2d_sweep.ex
log "BUILD_OK $(ls -la main2d_sweep.ex | awk '{print $5}') bytes"

./go_sweep.sh smoke
if ! tail -5 "$S" | grep -q "NPRSWEEP_SMOKE_END rc=0"; then
  log "SMOKE_FAIL — chain NOT started (see smoke/run.log)"
  exit 1
fi
./go_sweep.sh pace >> "$S" 2>&1

./go_sweep.sh chain
log "QUEUE_WATCHER_EXIT"
