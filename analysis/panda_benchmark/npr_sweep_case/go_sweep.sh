#!/bin/bash
# =============================================================================
# NPR sweep chain driver — 2D RZ Panda jet, staircase NPR 2.0 -> 5.746, D/256.
#
#   ./go_sweep.sh dryrun   print the leg plan + commands, execute nothing
#   ./go_sweep.sh smoke    0.5 ms fresh-start pace test (smoke/ dir), report s/step
#   ./go_sweep.sh pace     parse smoke log, report pace + projected chain wall time
#   ./go_sweep.sh chain    run the 16-leg sweep sequentially (resumable)
#
# Legs alternate flush (transient discarded) / stats (per-plateau window).
# Each leg restarts from the newest checkpoint of the previous leg; the reflux
# tree zeroes online statistics on restart, so every stats leg's final
# plt/chk carries exactly that plateau's [flush_end, plateau_end] window.
# NPR(t) lives in prob.h (npr_of_t) as a function of ABSOLUTE time, so legs
# need no per-leg parameters beyond stop_time.
# NO JOB IS LAUNCHED unless this script is invoked with smoke|chain.
# =============================================================================
set -u
cd "$(dirname "$0")"

EXE=./main2d_sweep.ex
INPUTS=inputs_sweep
NP=${NP:-8}
BIND=./bind8_a100.sh
STATUS=NPRSWEEP.status
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:/shared/cuda-12.6/lib64:${LD_LIBRARY_PATH:-}

# leg table: name  stop_time(s)   (plateau NPR in the name is nominal)
LEGS=(
  "00_p1.893_flush  2.0e-3"
  "01_p1.893_stats  3.3e-3"
  "02_p2.394_flush  4.8e-3"
  "03_p2.394_stats  6.1e-3"
  "04_p2.800_flush  7.6e-3"
  "05_p2.800_stats  8.9e-3"
  "06_p3.273_flush  10.4e-3"
  "07_p3.273_stats  11.7e-3"
  "08_p3.671_flush  13.2e-3"
  "09_p3.671_stats  14.5e-3"
  "10_p4.200_flush  16.0e-3"
  "11_p4.200_stats  17.3e-3"
  "12_p4.900_flush  18.8e-3"
  "13_p4.900_stats  20.1e-3"
  "14_p5.746_flush  21.6e-3"
  "15_p5.746_stats  22.9e-3"
)

log_status () { echo "[$(date -u +%FT%TZ)] $1" | tee -a "$STATUS"; }

newest_chk () {  # newest_chk <dir>
  ls -d "$1"/chk[0-9]* 2>/dev/null | grep -v temp | sort -V | tail -1
}

run_leg () {  # run_leg <legdir> <stop_time> <restart_chk|-> <logfile>
  local dir=$1 stop=$2 rst=$3 lg=$4
  mkdir -p "$dir"
  local extra=""
  [ "$rst" != "-" ] && extra="amr.restart=$rst amr.regrid_on_restart=1"
  mpirun -np "$NP" --bind-to none -x LD_LIBRARY_PATH "$BIND" "$EXE" "$INPUTS" \
      stop_time="$stop" \
      amr.plot_file="$dir/plt" amr.check_file="$dir/chk" \
      amr.data_log="$dir/probes.log" \
      $extra > "$lg" 2>&1
  return $?
}

mode=${1:-dryrun}

case "$mode" in
# -----------------------------------------------------------------------------
dryrun)
  echo "exe=$EXE  np=$NP  inputs=$INPUTS"
  echo "leg plan (name / stop_time / restart source):"
  prev="-"
  for entry in "${LEGS[@]}"; do
    set -- $entry; name=$1; stop=$2
    echo "  legs/$name   stop_time=$stop   restart=${prev}"
    prev="legs/$name/<newest chk>"
  done
  echo "smoke: fresh start, stop_time=0.5e-3 -> smoke/"
  echo "(dryrun: nothing executed)"
  ;;
# -----------------------------------------------------------------------------
smoke)
  log_status "NPRSWEEP_SMOKE_START np=$NP"
  mkdir -p smoke
  run_leg smoke 0.5e-3 - smoke/run.log
  rc=$?
  nplt=$(ls -d smoke/plt[0-9]* 2>/dev/null | wc -l)
  pace=$(grep -a "Coarse TimeStep time" smoke/run.log | tail -20 \
        | awk '{s+=$NF; n++} END{if(n>0) printf "%.3f", s/n}')
  log_status "NPRSWEEP_SMOKE_END rc=$rc pace=${pace:-NA}s/step nplt=$nplt"
  ;;
# -----------------------------------------------------------------------------
pace)
  pace=$(grep -a "Coarse TimeStep time" smoke/run.log 2>/dev/null | tail -20 \
        | awk '{s+=$NF; n++} END{if(n>0) printf "%.3f", s/n}')
  dt=$(grep -a "^STEP = " smoke/run.log 2>/dev/null | tail -1 | awk '{print $9}')
  if [ -n "$pace" ] && [ -n "$dt" ]; then
    awk -v p="$pace" -v dt="$dt" 'BEGIN{
      steps = 22.9e-3/dt;
      printf "pace = %s s/step, dt_L0 = %s s -> ~%d coarse steps total, ~%.1f h wall\n",
             p, dt, steps, steps*p/3600.0 }'
  else
    echo "no smoke data (run ./go_sweep.sh smoke first)"
  fi
  ;;
# -----------------------------------------------------------------------------
chain)
  log_status "NPRSWEEP_CHAIN_START np=$NP"
  prevdir=""
  for entry in "${LEGS[@]}"; do
    set -- $entry; name=$1; stop=$2
    dir="legs/$name"
    if [ -f "$dir/DONE" ]; then
      log_status "SKIP $name (DONE marker present)"
      prevdir="$dir"; continue
    fi
    rst="-"
    if [ -n "$prevdir" ]; then
      rst=$(newest_chk "$prevdir")
      if [ -z "$rst" ]; then
        log_status "ABORT $name: no checkpoint found in $prevdir"
        exit 1
      fi
    fi
    log_status "LEG_START $name stop=$stop restart=${rst}"
    run_leg "$dir" "$stop" "$rst" "$dir/run.log"
    rc=$?
    if [ $rc -ne 0 ]; then
      log_status "LEG_FAIL $name rc=$rc (chain halted; fix and re-run chain to resume)"
      exit $rc
    fi
    touch "$dir/DONE"
    tfin=$(grep -a "^STEP = " "$dir/run.log" | tail -1 | awk '{print $6}')
    log_status "LEG_DONE $name t=$tfin"
    prevdir="$dir"
  done
  log_status "NPRSWEEP_CHAIN_ALL_DONE"
  ;;
# -----------------------------------------------------------------------------
*)
  echo "usage: $0 {dryrun|smoke|pace|chain}"; exit 2 ;;
esac
