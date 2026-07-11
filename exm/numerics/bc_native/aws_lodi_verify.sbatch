#!/bin/bash
#SBATCH -J lodi_verify
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH -o lodi_verify_%j.out
#SBATCH -e lodi_verify_%j.err

# Independent verification of the multi-D LODI hook:
#   (1) PERFORMANCE: same binary, same IC/grid, vary only the BC (lodi-active
#       vs lodi-fallback=char vs foextrap) -> isolate the hook's per-step cost.
#   (2) PHYSICS: reproduce the reflection metric matrix (normal / oblique blob /
#       oblique_train) x (char / foextrap / lodi) at N=96,192.

set -uo pipefail
CASE=/shared/cerisse/exm/numerics/bc_native
cd "$CASE"
module purge; module load openmpi/4.1.7

echo "=== BUILD (WenoZ5, with bcnormal_lodi hook) $(date) ==="
rm -rf tmp_build_dir main2d.gnu.MPI.ex
make -j8 2>&1 | tail -3
ls -la main2d.gnu.MPI.ex || { echo "BUILD FAILED"; exit 1; }
EXE=$CASE/main2d.gnu.MPI.ex

# ---------- (1) PERFORMANCE ----------
# Identical grid + IC (mode=2 oblique carrier) for all three; only the BC differs.
#   lodi : hi_bc=7 + bc_mode=2  (hook ACTIVE, transverse derivatives)
#   char : hi_bc=7 + bc_mode=1  (hook returns false -> 1D characteristic fallback)
#   foe  : hi_bc=2              (foextrap, AMReX direct, bcnormal not called)
COMMON="max_step=3000 stop_time=1e9 amr.n_cell=288 192 geometry.prob_lo=0.0 0.0 \
geometry.prob_hi=1.5 1.0 geometry.is_periodic=0 1 amr.max_grid_size=512 \
amr.blocking_factor=2 amr.plot_files_output=0 amr.checkpoint_files_output=0 \
cns.nstep_screen_output=100000 prob.mode=2 prob.mach=0.30 prob.amp=1.0e-4 \
prob.theta_deg=30.0 cns.lo_bc=1 -1"

echo ""; echo "=== (1) PERFORMANCE  (3000 steps, 288x192, 4 ranks, 2 reps each) ==="
run_time () {  # $1=label  $2=hi_bc  $3=bc_mode
  for rep in 1 2; do
    t0=$(date +%s.%N)
    mpirun --oversubscribe -np 4 "$EXE" inputs $COMMON \
        cns.hi_bc="$2" prob.bc_mode=$3 > /dev/null 2>&1
    t1=$(date +%s.%N)
    echo "  $1  rep$rep  wall = $(echo "$t0 $t1" | awk '{printf "%.2f", $2-$1}') s"
  done
}
run_time "lodi (hook ACTIVE)   " "7 -1" 2
run_time "char (hook fallback) " "7 -1" 1
run_time "foe  (foextrap)      " "2 -1" 1

# ---------- (2) PHYSICS reproduction ----------
echo ""; echo "=== (2) PHYSICS metric matrix via runner $(date) ==="
cd /shared/cerisse
python3 tools/run_native_bc_validation.py \
    --outdir /shared/cerisse/temp/bc_lodi_verify \
    --cases normal,oblique,oblique_train \
    --bcs char,foextrap,lodi \
    --res 96,192 2>&1 | tail -40 || echo "(runner failed -- will analyze plotfiles locally)"

echo ""; echo "=== metrics CSV ==="
cat /shared/cerisse/temp/bc_lodi_verify/native_bc_metrics.csv 2>/dev/null || echo "no csv"
echo "=== done $(date) ==="
