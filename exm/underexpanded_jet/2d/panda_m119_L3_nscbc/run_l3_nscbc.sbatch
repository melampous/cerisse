#!/bin/bash
#SBATCH -J m119L3nscbc
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=96
#SBATCH --time=12:00:00
#SBATCH -o l3nscbc_%j.out
#SBATCH -e l3nscbc_%j.err
# m119 L3 NSCBC/LODI pilot + legacy baseline twin (A/B lip receptivity).
set -uo pipefail
CASE=/shared/cerisse_reflux/exm/underexpanded_jet/2d/panda_m119_L3_nscbc
cd "$CASE"
module purge; module load openmpi/4.1.7

echo "=== BUILD (CPU, reflux tree) $(date -u) ==="
rm -rf tmp_build_dir main2d.gnu.MPI.ex
make -j 48 2>&1 | tail -3
ls -la main2d.gnu.MPI.ex || { echo BUILD_FAILED; exit 1; }
EXE=$CASE/main2d.gnu.MPI.ex

echo "=== SMOKE 40 steps, both modes $(date -u) ==="
rm -rf smoke_on smoke_off run_lodi run_base
mkdir -p smoke_on smoke_off run_lodi/plot run_base/plot
SMK="inputs max_step=40 stop_time=1e9 amr.plot_int=40 amr.check_int=100000 amr.data_log=/dev/null"
mpirun -np 8 $EXE $SMK prob.use_lodi=1 amr.plot_file=smoke_on/plt  > smoke_on.log 2>&1 || { echo SMOKE_ON_FAILED; exit 1; }
mpirun -np 8 $EXE $SMK prob.use_lodi=0 cns.hi_bc="2 2" amr.plot_file=smoke_off/plt > smoke_off.log 2>&1 || { echo SMOKE_OFF_FAILED; exit 1; }
HA=$(find smoke_on  -name "Cell_D_*" -exec md5sum {} + | sort -k2 | md5sum | cut -d" " -f1)
HB=$(find smoke_off -name "Cell_D_*" -exec md5sum {} + | sort -k2 | md5sum | cut -d" " -f1)
[ "$HA" != "$HB" ] && echo "HOOK_ACTIVE (plotfiles differ)" || { echo "HOOK_NOT_FIRING"; exit 1; }

echo "=== RUN A: NSCBC/LODI to 10 ms $(date -u) ==="
mpirun -np 96 $EXE inputs prob.use_lodi=1 \
  amr.plot_file=run_lodi/plot/plt amr.check_file=run_lodi/plot/chk \
  amr.data_log=run_lodi/probes.log > run_lodi/run.log 2>&1
RCA=$?; echo "RUN_A rc=$RCA $(date -u)"; tail -2 run_lodi/run.log

echo "=== RUN B: legacy baseline twin to 10 ms $(date -u) ==="
mpirun -np 96 $EXE inputs prob.use_lodi=0 cns.hi_bc="2 2" \
  amr.plot_file=run_base/plot/plt amr.check_file=run_base/plot/chk \
  amr.data_log=run_base/probes.log > run_base/run.log 2>&1
RCB=$?; echo "RUN_B rc=$RCB $(date -u)"; tail -2 run_base/run.log

echo "=== DONE A=$RCA B=$RCB $(date -u) ==="
