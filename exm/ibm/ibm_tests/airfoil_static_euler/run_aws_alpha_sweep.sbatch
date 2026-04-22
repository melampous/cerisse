#!/bin/bash
#SBATCH --job-name=aer_euler
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler/slurm_%j.out
#SBATCH --error=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler/slurm_%j.err
#
# Canonical static diamond wedge (Euler) - full alpha sweep
#   2 geometries (mid, x018) x 6 AoA (0, 4, 8, 12, 16, 18) = 12 cases
#   512 x 544 base, max_level=3, CFL 0.25, 5000 steps each
#
# On an hpc8a node (96 ranks) each case takes ~2 min; the whole matrix
# finishes in well under 1 h.

set -euo pipefail

source /etc/profile.d/modules.sh 2>/dev/null || true
module load openmpi/4.1.7 2>/dev/null || true

CASE_DIR=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler
cd "${CASE_DIR}"

EXE=$(ls main2d.*.ex 2>/dev/null | head -1)
if [ -z "$EXE" ] || [ ! -x "$EXE" ]; then
  echo "ERROR: binary not built; run 'make -j' in ${CASE_DIR} first"
  exit 1
fi
echo "Binary: ${EXE}"

MAX_STEP=${MAX_STEP:-5000}
PLOT_INT=${PLOT_INT:-250}

for geom in mid x018; do
  for a in 00 04 08 12 16 18; do
    tag=${geom}_a${a}
    echo "================================================================"
    echo " case ${tag}   start $(date)"
    echo "================================================================"
    rm -rf plot_${tag} surf_${tag}
    srun --mpi=pmix ./"${EXE}" inputs_${tag} \
         max_step=${MAX_STEP} \
         amr.plot_int=${PLOT_INT} \
         amr.check_int=-1 \
         2>&1 | tee run_${tag}.log
    echo " case ${tag}   done  $(date)"
  done
done

echo "[$(date)] full alpha sweep finished"
