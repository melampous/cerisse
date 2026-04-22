#!/bin/bash
#SBATCH --job-name=aer_all
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler/slurm_all_%j.out
#SBATCH --error=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler/slurm_all_%j.err
#
# Full static-airfoil matrix on one 96-rank hpc8a node.
#
#   Block 1: Euler canonical         (build once, 12 runs share binary)
#             2 geometries x 6 AoA = 12 cases
#   Block 2: NS robustness           (rebuild per Re, 3 Re values)
#             mid_a00 only, Re = 1e4 / 1e5 / 1e6
#
# Expected wall-clock on 96 ranks: ~30-45 min total.

set -euo pipefail

source /etc/profile.d/modules.sh 2>/dev/null || true
module load openmpi/4.1.7 2>/dev/null || true

ROOT=/shared/cerisse
MAX_STEP=${MAX_STEP:-5000}
PLOT_INT=${PLOT_INT:-250}

echo "================================================================"
echo " Static-airfoil full matrix"
echo " Job ${SLURM_JOB_ID} on node(s) ${SLURM_JOB_NODELIST}"
echo " Start $(date)"
echo " MAX_STEP=${MAX_STEP}  PLOT_INT=${PLOT_INT}"
echo "================================================================"

########################
# Block 1: Euler sweep #
########################
EULER_DIR=${ROOT}/exm/ibm/ibm_tests/airfoil_static_euler
cd "${EULER_DIR}"
echo "=== Block 1: Euler canonical ==="
echo "[$(date)] building ${EULER_DIR}"
make -j 96 2>&1 | tail -20

EXE=$(ls main2d.*.ex 2>/dev/null | head -1)
if [ -z "$EXE" ] || [ ! -x "$EXE" ]; then
  echo "ERROR: Euler binary not built"; exit 1
fi
echo "Euler binary: ${EXE}"

for geom in mid x018; do
  for a in 00 04 08 12 16 18; do
    tag=${geom}_a${a}
    echo "---- Euler case ${tag}  start $(date)"
    rm -rf plot_${tag} surf_${tag}
    srun --mpi=pmix ./"${EXE}" inputs_${tag} \
         max_step=${MAX_STEP} \
         amr.plot_int=${PLOT_INT} \
         amr.check_int=-1 \
         2>&1 | tee run_${tag}.log
    echo "---- Euler case ${tag}  done  $(date)"
  done
done

#######################
# Block 2: NS Re sweep#
#######################
NS_DIR=${ROOT}/exm/ibm/ibm_tests/airfoil_static
cd "${NS_DIR}"
echo "=== Block 2: NS Reynolds sweep ==="

for re in 1e4 1e5 1e6; do
  echo "---- NS Re=${re}  start $(date)"
  echo "     rebuilding with REYNOLDS_USER=${re}"
  make realclean > /dev/null 2>&1
  make -j 96 REYNOLDS_USER=${re} 2>&1 | tail -8

  EXE_NS=$(ls main2d.*.ex 2>/dev/null | head -1)
  if [ -z "$EXE_NS" ] || [ ! -x "$EXE_NS" ]; then
    echo "ERROR: NS binary not built for Re=${re}"; exit 1
  fi
  # stage the same pre-rotated .dat used by the Euler block so both cases
  # reference identical geometry without duplicating files
  ln -sf ${EULER_DIR}/diamond_wedge_mid_a00.dat ./diamond_wedge.dat

  rm -rf plot_re${re} surf_re${re}
  srun --mpi=pmix ./"${EXE_NS}" inputs_re${re} \
       max_step=${MAX_STEP} \
       amr.plot_int=${PLOT_INT} \
       amr.check_int=-1 \
       2>&1 | tee run_re${re}.log
  echo "---- NS Re=${re}  done  $(date)"
done

echo "================================================================"
echo " All cases finished at $(date)"
echo "================================================================"
