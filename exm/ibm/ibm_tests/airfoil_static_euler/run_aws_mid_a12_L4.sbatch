#!/bin/bash
#SBATCH --job-name=mid_a12_L4
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --threads-per-core=1
#SBATCH --time=02:00:00
#SBATCH --output=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler/slurm_L4_%j.out
#SBATCH --error=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler/slurm_L4_%j.err

set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load openmpi/4.1.7 2>/dev/null || true

ROOT=/shared/cerisse
DIR=${ROOT}/exm/ibm/ibm_tests/airfoil_static_euler
cd "${DIR}"

MAX_STEP=${MAX_STEP:-5000}
PLOT_INT=${PLOT_INT:-250}

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
SRUN_OPTS="--mpi=pmix --ntasks=96 --cpus-per-task=1 --hint=nomultithread --threads-per-core=1"

echo "================================================================"
echo " mid_a12  L4 convergence  Job ${SLURM_JOB_ID}  $(date)"
echo "================================================================"

EXE=$(ls main2d.*.ex 2>/dev/null | head -1)
if [ -z "$EXE" ] || [ ! -x "$EXE" ]; then
  echo "building..."
  make -j 96 2>&1 | tail -10
  EXE=$(ls main2d.*.ex 2>/dev/null | head -1)
fi
echo "binary: ${EXE}"

rm -rf plot_mid_a12_L4 surf_mid_a12_L4
srun ${SRUN_OPTS} ./"${EXE}" inputs_mid_a12_L4 \
     max_step=${MAX_STEP} \
     amr.plot_int=${PLOT_INT} \
     amr.check_int=-1 \
     2>&1 | tee run_mid_a12_L4.log

echo "=== done $(date) ==="
