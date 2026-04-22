#!/bin/bash
#SBATCH --job-name=ab_check
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --threads-per-core=1
#SBATCH --time=01:00:00
#SBATCH --output=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler_abcheck/slurm_ab_%j.out
#SBATCH --error=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static_euler_abcheck/slurm_ab_%j.err

set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load openmpi/4.1.7 2>/dev/null || true

ROOT=/shared/cerisse
DIR=${ROOT}/exm/ibm/ibm_tests/airfoil_static_euler_abcheck
cd "${DIR}"

MAX_STEP=${MAX_STEP:-5000}
PLOT_INT=${PLOT_INT:-250}

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
SRUN_OPTS="--mpi=pmix --ntasks=96 --cpus-per-task=1 --hint=nomultithread --threads-per-core=1"

echo "================================================================"
echo " A/B equivalence check  Job ${SLURM_JOB_ID}  $(date)"
echo "================================================================"

echo "building..."
make -j 96 2>&1 | tail -10

EXE=$(ls main2d.*.ex 2>/dev/null | head -1)
if [ -z "$EXE" ] || [ ! -x "$EXE" ]; then
  echo "ERROR: build failed"; exit 1
fi
echo "binary: ${EXE}"

rm -rf plot_abcheck surf_abcheck
srun ${SRUN_OPTS} ./"${EXE}" inputs_abcheck \
     max_step=${MAX_STEP} \
     amr.plot_int=${PLOT_INT} \
     amr.check_int=-1 \
     2>&1 | tee run_abcheck.log

echo "=== done $(date) ==="
