#!/bin/bash
#SBATCH --job-name=nsmidR4
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --threads-per-core=1
#SBATCH --time=01:00:00
#SBATCH --output=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static/slurm_mid_a12_re1e4_%j.out
#SBATCH --error=/shared/cerisse/exm/ibm/ibm_tests/airfoil_static/slurm_mid_a12_re1e4_%j.err

set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load openmpi/4.1.7 2>/dev/null || true

ROOT=/shared/cerisse
DIR=${ROOT}/exm/ibm/ibm_tests/airfoil_static
EULER_DIR=${ROOT}/exm/ibm/ibm_tests/airfoil_static_euler
cd "${DIR}"

MAX_STEP=${MAX_STEP:-5000}
PLOT_INT=${PLOT_INT:-250}
RE=${RE:-1e4}

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
SRUN_OPTS="--mpi=pmix --ntasks=96 --cpus-per-task=1 --hint=nomultithread --threads-per-core=1"

echo "================================================================"
echo " NS mid_a12 Re=${RE}  Job ${SLURM_JOB_ID}  $(date)"
echo "================================================================"

# rebuild with REYNOLDS_USER baked in (exact same convention as the
# airfoil_static Re sweep in run_aws_all_airfoil.sbatch)
make realclean > /dev/null 2>&1
echo "building Re=${RE} ..."
make -j 96 REYNOLDS_USER=${RE} 2>&1 | tail -10

EXE=$(ls main2d.*.ex 2>/dev/null | head -1)
if [ -z "$EXE" ] || [ ! -x "$EXE" ]; then
  echo "ERROR: build failed"; exit 1
fi
echo "binary: ${EXE}"

# stage the pre-rotated .dat from the Euler case dir as a symlink so
# both cases reference identical geometry
ln -sf ${EULER_DIR}/diamond_wedge_mid_a12.dat ./diamond_wedge_mid_a12.dat

rm -rf plot_mid_a12_re1e4 surf_mid_a12_re1e4
srun ${SRUN_OPTS} ./"${EXE}" inputs_mid_a12_re1e4 \
     max_step=${MAX_STEP} \
     amr.plot_int=${PLOT_INT} \
     amr.check_int=-1 \
     2>&1 | tee run_mid_a12_re1e4.log

echo "=== done $(date) ==="
