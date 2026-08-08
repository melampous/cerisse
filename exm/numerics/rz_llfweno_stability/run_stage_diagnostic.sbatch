#!/bin/bash
# Build the current source snapshot and locate the first inadmissible SSPRK
# stage in the frozen all-fluid Panda RZ L4 case.
#
# Required submission environment:
#   AMR_SOLVER=/shared/cerisse_llfweno_rz_diag_<tag>/source
#   CASE_ROOT=/shared/cerisse_llfweno_rz_diag_<tag>/case
#   RESULT_ROOT=/shared/cerisse_llfweno_rz_diag_<tag>/stage_baseline

#SBATCH --job-name=rz_llfweno_stage
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --time=00:25:00
#SBATCH --output=slurm-%x-%j.out

set -u

AMR_SOLVER="${AMR_SOLVER:?set AMR_SOLVER}"
CASE_ROOT="${CASE_ROOT:?set CASE_ROOT}"
RESULT_ROOT="${RESULT_ROOT:?set RESULT_ROOT to a fresh directory}"

if [ ! -f "${AMR_SOLVER}/src/Make.CNS" ] ||
   [ ! -d "${AMR_SOLVER}/lib/amrex" ]; then
  echo "ERROR: invalid source snapshot ${AMR_SOLVER}" >&2
  exit 2
fi
if [ ! -f "${CASE_ROOT}/prob.h" ] ||
   [ ! -f "${CASE_ROOT}/inputs_L4" ] ||
   [ ! -f "${CASE_ROOT}/GNUmakefile.cpu" ]; then
  echo "ERROR: incomplete case snapshot ${CASE_ROOT}" >&2
  exit 3
fi
if [ -e "${RESULT_ROOT}" ]; then
  echo "ERROR: refusing to reuse result directory ${RESULT_ROOT}" >&2
  exit 4
fi

mkdir -p "${RESULT_ROOT}"
date -u +%FT%TZ > "${RESULT_ROOT}/STARTED_UTC"

module purge
module load openmpi/4.1.7
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${CASE_ROOT}" || exit 5
BUILD_NAME="rz_llfweno_stage_${SLURM_JOB_ID}"
BUILD_DIR="tmp_build_${BUILD_NAME}"
set +e
make -f GNUmakefile.cpu -j64 \
  AMR_SOLVER="${AMR_SOLVER}" \
  EBASE="${BUILD_NAME}" \
  TMP_BUILD_DIR="${BUILD_DIR}" \
  > "${RESULT_ROOT}/build.log" 2>&1
build_rc=$?
set -e
echo "${build_rc}" > "${RESULT_ROOT}/build_exit_code"
if [ "${build_rc}" -ne 0 ]; then
  exit "${build_rc}"
fi

EXE=$(find . -maxdepth 1 -type f \
  -name "${BUILD_NAME}2d.gnu.*MPI.ex" -perm -111 | sort | head -1)
if [ -z "${EXE}" ]; then
  echo "ERROR: build produced no executable" >&2
  exit 6
fi

{
  echo "job_id=${SLURM_JOB_ID}"
  echo "node_list=${SLURM_JOB_NODELIST}"
  echo "source=${AMR_SOLVER}"
  echo "case=${CASE_ROOT}"
  sha256sum "${EXE}" "${CASE_ROOT}/inputs_L4" "${CASE_ROOT}/prob.h" \
    "${AMR_SOLVER}/src/rhs/Weno.h" \
    "${AMR_SOLVER}/src/tim/advance.cpp" \
    "${AMR_SOLVER}/src/tim/compute_rhs.cpp" \
    "${AMR_SOLVER}/src/CNS.cpp"
} > "${RESULT_ROOT}/manifest.txt"

set +e
/opt/amazon/openmpi/bin/mpirun \
  --bind-to core --map-by ppr:96:node \
  --mca orte_abort_on_non_zero_status 1 \
  -np 96 "${EXE}" "${CASE_ROOT}/inputs_L4" \
    max_step=130 \
    stop_time=2.0e-3 \
    cfl=0.30 \
    cns.rz_weno_epsilon_relative=1.0e-6 \
    cns.stage_state_diagnostics=1 \
    cns.check_state=1 \
    cns.record_stats=0 \
    cns.nstep_screen_output=1 \
    cns.verbose=0 \
    amr.plot_files_output=0 \
    amr.checkpoint_files_output=0 \
    amr.run_log="${RESULT_ROOT}/amr_run.log" \
    amr.grid_log="${RESULT_ROOT}/amr_grid.log" \
    > "${RESULT_ROOT}/run.log" 2>&1
run_rc=$?
set -e
echo "${run_rc}" > "${RESULT_ROOT}/run_exit_code"
grep -E "SSPRK STAGE STATE FAILED|stage =|cell =|rho\*e_int|raw conservative|STEP =|AMReX .* finalized" \
  "${RESULT_ROOT}/run.log" | tail -80 > "${RESULT_ROOT}/summary.txt" || true

date -u +%FT%TZ > "${RESULT_ROOT}/COMPLETED_UTC"
echo "RZ_LLFWENO_STAGE_DIAGNOSTIC_COMPLETE rc=${run_rc}"
