#!/bin/bash
# Binary-level LLF-WENO-Z5 RZ stability matrix for the frozen Panda L4 case.
# This job never writes into the source case.  It runs each variant in a fresh
# subdirectory and records the executable, input, and problem-definition hashes.
#
# Required submission environment:
#   SOURCE_CASE=/shared/ujet_2d/b2d_baseline_newsrc_L4
#   CAMPAIGN_ROOT=/shared/ujet_2d/rz_llfweno_binary_ab_<tag>

#SBATCH --job-name=rz_llfweno_ab
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%x-%j.out

set -u

SOURCE_CASE="${SOURCE_CASE:-/shared/ujet_2d/b2d_baseline_newsrc_L4}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:?set CAMPAIGN_ROOT to a fresh directory}"
EXE="${SOURCE_CASE}/main2d.gnu.TPROF.MPI.ex"
INPUTS="${SOURCE_CASE}/inputs_L4"

if [ ! -x "${EXE}" ] || [ ! -f "${INPUTS}" ]; then
  echo "ERROR: frozen executable or inputs are missing from ${SOURCE_CASE}" >&2
  exit 2
fi

if [ -e "${CAMPAIGN_ROOT}/STARTED_UTC" ]; then
  echo "ERROR: refusing to reuse campaign ${CAMPAIGN_ROOT}" >&2
  exit 3
fi

mkdir -p "${CAMPAIGN_ROOT}"
date -u +%FT%TZ > "${CAMPAIGN_ROOT}/STARTED_UTC"
{
  echo "job_id=${SLURM_JOB_ID}"
  echo "node_list=${SLURM_JOB_NODELIST}"
  echo "source_case=${SOURCE_CASE}"
  sha256sum "${EXE}" "${INPUTS}" "${SOURCE_CASE}/prob.h"
} > "${CAMPAIGN_ROOT}/manifest.txt"

module purge
module load openmpi/4.1.7
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

run_variant () {
  local name="$1"
  local cfl_value="$2"
  local epsilon_value="$3"
  local regrid_interval="$4"
  local run_dir="${CAMPAIGN_ROOT}/${name}"
  local log_file="${run_dir}/run.log"

  mkdir -p "${run_dir}"
  {
    echo "name=${name}"
    echo "cfl=${cfl_value}"
    echo "cns.rz_weno_epsilon_relative=${epsilon_value}"
    echo "amr.regrid_int=${regrid_interval}"
  } > "${run_dir}/variant.txt"

  cd "${run_dir}" || exit 4
  set +e
  timeout 90 /opt/amazon/openmpi/bin/mpirun \
    --bind-to core --map-by ppr:96:node \
    --mca orte_abort_on_non_zero_status 1 \
    -np 96 "${EXE}" "${INPUTS}" \
      max_step=400 \
      stop_time=1.0e-4 \
      cfl="${cfl_value}" \
      amr.regrid_int="${regrid_interval}" \
      cns.rz_weno_epsilon_relative="${epsilon_value}" \
      cns.record_stats=0 \
      cns.nstep_screen_output=1 \
      cns.verbose=0 \
      amr.plot_files_output=0 \
      amr.checkpoint_files_output=0 \
      amr.run_log="${run_dir}/amr_run.log" \
      amr.grid_log="${run_dir}/amr_grid.log" \
      > "${log_file}" 2>&1
  local rc=$?
  set -e

  echo "${rc}" > "${run_dir}/exit_code"
  grep -E "STATE CHECK FAILED|min\(rho\)|argmin\(rho\*e_int\)|STEP =|AMReX .* finalized" \
    "${log_file}" | tail -40 > "${run_dir}/summary.txt" || true
  cd "${CAMPAIGN_ROOT}" || exit 4
}

set -e
run_variant cfl015_eps1e6 0.15 1.0e-6 80
run_variant cfl015_eps0    0.15 0.0    80
run_variant cfl030_eps0    0.30 0.0    40

date -u +%FT%TZ > "${CAMPAIGN_ROOT}/COMPLETED_UTC"
echo "RZ_LLFWENO_BINARY_AB_COMPLETE"
