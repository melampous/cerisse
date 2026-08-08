#!/usr/bin/env bash
set -euo pipefail

case_root=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_llfdiag_3way_20260808
build_case=/shared/cerisse_local_0808_admissible_prolongation/exm/test1853_admissible_prolongation
run_dir="${case_root}/testE_admissible_prolongation"
restart="${case_root}/testB_normal_regrid/chk03059"
geometry=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_chain_20260807T115518Z/geometry/test1853_run165_single_master_SI.stl
exe=main3d.gnu.TPROF.MPI.CUDA.ex

test -x "${build_case}/${exe}"
test -f "${build_case}/repair_build_manifest.txt"
test -d "${restart}"
test -f "${geometry}"
if [[ -e "${run_dir}" ]]; then
  printf 'ERROR: refusing to overwrite existing run directory: %s\n' "${run_dir}" >&2
  exit 2
fi

mkdir -p "${run_dir}/output/run165"
cp "${build_case}/${exe}" "${run_dir}/"
cp "${build_case}/repair_build_manifest.txt" "${run_dir}/"
cp "${case_root}/inputs.run165" "${run_dir}/"
cp "${case_root}/prob.h" "${run_dir}/"

cd "${run_dir}"
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'restart=%s\n' "${restart}"
  printf 'geometry=%s\n' "${geometry}"
  printf 'target_max_step=3061\n'
  printf 'time_step=normal_checkpoint_cfl_path\n'
  sha256sum "${exe}" inputs.run165 prob.h
} > run_manifest.txt

export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

set +e
mpirun -np 8 "./${exe}" inputs.run165 \
  amr.restart="${restart}" \
  ib.filename="${geometry}" \
  amr.regrid_on_restart=1 \
  max_step=3061 \
  stop_time=1.0 \
  cns.restart_first_dt_from_cfl=0 \
  cns.stage_positivity=0 \
  cns.ibm_positivity_flux_limiter=0 \
  cns.llf_fallback_mask_diagnostics=1 \
  cns.regrid_prolongation_diagnostics=1 \
  amr.plot_files_output=0 \
  amr.checkpoint_files_output=0 \
  cns.verbose=0 \
  > run.log 2>&1
status=$?
set -e

printf '%s\n' "${status}" > exit_status.txt
printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  >> run_manifest.txt
exit "${status}"
