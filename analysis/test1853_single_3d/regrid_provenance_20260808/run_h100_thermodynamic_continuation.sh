#!/usr/bin/env bash
set -euo pipefail

case_root=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_llfdiag_3way_20260808
build_case=/shared/cerisse_local_0808_thermodynamic_prolongation/exm/test1853_thermodynamic_prolongation
run_dir="${case_root}/testG_thermodynamic_continuation"
restart="${case_root}/testB_normal_regrid/chk03059"
geometry=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_chain_20260807T115518Z/geometry/test1853_run165_single_master_SI.stl
exe=main3d.gnu.TPROF.MPI.CUDA.ex
expected_exe_sha=90a733b431a9a27782862877ee29f3aa8025910c8571032ff0c4cc077f131a5d
target_max_step=4000
deadline_utc=2026-08-08T10:45:00Z

test -x "${build_case}/${exe}"
test -f "${build_case}/thermodynamic_build_manifest.txt"
test -d "${restart}"
test -f "${geometry}"
test "$(sha256sum "${build_case}/${exe}" | cut -d' ' -f1)" = "${expected_exe_sha}"
grep -Fqx 'method=pure_shared_gp_full_cartesian' \
  "${build_case}/thermodynamic_build_manifest.txt"
grep -Fqx 'cut_control=disabled' \
  "${build_case}/thermodynamic_build_manifest.txt"
grep -Fqx \
  'amr_solver=/shared/cerisse_local_0808_thermodynamic_prolongation' \
  "${build_case}/thermodynamic_build_manifest.txt"

if [[ -e "${run_dir}" ]]; then
  printf 'ERROR: refusing to overwrite existing run directory: %s\n' \
    "${run_dir}" >&2
  exit 2
fi

source /shared/gpu_kit/gpu_env.sh
mkdir -p "${run_dir}/output/run165"
cp "${build_case}/${exe}" "${run_dir}/"
cp "${build_case}/thermodynamic_build_manifest.txt" "${run_dir}/"
cp "${case_root}/inputs.run165" "${run_dir}/"
cp "${case_root}/prob.h" "${run_dir}/"

cd "${run_dir}"
printf '%s\n' running > run_state.txt
printf '%s\n' "$$" > launcher_pid.txt
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'restart=%s\n' "${restart}"
  printf 'geometry=%s\n' "${geometry}"
  printf 'target_max_step=%s\n' "${target_max_step}"
  printf 'capacity_block_deadline_guard=%s\n' "${deadline_utc}"
  printf 'cfl=0.01\n'
  printf 'regrid_int=4,4\n'
  printf 'check_int=100\n'
  printf 'plot_int=200\n'
  printf 'surface_int=100\n'
  printf 'strict_positivity=1\n'
  printf 'soft_positivity=0\n'
  printf 'prolongation_floor=max(pressure_rhoe,rho*cv*T_min)\n'
  sha256sum "${exe}" inputs.run165 prob.h
} > run_manifest.txt

deadline_epoch=$(date -ud "${deadline_utc}" +%s)
now_epoch=$(date -u +%s)
if (( now_epoch >= deadline_epoch )); then
  printf '%s\n' 'ERROR: continuation requested after safety deadline' \
    > run_state.txt
  exit 3
fi

guard_delay=$((deadline_epoch - now_epoch))
(
  sleep "${guard_delay}"
  touch dump_and_stop
  printf '[%s] capacity-block guard requested graceful checkpoint and stop\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> deadline_guard.log
) &
guard_pid=$!
printf '%s\n' "${guard_pid}" > deadline_guard_pid.txt

stop_guard() {
  if kill -0 "${guard_pid}" 2>/dev/null; then
    kill "${guard_pid}" 2>/dev/null || true
    wait "${guard_pid}" 2>/dev/null || true
  fi
}
trap stop_guard EXIT

export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

set +e
mpirun -np 8 "./${exe}" inputs.run165 \
  amr.restart="${restart}" \
  ib.filename="${geometry}" \
  amr.regrid_on_restart=1 \
  max_step="${target_max_step}" \
  stop_time=1.0 \
  cns.restart_first_dt_from_cfl=0 \
  cns.strict_positivity=1 \
  cns.stage_positivity=0 \
  cns.soft_positivity=0 \
  cns.ibm_positivity_flux_limiter=0 \
  cns.llf_fallback_mask_diagnostics=1 \
  cns.regrid_prolongation_diagnostics=1 \
  amr.plot_files_output=1 \
  amr.plot_file=output/run165/plt \
  amr.plot_int=200 \
  amr.checkpoint_files_output=1 \
  amr.check_file=output/run165/chk \
  amr.check_int=100 \
  ib.surf_file=output/run165/surface \
  ib.surf_int=100 \
  cns.verbose=0 \
  > run.log 2>&1
status=$?
set -e

printf '%s\n' "${status}" > exit_status.txt
final_step=$(awk '/^STEP = [0-9]+/ {step=$3} END {print step+0}' run.log)
printf 'final_step=%s\n' "${final_step}" >> run_manifest.txt
if [[ "${status}" -eq 0 ]]; then
  if (( final_step >= target_max_step )); then
    printf '%s\n' complete_target > run_state.txt
  else
    printf 'stopped_cleanly_at_step:%s\n' "${final_step}" > run_state.txt
  fi
else
  printf 'failed:%s\n' "${status}" > run_state.txt
fi
printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  >> run_manifest.txt
exit "${status}"
