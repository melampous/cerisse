#!/usr/bin/env bash
set -euo pipefail

case_root=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_llfdiag_3way_20260808
source_run="${case_root}/testH_soft_Tmin_continuation"
run_dir="${case_root}/testI_energy_history_replay"
restart="${case_root}/testG_thermodynamic_continuation/output/run165/chk03300"
geometry=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_chain_20260807T115518Z/geometry/test1853_run165_single_master_SI.stl
exe=main3d.gnu.TPROF.MPI.CUDA.ex
expected_exe_sha=90a733b431a9a27782862877ee29f3aa8025910c8571032ff0c4cc077f131a5d
latest_start_utc=2026-08-08T10:20:00Z
stop_guard_utc=2026-08-08T10:55:00Z

while [[ "$(cat "${source_run}/run_state.txt" 2>/dev/null || true)" == running ]]; do
  if (( $(date -u +%s) >= $(date -ud "${latest_start_utc}" +%s) )); then
    printf 'source continuation did not finish before %s\n' "${latest_start_utc}" >&2
    exit 3
  fi
  sleep 20
done

test -x "${source_run}/${exe}"
test -f "${source_run}/thermodynamic_build_manifest.txt"
test -f "${restart}/Header"
test -f "${geometry}"
test "$(sha256sum "${source_run}/${exe}" | cut -d' ' -f1)" = "${expected_exe_sha}"
grep -Fqx 'method=pure_shared_gp_full_cartesian' \
  "${source_run}/thermodynamic_build_manifest.txt"
grep -Fqx 'cut_control=disabled' \
  "${source_run}/thermodynamic_build_manifest.txt"

if [[ -e "${run_dir}" ]]; then
  printf 'refusing to overwrite %s\n' "${run_dir}" >&2
  exit 2
fi

source /shared/gpu_kit/gpu_env.sh
mkdir -p "${run_dir}/output/run165"
cp "${source_run}/${exe}" "${run_dir}/"
cp "${source_run}/thermodynamic_build_manifest.txt" "${run_dir}/"
cp "${source_run}/inputs.run165" "${run_dir}/"
cp "${source_run}/prob.h" "${run_dir}/"
cd "${run_dir}"

printf '%s\n' running > run_state.txt
printf '%s\n' "$$" > launcher_pid.txt
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'run_class=deterministic_pre_failure_energy_history\n'
  printf 'restart=%s\n' "${restart}"
  printf 'target_max_step=3366\n'
  printf 'plot_interval=10\n'
  printf 'checkpoint_interval=60\n'
  printf 'strict_positivity=1\n'
  printf 'soft_positivity=0\n'
  printf 'method=pure_shared_gp_full_cartesian\n'
  printf 'cut_control=disabled\n'
  sha256sum "${exe}" inputs.run165 prob.h "${restart}/Header"
} > run_manifest.txt

guard_delay=$(( $(date -ud "${stop_guard_utc}" +%s) - $(date -u +%s) ))
if (( guard_delay <= 0 )); then
  printf '%s\n' deadline_before_start > run_state.txt
  exit 4
fi
(
  sleep "${guard_delay}"
  touch dump_and_stop
) &
guard_pid=$!
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
  max_step=3366 \
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
  amr.plot_int=10 \
  amr.checkpoint_files_output=1 \
  amr.check_file=output/run165/chk \
  amr.check_int=60 \
  ib.surf_int=1000000 \
  cns.verbose=0 \
  > run.log 2>&1
status=$?
set -e

printf '%s\n' "${status}" > exit_status.txt
final_step=$(awk '/^STEP = [0-9]+/ {step=$3} END {print step+0}' run.log)
printf 'final_step=%s\n' "${final_step}" >> run_manifest.txt
printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> run_manifest.txt
if [[ "${status}" -eq 0 && "${final_step}" -eq 3366 ]]; then
  printf '%s\n' complete_target > run_state.txt
elif [[ "${status}" -eq 0 ]]; then
  printf 'stopped_cleanly_at_step:%s\n' "${final_step}" > run_state.txt
else
  printf 'failed:%s_at_step:%s\n' "${status}" "${final_step}" > run_state.txt
fi
exit "${status}"
