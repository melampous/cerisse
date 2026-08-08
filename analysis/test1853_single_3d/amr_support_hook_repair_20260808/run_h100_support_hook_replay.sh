#!/usr/bin/env bash
set -euo pipefail

old_root=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_llfdiag_3way_20260808
build_case=/shared/cerisse_local_0808_support_hook/exm/test1853_support_hook
run_root=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_support_hook_repair_20260808
run_dir=${run_root}/testB_chk03300_regrid_step3370
restart=${old_root}/testG_thermodynamic_continuation/output/run165/chk03300
geometry=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_chain_20260807T115518Z/geometry/test1853_run165_single_master_SI.stl
exe=main3d.gnu.TPROF.MPI.CUDA.ex

test -x "${build_case}/${exe}"
test -f "${build_case}/support_hook_build_manifest.txt"
test -f "${restart}/Header"
test -f "${geometry}"
grep -Fqx 'method=pure_shared_gp_full_cartesian' \
  "${build_case}/support_hook_build_manifest.txt"
grep -Fqx 'cut_control=disabled' \
  "${build_case}/support_hook_build_manifest.txt"
grep -Fqx 'support_hook_dispatch=automatic' \
  "${build_case}/support_hook_build_manifest.txt"
grep -Fqx 'legacy_support_cap_macro=absent' \
  "${build_case}/support_hook_build_manifest.txt"

if [[ -e "${run_dir}" ]]; then
  printf 'refusing to overwrite %s\n' "${run_dir}" >&2
  exit 2
fi

source /shared/gpu_kit/gpu_env.sh
mkdir -p "${run_dir}/output/run165"
cp "${build_case}/${exe}" "${run_dir}/"
cp "${build_case}/support_hook_build_manifest.txt" "${run_dir}/"
cp "${build_case}/inputs.run165" "${run_dir}/"
cp "${build_case}/prob.h" "${run_dir}/"
cd "${run_dir}"

printf '%s\n' running > run_state.txt
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'restart=%s\n' "${restart}"
  printf 'target_step=3370\n'
  printf 'regrid_on_restart=1\n'
  printf 'regrid_int=1000000 1000000\n'
  printf 'subsequent_regrid=outside_test_window\n'
  printf 'strict_positivity=1\n'
  printf 'soft_positivity=0\n'
  printf 'method=pure_shared_gp_full_cartesian\n'
  printf 'cut_control=disabled\n'
  sha256sum "${exe}" inputs.run165 prob.h "${restart}/Header"
} > run_manifest.txt

export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7
set +e
mpirun -np 8 "./${exe}" inputs.run165 \
  amr.restart="${restart}" \
  ib.filename="${geometry}" \
  amr.regrid_on_restart=1 \
  amr.regrid_int='1000000 1000000' \
  max_step=3370 \
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
  amr.check_int=70 \
  ib.surf_int=1000000 \
  cns.verbose=0 \
  > run.log 2>&1
status=$?
set -e

printf '%s\n' "${status}" > exit_status.txt
final_step=$(awk '/^STEP = [0-9]+/ {step=$3} END {print step+0}' run.log)
hook_line=$(grep -F '[IBM-AMR-Support-Hook]' run.log | head -n 1 || true)
level2_gp_line=$(grep -F '[IBM-AMR-Support] GP level=2' run.log | tail -n 1 || true)
failure_count=$(
  grep -Ec 'STATE CHECK FAILED|INADMISSIBLE_BEFORE|\[soft_positivity\]' \
    run.log || true
)
{
  printf 'final_step=%s\n' "${final_step}"
  printf 'exit_status=%s\n' "${status}"
  printf 'failure_count=%s\n' "${failure_count}"
  printf 'hook_line=%s\n' "${hook_line}"
  printf 'level2_gp_line=%s\n' "${level2_gp_line}"
  printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >> run_manifest.txt

if [[ "${status}" -eq 0 && "${final_step}" -ge 3370 && \
      "${failure_count}" -eq 0 && \
      "${hook_line}" == *'detected=1'* && \
      "${level2_gp_line}" == *'targets=0 '* ]]; then
  printf '%s\n' complete_pass > run_state.txt
else
  printf 'failed:%s_at_step:%s\n' "${status}" "${final_step}" > run_state.txt
  exit 4
fi
