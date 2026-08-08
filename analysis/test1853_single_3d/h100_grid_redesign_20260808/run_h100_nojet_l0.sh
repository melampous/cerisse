#!/usr/bin/env bash
set -euo pipefail

root=/shared/cerisse_test1853_3d_20260807
build_case=/shared/cerisse_local_0808_support_hook/exm/test1853_support_hook_nojet
config_dir=${root}/analysis/h100_grid_redesign_20260808
run_dir=${root}/runs/test1853_single_h100/run165_h100_largebox_nojet_20260808/stage0c_external_l0
geometry=${root}/runs/test1853_single_h100/run165_chain_20260807T115518Z/geometry/test1853_run165_single_master_SI.stl
binding=${root}/IBM/cases/bind_test1853_h100_rank.sh
exe=main3d.gnu.TPROF.MPI.CUDA.ex
inputs=inputs.fresh_256x192x192_h100
stop_time=5.0e-4

test -x "${build_case}/${exe}"
test -f "${build_case}/prob.h"
test -f "${build_case}/support_hook_build_manifest.txt"
test -f "${config_dir}/${inputs}"
test -f "${geometry}"
test -x "${binding}"
grep -Fqx 'method=pure_shared_gp_full_cartesian' \
  "${build_case}/support_hook_build_manifest.txt"
grep -Fqx 'cut_control=disabled' \
  "${build_case}/support_hook_build_manifest.txt"
grep -Fqx 'support_hook_dispatch=automatic' \
  "${build_case}/support_hook_build_manifest.txt"
grep -Fqx 'jet_enabled_compile_time=0' \
  "${build_case}/support_hook_build_manifest.txt"

if [[ -e "${run_dir}" ]]; then
  printf 'refusing to overwrite %s\n' "${run_dir}" >&2
  exit 2
fi

source /shared/gpu_kit/gpu_env.sh
mkdir -p "${run_dir}/output/run165"
cp "${build_case}/${exe}" "${run_dir}/"
cp "${build_case}/prob.h" "${run_dir}/"
cp "${build_case}/support_hook_build_manifest.txt" "${run_dir}/"
cp "${config_dir}/${inputs}" "${run_dir}/"
cd "${run_dir}"

printf '%s\n' running > run_state.txt
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'stage=external_nojet_l0\n'
  printf 'cold_start=1\n'
  printf 'jet_enabled=0\n'
  printf 'stop_time=%s\n' "${stop_time}"
  printf 'cfl=0.15\n'
  printf 'n_cell=256 192 192\n'
  printf 'domain=-0.192:0.384 -0.216:0.216 -0.216:0.216\n'
  printf 'max_level=0\n'
  printf 'blocking_factor=32\n'
  printf 'max_grid_size=128\n'
  printf 'distribution=SFC\n'
  printf 'expected_l0_boxes=8\n'
  printf 'strict_positivity=1\n'
  printf 'stage_positivity=0\n'
  printf 'positivity_flux_limiter=0_integration_gate_limitation\n'
  printf 'soft_positivity=0\n'
  printf 'temperature_floor_injection=disabled\n'
  printf 'method=pure_shared_gp_full_cartesian\n'
  printf 'cut_control=disabled\n'
  sha256sum "${exe}" "${inputs}" prob.h "${geometry}"
} > run_manifest.txt

export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7
set +e
mpirun -np 8 --map-by slot "${binding}" "./${exe}" "${inputs}" \
  max_step=2000000 \
  stop_time="${stop_time}" \
  cfl=0.15 \
  prob.jet_enabled=0 \
  prob.expected_run_id=165 \
  ib.filename="${geometry}" \
  amr.max_level=0 \
  amr.blocking_factor=32 \
  amr.max_grid_size=128 \
  amr.regrid_int=100000 \
  amr.grid_eff=0.75 \
  DistributionMapping.strategy=SFC \
  cns.strict_positivity=1 \
  cns.stage_positivity=0 \
  cns.soft_positivity=0 \
  cns.ibm_positivity_flux_limiter=0 \
  cns.llf_fallback_mask_diagnostics=0 \
  cns.regrid_prolongation_diagnostics=0 \
  amr.plot_files_output=1 \
  amr.plot_int=200 \
  amr.plot_file=output/run165/plt \
  amr.checkpoint_files_output=1 \
  amr.check_int=100 \
  amr.check_file=output/run165/chk \
  amr.run_log=output/run165/runlog \
  amr.grid_log=output/run165/gridlog \
  ib.plot_surf=1 \
  ib.surf_int=200 \
  ib.surf_file=output/run165/surface \
  amrex.use_gpu_aware_mpi=0 \
  > run.log 2>&1
status=$?
set -e

printf '%s\n' "${status}" > exit_status.txt
final_step=$(awk '/^STEP = [0-9]+/ {step=$3} END {print step+0}' run.log)
final_time=$(awk '/^STEP = [0-9]+/ {time=$6} END {print time+0}' run.log)
failure_count=$(
  grep -Eic 'STATE CHECK FAILED|INADMISSIBLE_BEFORE|\[soft_positivity\]|(^|[^[:alpha:]])nan([^[:alpha:]]|$)|(^|[^[:alpha:]])inf([^[:alpha:]]|$)' \
    run.log || true
)
support_failure_count=$(
  grep -Eic 'coarse-fine-missing-targets=[1-9]|coarse-fine-missing-supports=[1-9]' \
    run.log || true
)
{
  printf 'final_step=%s\n' "${final_step}"
  printf 'final_time=%s\n' "${final_time}"
  printf 'exit_status=%s\n' "${status}"
  printf 'failure_count=%s\n' "${failure_count}"
  printf 'support_failure_count=%s\n' "${support_failure_count}"
  printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >> run_manifest.txt

reached_stop=$(
  awk -v actual="${final_time}" -v target="${stop_time}" \
    'BEGIN {print (actual >= target * (1.0 - 1.0e-10)) ? 1 : 0}'
)
if [[ "${status}" -eq 0 && "${reached_stop}" -eq 1 && \
      "${failure_count}" -eq 0 && "${support_failure_count}" -eq 0 ]]; then
  printf '%s\n' complete_pass > run_state.txt
else
  printf 'failed:%s_at_step:%s_time:%s\n' \
    "${status}" "${final_step}" "${final_time}" > run_state.txt
  exit 4
fi
