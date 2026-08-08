#!/usr/bin/env bash
set -euo pipefail

root=/shared/cerisse_test1853_3d_20260807
chain="${root}/runs/test1853_single_h100/run165_chain_20260807T115518Z"
diagnostics="${root}/runs/test1853_single_h100/run165_llfdiag_3way_20260808"
analysis_dir="${root}/analysis/energy_history_20260808"
replay="${diagnostics}/testI_energy_history_replay"

while [[ "$(cat "${replay}/run_state.txt" 2>/dev/null || true)" == running ]] || \
      [[ ! -e "${replay}/run_state.txt" ]]; do
  sleep 20
done

state=$(cat "${replay}/run_state.txt")
if [[ "${state}" != complete_target ]]; then
  printf 'energy replay did not complete: %s\n' "${state}" >&2
  exit 2
fi

plots=(
  "${chain}/stage1_external_l0/output/plt00000"
  "${chain}/stage1_external_l0/output/plt00400"
  "${chain}/stage1_external_l0/output/plt00800"
  "${chain}/stage1_external_l0/output/plt01200"
  "${chain}/stage1_external_l0/output/plt01339"
  "${chain}/stage2_jet_l1/output/plt01600"
  "${chain}/stage2_jet_l1/output/plt02000"
  "${chain}/stage2_jet_l1/output/plt02400"
  "${chain}/stage2_jet_l1/output/plt02800"
  "${diagnostics}/testG_thermodynamic_continuation/output/run165/plt03200"
  "${replay}/output/run165/plt03310"
  "${replay}/output/run165/plt03320"
  "${replay}/output/run165/plt03330"
  "${replay}/output/run165/plt03340"
  "${replay}/output/run165/plt03350"
  "${replay}/output/run165/plt03360"
  "${diagnostics}/testH_soft_Tmin_continuation/output/run165/plt03400"
  "${diagnostics}/testH_soft_Tmin_continuation/output/run165/plt03600"
  "${diagnostics}/testH_soft_Tmin_continuation/output/run165/plt03800"
  "${diagnostics}/testH_soft_Tmin_continuation/output/run165/plt04000"
)
labels=(
  step0000
  step0400
  step0800
  step1200
  step1339
  step1600
  step2000
  step2400
  step2800
  step3200
  step3310
  step3320
  step3330
  step3340
  step3350
  step3360
  step3400_soft
  step3600_soft
  step3800_soft
  step4000_soft
)

for plot in "${plots[@]}"; do
  test -f "${plot}/Header"
done

label_args=()
for label in "${labels[@]}"; do
  label_args+=(--label "${label}")
done

printf '%s\n' running > "${analysis_dir}/full_audit_state.txt"
set +e
nice -n 10 /tmp/run165_yt/bin/python \
  "${analysis_dir}/analyze_run165_energy_history.py" \
  "${plots[@]}" \
  --output "${analysis_dir}/run165_energy_history_full.json" \
  --low-records 64 \
  "${label_args[@]}" \
  > "${analysis_dir}/full_audit.log" 2>&1
status=$?
set -e
printf '%s\n' "${status}" > "${analysis_dir}/full_audit_exit_status.txt"
if [[ "${status}" -eq 0 ]]; then
  printf '%s\n' complete > "${analysis_dir}/full_audit_state.txt"
else
  printf 'failed:%s\n' "${status}" > "${analysis_dir}/full_audit_state.txt"
fi
exit "${status}"
