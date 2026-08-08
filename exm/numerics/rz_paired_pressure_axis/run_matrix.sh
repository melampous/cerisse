#!/usr/bin/env bash

set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
scheme="${RZ_AXIS_SCHEME:-llf-wenoz5}"
exe=""
result_root=""

usage() {
  cat <<EOF
Usage: $0 [--scheme SCHEME] [--binary PATH] [--result-root DIR] [DIR]

SCHEME is one of llf-wenoz5, llf-teno5, llf-teno6, afd-hllc-wenoz5, or
afd-hllc-teno5.  A final positional DIR is retained as a backwards-compatible
alias for --result-root.
EOF
}

while (($#)); do
  case "$1" in
    --scheme)
      [[ $# -ge 2 ]] || { echo "--scheme requires a value" >&2; exit 2; }
      scheme="$2"
      shift 2
      ;;
    --binary)
      [[ $# -ge 2 ]] || { echo "--binary requires a value" >&2; exit 2; }
      exe="$2"
      shift 2
      ;;
    --result-root)
      [[ $# -ge 2 ]] || { echo "--result-root requires a value" >&2; exit 2; }
      result_root="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      [[ -z "${result_root}" ]] || {
        echo "Only one result directory may be supplied" >&2
        exit 2
      }
      result_root="$1"
      shift
      ;;
  esac
done

case "${scheme}" in
  llf-wenoz5)
    scheme_tag=llf_wenoz5
    numerics_scheme=characteristic_llf_weno_z5
    operator_manifest='rz_llf_weno_pressure_operator=paired_metric_advection_and_pressure_gradient_when_RZ'
    ;;
  llf-teno5)
    scheme_tag=llf_teno5
    numerics_scheme=characteristic_llf_teno5
    operator_manifest='rz_llf_weno_pressure_operator=paired_metric_advection_and_pressure_gradient_when_RZ'
    ;;
  llf-teno6)
    scheme_tag=llf_teno6
    numerics_scheme=characteristic_llf_teno6
    operator_manifest='rz_llf_weno_pressure_operator=paired_metric_advection_and_pressure_gradient_when_RZ'
    ;;
  afd-hllc-wenoz5)
    scheme_tag=afd_hllc_wenoz5
    numerics_scheme=afd_hllc_weno_z5
    operator_manifest=''
    ;;
  afd-hllc-teno5)
    scheme_tag=afd_hllc_teno5
    numerics_scheme=afd_hllc_teno5
    operator_manifest=''
    ;;
  *)
    echo "Unsupported scheme: ${scheme}" >&2
    usage >&2
    exit 2
    ;;
esac

exe="${exe:-${case_dir}/main2d.gnu.${scheme_tag}.ex}"
result_root="${result_root:-${case_dir}/results/axis_pressure_${scheme_tag}_$(date -u +%Y%m%dT%H%M%SZ)}"
test_manifest="[TestManifest] case=rz_paired_pressure_axis scheme=${scheme}"
scheme_manifest="[Numerics] inviscid_scheme=${numerics_scheme}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/cerisse-rz-axis-mpl}"

if [[ -e "${result_root}" ]]; then
  echo "Refusing to overwrite existing result directory: ${result_root}" >&2
  exit 1
fi
mkdir -p "${MPLCONFIGDIR}"

if [[ ! -x "${exe}" ]]; then
  echo "Missing ${exe}; build it with:" >&2
  echo "  make -C ${case_dir} RZ_AXIS_SCHEME=${scheme} -j" >&2
  exit 1
fi

run_case() {
  local name="$1"
  local nr="$2"
  local nz="$3"
  local profile="$4"
  local amplitude="$5"
  local quartic="$6"
  local out="${result_root}/${name}_nr${nr}"

  mkdir -p "${out}"
  "${exe}" "${case_dir}/inputs" \
    amr.n_cell="${nr} ${nz}" \
    prob.profile="${profile}" \
    prob.amplitude="${amplitude}" \
    prob.quartic="${quartic}" \
    amr.plot_file="${out}/plt" \
    >"${out}/run.log" 2>&1

  if ! grep -Fq "${scheme_manifest}" "${out}/run.log"; then
    echo "Production scheme manifest mismatch in ${out}/run.log; refusing stale binary" >&2
    exit 2
  fi
  if [[ -n "${operator_manifest}" ]] &&
     ! grep -Fq "${operator_manifest}" "${out}/run.log"; then
    echo "Paired-pressure manifest missing in ${out}/run.log" >&2
    exit 2
  fi
  if ! grep -Fq "${test_manifest}" "${out}/run.log"; then
    echo "Scheme manifest mismatch in ${out}/run.log; refusing stale binary" >&2
    exit 2
  fi

  python3 "${case_dir}/analyze.py" "${out}/plt00001" \
    --profile "${profile}" \
    --amplitude "${amplitude}" \
    --quartic "${quartic}" \
    --json "${out}/metrics.json" \
    --figure "${out}/rhs_error.png" \
    >"${out}/analysis.log"
}

for nr in 16 32 64 128; do
  run_case constant "${nr}" 16 0 0.0 0.0
  run_case quadratic "${nr}" 16 0 0.25 0.0
  run_case quartic "${nr}" 16 2 0.25 0.10
  run_case curved "${nr}" "$((4 * nr))" 1 0.25 0.10
done

python3 "${case_dir}/summarize_matrix.py" "${result_root}" \
  --scheme "${scheme}" \
  --output "${result_root}/summary.json"

echo "Axis-pressure matrix complete: ${result_root}"
