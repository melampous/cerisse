#!/usr/bin/env bash

set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
scheme="${RZ_PLANAR_SHOCK_SCHEME:-llf-wenoz5}"
binary=""
result_root=""

usage() {
  cat <<EOF
Usage: $0 [--scheme SCHEME] [--binary PATH] [--result-root DIR] [DIR]

SCHEME is one of llf-wenoz5, llf-teno5, llf-teno6, afd-hllc-wenoz5, or
afd-hllc-teno5.  The default binary is the uniquely named non-MPI build.  Use
--binary for an MPI or custom build.  A final positional DIR is retained as an
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
      binary="$2"
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

binary="${binary:-${case_dir}/main2d.gnu.${scheme_tag}.ex}"
result_root="${result_root:-${case_dir}/results/dynamic_axis_gate_${scheme_tag}_$(date -u +%Y%m%dT%H%M%SZ)}"
test_manifest="[TestManifest] case=rz_planar_shock scheme=${scheme}"
scheme_manifest="[Numerics] inviscid_scheme=${numerics_scheme}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/cerisse-rz-axis-mpl}"

if [[ ! -x "${binary}" ]]; then
  echo "Missing ${binary}; build it with:" >&2
  echo "  make -C ${case_dir} USE_MPI=FALSE RZ_PLANAR_SHOCK_SCHEME=${scheme} -j" >&2
  exit 1
fi
if [[ -e "${result_root}" ]]; then
  echo "Refusing to overwrite existing result directory: ${result_root}" >&2
  exit 1
fi
mkdir -p "${result_root}" "${MPLCONFIGDIR}"

run_case() {
  local name="$1"
  local mach="$2"
  local flow_sign="$3"
  local out="${result_root}/${name}"
  mkdir -p "${out}"

  "${binary}" "${case_dir}/inputs" \
    max_step=80 stop_time=2.0e-2 \
    cns.order_rk=3 cns.stages_rk=4 \
    amr.max_level=0 amr.n_cell="16 64" \
    amr.plot_file="${out}/plt" amr.plot_int=20 \
    prob.mach_up="${mach}" prob.flow_sign="${flow_sign}" \
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
  python3 "${case_dir}/analyze_current_axis_gate.py" "${out}" \
    --output "${out}/metrics.json" >"${out}/analysis.log"
}

run_case m3_plus 3.0 1
run_case m10_plus 10.0 1
run_case m10_minus 10.0 -1

python3 - "${result_root}" "${scheme}" <<'PY'
import json
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
scheme = sys.argv[2]
cases = {}
threshold = 1.0e-10
for path in sorted(root.glob("*/metrics.json")):
    data = json.loads(path.read_text(encoding="utf-8"))
    maxima = data["history_maxima"]
    cases[path.parent.name] = {
        "max_radial_pressure_spread_rel_local":
            maxima["max_radial_pressure_spread_rel_local"],
        "max_axis_to_second_ring_pressure_delta_rel_local":
            maxima["max_axis_to_second_ring_pressure_delta_rel_local"],
        "max_abs_radial_velocity_over_sound":
            maxima["max_abs_radial_velocity_over_sound"],
    }
passed = len(cases) == 3 and all(
    all(math.isfinite(value) and value < threshold for value in values.values())
    for values in cases.values()
)
summary = {
    "scheme": scheme,
    "threshold": threshold,
    "passed": passed,
    "cases": cases,
}
(root / "summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(summary, indent=2, sort_keys=True))
if not passed:
    raise SystemExit(3)
PY

python3 "${case_dir}/plot_dynamic_axis_gate.py" "${result_root}" \
  --scheme "${scheme}" \
  --output "${result_root}/dynamic_axis_pressure_gate_summary.png"

echo "Dynamic axis-pressure gate complete: ${result_root}"
