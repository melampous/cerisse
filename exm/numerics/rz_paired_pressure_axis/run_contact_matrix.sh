#!/usr/bin/env bash
set -euo pipefail

case_dir=$(cd -- "$(dirname -- "$0")" && pwd)
cd "$case_dir"

scheme="${RZ_AXIS_SCHEME:-llf-wenoz5}"
binary=""
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
      echo "unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      [[ -z "$result_root" ]] || {
        echo "only one result directory may be supplied" >&2
        exit 2
      }
      result_root="$1"
      shift
      ;;
  esac
done

case "$scheme" in
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
    echo "unsupported scheme: $scheme" >&2
    usage >&2
    exit 2
    ;;
esac

binary="${binary:-./main2d.gnu.${scheme_tag}.ex}"
result_root="${result_root:-results/contact_${scheme_tag}_$(date -u +%Y%m%dT%H%M%SZ)}"
test_manifest="[TestManifest] case=rz_paired_pressure_axis scheme=${scheme}"
scheme_manifest="[Numerics] inviscid_scheme=${numerics_scheme}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/cerisse-rz-axis-mpl}"
mkdir -p "$MPLCONFIGDIR"

if [[ ! -x "$binary" ]]; then
  echo "missing $binary; build it with:" >&2
  echo "  make -C $case_dir RZ_AXIS_SCHEME=$scheme -j" >&2
  exit 2
fi

run_case() {
  local name=$1
  local orientation=$2
  local location=$3
  local out="$result_root/$name"
  if [[ -e "$out" ]]; then
    echo "refusing to overwrite existing result directory: $out" >&2
    exit 2
  fi
  mkdir -p "$out"
  "$binary" inputs_contact \
    prob.contact_orientation="$orientation" \
    prob.contact_location="$location" \
    amr.plot_file="$out/plt" >"$out/run.log" 2>&1
  if ! grep -Fq "$scheme_manifest" "$out/run.log"; then
    echo "production scheme manifest mismatch in $out/run.log; refusing stale binary" >&2
    exit 3
  fi
  if [[ -n "$operator_manifest" ]] &&
     ! grep -Fq "$operator_manifest" "$out/run.log"; then
    echo "paired-pressure manifest missing from $out/run.log" >&2
    exit 3
  fi
  if ! grep -Fq "$test_manifest" "$out/run.log"; then
    echo "scheme manifest mismatch in $out/run.log; refusing stale binary" >&2
    exit 3
  fi
  python3 analyze_contact.py "$out/plt00100" \
    --json "$out/metrics.json" >"$out/analysis.log"
}

# Strongest near-axis gate: only the first radial ring is on the dense side.
run_case radial_face1 0 0.015625
# A slightly wider axis collar distinguishes first-face effects from a generic
# radial contact response.
run_case radial_face4 0 0.0625
# Orthogonal control: a radially uniform planar contact, normal to z.
run_case axial_plane 1 0.5

python3 - "$result_root" "$scheme" <<'PY'
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
scheme = sys.argv[2]
cases = {}
for path in sorted(root.glob("*/metrics.json")):
    cases[path.parent.name] = json.loads(path.read_text(encoding="utf-8"))
threshold = 1.0e-10
passed = len(cases) == 3 and all(
    row["pressure_linf_rel"] < threshold
    and row["near_axis_pressure_linf_rel"] < threshold
    and row["radial_mach_linf"] < threshold
    and row["near_axis_radial_mach_linf"] < threshold
    for row in cases.values()
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
    raise SystemExit(4)
PY
