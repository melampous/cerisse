#!/usr/bin/env bash

# Short, all-fluid R-Z strong-gradient gate.  This script deliberately stays
# case-local: it selects already-implemented production operators but does not
# alter any flux, IBM, cut-cell, EB, or update code.
set -u
set -o pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
result_root="${1:-${case_dir}/results/nonlinear_axis_matrix_$(date -u +%Y%m%dT%H%M%SZ)}"
seed_cells="${SEED_CELLS:-0.0}"
seed_rings="${SEED_RINGS:-0}"
schemes=(llf-wenoz5 afd-hllc-wenoz5 skew4-central skew4-jst)

if [[ -e "${result_root}" ]]; then
  echo "Refusing to overwrite ${result_root}" >&2
  exit 2
fi
mkdir -p "${result_root}"

scheme_tag() {
  case "$1" in
    llf-wenoz5) echo llf_wenoz5 ;;
    afd-hllc-wenoz5) echo afd_hllc_wenoz5 ;;
    skew4-central) echo skew4_central ;;
    skew4-jst) echo skew4_jst ;;
    *) return 2 ;;
  esac
}

for scheme in "${schemes[@]}"; do
  tag="$(scheme_tag "${scheme}")"
  binary="${case_dir}/main2d.gnu.${tag}.ex"
  if [[ ! -x "${binary}" ]]; then
    echo "Missing ${binary}; build RZ_PLANAR_SHOCK_SCHEME=${scheme}" >&2
    exit 2
  fi
  for fallback in 0 1; do
    out="${result_root}/${tag}_fallback${fallback}"
    mkdir -p "${out}"
    fallback_args=(cns.llf_first_order_fallback="${fallback}")
    if [[ "${scheme}" == skew4-* ]]; then
      fallback_args=(cns.skew_first_order_fallback="${fallback}")
    fi

    set +e
    "${binary}" "${case_dir}/inputs" \
      max_step=80 stop_time=2.0e-2 cfl=0.25 \
      cns.order_rk=3 cns.stages_rk=4 \
      amr.max_level=0 amr.n_cell="16 64" \
      amr.plot_file="${out}/plt" amr.plot_int=10 \
      prob.mach_up=10.0 prob.flow_sign=1 \
      prob.shock_checkerboard_cells="${seed_cells}" \
      prob.shock_checkerboard_radial_cells="${seed_rings}" \
      "${fallback_args[@]}" \
      >"${out}/run.log" 2>&1
    status=$?
    set -e

    printf '%s\n' "${status}" >"${out}/exit_status.txt"
    if compgen -G "${out}/plt[0-9]*/Header" >/dev/null; then
      python3 "${case_dir}/analyze_nonlinear_axis_gate.py" "${out}" \
        --output "${out}/metrics.json" >"${out}/analysis.log" 2>&1
    fi
  done
done

python3 - "${result_root}" "${seed_cells}" "${seed_rings}" <<'PY'
import csv
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
rows = []
for case in sorted(path for path in root.iterdir() if path.is_dir()):
    status_path = case / "exit_status.txt"
    metrics_path = case / "metrics.json"
    row = {
        "case": case.name,
        "exit_status": int(status_path.read_text().strip()),
        "last_step": -1,
        "max_dp_radial_rel": float("nan"),
        "max_drho_radial_rel": float("nan"),
        "max_abs_ur_over_a": float("nan"),
        "max_abs_axis_ur_over_a": float("nan"),
        "max_abs_even_p_axis_proxy": float("nan"),
        "max_abs_even_rho_axis_proxy": float("nan"),
        "max_abs_odd_ur_axis_proxy": float("nan"),
        "shock_corrugation_cells": float("nan"),
    }
    if metrics_path.is_file():
        data = json.loads(metrics_path.read_text())
        hist = data["history_maxima"]
        row.update({
            "last_step": data["final"]["step"],
            "max_dp_radial_rel": hist["max_radial_pressure_spread_rel_local"],
            "max_drho_radial_rel": hist["max_radial_density_spread_rel_local"],
            "max_abs_ur_over_a": hist["max_abs_radial_velocity_over_sound"],
            "max_abs_axis_ur_over_a": hist["max_axis_abs_radial_velocity_over_sound"],
            "max_abs_even_p_axis_proxy": hist[
                "max_abs_signed_axis_to_second_ring_pressure_rel"
            ],
            "max_abs_even_rho_axis_proxy": hist[
                "max_abs_signed_axis_to_second_ring_density_rel"
            ],
            "max_abs_odd_ur_axis_proxy": hist[
                "max_abs_signed_odd_axis_extrapolation_ur_over_sound"
            ],
            "shock_corrugation_cells": hist["shock_front_corrugation_cells"],
        })
    rows.append(row)

summary = {
    "case": "rz_planar_shock_short_nonlinear_axis_matrix",
    "comparison_scope": (
        "whole spatial-operator A/B; LLF-WENO-Z5 versus AFD-HLLC-WENO-Z5 "
        "is not a Riemann-solver-only comparison"
    ),
    "grid": [16, 64],
    "cfl": 0.25,
    "max_step": 80,
    "mach": 10.0,
    "shock_checkerboard_cells": float(sys.argv[2]),
    "shock_checkerboard_radial_cells": int(sys.argv[3]),
    "rows": rows,
}
(root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
with (root / "summary.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(summary, indent=2))
PY

