#!/usr/bin/env python3
"""Collect axis-pressure matrix metrics and observed refinement rates."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--scheme", required=True)
    args = parser.parse_args()

    rows = []
    for path in sorted(args.result_root.glob("*_nr*/metrics.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        row["case"] = path.parent.name.rsplit("_nr", 1)[0]
        rows.append(row)

    summary: dict[str, object] = {
        "scheme": args.scheme,
        "cases": {},
        "checks": {},
    }
    passed = True
    for name in ("constant", "quadratic", "quartic", "curved"):
        selected = sorted(
            (row for row in rows if row["case"] == name),
            key=lambda row: row["nr"],
        )
        for index, row in enumerate(selected):
            row["axis_order"] = None
            row["global_linf_order"] = None
            if index:
                coarse = selected[index - 1]
                ratio = coarse["dr"] / row["dr"]
                if coarse["axis_error_linf"] > 0.0 and row["axis_error_linf"] > 0.0:
                    row["axis_order"] = math.log(
                        coarse["axis_error_linf"] / row["axis_error_linf"]
                    ) / math.log(ratio)
                if coarse["error_linf"] > 0.0 and row["error_linf"] > 0.0:
                    row["global_linf_order"] = math.log(
                        coarse["error_linf"] / row["error_linf"]
                    ) / math.log(ratio)
        summary["cases"][name] = selected
        complete = len(selected) == 4
        finite = complete and all(bool(row["finite"]) for row in selected)
        check: dict[str, object] = {
            "complete_four_grids": complete,
            "finite": finite,
        }
        if name == "constant" and complete:
            check["max_rhs_over_p_by_dr"] = max(
                float(row["rhs_over_p_by_dr"]) for row in selected
            )
            check["roundoff_balance"] = check["max_rhs_over_p_by_dr"] < 1.0e-12
        elif complete:
            finest = selected[-1]
            pressure_gradient_scale = max(2.0 / float(finest["dr"]), 1.0)
            check["finest_axis_error_over_p_by_dr"] = float(
                finest["axis_error_over_p_by_dr"]
            )
            check["finest_first4_error_over_p_by_dr"] = float(
                finest["first4_error_linf"]
            ) / pressure_gradient_scale
            check["axis_accurate"] = (
                check["finest_axis_error_over_p_by_dr"] < 1.0e-10
                and check["finest_first4_error_over_p_by_dr"] < 1.0e-9
            )
            if name == "curved":
                check["finest_global_relative_error"] = float(
                    finest["error_linf"]
                ) / max(float(finest["exact_absmax"]), 1.0)
                check["smooth_field_accurate"] = (
                    check["finest_global_relative_error"] < 1.0e-4
                )
        check_pass = complete and finite and all(
            value for key, value in check.items()
            if key in ("roundoff_balance", "axis_accurate", "smooth_field_accurate")
        )
        check["passed"] = check_pass
        summary["checks"][name] = check
        passed = passed and check_pass

    summary["passed"] = passed

    encoded = json.dumps(summary, indent=2)
    print(encoded)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded + "\n", encoding="utf-8")
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
