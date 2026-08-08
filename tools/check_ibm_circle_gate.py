#!/usr/bin/env python3
"""Enforce the steady Mach-4 circle bow-shock regression gate.

This is a software-regression gate for one frozen N=320 GPIBM protocol.  It is
not a substitute for a matched body-fitted validation solution.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path
from typing import Any


DEFAULT_N320_REFERENCE = 0.5643135034711624
DEFAULT_CELLS_PER_DIAMETER = 64.0


def load_history(directory: Path) -> list[dict[str, Any]]:
    reports: list[dict[str, Any]] = []
    for filename in glob.glob(str(directory / "plt*_standoff.json")):
        with open(filename, encoding="utf-8") as stream:
            report = json.load(stream)
        report["_filename"] = filename
        reports.append(report)
    reports.sort(
        key=lambda report: (
            float(report.get("time", -math.inf)),
            report["_filename"],
        )
    )
    return reports


def finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def assess_history(
    directory: Path,
    u_inf: float,
    diameter: float,
    expected_cells_per_diameter: float,
    cells_per_diameter_tolerance: float,
    minimum_t_star: float,
    plateau_points: int,
    maximum_plateau_span: float,
    maximum_method_span_cells: float,
    minimum_standoff_cells: float,
) -> tuple[dict[str, Any], list[str]]:
    failures: list[str] = []
    reports = load_history(directory)
    result: dict[str, Any] = {
        "directory": str(directory),
        "reports": len(reports),
    }
    if len(reports) < plateau_points:
        failures.append(f"only {len(reports)} reports; need {plateau_points}")
        return result, failures

    tail = reports[-plateau_points:]
    invalid_tail = [
        Path(report["_filename"]).name
        for report in tail
        if report.get("recommended") is None
    ]
    if invalid_tail:
        failures.append(
            "detached bow shock was not identified in the final steady window: "
            + ", ".join(invalid_tail)
        )
        return result, failures

    final = tail[-1]
    estimate = final["recommended"]
    final_time = finite_float(final.get("time"))
    cells_per_diameter = finite_float(final.get("cells_per_diameter"))
    standoff_over_radius = finite_float(estimate.get("standoff_over_radius"))
    standoff_cells = finite_float(estimate.get("standoff_cells"))
    method_span_cells = finite_float(estimate.get("method_span_cells"))
    values = [
        finite_float(report["recommended"].get("standoff_over_radius"))
        for report in tail
    ]

    required_values = {
        "final time": final_time,
        "cells per diameter": cells_per_diameter,
        "stand-off distance": standoff_over_radius,
        "stand-off cells": standoff_cells,
        "density/pressure method span": method_span_cells,
    }
    for name, value in required_values.items():
        if value is None:
            failures.append(f"{name} is missing or non-finite")
    if any(value is None for value in values):
        failures.append("steady-window stand-off history contains a non-finite value")
    if failures:
        return result, failures

    assert final_time is not None
    assert cells_per_diameter is not None
    assert standoff_over_radius is not None
    assert standoff_cells is not None
    assert method_span_cells is not None
    finite_values = [float(value) for value in values if value is not None]

    final_t_star = final_time * u_inf / diameter
    mean_value = sum(finite_values) / len(finite_values)
    relative_plateau_span = (max(finite_values) - min(finite_values)) / abs(mean_value)

    if (
        abs(cells_per_diameter - expected_cells_per_diameter)
        > cells_per_diameter_tolerance
    ):
        failures.append(
            "wrong frozen grid: cells/diameter "
            f"{cells_per_diameter:.12g} differs from {expected_cells_per_diameter:.12g}"
        )
    if final_t_star < minimum_t_star:
        failures.append(
            f"final convective time {final_t_star:.6g} < {minimum_t_star:.6g}"
        )
    if relative_plateau_span > maximum_plateau_span:
        failures.append(
            "stand-off history is not stationary: relative span "
            f"{relative_plateau_span:.6g} > {maximum_plateau_span:.6g}"
        )
    if standoff_cells < minimum_standoff_cells:
        failures.append(
            "detached shock is unresolved or absent: stand-off "
            f"{standoff_cells:.6g} cells < {minimum_standoff_cells:.6g}"
        )
    if method_span_cells > maximum_method_span_cells:
        failures.append(
            "density/pressure shock locations disagree by "
            f"{method_span_cells:.6g} cells > {maximum_method_span_cells:.6g}"
        )

    result.update(
        {
            "final_time": final_time,
            "final_convective_time": final_t_star,
            "cells_per_diameter": cells_per_diameter,
            "standoff_over_radius": standoff_over_radius,
            "standoff_cells": standoff_cells,
            "method_span_cells": method_span_cells,
            "plateau_points": plateau_points,
            "relative_plateau_span": relative_plateau_span,
        }
    )
    return result, failures


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_directory", type=Path)
    parser.add_argument(
        "--reference-value",
        type=float,
        default=DEFAULT_N320_REFERENCE,
        help="frozen N=320 steady software reference for Delta/R",
    )
    parser.add_argument("--u-inf", type=float, default=1389.00065)
    parser.add_argument("--diameter", type=float, default=0.2)
    parser.add_argument(
        "--expected-cells-per-diameter",
        type=float,
        default=DEFAULT_CELLS_PER_DIAMETER,
    )
    parser.add_argument("--cells-per-diameter-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--minimum-t-star", type=float, default=40.0)
    parser.add_argument("--plateau-points", type=int, default=4)
    parser.add_argument("--maximum-plateau-span", type=float, default=0.002)
    parser.add_argument("--maximum-method-span-cells", type=float, default=0.5)
    parser.add_argument("--minimum-standoff-cells", type=float, default=4.0)
    parser.add_argument("--reference-relative-tolerance", type=float, default=0.02)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.u_inf <= 0.0 or args.diameter <= 0.0:
        parser.error("--u-inf and --diameter must be positive")
    if args.plateau_points < 2:
        parser.error("--plateau-points must be at least two")
    if args.reference_value <= 0.0:
        parser.error("--reference-value must be positive")

    candidate, failures = assess_history(
        args.candidate_directory,
        args.u_inf,
        args.diameter,
        args.expected_cells_per_diameter,
        args.cells_per_diameter_tolerance,
        args.minimum_t_star,
        args.plateau_points,
        args.maximum_plateau_span,
        args.maximum_method_span_cells,
        args.minimum_standoff_cells,
    )

    if "standoff_over_radius" in candidate:
        relative_difference = (
            candidate["standoff_over_radius"] / args.reference_value - 1.0
        )
        candidate["relative_difference_from_frozen_reference"] = relative_difference
        if abs(relative_difference) > args.reference_relative_tolerance:
            failures.append(
                "stand-off differs from the frozen N=320 reference by "
                f"{relative_difference:.6%}; tolerance is "
                f"{args.reference_relative_tolerance:.6%}"
            )

        billig = 0.386 * math.exp(4.67 / 4.0**2)
        gamma = 1.4
        epsilon = (gamma - 1.0 + 2.0 / 4.0**2) / (gamma + 1.0)
        hornung = 2.14 * epsilon * (1.0 + 0.5 * epsilon)
        candidate["relative_difference_from_billig"] = (
            candidate["standoff_over_radius"] / billig - 1.0
        )
        candidate["relative_difference_from_hornung"] = (
            candidate["standoff_over_radius"] / hornung - 1.0
        )

    report = {
        "status": "PASS" if not failures else "FAIL",
        "candidate": candidate,
        "frozen_n320_reference": args.reference_value,
        "criteria": {
            "expected_cells_per_diameter": args.expected_cells_per_diameter,
            "minimum_convective_time": args.minimum_t_star,
            "plateau_points": args.plateau_points,
            "maximum_relative_plateau_span": args.maximum_plateau_span,
            "minimum_standoff_cells": args.minimum_standoff_cells,
            "maximum_method_span_cells": args.maximum_method_span_cells,
            "reference_relative_tolerance": args.reference_relative_tolerance,
        },
        "failures": failures,
        "stop_policy": (
            "On FAIL, stop unrelated IBM development and diagnose the Mach-4 "
            "circle bow shock before proceeding."
        ),
        "note": (
            "This is a frozen-grid software-regression gate. Billig and Hornung "
            "are diagnostics, not exact matched Euler acceptance references."
        ),
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    raise SystemExit(0 if not failures else 1)


if __name__ == "__main__":
    main()
