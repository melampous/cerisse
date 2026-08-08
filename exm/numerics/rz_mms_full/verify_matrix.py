#!/usr/bin/env python3
"""Analyze and gate the R-Z full-Euler MMS convergence matrix."""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path


RESOLUTIONS = (16, 32, 64, 128, 256)
DYNAMIC_FIELDS = ("density", "radial_momentum", "axial_momentum", "energy")


def parse_json(text: str) -> dict[str, object]:
    start = text.find("{")
    if start < 0:
        raise ValueError(f"analyzer did not emit JSON:\n{text}")
    return json.loads(text[start:])


def analyze(
    analyzer: Path,
    plotfile: Path,
    quantity: str,
    exact_time: float | None = None,
    extra_arguments: tuple[str, ...] = (),
) -> dict[str, object]:
    command = [
        sys.executable,
        str(analyzer),
        str(plotfile),
        "--quantity",
        quantity,
    ]
    if exact_time is not None:
        command.extend(("--exact-time", str(exact_time)))
    command.extend(extra_arguments)
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
    )
    return parse_json(completed.stdout)


def plotfiles(case_root: Path) -> list[Path]:
    paths = sorted((case_root / "output").glob("plt*"))
    if not paths:
        raise SystemExit(f"no plotfiles found under {case_root}")
    return paths


def orders(values: list[float]) -> list[float]:
    return [math.log(coarse / fine, 2.0) for coarse, fine in zip(values, values[1:])]


def strictly_decreasing(values: list[float]) -> bool:
    return all(coarse > fine > 0.0 for coarse, fine in zip(values, values[1:]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix_root", type=Path)
    parser.add_argument(
        "--analyzer", type=Path, default=Path(__file__).with_name("analyze.py")
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rhs: dict[str, dict[str, object]] = {}
    pressure_rhs: dict[str, dict[str, object]] = {}
    evolution: dict[str, dict[str, object]] = {}
    for nr in RESOLUTIONS:
        rhs_paths = plotfiles(args.matrix_root / "runs" / f"rhs_n{nr}")
        evolution_paths = plotfiles(
            args.matrix_root / "runs" / f"evolve_n{nr}"
        )
        pressure_paths = plotfiles(
            args.matrix_root / "runs" / f"pressure_rhs_n{nr}"
        )
        rhs[str(nr)] = analyze(args.analyzer, rhs_paths[-1], "rhs", 0.0)
        pressure_rhs[str(nr)] = analyze(
            args.analyzer,
            pressure_paths[-1],
            "rhs",
            0.0,
            (
                "--rho-amp",
                "0",
                "--ur-amp",
                "0",
                "--uz-amp",
                "0",
                "--uz0",
                "0",
            ),
        )
        evolution[str(nr)] = analyze(
            args.analyzer, evolution_paths[-1], "state"
        )

    initial_path = plotfiles(args.matrix_root / "runs" / "evolve_n16")[0]
    initial = analyze(args.analyzer, initial_path, "state", 0.0)

    convergence: dict[str, dict[str, dict[str, list[float]]]] = {
        "rhs": {},
        "pressure_rhs": {},
        "evolution": {},
    }
    for family, metrics in (
        ("rhs", rhs),
        ("pressure_rhs", pressure_rhs),
        ("evolution", evolution),
    ):
        for field in DYNAMIC_FIELDS:
            convergence[family][field] = {}
            for norm in ("l1", "l2", "linf", "axis_linf", "first_three_rings_linf"):
                values = [
                    float(metrics[str(nr)]["errors"][field][norm])
                    for nr in RESOLUTIONS
                ]
                convergence[family][field][norm] = values
                convergence[family][field][f"{norm}_orders"] = orders(values)

    all_numeric = []
    for family in (rhs, pressure_rhs, evolution):
        for case in family.values():
            for field_metrics in case["errors"].values():
                all_numeric.extend(float(value) for value in field_metrics.values())

    initial_linf = max(
        float(initial["errors"][field]["linf"]) for field in DYNAMIC_FIELDS
    )
    final_times = [float(evolution[str(nr)]["plot_time"]) for nr in RESOLUTIONS]
    minimum_density = min(
        float(evolution[str(nr)]["admissibility"]["minimum_density"])
        for nr in RESOLUTIONS
    )
    minimum_pressure = min(
        float(evolution[str(nr)]["admissibility"]["minimum_cell_pressure"])
        for nr in RESOLUTIONS
    )
    full_axis_errors = convergence["rhs"]["radial_momentum"]["axis_linf"]
    pressure_axis_errors = convergence["pressure_rhs"]["radial_momentum"][
        "axis_linf"
    ]
    pressure_axis_orders = convergence["pressure_rhs"]["radial_momentum"][
        "axis_linf_orders"
    ]
    full_axis_orders = convergence["rhs"]["radial_momentum"][
        "axis_linf_orders"
    ]
    full_first_three_orders = convergence["rhs"]["radial_momentum"][
        "first_three_rings_linf_orders"
    ]
    full_linf_orders = convergence["rhs"]["radial_momentum"]["linf_orders"]
    evolution_first_three_orders = convergence["evolution"]["radial_momentum"][
        "first_three_rings_linf_orders"
    ]
    pressure_linf_errors = convergence["pressure_rhs"]["radial_momentum"][
        "linf"
    ]
    fine_pressure_fraction = pressure_axis_errors[-1] / full_axis_errors[-1]

    checks = {
        "all_metrics_finite": all(math.isfinite(value) for value in all_numeric),
        "initial_cell_averages_match_oracle": initial_linf < 5.0e-13,
        "all_evolution_runs_reach_target_time": all(
            abs(time - 2.0e-2) < 5.0e-13 for time in final_times
        ),
        "density_and_pressure_remain_positive": (
            minimum_density > 0.5 and minimum_pressure > 0.5
        ),
        "azimuthal_momentum_remains_zero": all(
            float(case["errors"]["azimuthal_momentum"]["linf"]) < 1.0e-13
            for family in (rhs, pressure_rhs, evolution)
            for case in family.values()
        ),
        "rhs_errors_decrease": all(
            strictly_decreasing(convergence["rhs"][field]["l1"])
            for field in DYNAMIC_FIELDS
        ),
        "rhs_asymptotic_order_is_at_least_1p7": all(
            min(convergence["rhs"][field]["l1_orders"][-2:]) > 1.7
            for field in DYNAMIC_FIELDS
        ),
        "rhs_axis_radial_momentum_is_second_order": min(
            full_axis_orders[-2:]
        )
        > 1.7,
        "rhs_first_three_radial_momentum_rings_are_second_order": min(
            full_first_three_orders[-2:]
        ) > 1.7,
        "rhs_radial_momentum_linf_is_second_order": min(
            full_linf_orders[-2:]
        ) > 1.7,
        "pressure_operator_preserves_regular_quadratic": max(
            pressure_linf_errors
        ) < 5.0e-11,
        "evolution_errors_decrease": all(
            strictly_decreasing(convergence["evolution"][field]["l1"])
            for field in DYNAMIC_FIELDS
        ),
        "evolution_asymptotic_order_is_at_least_1p7": all(
            min(convergence["evolution"][field]["l1_orders"][-2:]) > 1.7
            for field in DYNAMIC_FIELDS
        ),
        "evolution_first_three_radial_momentum_rings_converge": min(
            evolution_first_three_orders[-2:]
        ) > 1.7,
    }
    passed = all(checks.values())
    report = {
        "schema_version": 1,
        "matrix_root": str(args.matrix_root),
        "resolutions": list(RESOLUTIONS),
        "rhs": rhs,
        "pressure_rhs": pressure_rhs,
        "evolution": evolution,
        "initial_oracle": initial,
        "convergence": convergence,
        "summary": {
            "initial_max_linf": initial_linf,
            "final_times": final_times,
            "minimum_density": minimum_density,
            "minimum_cell_pressure": minimum_pressure,
            "full_radial_momentum_axis_orders": full_axis_orders,
            "full_radial_momentum_first_three_ring_orders": (
                full_first_three_orders
            ),
            "pressure_only_radial_momentum_axis_orders": pressure_axis_orders,
            "pressure_only_radial_momentum_max_linf": max(
                pressure_linf_errors
            ),
            "fine_pressure_only_fraction_of_full_axis_error": (
                fine_pressure_fraction
            ),
        },
        "checks": checks,
        "passed": passed,
    }

    output = args.output or args.matrix_root / "verification.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    marker = args.matrix_root / (
        "VERIFICATION_PASS" if passed else "VERIFICATION_FAIL"
    )
    marker.write_text(json.dumps(checks, indent=2) + "\n", encoding="ascii")
    print(json.dumps(report, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
