#!/usr/bin/env python3
"""Assemble the no-solid R-Z Euler (RZ-1) verification record."""

from __future__ import annotations

import argparse
import csv
import filecmp
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")
import analyze as rz_analyze


RESOLUTIONS = (16, 32, 64, 128, 256)
FIELDS = ("density", "radial_momentum", "axial_momentum", "energy")


def parse_analyzer_output(text: str) -> dict[str, object]:
    start = text.find("{")
    if start < 0:
        raise RuntimeError(f"analyze.py did not emit JSON:\n{text}")
    return json.loads(text[start:])


def analyze_plot(
    analyzer: Path,
    plotfile: Path,
    quantity: str,
    extra_arguments: tuple[str, ...] = (),
) -> dict[str, object]:
    command = [
        sys.executable,
        str(analyzer),
        str(plotfile),
        "--quantity",
        quantity,
    ]
    if quantity == "rhs":
        command.extend(("--exact-time", "0"))
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
    return parse_analyzer_output(completed.stdout)


def latest_plotfile(case_root: Path) -> Path:
    paths = sorted(case_root.glob("plt*"))
    if not paths:
        raise RuntimeError(f"no plotfile found under {case_root}")
    return paths[-1]


def observed_orders(values: list[float]) -> list[float]:
    return [
        math.log(coarse / fine, 2.0)
        for coarse, fine in zip(values, values[1:])
    ]


def convergence_table(
    cases: dict[str, dict[str, object]],
) -> dict[str, dict[str, list[float]]]:
    result: dict[str, dict[str, list[float]]] = {}
    for field in FIELDS:
        result[field] = {}
        for norm in (
            "l1",
            "l2",
            "linf",
            "axis_linf",
            "first_three_rings_linf",
        ):
            values = [
                float(cases[str(n)]["errors"][field][norm])
                for n in RESOLUTIONS
            ]
            result[field][norm] = values
            result[field][f"{norm}_orders"] = observed_orders(values)
    return result


def audit_fe_csv(path: Path) -> dict[str, object]:
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    per_step: dict[int, int] = {}
    for row in rows:
        step = int(row["step"])
        per_step[step] = per_step.get(step, 0) + 1
    return {
        "path": str(path),
        "records": len(rows),
        "steps": len(per_step),
        "all_steps_have_four_brackets": all(v == 4 for v in per_step.values()),
        "minimum_density": min(float(row["min_rho"]) for row in rows),
        "minimum_internal_energy_density": min(
            float(row["min_rhoe"]) for row in rows
        ),
        "invalid_cells": sum(int(row["invalid_cells"]) for row in rows),
    }


def load_fields(plotfile: Path) -> dict[str, np.ndarray]:
    rz_analyze.yt.funcs.mylog.setLevel(50)
    dataset = rz_analyze.load_boxlib_rz(plotfile)
    dimensions = np.asarray(dataset.domain_dimensions, dtype=int)
    grid = dataset.covering_grid(0, dataset.domain_left_edge, dimensions)
    return {
        name: rz_analyze.raw_field(grid, plot_name)
        for name, plot_name in rz_analyze.FIELD_NAMES.items()
    }


def compare_plotfiles(reference: Path, candidate: Path) -> dict[str, object]:
    reference_fields = load_fields(reference)
    candidate_fields = load_fields(candidate)
    fields: dict[str, dict[str, float | bool]] = {}
    for name in rz_analyze.FIELD_NAMES:
        lhs = reference_fields[name]
        rhs = candidate_fields[name]
        difference = np.abs(lhs - rhs)
        fields[name] = {
            "array_equal": bool(np.array_equal(lhs, rhs)),
            "maximum_absolute_difference": float(np.max(difference)),
            "relative_to_reference_linf": float(
                np.max(difference) / max(float(np.max(np.abs(lhs))), 1.0e-300)
            ),
        }
    return {
        "reference": str(reference),
        "candidate": str(candidate),
        "fields": fields,
        "all_fields_bitwise_equal": all(
            bool(value["array_equal"]) for value in fields.values()
        ),
        "maximum_absolute_difference": max(
            float(value["maximum_absolute_difference"])
            for value in fields.values()
        ),
    }


def compare_plot_trees(reference: Path, candidate: Path) -> bool:
    relative_files = sorted(
        path.relative_to(reference) for path in reference.rglob("*") if path.is_file()
    )
    candidate_files = sorted(
        path.relative_to(candidate) for path in candidate.rglob("*") if path.is_file()
    )
    return relative_files == candidate_files and all(
        filecmp.cmp(reference / path, candidate / path, shallow=False)
        for path in relative_files
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parent
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    diagnostics = root / "diagnostics_deconv"
    analyzer = root / "analyze.py"
    rhs: dict[str, dict[str, object]] = {}
    pressure_rhs: dict[str, dict[str, object]] = {}
    evolution: dict[str, dict[str, object]] = {}
    pressure_arguments = (
        "--rho-amp",
        "0",
        "--ur-amp",
        "0",
        "--uz-amp",
        "0",
        "--uz0",
        "0",
    )
    fe_audits: list[dict[str, object]] = []

    for resolution in RESOLUTIONS:
        key = str(resolution)
        rhs[key] = analyze_plot(
            analyzer,
            latest_plotfile(diagnostics / f"full_n{resolution}"),
            "rhs",
        )
        pressure_rhs[key] = analyze_plot(
            analyzer,
            latest_plotfile(
                root
                / "preflight_matrix_deconv_r1_serial"
                / "runs"
                / f"pressure_rhs_n{resolution}"
                / "output"
            ),
            "rhs",
            pressure_arguments,
        )
        evolution[key] = analyze_plot(
            analyzer,
            latest_plotfile(diagnostics / f"evolve_n{resolution}"),
            "state",
        )
        fe_audits.append(audit_fe_csv(diagnostics / f"evolve_n{resolution}" / "fe.csv"))

    rhs_convergence = convergence_table(rhs)
    pressure_convergence = convergence_table(pressure_rhs)
    evolution_convergence = convergence_table(evolution)
    cpu_fixed_box = latest_plotfile(diagnostics / "fab_np1_n64")
    mpi_fixed_box = latest_plotfile(diagnostics / "fab_np2_n64")
    gpu_fixed_box = latest_plotfile(diagnostics / "gpu_n64")
    cpu_single_box = latest_plotfile(diagnostics / "evolve_n64")

    mpi_comparison = compare_plotfiles(cpu_fixed_box, mpi_fixed_box)
    box_comparison = compare_plotfiles(cpu_single_box, cpu_fixed_box)
    gpu_comparison = compare_plotfiles(cpu_fixed_box, gpu_fixed_box)
    gpu_fe_audit = audit_fe_csv(diagnostics / "gpu_n64" / "fe.csv")

    default_off_reference = (
        root
        / "preflight_matrix_fix_r4_serial"
        / "runs"
        / "evolve_n64"
        / "output"
        / "plt00021"
    )
    default_off_candidate = root / "diagnostics_default_off_n64" / "plt00021"
    default_off_bitwise = compare_plot_trees(
        default_off_reference, default_off_candidate
    )

    fe_records = sum(int(audit["records"]) for audit in fe_audits)
    fe_invalid = sum(int(audit["invalid_cells"]) for audit in fe_audits)
    minimum_density = min(float(audit["minimum_density"]) for audit in fe_audits)
    minimum_rhoe = min(
        float(audit["minimum_internal_energy_density"]) for audit in fe_audits
    )
    pressure_radial_linf = pressure_convergence["radial_momentum"]["linf"]

    checks = {
        "semidiscrete_global_l1_asymptotic_order_at_least_1p8": all(
            min(rhs_convergence[field]["l1_orders"][-2:]) > 1.8
            for field in FIELDS
        ),
        "semidiscrete_axis_radial_momentum_order_at_least_2": min(
            rhs_convergence["radial_momentum"]["axis_linf_orders"][-2:]
        ) > 2.0,
        "pressure_flux_source_balance_is_roundoff_limited": max(
            pressure_radial_linf
        ) < 5.0e-11,
        "ssprk43_global_l1_asymptotic_order_at_least_1p8": all(
            min(evolution_convergence[field]["l1_orders"][-2:]) > 1.8
            for field in FIELDS
        ),
        "ssprk43_axis_radial_momentum_order_at_least_1p8": min(
            evolution_convergence["radial_momentum"]["axis_linf_orders"][-2:]
        ) > 1.8,
        "all_632_forward_euler_brackets_are_admissible": (
            fe_records == 632
            and fe_invalid == 0
            and minimum_density > 0.0
            and minimum_rhoe > 0.0
            and all(
                bool(audit["all_steps_have_four_brackets"])
                for audit in fe_audits
            )
        ),
        "default_off_path_is_bitwise_unchanged": default_off_bitwise,
        "fixed_boxarray_is_mpi_bitwise_invariant": bool(
            mpi_comparison["all_fields_bitwise_equal"]
        ),
        "boxarray_difference_is_below_5e_minus_9": (
            float(box_comparison["maximum_absolute_difference"]) < 5.0e-9
        ),
        "cpu_gpu_difference_is_below_1e_minus_12": (
            float(gpu_comparison["maximum_absolute_difference"]) < 1.0e-12
            and int(gpu_fe_audit["invalid_cells"]) == 0
        ),
    }

    report = {
        "schema_version": 1,
        "scope": (
            "no-solid uniform level-0 ideal-gas R-Z Euler, LLF-WENO-Z5, "
            "SSPRK(4,3)"
        ),
        "resolutions": list(RESOLUTIONS),
        "rhs": rhs,
        "pressure_rhs": pressure_rhs,
        "evolution": evolution,
        "convergence": {
            "rhs": rhs_convergence,
            "pressure_rhs": pressure_convergence,
            "evolution": evolution_convergence,
        },
        "forward_euler_audits": fe_audits,
        "forward_euler_summary": {
            "records": fe_records,
            "invalid_cells": fe_invalid,
            "minimum_density": minimum_density,
            "minimum_internal_energy_density": minimum_rhoe,
        },
        "default_off_bitwise": default_off_bitwise,
        "fixed_boxarray_mpi_comparison": mpi_comparison,
        "boxarray_comparison": box_comparison,
        "cpu_gpu_comparison": gpu_comparison,
        "gpu_forward_euler_audit": gpu_fe_audit,
        "checks": checks,
        "passed": all(checks.values()),
    }

    output = args.output or root / "rz1_verification.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
