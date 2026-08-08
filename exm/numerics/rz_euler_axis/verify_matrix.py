#!/usr/bin/env python3
"""Analyze and gate the complete R-Z Euler-axis regression matrix."""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path


CASES = {
    **{f"m0_legacy_n{nr}": 0 for nr in (16, 32, 64, 128)},
    "m0_annular_n16": 0,
    "m0_corrected_n16": 0,
    **{f"m1_corrected_n{nr}": 1 for nr in (16, 32, 64, 128)},
}


def parse_json(text: str) -> dict[str, object]:
    start = text.find("{")
    if start < 0:
        raise ValueError(f"analyzer did not emit JSON:\n{text}")
    return json.loads(text[start:])


def analyze_case(analyzer: Path, plotfile: Path, mode: int) -> dict[str, object]:
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")
    completed = subprocess.run(
        [sys.executable, str(analyzer), str(plotfile), "--mode", str(mode)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
    )
    return parse_json(completed.stdout)


def orders(values: list[float]) -> list[float]:
    return [math.log(a / b, 2.0) for a, b in zip(values[:-1], values[1:])]


def decreasing(values: list[float]) -> bool:
    return all(a > b > 0.0 for a, b in zip(values[:-1], values[1:]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix_root", type=Path)
    parser.add_argument(
        "--analyzer", type=Path, default=Path(__file__).with_name("analyze.py")
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    metrics: dict[str, dict[str, object]] = {}
    for label, mode in CASES.items():
        plotfiles = sorted((args.matrix_root / "runs" / label / "output").glob("plt*"))
        if not plotfiles:
            raise SystemExit(f"no plotfile found for {label}")
        metrics[label] = analyze_case(args.analyzer, plotfiles[-1], mode)

    resolutions = [16, 32, 64, 128]
    legacy_density = [
        float(metrics[f"m0_legacy_n{nr}"]["density_error_linf"])
        for nr in resolutions
    ]
    mode1_density = [
        float(metrics[f"m1_corrected_n{nr}"]["density_error_linf"])
        for nr in resolutions
    ]
    mode1_momentum = [
        float(metrics[f"m1_corrected_n{nr}"]["momentum_error_linf"])
        for nr in resolutions
    ]
    legacy_density_orders = orders(legacy_density)
    density_orders = orders(mode1_density)
    momentum_orders = orders(mode1_momentum)

    annular = metrics["m0_annular_n16"]
    corrected = metrics["m0_corrected_n16"]
    checks = {
        "all_metrics_finite": all(
            math.isfinite(float(value))
            for case in metrics.values()
            for key, value in case.items()
            if key.endswith("_linf")
        ),
        "legacy_density_defect_is_nonconvergent": (
            min(legacy_density) > 5.0e-2
            and max(abs(value) for value in legacy_density_orders) < 5.0e-4
        ),
        "annular_only_fixes_continuity": float(annular["density_error_linf"])
        < 1.0e-10,
        "annular_only_leaves_momentum_defect": float(annular["momentum_error_linf"])
        > 1.0e-5,
        "corrected_mode0_is_roundoff": (
            float(corrected["density_error_linf"]) < 1.0e-10
            and float(corrected["momentum_error_linf"]) < 1.0e-10
        ),
        "corrected_mode1_errors_decrease": decreasing(mode1_density)
        and decreasing(mode1_momentum),
        "corrected_mode1_density_is_second_order": min(density_orders) > 1.95,
        "corrected_mode1_momentum_is_second_order": min(momentum_orders) > 1.90,
    }
    passed = all(checks.values())
    report = {
        "schema_version": 1,
        "matrix_root": str(args.matrix_root),
        "resolutions": resolutions,
        "metrics": metrics,
        "convergence": {
            "mode1_density_linf": mode1_density,
            "mode1_density_orders": density_orders,
            "mode1_momentum_linf": mode1_momentum,
            "mode1_momentum_orders": momentum_orders,
            "legacy_mode0_density_linf": legacy_density,
            "legacy_mode0_density_orders": legacy_density_orders,
        },
        "checks": checks,
        "passed": passed,
    }

    output = args.output or args.matrix_root / "verification.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    marker = args.matrix_root / ("VERIFICATION_PASS" if passed else "VERIFICATION_FAIL")
    marker.write_text(json.dumps(checks, indent=2) + "\n", encoding="ascii")
    print(json.dumps(report, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
