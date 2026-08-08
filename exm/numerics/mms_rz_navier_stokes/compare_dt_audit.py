#!/usr/bin/env python3
"""Compare exact-solution errors for the N=64 time-step audit."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


VARIABLES = (
    "density",
    "radial_momentum",
    "axial_momentum",
    "energy",
)


def select_row(path: Path, grid: int):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    matching = [row for row in rows if int(row["base_n"]) == grid]
    if len(matching) != 1:
        raise ValueError(
            f"{path} contains {len(matching)} rows for N={grid}"
        )
    return matching[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline_csv", type=Path)
    parser.add_argument("half_dt_csv", type=Path)
    parser.add_argument("--grid", type=int, default=64)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    baseline = select_row(args.baseline_csv, args.grid)
    half_dt = select_row(args.half_dt_csv, args.grid)
    baseline_time = float(baseline["plot_time"])
    half_dt_time = float(half_dt["plot_time"])
    if abs(baseline_time - half_dt_time) > 1.0e-14:
        raise ValueError(
            f"final-time mismatch: {baseline_time} versus {half_dt_time}"
        )

    rows = []
    for variable in VARIABLES:
        baseline_error = float(baseline[f"{variable}_L2"])
        half_dt_error = float(half_dt[f"{variable}_L2"])
        change = half_dt_error - baseline_error
        relative_change = abs(change) / baseline_error
        rows.append(
            {
                "base_n": args.grid,
                "final_time": baseline_time,
                "variable": variable,
                "baseline_L2_error": baseline_error,
                "half_dt_L2_error": half_dt_error,
                "half_to_baseline_ratio": half_dt_error
                / baseline_error,
                "absolute_error_change": change,
                "relative_error_change": relative_change,
            }
        )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"N={row['base_n']} {row['variable']:17s} "
            f"baseline={row['baseline_L2_error']:.6e} "
            f"half_dt={row['half_dt_L2_error']:.6e} "
            f"relative_change={row['relative_error_change']:.3e}"
        )
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    main()
