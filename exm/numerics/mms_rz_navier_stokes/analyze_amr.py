#!/usr/bin/env python3
"""Leaf-cell errors for the fixed-band R-Z MMS AMR verification."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

import numpy as np

# yt imports Matplotlib even though this script does not create figures.  Use
# a writable cache location on restricted compute nodes.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib-cache")

import analyze as uniform_analysis


VARIABLE_FIELDS = {
    "density": "Density",
    "radial_momentum": "Xmom",
    "axial_momentum": "Ymom",
    "energy": "Energy",
}


def masked_error_metrics(error, volume, mask):
    """Return normalized norms over one nonempty leaf-cell selection."""
    if not np.any(mask):
        return {"L1": float("nan"), "L2": float("nan"), "Linf": float("nan")}
    selected_volume = volume[mask]
    selected_error = error[mask]
    total_volume = float(np.sum(selected_volume))
    return {
        "L1": float(np.sum(np.abs(selected_error) * selected_volume) / total_volume),
        "L2": float(
            np.sqrt(
                np.sum(selected_error * selected_error * selected_volume)
                / total_volume
            )
        ),
        "Linf": float(np.max(np.abs(selected_error))),
    }


def add_prefixed_metrics(row, prefix, metrics):
    for name, value in metrics.items():
        row[f"{prefix}_{name}"] = value


def analyze_case(path, quantity, physics, exact_time_override, args):
    plotfile = uniform_analysis.latest_plotfile(path)
    dataset = uniform_analysis.load_boxlib_rz(plotfile)
    if int(dataset.index.max_level) != 1:
        raise ValueError(
            f"{plotfile} must contain exactly one refined AMR level"
        )

    data = dataset.all_data()
    radial_coordinate = np.asarray(data[("index", "r")].d, dtype=float)
    axial_coordinate = np.asarray(data[("index", "z")].d, dtype=float)
    radial_spacing = np.asarray(data[("index", "dr")].d, dtype=float)
    cell_volume = np.asarray(data[("index", "cell_volume")].d, dtype=float)
    grid_level = np.rint(
        np.asarray(data[("index", "grid_level")].d, dtype=float)
    ).astype(int)

    numerical = {
        name: np.asarray(
            data[uniform_analysis.field_key(dataset, plot_name)].d,
            dtype=float,
        )
        for name, plot_name in VARIABLE_FIELDS.items()
    }

    plot_time = float(dataset.current_time)
    exact_time = (
        plot_time if exact_time_override is None else exact_time_override
    )
    exact = (
        uniform_analysis.exact_state(radial_coordinate, axial_coordinate, exact_time)
        if quantity == "state"
        else uniform_analysis.exact_time_derivative(
            radial_coordinate, axial_coordinate, exact_time
        )
    )

    radial_lo = float(dataset.domain_left_edge[0])
    radial_hi = float(dataset.domain_right_edge[0])
    if not (
        radial_lo < args.refine_radial_lo < args.refine_radial_hi < radial_hi
    ):
        raise ValueError(
            "the requested radial refinement band must lie strictly inside "
            f"the domain [{radial_lo}, {radial_hi}]"
        )

    distance_to_interface = np.minimum(
        np.abs(radial_coordinate - args.refine_radial_lo),
        np.abs(radial_coordinate - args.refine_radial_hi),
    )
    interface_mask = (
        distance_to_interface
        <= args.interface_width_cells * radial_spacing
    )
    maximum_level = int(np.max(grid_level))
    coarse_interior_mask = (grid_level == 0) & ~interface_mask
    fine_interior_mask = (grid_level == maximum_level) & ~interface_mask

    expected_first_ring_centre = radial_lo + 0.5 * radial_spacing
    first_ring_tolerance = np.maximum(
        64.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(radial_coordinate)),
        1.0e-10 * radial_spacing,
    )
    first_ring_mask = (
        np.abs(radial_coordinate - expected_first_ring_centre)
        <= first_ring_tolerance
    )

    if not np.any(interface_mask):
        raise ValueError(f"no leaf cells selected near the interfaces in {plotfile}")
    if not np.any(coarse_interior_mask) or not np.any(fine_interior_mask):
        raise ValueError(
            f"{plotfile} does not contain both coarse and fine cells away from "
            "the requested interface band"
        )
    if not np.any(first_ring_mask):
        raise ValueError(f"no first-axis-ring leaf cells found in {plotfile}")

    base_dimensions = np.asarray(dataset.domain_dimensions, dtype=int)
    row = {
        "physics": physics,
        "quantity": quantity,
        "state_semantics": "point_sample_amr_leaf",
        "plotfile": str(plotfile),
        "base_n": int(base_dimensions[0]),
        "base_nz": int(base_dimensions[1]),
        "max_level": maximum_level,
        "active_leaf_cells": int(radial_coordinate.size),
        "interface_leaf_cells": int(np.count_nonzero(interface_mask)),
        "coarse_interior_leaf_cells": int(
            np.count_nonzero(coarse_interior_mask)
        ),
        "fine_interior_leaf_cells": int(np.count_nonzero(fine_interior_mask)),
        "plot_time": plot_time,
        "exact_time": exact_time,
        "refine_radial_lo": args.refine_radial_lo,
        "refine_radial_hi": args.refine_radial_hi,
        "interface_width_cells": args.interface_width_cells,
    }

    all_leaf_cells = np.ones(radial_coordinate.shape, dtype=bool)
    for variable in VARIABLE_FIELDS:
        error = numerical[variable] - exact[variable]
        add_prefixed_metrics(
            row,
            f"{variable}_global",
            masked_error_metrics(error, cell_volume, all_leaf_cells),
        )
        add_prefixed_metrics(
            row,
            f"{variable}_interface",
            masked_error_metrics(error, cell_volume, interface_mask),
        )
        add_prefixed_metrics(
            row,
            f"{variable}_coarse_interior",
            masked_error_metrics(error, cell_volume, coarse_interior_mask),
        )
        add_prefixed_metrics(
            row,
            f"{variable}_fine_interior",
            masked_error_metrics(error, cell_volume, fine_interior_mask),
        )
        row[f"{variable}_first_ring_Linf"] = float(
            np.max(np.abs(error[first_ring_mask]))
        )

        numerical_integral = float(np.sum(numerical[variable] * cell_volume))
        exact_integral = float(np.sum(exact[variable] * cell_volume))
        integral_scale = float(np.sum(np.abs(exact[variable]) * cell_volume))
        integral_error = abs(numerical_integral - exact_integral)
        row[f"{variable}_numerical_integral"] = numerical_integral
        row[f"{variable}_exact_integral"] = exact_integral
        row[f"{variable}_integral_abs_error"] = integral_error
        row[f"{variable}_integral_relative_error"] = (
            integral_error / integral_scale
            if integral_scale > np.finfo(float).tiny
            else float("nan")
        )

    return row


def add_convergence_rates(rows):
    rate_metrics = []
    for variable in VARIABLE_FIELDS:
        rate_metrics.extend(
            (
                f"{variable}_global_L1",
                f"{variable}_global_L2",
                f"{variable}_interface_L1",
                f"{variable}_interface_L2",
                f"{variable}_coarse_interior_L2",
                f"{variable}_fine_interior_L2",
                f"{variable}_integral_relative_error",
            )
        )

    for row in rows:
        for metric in rate_metrics:
            row[f"{metric}_rate"] = float("nan")

    for coarse, fine in zip(rows, rows[1:]):
        refinement = fine["base_n"] / coarse["base_n"]
        for metric in rate_metrics:
            coarse_error = coarse[metric]
            fine_error = fine[metric]
            if (
                math.isfinite(coarse_error)
                and math.isfinite(fine_error)
                and coarse_error > 0.0
                and fine_error > 0.0
            ):
                fine[f"{metric}_rate"] = math.log(
                    coarse_error / fine_error
                ) / math.log(refinement)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument(
        "--quantity", choices=("state", "rhs"), default="state"
    )
    parser.add_argument(
        "--physics", choices=("euler", "navier-stokes"), default="euler"
    )
    parser.add_argument("--exact-time", type=float)
    parser.add_argument("--refine-radial-lo", type=float, default=0.25)
    parser.add_argument("--refine-radial-hi", type=float, default=0.75)
    parser.add_argument("--interface-width-cells", type=float, default=3.0)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    if not args.interface_width_cells > 0.0:
        raise ValueError("--interface-width-cells must be positive")

    rows = sorted(
        (
            analyze_case(
                case,
                args.quantity,
                args.physics,
                args.exact_time,
                args,
            )
            for case in args.cases
        ),
        key=lambda row: row["base_n"],
    )
    if any(row["base_n"] != row["base_nz"] for row in rows):
        raise ValueError("MMS convergence requires square base grids")
    if len({row["base_n"] for row in rows}) != len(rows):
        raise ValueError("the AMR convergence sequence contains duplicate grids")
    if len({row["max_level"] for row in rows}) != 1:
        raise ValueError("all AMR cases must use the same maximum level")

    add_convergence_rates(rows)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{args.physics:14s} {args.quantity:5s} "
            f"N0={row['base_n']:4d} L={row['max_level']} "
            f"rho_global_L2={row['density_global_L2']:.6e} "
            f"rho_interface_L2={row['density_interface_L2']:.6e} "
            f"mr_global_L2={row['radial_momentum_global_L2']:.6e} "
            f"mr_interface_L2={row['radial_momentum_interface_L2']:.6e}"
        )
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    uniform_analysis.yt.set_log_level(40)
    main()
