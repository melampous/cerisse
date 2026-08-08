#!/usr/bin/env python3
"""Measure radial-uniformity errors in an all-fluid R-Z plotfile series."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yt
import yt.frontends.boxlib.data_structures as boxlib_data_structures


# yt 4.3 plus NumPy 2 can expose a read-only padded domain edge while loading
# a 2-D cylindrical plotfile.  Treating the plot as Cartesian is sufficient for
# reading these cell-centered fields and avoids changing the installed package.
_geometry = boxlib_data_structures.Geometry
boxlib_data_structures.Geometry = lambda name: _geometry(
    "cartesian" if name == "cylindrical" else name
)


def field(grid, name: str) -> np.ndarray:
    return np.asarray(grid[("boxlib", name)].d, dtype=np.float64).squeeze(axis=2)


def analyze_plot(plotfile: Path, gamma: float) -> dict[str, float | int | str]:
    dataset = yt.load(str(plotfile))
    grid = dataset.covering_grid(
        level=0,
        left_edge=dataset.domain_left_edge,
        dims=dataset.domain_dimensions,
    )
    pressure = field(grid, "pressure")
    density = field(grid, "density")
    radial_velocity = field(grid, "x_velocity")
    sound_speed = np.sqrt(gamma * pressure / density)

    radial_spread = np.ptp(pressure, axis=0)
    local_scale = np.maximum(np.max(np.abs(pressure), axis=0), np.finfo(float).tiny)
    axis_ring_delta = np.abs(pressure[0, :] - pressure[1, :])
    axis_mean_delta = np.abs(pressure[0, :] - np.mean(pressure, axis=0))

    dimensions = np.asarray(dataset.domain_dimensions, dtype=int)
    left = np.asarray(dataset.domain_left_edge.d, dtype=float)
    right = np.asarray(dataset.domain_right_edge.d, dtype=float)
    return {
        "plotfile": str(plotfile),
        "step": int(plotfile.name.removeprefix("plt")),
        "time": float(dataset.current_time),
        "nr": int(dimensions[0]),
        "nz": int(dimensions[1]),
        "dr": float((right[0] - left[0]) / dimensions[0]),
        "dz": float((right[1] - left[1]) / dimensions[1]),
        "pressure_min": float(np.min(pressure)),
        "pressure_max": float(np.max(pressure)),
        "max_radial_pressure_spread_abs": float(np.max(radial_spread)),
        "max_radial_pressure_spread_rel_local": float(
            np.max(radial_spread / local_scale)
        ),
        "max_axis_to_second_ring_pressure_delta_abs": float(
            np.max(axis_ring_delta)
        ),
        "max_axis_to_second_ring_pressure_delta_rel_local": float(
            np.max(axis_ring_delta / local_scale)
        ),
        "max_axis_to_radial_mean_pressure_delta_abs": float(
            np.max(axis_mean_delta)
        ),
        "max_axis_to_radial_mean_pressure_delta_rel_local": float(
            np.max(axis_mean_delta / local_scale)
        ),
        "max_abs_radial_velocity": float(np.max(np.abs(radial_velocity))),
        "max_abs_radial_velocity_over_sound": float(
            np.max(np.abs(radial_velocity) / sound_speed)
        ),
        "max_axis_abs_radial_velocity_over_sound": float(
            np.max(np.abs(radial_velocity[0, :]) / sound_speed[0, :])
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    plots = sorted(
        path for path in args.case_directory.glob("plt[0-9]*")
        if (path / "Header").is_file()
    )
    if not plots:
        raise SystemExit(f"no plotfiles found in {args.case_directory}")
    records = [analyze_plot(plot, args.gamma) for plot in plots]
    keys = (
        "max_radial_pressure_spread_abs",
        "max_radial_pressure_spread_rel_local",
        "max_axis_to_second_ring_pressure_delta_abs",
        "max_axis_to_second_ring_pressure_delta_rel_local",
        "max_axis_to_radial_mean_pressure_delta_abs",
        "max_axis_to_radial_mean_pressure_delta_rel_local",
        "max_abs_radial_velocity",
        "max_abs_radial_velocity_over_sound",
        "max_axis_abs_radial_velocity_over_sound",
    )
    result = {
        "case_directory": str(args.case_directory),
        "number_of_plotfiles": len(records),
        "history_maxima": {key: max(float(row[key]) for row in records) for key in keys},
        "final": records[-1],
        "plotfiles": records,
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
