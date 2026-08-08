#!/usr/bin/env python3
"""Quantify near-axis regularity for the short nonlinear R-Z shock gate.

The first plotted radial row is the annular cell centred at r=dr/2, not a
point value at r=0.  Consequently this script reports signed one-sided
regularity proxies instead of incorrectly demanding that the first-ring
radial velocity itself vanish.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yt
import yt.frontends.boxlib.data_structures as boxlib_data_structures


# yt 4.3 plus NumPy 2 can expose a read-only padded edge for 2-D cylindrical
# plotfiles.  Cartesian geometry is sufficient for reading cell-centred data.
_geometry = boxlib_data_structures.Geometry
boxlib_data_structures.Geometry = lambda name: _geometry(
    "cartesian" if name == "cylindrical" else name
)


def field(grid, name: str) -> np.ndarray:
    return np.asarray(grid[("boxlib", name)].d, dtype=np.float64).squeeze(axis=2)


def signed_extrema(values: np.ndarray) -> dict[str, float]:
    return {"min": float(np.min(values)), "max": float(np.max(values))}


def analyze_plot(plotfile: Path, gamma: float, rings: int) -> dict:
    dataset = yt.load(str(plotfile))
    grid = dataset.covering_grid(
        level=0,
        left_edge=dataset.domain_left_edge,
        dims=dataset.domain_dimensions,
    )
    pressure = field(grid, "pressure")
    density = field(grid, "density")
    radial_velocity = field(grid, "x_velocity")
    axial_velocity = field(grid, "y_velocity")
    sound = np.sqrt(gamma * pressure / density)
    radial_mach = radial_velocity / sound

    dimensions = np.asarray(dataset.domain_dimensions, dtype=int)
    left = np.asarray(dataset.domain_left_edge.d, dtype=float)
    right = np.asarray(dataset.domain_right_edge.d, dtype=float)
    dr = float((right[0] - left[0]) / dimensions[0])
    dz = float((right[1] - left[1]) / dimensions[1])
    nrings = min(rings, pressure.shape[0])

    p_scale = np.maximum(np.max(np.abs(pressure), axis=0), np.finfo(float).tiny)
    rho_scale = np.maximum(np.max(np.abs(density), axis=0), np.finfo(float).tiny)
    p_axis_second = (pressure[0] - pressure[1]) / p_scale
    rho_axis_second = (density[0] - density[1]) / rho_scale

    # For a smooth odd u_r=c_1 r+c_3 r^3+..., the linear extrapolation from
    # r=dr/2 and 3dr/2 to r=0 is (3*u_0-u_1)/2 = O(dr^3).  It is a signed
    # near-axis regularity proxy; it is not a sampled value at the axis.
    axis_sound_scale = np.maximum(sound[0], np.finfo(float).tiny)
    odd_axis_proxy = (3.0 * radial_velocity[0] - radial_velocity[1]) / (
        2.0 * axis_sound_scale
    )

    # Track the strongest axial pressure jump in each radial ring.  A truly
    # planar front has the same face index for all rings.
    dp = np.abs(np.diff(pressure, axis=1))
    shock_face_index = np.argmax(dp, axis=1)
    shock_z = left[1] + (shock_face_index.astype(float) + 1.0) * dz
    shock_offset_cells = (shock_z - np.median(shock_z)) / dz

    ring_records = []
    for ring in range(nrings):
        p_corr = (pressure[ring] - np.mean(pressure, axis=0)) / p_scale
        rho_corr = (density[ring] - np.mean(density, axis=0)) / rho_scale
        ring_records.append(
            {
                "ring": ring,
                "r_over_dr": ring + 0.5,
                "pressure": signed_extrema(pressure[ring]),
                "density": signed_extrema(density[ring]),
                "radial_velocity_over_sound": signed_extrema(radial_mach[ring]),
                "axial_velocity": signed_extrema(axial_velocity[ring]),
                "signed_pressure_corrugation_rel": signed_extrema(p_corr),
                "signed_density_corrugation_rel": signed_extrema(rho_corr),
                "shock_face_index": int(shock_face_index[ring]),
                "shock_z": float(shock_z[ring]),
                "shock_offset_cells_from_radial_median": float(
                    shock_offset_cells[ring]
                ),
            }
        )

    return {
        "plotfile": str(plotfile),
        "step": int(plotfile.name.removeprefix("plt")),
        "time": float(dataset.current_time),
        "nr": int(dimensions[0]),
        "nz": int(dimensions[1]),
        "dr": dr,
        "dz": dz,
        "axis_note": "ring 0 is centred at r=dr/2; no point sample exists at r=0",
        "signed_axis_to_second_ring_pressure_rel": signed_extrema(p_axis_second),
        "signed_axis_to_second_ring_density_rel": signed_extrema(rho_axis_second),
        "signed_odd_axis_extrapolation_ur_over_sound": signed_extrema(
            odd_axis_proxy
        ),
        "max_radial_pressure_spread_rel_local": float(
            np.max(np.ptp(pressure, axis=0) / p_scale)
        ),
        "max_radial_density_spread_rel_local": float(
            np.max(np.ptp(density, axis=0) / rho_scale)
        ),
        "max_abs_radial_velocity_over_sound": float(np.max(np.abs(radial_mach))),
        "max_axis_abs_radial_velocity_over_sound": float(
            np.max(np.abs(radial_mach[0]))
        ),
        "shock_front_corrugation_cells": float(np.ptp(shock_z) / dz),
        "first_rings": ring_records,
    }


def max_abs_signed(record: dict[str, float]) -> float:
    return max(abs(float(record["min"])), abs(float(record["max"])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--rings", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    plots = sorted(
        path
        for path in args.case_directory.glob("plt[0-9]*")
        if (path / "Header").is_file()
    )
    if not plots:
        raise SystemExit(f"no plotfiles found in {args.case_directory}")
    records = [analyze_plot(path, args.gamma, args.rings) for path in plots]

    scalar_keys = (
        "max_radial_pressure_spread_rel_local",
        "max_radial_density_spread_rel_local",
        "max_abs_radial_velocity_over_sound",
        "max_axis_abs_radial_velocity_over_sound",
        "shock_front_corrugation_cells",
    )
    signed_keys = (
        "signed_axis_to_second_ring_pressure_rel",
        "signed_axis_to_second_ring_density_rel",
        "signed_odd_axis_extrapolation_ur_over_sound",
    )
    history = {key: max(float(row[key]) for row in records) for key in scalar_keys}
    history.update(
        {f"max_abs_{key}": max(max_abs_signed(row[key]) for row in records)
         for key in signed_keys}
    )
    result = {
        "case_directory": str(args.case_directory),
        "number_of_plotfiles": len(records),
        "history_maxima": history,
        "final": records[-1],
        "plotfiles": records,
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
