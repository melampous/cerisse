#!/usr/bin/env python3
"""Measure odd-even growth and admissibility in Cartesian shock plotfiles."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import yt


def joint_sensor_metrics(
    density: np.ndarray,
    pressure: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    dx: float,
    dy: float,
    gamma: float,
    smoothness_threshold: float,
    pressure_jump_threshold: float,
    compression_threshold: float,
) -> dict[str, int | float]:
    """Mirror the bulk Cartesian AFD shock fallback sensor.

    A face falls back to high-order LLF only when the broad six-cell
    smoothness test, the pressure-jump gate, and the multidimensional
    compression gate all trigger.  Periodic indexing is used in y.  Faces
    whose diagnostic stencil would cross a physical x boundary are omitted.
    """

    nx, ny = density.shape
    tiny = np.finfo(float).tiny
    sound = np.sqrt(gamma * pressure / density)
    zero_velocity = np.zeros_like(ux)
    fields = (density, pressure, ux, uy, zero_velocity)
    minimum_cell_size = min(dx, dy)

    def value(values: np.ndarray, i: int, j: int) -> float:
        return float(values[i, j % ny])

    def offset(i: int, j: int, direction: int, amount: int) -> tuple[int, int]:
        if direction == 0:
            return i + amount, j
        return i, (j + amount) % ny

    def classify_face(i: int, j: int, direction: int) -> tuple[
        bool, bool, bool, float, float
    ]:
        stencil_coordinates = [
            offset(i, j, direction, m - 3) for m in range(6)
        ]
        sound_scale = max(
            abs(value(sound, ii, jj)) for ii, jj in stencil_coordinates
        )

        nonsmooth = False
        for field_index, values in enumerate(fields):
            samples = [
                value(values, ii, jj) for ii, jj in stencil_coordinates
            ]
            field_scale = max(abs(sample) for sample in samples)
            floor = sound_scale if field_index >= 2 else field_scale
            for m in range(4):
                a, b, c = samples[m:m + 3]
                denominator = (
                    abs(a) + 2.0 * abs(b) + abs(c) + floor + tiny
                )
                if abs(a - 2.0 * b + c) / denominator > smoothness_threshold:
                    nonsmooth = True
                    break
            if nonsmooth:
                break

        maximum_pressure_jump = 0.0
        for left, right in zip(
            stencil_coordinates[:-1], stencil_coordinates[1:]
        ):
            pressure_left = value(pressure, *left)
            pressure_right = value(pressure, *right)
            pressure_scale = max(pressure_left, pressure_right)
            relative_jump = (
                abs(pressure_right - pressure_left) /
                max(pressure_scale, tiny)
            )
            maximum_pressure_jump = max(
                maximum_pressure_jump, relative_jump
            )
        pressure_gate = maximum_pressure_jump > pressure_jump_threshold

        maximum_compression = 0.0
        for m in range(1, 5):
            center_i, center_j = offset(i, j, direction, m - 3)
            divergence = (
                (
                    value(ux, center_i + 1, center_j) -
                    value(ux, center_i - 1, center_j)
                ) / (2.0 * dx) +
                (
                    value(uy, center_i, center_j + 1) -
                    value(uy, center_i, center_j - 1)
                ) / (2.0 * dy)
            )
            local_sound_scale = value(sound, center_i, center_j)
            for neighbour_i, neighbour_j in (
                (center_i - 1, center_j),
                (center_i + 1, center_j),
                (center_i, center_j - 1),
                (center_i, center_j + 1),
            ):
                local_sound_scale = max(
                    local_sound_scale,
                    value(sound, neighbour_i, neighbour_j),
                )
            normalized_compression = (
                max(-divergence, 0.0) * minimum_cell_size /
                max(local_sound_scale, tiny)
            )
            maximum_compression = max(
                maximum_compression, normalized_compression
            )
        compression_gate = maximum_compression > compression_threshold

        return (
            nonsmooth,
            pressure_gate,
            compression_gate,
            maximum_pressure_jump,
            maximum_compression,
        )

    metrics: dict[str, int | float] = {
        "sensor_smoothness_threshold": smoothness_threshold,
        "sensor_pressure_jump_threshold": pressure_jump_threshold,
        "sensor_compression_threshold": compression_threshold,
    }
    for direction, name in ((0, "x"), (1, "y")):
        faces = (
            ((i, j) for i in range(3, nx - 2) for j in range(ny))
            if direction == 0 else
            ((i, j) for i in range(1, nx - 1) for j in range(ny))
        )
        counts = {
            "considered": 0,
            "nonsmooth": 0,
            "pressure": 0,
            "compression": 0,
            "fallback": 0,
        }
        fallback_rows: set[int] = set()
        fallback_columns: set[int] = set()
        maximum_pressure_jump = 0.0
        maximum_compression = 0.0
        for i, j in faces:
            (
                nonsmooth,
                pressure_gate,
                compression_gate,
                local_pressure_jump,
                local_compression,
            ) = classify_face(i, j, direction)
            fallback = nonsmooth and pressure_gate and compression_gate
            counts["considered"] += 1
            counts["nonsmooth"] += int(nonsmooth)
            counts["pressure"] += int(pressure_gate)
            counts["compression"] += int(compression_gate)
            counts["fallback"] += int(fallback)
            maximum_pressure_jump = max(
                maximum_pressure_jump, local_pressure_jump
            )
            maximum_compression = max(
                maximum_compression, local_compression
            )
            if fallback:
                fallback_rows.add(j)
                fallback_columns.add(i)

        prefix = f"sensor_{name}"
        for label, count in counts.items():
            metrics[f"{prefix}_{label}_faces"] = count
        metrics[f"{prefix}_fallback_rows"] = len(fallback_rows)
        metrics[f"{prefix}_fallback_columns"] = len(fallback_columns)
        metrics[f"{prefix}_maximum_pressure_jump"] = maximum_pressure_jump
        metrics[f"{prefix}_maximum_compression"] = maximum_compression

    return metrics


def field(grid, name: str) -> np.ndarray:
    for candidate in (("boxlib", name), ("gas", name), ("stream", name)):
        try:
            return np.squeeze(np.asarray(grid[candidate], dtype=float))
        except Exception:
            pass
    raise KeyError(f"plotfile does not contain {name!r}")


def orient(values: np.ndarray, nx: int, ny: int) -> np.ndarray:
    values = np.squeeze(values)
    if values.shape == (nx, ny):
        return values
    if values.shape == (ny, nx):
        return values.T
    raise ValueError(f"unexpected field shape {values.shape}, expected {(nx, ny)}")


def shock_locus(pressure: np.ndarray, x: np.ndarray, midpoint: float,
                shock_x: float) -> np.ndarray:
    locus = np.full(pressure.shape[1], np.nan)
    for j, profile in enumerate(pressure.T):
        offset = profile - midpoint
        crossings = np.flatnonzero(offset[:-1] * offset[1:] <= 0.0)
        if crossings.size == 0:
            continue
        crossing_x = 0.5 * (x[crossings] + x[crossings + 1])
        i = int(crossings[np.argmin(np.abs(crossing_x - shock_x))])
        p0 = profile[i]
        p1 = profile[i + 1]
        fraction = 0.5 if p1 == p0 else (midpoint - p0) / (p1 - p0)
        locus[j] = x[i] + fraction * (x[i + 1] - x[i])
    return locus


def read_plotfile(
    path: Path,
    mach: float,
    shock_x: float,
    gamma: float,
    smoothness_threshold: float,
    pressure_jump_threshold: float,
    compression_threshold: float,
) -> dict[str, object]:
    yt.funcs.mylog.setLevel(50)
    ds = yt.load(str(path))
    if int(ds.index.max_level) != 0:
        raise ValueError("this uniform-grid diagnostic requires max_level=0")
    dims = np.asarray(ds.domain_dimensions, dtype=int)
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    nx, ny = int(dims[0]), int(dims[1])

    density = orient(field(grid, "Density"), nx, ny)
    pressure = orient(field(grid, "pressure"), nx, ny)
    try:
        ux = orient(field(grid, "x_velocity"), nx, ny)
        uy = orient(field(grid, "y_velocity"), nx, ny)
    except KeyError:
        ux = orient(field(grid, "Xmom"), nx, ny) / density
        uy = orient(field(grid, "Ymom"), nx, ny) / density

    lo = np.asarray(ds.domain_left_edge, dtype=float)
    hi = np.asarray(ds.domain_right_edge, dtype=float)
    dx = (hi[0] - lo[0]) / nx
    dy = (hi[1] - lo[1]) / ny
    x = lo[0] + (np.arange(nx) + 0.5) * dx
    pressure_ratio = (
        1.0 + 2.0 * gamma / (gamma + 1.0) * (mach**2 - 1.0)
    )
    midpoint = 0.5 * (1.0 + pressure_ratio)
    locus = shock_locus(pressure, x, midpoint, shock_x)
    if not np.all(np.isfinite(locus)):
        raise ValueError(
            f"shock crossing missing in {np.count_nonzero(~np.isfinite(locus))} rows"
        )

    displacement_cells = (locus - np.mean(locus)) / dx
    parity = np.where((np.arange(ny) & 1) == 0, -1.0, 1.0)
    odd_even = float(np.mean(parity * displacement_cells))
    window = np.abs(x - shock_x) <= 0.08
    uy_window = uy[window, :]
    density_window = density[window, :]
    pressure_window = pressure[window, :]

    metrics = {
        "plotfile": str(path),
        "time": float(ds.current_time),
        "nx": nx,
        "ny": ny,
        "dx": float(dx),
        "density_min_global": float(np.min(density)),
        "pressure_min_global": float(np.min(pressure)),
        "density_min_shock_window": float(np.min(density_window)),
        "pressure_min_shock_window": float(np.min(pressure_window)),
        "shock_mean": float(np.mean(locus)),
        "shock_mean_offset_cells": float((np.mean(locus) - shock_x) / dx),
        "shock_peak_to_peak_cells": float(np.ptp(locus) / dx),
        "shock_rms_cells": float(np.sqrt(np.mean(displacement_cells**2))),
        "shock_odd_even_amplitude_cells": odd_even,
        "shock_odd_even_abs_cells": abs(odd_even),
        "uy_absmax_global": float(np.max(np.abs(uy))),
        "uy_absmax_shock_window": float(np.max(np.abs(uy_window))),
        "uy_rms_shock_window": float(np.sqrt(np.mean(uy_window**2))),
        "ux_transverse_range_max": float(np.max(np.ptp(ux, axis=1))),
        "pressure_transverse_range_max": float(
            np.max(np.ptp(pressure, axis=1))
        ),
        "shock_locus": [float(value) for value in locus],
    }
    metrics.update(
        joint_sensor_metrics(
            density,
            pressure,
            ux,
            uy,
            float(dx),
            float(dy),
            gamma,
            smoothness_threshold,
            pressure_jump_threshold,
            compression_threshold,
        )
    )
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfiles", nargs="+", type=Path)
    parser.add_argument("--mach", type=float, default=10.0)
    parser.add_argument("--shock-x", type=float, default=0.5)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--smoothness-threshold", type=float, default=0.08)
    parser.add_argument("--pressure-jump-threshold", type=float, default=0.01)
    parser.add_argument("--compression-threshold", type=float, default=0.001)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()

    rows = [
        read_plotfile(
            path,
            args.mach,
            args.shock_x,
            args.gamma,
            args.smoothness_threshold,
            args.pressure_jump_threshold,
            args.compression_threshold,
        )
        for path in args.plotfiles
    ]
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows, indent=2) + "\n")
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        scalar_keys = [key for key in rows[0] if key != "shock_locus"]
        with args.csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=scalar_keys)
            writer.writeheader()
            writer.writerows(
                {key: row[key] for key in scalar_keys} for row in rows
            )
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
