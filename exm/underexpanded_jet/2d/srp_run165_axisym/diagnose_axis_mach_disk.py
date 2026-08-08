#!/usr/bin/env python3
"""Quantify the instantaneous and running-mean Run-165 Mach disk near r=0."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yt


def load_boxlib_rz(plotfile: Path):
    """Load old RZ plotfiles with yt versions that mark geometry attrs read-only."""
    from yt.data_objects.static_output import MutableAttribute

    original_set = MutableAttribute.__set__

    def writable_set(self, instance, value):
        self.data[instance] = value

    MutableAttribute.__set__ = writable_set
    try:
        return yt.load(str(plotfile))
    finally:
        MutableAttribute.__set__ = original_set


def field(grid, name: str) -> np.ndarray:
    for namespace in ("boxlib", "gas", "stream", "index"):
        try:
            values = np.asarray(grid[(namespace, name)].d, dtype=float)
            while values.ndim > 2:
                values = values[..., 0]
            return values
        except Exception:
            pass
    raise KeyError(name)


def local_peak_indices(signal: np.ndarray, spacing: float, count: int = 8) -> list[int]:
    order = np.argsort(signal)[::-1]
    separation = max(1, int(round(1.0e-3 / spacing)))
    selected: list[int] = []
    for candidate in order:
        index = int(candidate)
        if all(abs(index - other) >= separation for other in selected):
            selected.append(index)
        if len(selected) == count:
            break
    return selected


def shock_curve(
    pressure: np.ndarray,
    z: np.ndarray,
    z_bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = (z >= z_bounds[0]) & (z <= z_bounds[1])
    window_indices = np.flatnonzero(mask)
    if window_indices.size == 0:
        raise ValueError(f"empty shock window {z_bounds}")
    gradient = np.abs(
        np.gradient(np.log(np.maximum(pressure, 1.0e-30)), z, axis=1)
    )
    local = np.argmax(gradient[:, mask], axis=1)
    indices = window_indices[local]
    radial = np.arange(pressure.shape[0])
    return z[indices], gradient[radial, indices], indices


def second_moment_metrics(mean: np.ndarray, second: np.ndarray) -> dict[str, float]:
    variance = second - mean * mean
    scale = np.maximum(second, mean * mean)
    relative = variance / np.maximum(scale, 1.0e-300)
    return {
        "negative_fraction": float(np.mean(variance < 0.0)),
        "relative_min": float(np.nanmin(relative)),
        "relative_p001": float(np.nanpercentile(relative, 0.1)),
        "relative_p50": float(np.nanpercentile(relative, 50.0)),
        "relative_max": float(np.nanmax(relative)),
    }


def analyze(
    plotfile: Path,
    r_max: float,
    z_bounds: tuple[float, float],
) -> dict[str, object]:
    dataset = load_boxlib_rz(plotfile)
    level = int(dataset.index.max_level)
    base_dx = (
        np.asarray(dataset.domain_right_edge.d)
        - np.asarray(dataset.domain_left_edge.d)
    ) / np.asarray(dataset.domain_dimensions)
    dr = float(base_dx[0] / dataset.refine_by**level)
    dz = float(base_dx[1] / dataset.refine_by**level)
    nr = min(
        int(np.ceil(r_max / dr)),
        int(dataset.domain_dimensions[0] * dataset.refine_by**level),
    )
    nz = int(dataset.domain_dimensions[1] * dataset.refine_by**level)
    grid = dataset.covering_grid(
        level,
        dataset.domain_left_edge,
        [nr, nz, 1],
    )

    density = field(grid, "Density")
    pressure = field(grid, "pressure")
    radial_velocity = field(grid, "Xmom") / density
    axial_velocity = field(grid, "Ymom") / density
    pressure_mean = field(grid, "pressureMEAN")
    radial_velocity_mean = field(grid, "x_velocityMEAN")
    axial_velocity_mean = field(grid, "y_velocityMEAN")
    grid_level = field(grid, "grid_level")

    r = float(dataset.domain_left_edge[0]) + (np.arange(nr) + 0.5) * dr
    z = float(dataset.domain_left_edge[1]) + (np.arange(nz) + 0.5) * dz
    instant_z, instant_strength, instant_index = shock_curve(
        pressure, z, z_bounds
    )
    mean_z, mean_strength, mean_index = shock_curve(
        pressure_mean, z, z_bounds
    )

    axis_gradient = np.abs(
        np.gradient(np.log(np.maximum(pressure[0], 1.0e-30)), z)
    )
    axis_mean_gradient = np.abs(
        np.gradient(np.log(np.maximum(pressure_mean[0], 1.0e-30)), z)
    )
    broad_mask = (z >= 0.005) & (z <= 0.155)
    broad_indices = np.flatnonzero(broad_mask)
    instant_candidates = local_peak_indices(axis_gradient[broad_mask], dz)
    mean_candidates = local_peak_indices(axis_mean_gradient[broad_mask], dz)

    first_count = min(32, nr)
    first_rows = []
    for i in range(first_count):
        ji = int(instant_index[i])
        jm = int(mean_index[i])
        first_rows.append(
            {
                "i": i,
                "r_mm": float(r[i] * 1.0e3),
                "instant_z_mm": float(instant_z[i] * 1.0e3),
                "mean_z_mm": float(mean_z[i] * 1.0e3),
                "instant_grad_logp": float(instant_strength[i]),
                "mean_grad_logp": float(mean_strength[i]),
                "instant_source_level": int(round(grid_level[i, ji])),
                "mean_source_level": int(round(grid_level[i, jm])),
                "ur_at_instant_m_s": float(radial_velocity[i, ji]),
                "ur_mean_at_mean_m_s": float(radial_velocity_mean[i, jm]),
                "uz_at_instant_m_s": float(axial_velocity[i, ji]),
                "uz_mean_at_mean_m_s": float(axial_velocity_mean[i, jm]),
            }
        )

    return {
        "plotfile": str(plotfile),
        "time_s": float(dataset.current_time),
        "finest_level": level,
        "dr_mm": dr * 1.0e3,
        "dz_mm": dz * 1.0e3,
        "shock_window_mm": [z_bounds[0] * 1.0e3, z_bounds[1] * 1.0e3],
        "axis_instant_peak_candidates": [
            {
                "z_mm": float(z[broad_indices[index]] * 1.0e3),
                "grad_logp": float(axis_gradient[broad_indices[index]]),
            }
            for index in instant_candidates
        ],
        "axis_mean_peak_candidates": [
            {
                "z_mm": float(z[broad_indices[index]] * 1.0e3),
                "grad_logp": float(axis_mean_gradient[broad_indices[index]]),
            }
            for index in mean_candidates
        ],
        "instant_shape_first32": {
            "z_min_mm": float(np.min(instant_z[:first_count]) * 1.0e3),
            "z_max_mm": float(np.max(instant_z[:first_count]) * 1.0e3),
            "range_cells": float(np.ptp(instant_z[:first_count]) / dz),
        },
        "mean_shape_first32": {
            "z_min_mm": float(np.min(mean_z[:first_count]) * 1.0e3),
            "z_max_mm": float(np.max(mean_z[:first_count]) * 1.0e3),
            "range_cells": float(np.ptp(mean_z[:first_count]) / dz),
        },
        "second_moments": {
            "pressure": second_moment_metrics(
                pressure_mean, field(grid, "pressureSQR")
            ),
            "density": second_moment_metrics(
                field(grid, "DensityMEAN"), field(grid, "DensitySQR")
            ),
            "radial_velocity": second_moment_metrics(
                radial_velocity_mean, field(grid, "x_velocitySQR")
            ),
            "axial_velocity": second_moment_metrics(
                axial_velocity_mean, field(grid, "y_velocitySQR")
            ),
        },
        "first_radial_cells": first_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfiles", nargs="+", type=Path)
    parser.add_argument("--r-max", type=float, default=0.02)
    parser.add_argument("--shock-z-min", type=float, default=0.09)
    parser.add_argument("--shock-z-max", type=float, default=0.135)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    reports = [
        analyze(
            plotfile,
            args.r_max,
            (args.shock_z_min, args.shock_z_max),
        )
        for plotfile in args.plotfiles
    ]
    text = json.dumps(reports, indent=2)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
