#!/usr/bin/env python3
"""Audit Mach-disk shape, motion, and AMR coverage in the NPR=15 R-Z case.

This is a read-only plotfile diagnostic. It deliberately tracks the first
compression in a bounded axial window instead of taking the largest gradient
anywhere downstream, where a recirculation-generated compression can be
stronger than the Mach disk.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt


GAMMA = 1.4
RGAS = 287.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--diameter", type=float, default=0.0254)
    parser.add_argument("--r-max", type=float, default=0.8,
                        help="extracted radial extent in nozzle diameters")
    parser.add_argument("--extract-z-min", type=float, default=1.8)
    parser.add_argument("--extract-z-max", type=float, default=3.4)
    parser.add_argument("--track-z-min", type=float, default=2.0)
    parser.add_argument("--track-z-max", type=float, default=2.65)
    parser.add_argument("--shape-r-max", type=float, default=0.5)
    return parser.parse_args()


def as_numpy(value) -> np.ndarray:
    return np.asarray(getattr(value, "d", value))


def orient_2d(array: np.ndarray, dims: np.ndarray) -> np.ndarray:
    result = np.squeeze(as_numpy(array))
    expected = (int(dims[0]), int(dims[1]))
    if result.shape == expected:
        return result
    if result.T.shape == expected:
        return result.T
    raise ValueError(f"unexpected array shape {result.shape}, expected {expected}")


def composite_region(dataset, field_names: list[str], diameter: float,
                     r_max: float, z_min: float, z_max: float):
    level = int(dataset.max_level)
    refine = int(dataset.refine_by) ** level
    domain_dims = np.asarray(dataset.domain_dimensions, dtype=np.int64)
    fine_dims = domain_dims.copy()
    fine_dims[:int(dataset.dimensionality)] *= refine
    domain_left = as_numpy(dataset.domain_left_edge).astype(float)
    domain_right = as_numpy(dataset.domain_right_edge).astype(float)
    spacing = (domain_right - domain_left) / fine_dims

    requested_left = domain_left.copy()
    requested_right = domain_right.copy()
    requested_left[0] = 0.0
    requested_right[0] = r_max * diameter
    requested_left[1] = z_min * diameter
    requested_right[1] = z_max * diameter

    starts = np.floor(
        (requested_left - domain_left) / spacing + 1.0e-10
    ).astype(np.int64)
    stops = np.ceil(
        (requested_right - domain_left) / spacing - 1.0e-10
    ).astype(np.int64)
    starts = np.maximum(starts, 0)
    stops = np.minimum(stops, fine_dims)
    dims = stops - starts
    left_edge = domain_left + starts * spacing

    fields = [("boxlib", name) for name in field_names]
    grid = dataset.covering_grid(
        level, left_edge=left_edge, dims=dims, fields=fields,
        num_ghost_zones=0,
    )
    arrays = {
        name: orient_2d(grid[("boxlib", name)], dims).astype(np.float64)
        for name in field_names
    }
    r = left_edge[0] + (np.arange(dims[0]) + 0.5) * spacing[0]
    z = left_edge[1] + (np.arange(dims[1]) + 0.5) * spacing[1]

    level_map = np.zeros((dims[0], dims[1]), dtype=np.int8)
    for amr_grid in sorted(dataset.index.grids, key=lambda item: int(item.Level)):
        grid_level = int(amr_grid.Level)
        grid_left = as_numpy(amr_grid.LeftEdge).astype(float)
        grid_right = as_numpy(amr_grid.RightEdge).astype(float)
        local_lo = np.floor(
            (grid_left - left_edge) / spacing + 1.0e-8
        ).astype(np.int64)
        local_hi = np.ceil(
            (grid_right - left_edge) / spacing - 1.0e-8
        ).astype(np.int64)
        lo = np.maximum(local_lo, 0)
        hi = np.minimum(local_hi, dims)
        if np.all(hi > lo):
            level_map[lo[0]:hi[0], lo[1]:hi[1]] = np.maximum(
                level_map[lo[0]:hi[0], lo[1]:hi[1]], grid_level
            )

    return arrays, r, z, level_map, spacing, level


def smooth_axial(field: np.ndarray) -> np.ndarray:
    result = field.copy()
    for _ in range(2):
        old = result.copy()
        result[:, 1:-1] = (
            old[:, :-2] + 2.0 * old[:, 1:-1] + old[:, 2:]
        ) * 0.25
    return result


def positive_gradient(field: np.ndarray, z: np.ndarray) -> np.ndarray:
    return np.gradient(smooth_axial(field), z, axis=1, edge_order=2)


def peak_index(row: np.ndarray, z_d: np.ndarray, z_min: float,
               z_max: float) -> int:
    mask = (z_d >= z_min) & (z_d <= z_max)
    indices = np.flatnonzero(mask)
    if not len(indices):
        raise ValueError("empty shock tracking window")
    return int(indices[np.argmax(row[indices])])


def track_radial_front(gradient: np.ndarray, z_d: np.ndarray,
                       center_index: int, radial_count: int,
                       half_window: float = 0.14) -> tuple[np.ndarray, np.ndarray]:
    front = np.empty(radial_count, dtype=float)
    peaks = np.empty(radial_count, dtype=float)
    expected = float(z_d[center_index])
    for radial_index in range(radial_count):
        mask = np.abs(z_d - expected) <= half_window
        candidates = np.flatnonzero(mask)
        if not len(candidates):
            front[radial_index] = np.nan
            peaks[radial_index] = np.nan
            continue
        index = int(candidates[np.argmax(gradient[radial_index, candidates])])
        front[radial_index] = z_d[index]
        peaks[radial_index] = gradient[radial_index, index]
        expected = 0.75 * expected + 0.25 * float(z_d[index])
    return front, peaks


def gradient_fwhm(row: np.ndarray, index: int, spacing_d: float) -> float:
    peak = float(row[index])
    if not np.isfinite(peak) or peak <= 0.0:
        return float("nan")
    threshold = 0.5 * peak
    left = index
    right = index
    while left > 0 and row[left - 1] >= threshold:
        left -= 1
    while right + 1 < row.size and row[right + 1] >= threshold:
        right += 1
    return float((right - left + 1) * spacing_d)


def weighted_frame_mean(frames: np.ndarray, times: np.ndarray) -> np.ndarray:
    if len(times) == 1:
        return frames[0]
    weights = np.empty_like(times)
    weights[0] = 0.5 * (times[1] - times[0])
    weights[-1] = 0.5 * (times[-1] - times[-2])
    weights[1:-1] = 0.5 * (times[2:] - times[:-2])
    weights /= weights.sum()
    return np.tensordot(weights, frames, axes=(0, 0))


def linear_drift(times_ms: np.ndarray, positions: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(positions)
    slope, intercept = np.polyfit(times_ms[valid], positions[valid], 1)
    del intercept
    duration = times_ms[valid][-1] - times_ms[valid][0]
    return float(slope), float(slope * duration)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plotfiles = sorted((args.case_dir / "stat").glob("plt[0-9]*"))
    if not plotfiles:
        raise FileNotFoundError(f"no stat plotfiles under {args.case_dir}")

    yt.funcs.mylog.setLevel(50)
    frames = []
    final_means = None
    final_level_map = None
    final_level = None
    r = z = spacing = None

    for plot_index, plotfile in enumerate(plotfiles):
        dataset = yt.load(str(plotfile))
        field_names = [
            "Density", "pressure", "temperature", "x_velocity", "y_velocity"
        ]
        if plot_index == len(plotfiles) - 1:
            field_names += [
                "DensityMEAN", "temperatureMEAN",
                "x_velocityMEAN", "y_velocityMEAN",
            ]
        arrays, frame_r, frame_z, level_map, frame_spacing, level = composite_region(
            dataset, field_names, args.diameter, args.r_max,
            args.extract_z_min, args.extract_z_max,
        )
        if r is None:
            r, z, spacing = frame_r, frame_z, frame_spacing
        elif not (np.allclose(r, frame_r) and np.allclose(z, frame_z)):
            raise RuntimeError("composite extraction grid changed between frames")

        r_d = r / args.diameter
        z_d = z / args.diameter
        shape_count = int(np.count_nonzero(r_d <= args.shape_r_max))
        rho_gradient = positive_gradient(arrays["Density"], z)
        pressure_gradient = positive_gradient(arrays["pressure"], z)

        density_center_index = peak_index(
            rho_gradient[0], z_d, args.track_z_min, args.track_z_max
        )
        pressure_center_index = peak_index(
            pressure_gradient[0], z_d, args.track_z_min, args.track_z_max
        )
        misleading_index = peak_index(pressure_gradient[0], z_d, 1.8, 3.2)
        density_front, density_peaks = track_radial_front(
            rho_gradient, z_d, density_center_index, shape_count
        )
        pressure_front, pressure_peaks = track_radial_front(
            pressure_gradient, z_d, pressure_center_index, shape_count
        )
        sound_speed = np.sqrt(GAMMA * RGAS * arrays["temperature"])
        mach = np.hypot(arrays["x_velocity"], arrays["y_velocity"]) / sound_speed
        vorticity = (
            np.gradient(arrays["x_velocity"], z, axis=1, edge_order=2)
            - np.gradient(arrays["y_velocity"], r, axis=0, edge_order=2)
        )

        annulus = (r_d[:shape_count] >= 0.3) & (r_d[:shape_count] <= 0.5)
        level_roi = level_map[:shape_count, :]
        z_roi = (z_d >= args.track_z_min) & (z_d <= args.track_z_max)
        finest_fraction = float(np.mean(level_roi[:, z_roi] == level))
        disk_band = np.abs(z_d - density_front[0]) <= 0.18
        radial_velocity_absmax = float(np.max(np.abs(
            arrays["x_velocity"][:shape_count, :][:, disk_band]
        )))
        vorticity_absmax = float(np.max(np.abs(
            vorticity[:shape_count, :][:, disk_band]
        )))

        frames.append({
            "plotfile": plotfile.name,
            "time_s": float(dataset.current_time),
            "rho": arrays["Density"],
            "temperature": arrays["temperature"],
            "ur": arrays["x_velocity"],
            "uz": arrays["y_velocity"],
            "mach": mach,
            "density_front": density_front,
            "pressure_front": pressure_front,
            "density_peak": density_peaks,
            "pressure_peak": pressure_peaks,
            "density_center": float(density_front[0]),
            "pressure_center": float(pressure_front[0]),
            "pressure_global_max": float(z_d[misleading_index]),
            "density_shape_range": float(np.nanmax(density_front) - np.nanmin(density_front)),
            "pressure_shape_range": float(np.nanmax(pressure_front) - np.nanmin(pressure_front)),
            "density_center_minus_annulus": float(
                density_front[0] - np.nanmedian(density_front[annulus])
            ),
            "density_fwhm": gradient_fwhm(
                rho_gradient[0], density_center_index, spacing[1] / args.diameter
            ),
            "pressure_fwhm": gradient_fwhm(
                pressure_gradient[0], pressure_center_index,
                spacing[1] / args.diameter,
            ),
            "finest_fraction": finest_fraction,
            "axis_radial_velocity_at_front": float(
                arrays["x_velocity"][0, density_center_index]
            ),
            "disk_band_radial_velocity_absmax": radial_velocity_absmax,
            "disk_band_vorticity_absmax": vorticity_absmax,
        })
        if plot_index == len(plotfiles) - 1:
            final_means = {
                key: arrays[key] for key in (
                    "DensityMEAN", "temperatureMEAN",
                    "x_velocityMEAN", "y_velocityMEAN"
                )
            }
            final_level_map = level_map
            final_level = level
        print(
            f"{plotfile.name} t={float(dataset.current_time)*1e3:.4f} ms "
            f"x_disk={density_front[0]:.4f}D "
            f"shape={frames[-1]['density_shape_range']:.4f}D "
            f"L{level}={finest_fraction:.6f}",
            flush=True,
        )

    times = np.asarray([frame["time_s"] for frame in frames])
    times_ms = times * 1.0e3
    density_center = np.asarray([frame["density_center"] for frame in frames])
    pressure_center = np.asarray([frame["pressure_center"] for frame in frames])
    pressure_global = np.asarray([frame["pressure_global_max"] for frame in frames])
    density_shape = np.asarray([frame["density_shape_range"] for frame in frames])
    density_indent = np.asarray([
        frame["density_center_minus_annulus"] for frame in frames
    ])
    density_fwhm = np.asarray([frame["density_fwhm"] for frame in frames])
    finest_fraction = np.asarray([frame["finest_fraction"] for frame in frames])
    radial_velocity_absmax = np.asarray([
        frame["disk_band_radial_velocity_absmax"] for frame in frames
    ])
    vorticity_absmax = np.asarray([
        frame["disk_band_vorticity_absmax"] for frame in frames
    ])

    rho_frames = np.stack([frame["rho"] for frame in frames])
    temperature_frames = np.stack([frame["temperature"] for frame in frames])
    ur_frames = np.stack([frame["ur"] for frame in frames])
    uz_frames = np.stack([frame["uz"] for frame in frames])
    mach_frames = np.stack([frame["mach"] for frame in frames])
    offline_rho = weighted_frame_mean(rho_frames, times)
    offline_temperature = weighted_frame_mean(temperature_frames, times)
    offline_ur = weighted_frame_mean(ur_frames, times)
    offline_uz = weighted_frame_mean(uz_frames, times)
    offline_mach_frame_mean = weighted_frame_mean(mach_frames, times)
    offline_mach_ratio = np.hypot(offline_ur, offline_uz) / np.sqrt(
        GAMMA * RGAS * offline_temperature
    )
    online_mach_ratio = np.hypot(
        final_means["x_velocityMEAN"], final_means["y_velocityMEAN"]
    ) / np.sqrt(GAMMA * RGAS * final_means["temperatureMEAN"])

    density_slope, density_total_drift = linear_drift(times_ms, density_center)
    pressure_slope, pressure_total_drift = linear_drift(times_ms, pressure_center)
    metrics = {
        "case_dir": str(args.case_dir),
        "frame_count": len(frames),
        "time_start_ms": float(times_ms[0]),
        "time_end_ms": float(times_ms[-1]),
        "finest_level": int(final_level),
        "finest_spacing_over_D": float(spacing[1] / args.diameter),
        "disk_roi_finest_fraction_min": float(finest_fraction.min()),
        "disk_roi_finest_fraction_max": float(finest_fraction.max()),
        "density_front": {
            "mean_x_over_D": float(density_center.mean()),
            "std_x_over_D": float(density_center.std()),
            "min_x_over_D": float(density_center.min()),
            "max_x_over_D": float(density_center.max()),
            "linear_slope_D_per_ms": density_slope,
            "linear_total_drift_D": density_total_drift,
        },
        "pressure_front": {
            "mean_x_over_D": float(pressure_center.mean()),
            "std_x_over_D": float(pressure_center.std()),
            "min_x_over_D": float(pressure_center.min()),
            "max_x_over_D": float(pressure_center.max()),
            "linear_slope_D_per_ms": pressure_slope,
            "linear_total_drift_D": pressure_total_drift,
        },
        "unbounded_pressure_gradient_detector": {
            "mean_x_over_D": float(pressure_global.mean()),
            "std_x_over_D": float(pressure_global.std()),
            "min_x_over_D": float(pressure_global.min()),
            "max_x_over_D": float(pressure_global.max()),
        },
        "instantaneous_shape": {
            "median_radial_range_D": float(np.median(density_shape)),
            "max_radial_range_D": float(density_shape.max()),
            "median_center_minus_annulus_D": float(np.median(density_indent)),
            "max_abs_center_minus_annulus_D": float(np.max(np.abs(density_indent))),
        },
        "instantaneous_axis_shock_thickness": {
            "median_gradient_fwhm_D": float(np.nanmedian(density_fwhm)),
            "max_gradient_fwhm_D": float(np.nanmax(density_fwhm)),
            "median_gradient_fwhm_cells": float(
                np.nanmedian(density_fwhm) / (spacing[1] / args.diameter)
            ),
        },
        "transverse_mode": {
            "median_disk_band_abs_ur_m_per_s": float(
                np.median(radial_velocity_absmax)
            ),
            "max_disk_band_abs_ur_m_per_s": float(radial_velocity_absmax.max()),
            "median_disk_band_abs_vorticity_per_s": float(
                np.median(vorticity_absmax)
            ),
            "max_disk_band_abs_vorticity_per_s": float(vorticity_absmax.max()),
            "corr_shape_range_with_abs_ur": float(
                np.corrcoef(density_shape, radial_velocity_absmax)[0, 1]
            ),
            "corr_shape_range_with_abs_vorticity": float(
                np.corrcoef(density_shape, vorticity_absmax)[0, 1]
            ),
            "corr_abs_indent_with_abs_ur": float(
                np.corrcoef(np.abs(density_indent), radial_velocity_absmax)[0, 1]
            ),
            "corr_abs_indent_with_abs_vorticity": float(
                np.corrcoef(np.abs(density_indent), vorticity_absmax)[0, 1]
            ),
        },
        "mean_field_consistency": {
            "density_relative_l2_online_vs_frame_mean": float(
                np.linalg.norm(final_means["DensityMEAN"] - offline_rho)
                / np.linalg.norm(final_means["DensityMEAN"])
            ),
            "mach_ratio_relative_l2_online_vs_frame_ratio": float(
                np.linalg.norm(online_mach_ratio - offline_mach_ratio)
                / np.linalg.norm(online_mach_ratio)
            ),
            "mach_definition_relative_l2_frame_mean_vs_ratio_of_means": float(
                np.linalg.norm(offline_mach_frame_mean - offline_mach_ratio)
                / np.linalg.norm(offline_mach_frame_mean)
            ),
        },
    }

    with (args.output_dir / "metrics.json").open("w", encoding="utf-8") as stream:
        json.dump(metrics, stream, indent=2)
        stream.write("\n")
    with (args.output_dir / "frame_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "plotfile", "time_ms", "density_center_x_over_D",
            "pressure_center_x_over_D", "unbounded_pressure_peak_x_over_D",
            "density_shape_range_D", "density_center_minus_annulus_D",
            "density_gradient_fwhm_D", "finest_fraction",
            "axis_radial_velocity_at_front_m_per_s",
            "disk_band_radial_velocity_absmax_m_per_s",
            "disk_band_vorticity_absmax_per_s",
        ])
        for frame in frames:
            writer.writerow([
                frame["plotfile"], frame["time_s"] * 1.0e3,
                frame["density_center"], frame["pressure_center"],
                frame["pressure_global_max"], frame["density_shape_range"],
                frame["density_center_minus_annulus"], frame["density_fwhm"],
                frame["finest_fraction"],
                frame["axis_radial_velocity_at_front"],
                frame["disk_band_radial_velocity_absmax"],
                frame["disk_band_vorticity_absmax"],
            ])

    r_d = r / args.diameter
    z_d = z / args.diameter
    shape_count = int(np.count_nonzero(r_d <= args.shape_r_max))
    np.savez_compressed(
        args.output_dir / "compact_fields.npz",
        r_over_D=r_d,
        z_over_D=z_d,
        times_s=times,
        density_fronts=np.stack([frame["density_front"] for frame in frames]),
        pressure_fronts=np.stack([frame["pressure_front"] for frame in frames]),
        final_instantaneous_density=frames[-1]["rho"],
        final_instantaneous_mach=frames[-1]["mach"],
        online_density_mean=final_means["DensityMEAN"],
        offline_density_frame_mean=offline_rho,
        online_mach_ratio_of_means=online_mach_ratio,
        offline_mach_ratio_of_means=offline_mach_ratio,
        offline_mach_frame_mean=offline_mach_frame_mean,
        final_level_map=final_level_map,
    )
    figure, axes = plt.subplots(3, 2, figsize=(12.0, 12.0), constrained_layout=True)

    axis = axes[0, 0]
    axis.plot(times_ms, density_center, "o-", ms=3, label="first density jump")
    axis.plot(times_ms, pressure_center, "s-", ms=2.5, label="first pressure jump")
    axis.plot(times_ms, pressure_global, ".--", lw=0.8,
              label="largest pressure gradient in 1.8-3.2D")
    axis.set(xlabel="time (ms)", ylabel="axis front x/D", title="Disk tracking")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[0, 1]
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(frames)))
    for frame, color in zip(frames, colors):
        axis.plot(frame["density_front"], r_d[:shape_count], color=color,
                  alpha=0.55, lw=0.8)
    axis.set(xlabel="instantaneous density front x/D", ylabel="r/D",
             title="All instantaneous disk fronts")
    axis.grid(alpha=0.25)

    axis = axes[1, 0]
    axis.plot(times_ms, density_shape, "o-", ms=3, label="radial range")
    axis.plot(times_ms, np.abs(density_indent), "s-", ms=3,
              label="abs(center - r/D=0.3-0.5 median)")
    axis.plot(times_ms, density_fwhm, "^-", ms=3,
              label="axis gradient FWHM")
    axis.axhline(spacing[1] / args.diameter, color="black", ls=":", lw=1,
                 label="one L4 cell")
    axis.set(xlabel="time (ms)", ylabel="D", title="Shape and shock thickness")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[1, 1]
    level_image = axis.pcolormesh(
        z_d, r_d, final_level_map, shading="nearest", cmap="viridis",
        vmin=0, vmax=final_level,
    )
    figure.colorbar(level_image, ax=axis, label="AMR level")
    axis.axvspan(args.track_z_min, args.track_z_max, color="white", alpha=0.12)
    axis.axhline(args.shape_r_max, color="white", ls="--", lw=0.8)
    axis.set(xlabel="x/D", ylabel="r/D", title="Final composite level map")

    axis = axes[2, 0]
    mean_image = axis.pcolormesh(
        z_d, r_d, online_mach_ratio, shading="nearest", cmap="turbo",
        vmin=0.0, vmax=4.8,
    )
    axis.contour(z_d, r_d, online_mach_ratio, levels=[1.0], colors="white",
                 linewidths=0.8)
    selected = [0, len(frames) // 2, len(frames) - 1]
    for index, line_style in zip(selected, [":", "--", "-"]):
        axis.plot(frames[index]["density_front"], r_d[:shape_count],
                  color="black", ls=line_style, lw=1.0,
                  label=f"rho jump {times_ms[index]:.2f} ms")
    figure.colorbar(mean_image, ax=axis, label="online ratio-of-means Mach")
    axis.set(xlim=(1.9, 3.3), ylim=(0.0, args.r_max), xlabel="x/D", ylabel="r/D",
             title="Mean field with instantaneous fronts")
    axis.legend(fontsize=7, loc="upper right")

    axis = axes[2, 1]
    axis.plot(z_d, online_mach_ratio[0], label="online ratio of means")
    axis.plot(z_d, offline_mach_ratio[0], "--", label="frame ratio of means")
    axis.plot(z_d, offline_mach_frame_mean[0], ":", label="mean of frame Mach")
    axis.plot(z_d, frames[-1]["mach"][0], color="black", alpha=0.55,
              label="final instantaneous Mach")
    axis.axhline(1.0, color="gray", lw=0.8)
    axis.set(xlim=(1.9, 3.3), xlabel="x/D", ylabel="axis Mach",
             title="Statistics and Mach-definition audit")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    figure.suptitle("NPR=15 Mach-disk numerical audit", fontsize=15)
    figure.savefig(args.output_dir / "npr15_disk_audit.png", dpi=220)
    plt.close(figure)

    print(json.dumps(metrics, indent=2), flush=True)
    print(f"saved diagnostics to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
