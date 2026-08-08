#!/usr/bin/env python3
"""Extract conservative, AMR-aware SRP tri-nozzle axis profiles.

The primitive fields are sampled from a thin finest-resolution yt
``smoothed_covering_grid``.  The IBM flags are sampled independently with a
categorical ``covering_grid``.  A requested axis value generally lies between
four y-z cell centres, so a value is reported only when *all four* interpolation
corners are fluid.  This deliberately prevents a bilinear stencil from
crossing the solid/ghost region.  Invalid runs remain NaN and are never joined
by the profile plot or by the density-gradient calculation.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
import numpy as np


GAMMA = 1.4
DEFAULT_PINF = 574.56
DEFAULT_DB = 0.127
DEFAULT_BODY_CENTRE_X = 0.300
DEFAULT_NOZZLE_PLANE_X = 0.299
AXES = (
    ("body_axis", 0.275, 0.275, "body axis"),
    ("nozzle_0", 0.307, 0.275, "nozzle 0 axis"),
    ("nozzle_1", 0.259, 0.302, "nozzle 1 axis"),
    ("nozzle_2", 0.259, 0.248, "nozzle 2 axis"),
)
CONTINUOUS_FIELDS = (
    "Density", "pressure", "temperature",
    "x_velocity", "y_velocity", "z_velocity",
)
GEOMETRY_FIELDS = ("sld", "ghs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plotfile", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--x-min", type=float, default=0.009)
    parser.add_argument("--x-max", type=float, default=0.631)
    parser.add_argument("--p-inf", type=float, default=DEFAULT_PINF)
    parser.add_argument("--body-diameter", type=float, default=DEFAULT_DB)
    parser.add_argument("--body-centre-x", type=float, default=DEFAULT_BODY_CENTRE_X)
    parser.add_argument("--nozzle-plane-x", type=float, default=DEFAULT_NOZZLE_PLANE_X)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def yt_numpy(values) -> np.ndarray:
    if hasattr(values, "to_value"):
        return np.asarray(values.to_value(), dtype=np.float64)
    if hasattr(values, "d"):
        return np.asarray(values.d, dtype=np.float64)
    return np.asarray(values, dtype=np.float64)


def bracket(coords: np.ndarray, value: float) -> tuple[int, int, float]:
    """Return cell-centre interpolation indices and a bounded linear weight."""
    upper = int(np.searchsorted(coords, value, side="left"))
    upper = max(1, min(upper, len(coords) - 1))
    lower = upper - 1
    denominator = coords[upper] - coords[lower]
    weight = 0.0 if denominator == 0.0 else (value - coords[lower]) / denominator
    return lower, upper, float(np.clip(weight, 0.0, 1.0))


def bilinear_corners(
    array: np.ndarray,
    jy0: int,
    jy1: int,
    wy: float,
    kz0: int,
    kz1: int,
    wz: float,
) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    """Interpolate an x-line and return the four unblended corner lines."""
    corners = (
        array[:, jy0, kz0], array[:, jy1, kz0],
        array[:, jy0, kz1], array[:, jy1, kz1],
    )
    value_z0 = (1.0 - wy) * corners[0] + wy * corners[1]
    value_z1 = (1.0 - wy) * corners[2] + wy * corners[3]
    return (1.0 - wz) * value_z0 + wz * value_z1, corners


def segment_ids(valid: np.ndarray) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Label contiguous valid runs; invalid samples receive -1."""
    labels = np.full(valid.shape, -1, dtype=np.int32)
    padded = np.r_[False, valid, False]
    starts = np.flatnonzero(padded[1:] & ~padded[:-1])
    stops = np.flatnonzero(~padded[1:] & padded[:-1])
    segments = []
    for identifier, (start, stop) in enumerate(zip(starts, stops)):
        labels[start:stop] = identifier
        segments.append((int(start), int(stop)))  # stop is exclusive
    return labels, segments


def segmented_abs_gradient(
    values: np.ndarray,
    coords: np.ndarray,
    segments: list[tuple[int, int]],
) -> np.ndarray:
    result = np.full(values.shape, np.nan, dtype=np.float64)
    for start, stop in segments:
        if stop - start >= 2:
            result[start:stop] = np.abs(
                np.gradient(values[start:stop], coords[start:stop], edge_order=1)
            )
    return result


def source_level_line(ds, x: np.ndarray, y: float, z: float) -> np.ndarray:
    """Map the finest native AMR grid covering each requested cell centre."""
    levels = np.zeros(x.shape, dtype=np.uint8)
    tolerance = 1.0e-12 * float(np.max(yt_numpy(ds.domain_width)))
    for grid in sorted(ds.index.grids, key=lambda item: int(item.Level)):
        left = yt_numpy(grid.LeftEdge)
        right = yt_numpy(grid.RightEdge)
        if not (left[1] - tolerance <= y < right[1] - tolerance):
            continue
        if not (left[2] - tolerance <= z < right[2] - tolerance):
            continue
        first = int(np.searchsorted(x, left[0] - tolerance, side="left"))
        last = int(np.searchsorted(x, right[0] - tolerance, side="left"))
        first, last = max(0, first), min(len(x), last)
        if first < last:
            levels[first:last] = np.maximum(levels[first:last], np.uint8(grid.Level))
    return levels


def extract_profile(
    ds,
    arrays: dict[str, np.ndarray],
    local_coords: tuple[np.ndarray, np.ndarray, np.ndarray],
    global_coords: tuple[np.ndarray, np.ndarray, np.ndarray],
    y_value: float,
    z_value: float,
    p_inf: float,
    body_diameter: float,
    body_centre_x: float,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    x, y, z = local_coords
    _, global_y, global_z = global_coords
    jy0, jy1, wy = bracket(y, y_value)
    kz0, kz1, wz = bracket(z, z_value)

    interpolated: dict[str, np.ndarray] = {}
    corners_by_field: dict[str, tuple[np.ndarray, ...]] = {}
    for field in CONTINUOUS_FIELDS:
        interpolated[field], corners_by_field[field] = bilinear_corners(
            arrays[field], jy0, jy1, wy, kz0, kz1, wz
        )

    _, sld_corners = bilinear_corners(arrays["sld"], jy0, jy1, wy, kz0, kz1, wz)
    _, ghs_corners = bilinear_corners(arrays["ghs"], jy0, jy1, wy, kz0, kz1, wz)
    # Conservative categorical values: one solid/ghost stencil corner is enough
    # to invalidate the bilinear primitive sample.
    sld = np.max(np.stack(sld_corners), axis=0)
    ghs = np.max(np.stack(ghs_corners), axis=0)
    valid = np.ones(x.shape, dtype=bool)
    for corner_index in range(4):
        rho_corner = corners_by_field["Density"][corner_index]
        p_corner = corners_by_field["pressure"][corner_index]
        valid &= np.isfinite(rho_corner) & (rho_corner > 0.0)
        valid &= np.isfinite(p_corner) & (p_corner > 0.0)
        valid &= np.isfinite(sld_corners[corner_index]) & (sld_corners[corner_index] < 0.5)
        valid &= np.isfinite(ghs_corners[corner_index]) & (ghs_corners[corner_index] < 0.5)
        for field in ("temperature", "x_velocity", "y_velocity", "z_velocity"):
            valid &= np.isfinite(corners_by_field[field][corner_index])

    rho = interpolated["Density"]
    pressure = interpolated["pressure"]
    temperature = interpolated["temperature"]
    u = interpolated["x_velocity"]
    v = interpolated["y_velocity"]
    w = interpolated["z_velocity"]
    sound_speed = np.sqrt(
        GAMMA * np.maximum(pressure, 0.0) /
        np.maximum(rho, np.finfo(np.float64).tiny)
    )
    mach = np.sqrt(u * u + v * v + w * w) / np.maximum(
        sound_speed, np.finfo(np.float64).tiny
    )

    # Native source level is the least-resolved one among the four bilinear
    # corners.  It is diagnostic only; the primitives themselves are the yt
    # ghost-aware AMR composite values.
    corner_levels = []
    for j in (jy0, jy1):
        for k in (kz0, kz1):
            corner_levels.append(source_level_line(ds, x, y[j], z[k]))
    source_level = np.min(np.stack(corner_levels), axis=0).astype(np.uint8)

    labels, segments = segment_ids(valid)
    for values in (rho, pressure, temperature, u, v, w, mach):
        values[~valid] = np.nan
    density_gradient = segmented_abs_gradient(rho, x, segments)

    result = {
        "x_m": x.copy(),
        "x_over_Db": x / body_diameter,
        "X_body": (x - body_centre_x) / body_diameter,
        "Mach": mach,
        "p_Pa": pressure,
        "p_over_pinf": pressure / p_inf,
        "rho_kg_m3": rho,
        "T_K": temperature,
        "u_m_s": u,
        "v_m_s": v,
        "w_m_s": w,
        "axis_velocity_m_s": u.copy(),
        "abs_drhodx_kg_m4": density_gradient,
        "sld": sld,
        "ghs": ghs,
        "fluid_valid": valid.astype(np.uint8),
        "fluid_segment_id": labels,
        "source_level": source_level,
    }
    stencil = {
        "requested_y_m": float(y_value),
        "requested_z_m": float(z_value),
        "y_cell_centres_m": [float(y[jy0]), float(y[jy1])],
        "z_cell_centres_m": [float(z[kz0]), float(z[kz1])],
        "y_weight_upper": float(wy),
        "z_weight_upper": float(wz),
        "global_y_indices": [
            int(np.searchsorted(global_y, y[jy0])),
            int(np.searchsorted(global_y, y[jy1])),
        ],
        "global_z_indices": [
            int(np.searchsorted(global_z, z[kz0])),
            int(np.searchsorted(global_z, z[kz1])),
        ],
        "fluid_segments_index": [[start, stop - 1] for start, stop in segments],
        "fluid_segments_x_m": [
            [float(x[start]), float(x[stop - 1])] for start, stop in segments
        ],
    }
    return result, stencil


def finite_range(values: np.ndarray) -> list[float | None]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return [None, None]
    return [float(np.min(finite)), float(np.max(finite))]


def write_profile_csv(path: Path, profile: dict[str, np.ndarray]) -> None:
    columns = list(profile)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row_index in range(len(profile["x_m"])):
            row = []
            for column in columns:
                value = profile[column][row_index]
                if np.issubdtype(profile[column].dtype, np.integer):
                    row.append(str(int(value)))
                else:
                    row.append("{:.12g}".format(float(value)))
            writer.writerow(row)


def configure_matplotlib() -> None:
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8.5,
        "axes.labelsize": 9.0,
        "axes.titlesize": 9.0,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def render_figure(
    profiles: dict[str, dict[str, np.ndarray]],
    output_stem: Path,
    time_ms: float,
    nozzle_plane_x: float,
    p_inf: float,
    body_diameter: float,
    dpi: int,
) -> None:
    colors = {
        "body_axis": "#111111",
        "nozzle_0": "#d62728",
        "nozzle_1": "#2474b7",
        "nozzle_2": "#269c3c",
    }
    labels = {
        name: "{} ($y={:.3f}$ m, $z={:.3f}$ m)".format(label, y, z)
        for name, y, z, label in AXES
    }
    configure_matplotlib()
    fig, axes = plt.subplots(
        3, 1, figsize=(7.08, 7.5), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.08, 1.08]},
        constrained_layout=True,
    )
    for name, profile in profiles.items():
        style = "-" if name == "body_axis" else "--"
        linewidth = 1.35 if name == "body_axis" else 1.05
        axes[0].plot(profile["x_m"], profile["Mach"], style,
                     color=colors[name], lw=linewidth, label=labels[name])
        axes[1].plot(profile["x_m"], profile["p_over_pinf"], style,
                     color=colors[name], lw=linewidth)
        axes[2].plot(profile["x_m"], profile["abs_drhodx_kg_m4"], style,
                     color=colors[name], lw=linewidth)

    for index, axis in enumerate(axes):
        axis.axvline(
            nozzle_plane_x, color="#666666", ls=":", lw=0.9,
            label=r"nozzle exit plane, $x=0.299$ m" if index == 0 else None,
        )
        axis.grid(True, which="both", color="#b8b8b8", alpha=0.24, lw=0.45)
    axes[0].set_ylabel(r"Mach number, $M$")
    axes[0].set_ylim(bottom=0.0)
    axes[0].legend(loc="upper right", fontsize=6.9, ncol=2, handlelength=2.4)

    axes[1].set_ylabel(r"Pressure ratio, $p/p_\infty$")
    axes[1].set_yscale("log")
    pressure_values = np.concatenate([
        profile["p_over_pinf"][np.isfinite(profile["p_over_pinf"])]
        for profile in profiles.values()
    ])
    if pressure_values.size:
        lower = 10.0 ** np.floor(np.log10(max(np.min(pressure_values), 1.0e-12)))
        upper = 10.0 ** np.ceil(np.log10(np.max(pressure_values)))
        axes[1].set_ylim(lower, upper)
    axes[1].yaxis.set_major_formatter(LogFormatterMathtext())
    pressure_axis = axes[1].secondary_yaxis(
        "right", functions=(lambda ratio: ratio * p_inf, lambda p: p / p_inf)
    )
    pressure_axis.set_ylabel(r"Pressure, $p$ [Pa]")

    axes[2].set_ylabel(r"$|\partial \rho/\partial x|$ [kg m$^{-4}$]")
    axes[2].set_yscale("log")
    axes[2].set_ylim(bottom=1.0e-3)
    axes[2].set_xlabel(r"Axial coordinate, $x$ [m]")
    top_axis = axes[0].secondary_xaxis(
        "top", functions=(lambda x: x / body_diameter,
                          lambda x_over_db: x_over_db * body_diameter)
    )
    top_axis.set_xlabel(r"$x/D_b$")
    fig.suptitle(
        "SRP tri-nozzle axis profiles, $t={:.3f}$ ms".format(time_ms),
        fontsize=10.0,
    )
    fig.savefig(output_stem.with_suffix(".png"), dpi=dpi)
    fig.savefig(output_stem.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not (args.x_min < args.x_max):
        raise ValueError("x-min must be less than x-max")
    if args.p_inf <= 0.0 or args.body_diameter <= 0.0:
        raise ValueError("p-inf and body-diameter must be positive")

    try:
        import yt
    except ImportError as error:
        raise RuntimeError("Load the CX3 yt/4.3 module before running this script") from error
    yt.funcs.mylog.setLevel(30)

    plotfile = Path(args.plotfile).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    dataset = yt.load(str(plotfile))
    available = set(dataset.field_list)
    missing = [name for name in CONTINUOUS_FIELDS + GEOMETRY_FIELDS
               if ("boxlib", name) not in available]
    if missing:
        raise RuntimeError("Plotfile is missing fields: {}".format(", ".join(missing)))

    refine_factor = int(dataset.refine_by) ** int(dataset.max_level)
    fine_dims = np.asarray(dataset.domain_dimensions, dtype=np.int64) * refine_factor
    domain_left = yt_numpy(dataset.domain_left_edge)
    domain_right = yt_numpy(dataset.domain_right_edge)
    spacing = (domain_right - domain_left) / fine_dims
    global_coords = tuple(
        domain_left[axis] + (np.arange(fine_dims[axis]) + 0.5) * spacing[axis]
        for axis in range(3)
    )

    # smoothed_covering_grid recursively requests two parent cells.  Keep the
    # level-0-equivalent halo inside each non-periodic boundary.
    interpolation_halo = 2 * refine_factor
    x_selected = np.flatnonzero(
        (global_coords[0] >= args.x_min) & (global_coords[0] <= args.x_max)
    )
    if x_selected.size < 2:
        raise RuntimeError("Requested x range contains fewer than two fine cells")
    x_start = max(int(x_selected[0]), interpolation_halo)
    x_stop = min(int(x_selected[-1]) + 1, int(fine_dims[0]) - interpolation_halo)

    y_indices = []
    z_indices = []
    for _, y_value, z_value, _ in AXES:
        y0, y1, _ = bracket(global_coords[1], y_value)
        z0, z1, _ = bracket(global_coords[2], z_value)
        y_indices.extend([y0, y1])
        z_indices.extend([z0, z1])
    starts = np.asarray([x_start, min(y_indices), min(z_indices)], dtype=np.int64)
    stops = np.asarray([x_stop, max(y_indices) + 1, max(z_indices) + 1], dtype=np.int64)
    dims = stops - starts
    left_edge = domain_left + starts * spacing

    print(
        "Loading AMR composite prism dims={} from {}".format(
            tuple(int(value) for value in dims), plotfile
        ),
        flush=True,
    )
    continuous_fields = [("boxlib", name) for name in CONTINUOUS_FIELDS]
    geometry_fields = [("boxlib", name) for name in GEOMETRY_FIELDS]
    smooth = dataset.smoothed_covering_grid(
        int(dataset.max_level), left_edge=left_edge, dims=dims,
        fields=continuous_fields, num_ghost_zones=0, use_pbar=False,
    )
    categorical = dataset.covering_grid(
        int(dataset.max_level), left_edge=left_edge, dims=dims,
        fields=geometry_fields, num_ghost_zones=0, use_pbar=False,
    )
    arrays = {
        name: yt_numpy(smooth[("boxlib", name)]) for name in CONTINUOUS_FIELDS
    }
    arrays.update({
        name: yt_numpy(categorical[("boxlib", name)]) for name in GEOMETRY_FIELDS
    })
    local_coords = tuple(
        global_coords[axis][starts[axis]:stops[axis]] for axis in range(3)
    )

    profiles: dict[str, dict[str, np.ndarray]] = {}
    stencils: dict[str, dict[str, object]] = {}
    axis_metadata: dict[str, dict[str, object]] = {}
    for name, y_value, z_value, label in AXES:
        profile, stencil = extract_profile(
            dataset, arrays, local_coords, global_coords,
            y_value, z_value, args.p_inf, args.body_diameter,
            args.body_centre_x,
        )
        profiles[name] = profile
        stencils[name] = stencil
        csv_name = "{}_profile.csv".format(name)
        write_profile_csv(output_dir / csv_name, profile)
        valid = profile["fluid_valid"].astype(bool)
        levels, counts = np.unique(profile["source_level"], return_counts=True)
        axis_metadata[name] = {
            "label": label,
            "y_m": float(y_value),
            "z_m": float(z_value),
            "csv": csv_name,
            "sample_count": int(len(valid)),
            "valid_count": int(np.count_nonzero(valid)),
            "invalid_count": int(np.count_nonzero(~valid)),
            "fluid_segment_count": int(np.max(profile["fluid_segment_id"]) + 1)
            if np.any(valid) else 0,
            "mach_range_valid": finite_range(profile["Mach"]),
            "pressure_Pa_range_valid": finite_range(profile["p_Pa"]),
            "source_level_counts_all_samples": {
                str(int(level)): int(count) for level, count in zip(levels, counts)
            },
            "interpolation_stencil": stencil,
        }
        print(
            "{}: {}/{} valid samples, {} fluid segments".format(
                name, axis_metadata[name]["valid_count"], len(valid),
                axis_metadata[name]["fluid_segment_count"],
            ),
            flush=True,
        )

    combined_payload = {
        "plotfile": np.asarray(str(plotfile)),
        "simulation_time_s": np.asarray(float(yt_numpy(dataset.current_time))),
        "p_inf_Pa": np.asarray(args.p_inf),
        "body_diameter_m": np.asarray(args.body_diameter),
        "body_centre_x_m": np.asarray(args.body_centre_x),
        "nozzle_exit_plane_x_m": np.asarray(args.nozzle_plane_x),
    }
    for name, profile in profiles.items():
        for field, values in profile.items():
            combined_payload["{}__{}".format(name, field)] = values
    np.savez_compressed(output_dir / "axis_profiles_combined.npz", **combined_payload)

    time_s = float(yt_numpy(dataset.current_time))
    figure_stem = output_dir / "axis_profiles_Mach_pressure_schlieren"
    render_figure(
        profiles, figure_stem, time_s * 1.0e3,
        args.nozzle_plane_x, args.p_inf, args.body_diameter, args.dpi,
    )

    nearest_nozzle_index = int(np.argmin(np.abs(local_coords[0] - args.nozzle_plane_x)))
    summary = {
        "schema_version": 1,
        "plotfile": str(plotfile),
        "simulation_time_s": time_s,
        "simulation_time_ms": time_s * 1.0e3,
        "domain_left_m": domain_left.tolist(),
        "domain_right_m": domain_right.tolist(),
        "finest_level": int(dataset.max_level),
        "refine_by": int(dataset.refine_by),
        "fine_grid_shape": [int(value) for value in fine_dims],
        "fine_spacing_m": spacing.tolist(),
        "sampled_x_range_m": [float(local_coords[0][0]), float(local_coords[0][-1])],
        "sample_count_per_axis": int(len(local_coords[0])),
        "body_diameter_m": float(args.body_diameter),
        "body_centre_x_m": float(args.body_centre_x),
        "p_inf_Pa": float(args.p_inf),
        "nozzle_exit_plane_x_m": float(args.nozzle_plane_x),
        "nearest_sample_to_nozzle_exit": {
            "index": nearest_nozzle_index,
            "x_m": float(local_coords[0][nearest_nozzle_index]),
            "offset_m": float(local_coords[0][nearest_nozzle_index] - args.nozzle_plane_x),
        },
        "sampler": {
            "continuous_fields": "yt smoothed_covering_grid at max_level",
            "geometry_fields": "yt covering_grid at max_level (categorical)",
            "axis_interpolation": "bilinear in y-z cell centres",
            "validity_rule": (
                "all four y-z interpolation corners must have sld<0.5, "
                "ghs<0.5, finite primitives, rho>0 and p>0"
            ),
            "invalid_handling": (
                "continuous outputs are NaN; plot and d(rho)/dx are split into "
                "contiguous fluid segments and never bridge invalid cells"
            ),
            "source_level": "minimum native AMR level among four interpolation corners",
            "axis_velocity": "signed x_velocity because every nozzle axis is x-directed",
            "schlieren_proxy": (
                "absolute axial density derivative within each continuous fluid segment"
            ),
            "interpolation_halo_fine_cells": int(interpolation_halo),
        },
        "columns": {
            "x_over_Db": "x / D_b",
            "X_body": "(x - x_body_centre) / D_b",
            "p_over_pinf": "p / p_inf",
            "axis_velocity_m_s": "signed u = x_velocity",
            "fluid_segment_id": "zero-based within each axis; -1 is invalid",
        },
        "axes": axis_metadata,
        "outputs": {
            "combined_npz": "axis_profiles_combined.npz",
            "figure_png": figure_stem.with_suffix(".png").name,
            "figure_pdf": figure_stem.with_suffix(".pdf").name,
        },
    }
    with (output_dir / "axis_profiles_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print("Wrote axis-profile products to {}".format(output_dir), flush=True)


if __name__ == "__main__":
    main()
