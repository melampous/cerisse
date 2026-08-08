#!/usr/bin/env python3
"""Post-process the latest 3-D SRP tri-nozzle AMReX plotfile.

The script intentionally reads only selected AMR slices.  It produces
centre-plane Mach, pressure, and numerical-schlieren images together with
quantitative profiles along the body axis and all three nozzle axes.
"""

import argparse
import csv
import importlib.util
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


GAMMA = 1.4
BODY_AXIS = (0.275, 0.275)
NOZZLE_X = 0.299
NOZZLE_12_PLANE_Y = 0.259
NOZZLE_AXES = (
    ("nozzle_0", 0.307, 0.275),
    ("nozzle_1", 0.259, 0.302),
    ("nozzle_2", 0.259, 0.248),
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Mach, pressure, schlieren, and axis profiles for SRP tri-nozzle."
    )
    parser.add_argument("--plotfile", required=True)
    parser.add_argument("--reader", required=True, help="Path to the existing schlieren.py reader")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=0.56)
    parser.add_argument("--y-min", type=float, default=0.08)
    parser.add_argument("--y-max", type=float, default=0.47)
    parser.add_argument("--mach-max", type=float, default=13.0)
    parser.add_argument("--pressure-min", type=float, default=1.0e1)
    parser.add_argument("--pressure-max", type=float, default=5.0e6)
    parser.add_argument("--schlieren-k", type=float, default=18.0)
    return parser.parse_args()


def load_reader(path):
    spec = importlib.util.spec_from_file_location("amrex_schlieren_reader", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load reader: {}".format(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def nearest_index(coords, value):
    return int(np.argmin(np.abs(coords - value)))


def bracket(coords, value):
    """Return two bracketing cell-centre indices and interpolation weight."""
    hi = int(np.searchsorted(coords, value, side="left"))
    hi = max(1, min(hi, len(coords) - 1))
    lo = hi - 1
    denom = coords[hi] - coords[lo]
    weight = 0.0 if denom == 0.0 else float((value - coords[lo]) / denom)
    weight = min(1.0, max(0.0, weight))
    return lo, hi, weight


def interp_xy_plane(z_planes, y_coords, z_coords, y_value, z_value):
    """Bilinearly interpolate a line in y-z from XY planes; result varies in x."""
    jy0, jy1, wy = bracket(y_coords, y_value)
    kz0, kz1, wz = bracket(z_coords, z_value)
    line_z0 = (1.0 - wy) * z_planes[kz0][:, jy0] + wy * z_planes[kz0][:, jy1]
    line_z1 = (1.0 - wy) * z_planes[kz1][:, jy0] + wy * z_planes[kz1][:, jy1]
    return (1.0 - wz) * line_z0 + wz * line_z1


def interp_z_plane(z_planes, z_coords, z_value):
    """Linearly interpolate an XY plane at a requested z coordinate."""
    kz0, kz1, wz = bracket(z_coords, z_value)
    return (1.0 - wz) * z_planes[kz0] + wz * z_planes[kz1]


def mach_from_primitives(rho, pressure, ux, uy, uz):
    eps = 1.0e-30
    speed = np.sqrt(ux * ux + uy * uy + uz * uz)
    sound = np.sqrt(GAMMA * np.maximum(pressure, 0.0) / np.maximum(rho, eps))
    return speed / np.maximum(sound, eps)


def fluid_mask(sld, rho, pressure):
    return (sld <= 0.5) & np.isfinite(rho) & np.isfinite(pressure) & (rho > 0.0) & (pressure > 0.0)


def axis_profile(fields, y_coords, z_coords, y_value, z_value):
    values = {}
    for name, planes in fields.items():
        values[name] = interp_xy_plane(planes, y_coords, z_coords, y_value, z_value)

    valid = fluid_mask(values["sld"], values["Density"], values["pressure"])
    mach = mach_from_primitives(
        values["Density"], values["pressure"], values["x_velocity"],
        values["y_velocity"], values["z_velocity"],
    )
    mach[~valid] = np.nan
    pressure = values["pressure"].copy()
    pressure[~valid] = np.nan
    density = values["Density"].copy()
    density[~valid] = np.nan
    return {
        "mach": mach,
        "pressure": pressure,
        "density": density,
        "sld": values["sld"],
        "valid": valid,
    }


def segmented_abs_gradient(values, coords):
    """Gradient without bridging solid/NaN gaps."""
    out = np.full(values.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(values)
    starts = np.where(finite & np.r_[True, ~finite[:-1]])[0]
    ends = np.where(finite & np.r_[~finite[1:], True])[0]
    for start, end in zip(starts, ends):
        if end - start + 1 >= 2:
            out[start:end + 1] = np.abs(
                np.gradient(values[start:end + 1], coords[start:end + 1])
            )
    return out


def numerical_schlieren(field, valid, d0, d1, k):
    """Return a robust fluid-only density-gradient schlieren image and scale.

    The 99.5th-percentile reference keeps a single extreme cell from washing
    out the rest of the shocks.  Solid/non-finite cells are excluded from the
    normalisation and masked in the returned image.
    """
    clean = np.asarray(field, dtype=np.float64).copy()
    good = valid & np.isfinite(clean)
    fill = np.nanmedian(clean[good]) if np.any(good) else 0.0
    clean[~good] = fill
    g0, g1 = np.gradient(clean, d0, d1)
    grad = np.hypot(g0, g1)

    # Do not let the artificial density jump introduced when filling an IBM
    # cell determine the normalisation scale.
    reference_cells = good.copy()
    reference_cells[1:, :] &= good[:-1, :]
    reference_cells[:-1, :] &= good[1:, :]
    reference_cells[:, 1:] &= good[:, :-1]
    reference_cells[:, :-1] &= good[:, 1:]
    samples = grad[reference_cells & np.isfinite(grad)]
    scale = float(np.nanpercentile(samples, 99.5)) if samples.size else 0.0
    if not np.isfinite(scale) or scale <= 0.0:
        scale = float(np.nanmax(grad[good])) if np.any(good) else 1.0
    if not np.isfinite(scale) or scale <= 0.0:
        scale = 1.0

    result = np.exp(-k * grad / scale)
    result[~good] = np.nan
    return result, scale


def add_solid_overlay(ax, x, y, solid):
    if not np.any(solid):
        return
    overlay = np.ma.masked_where(~solid, np.ones_like(solid, dtype=float))
    ax.pcolormesh(x, y, overlay, shading="auto", cmap="Greys", vmin=0.0, vmax=1.0, alpha=0.72)
    try:
        ax.contour(x, y, solid.astype(float), levels=[0.5], colors="white", linewidths=0.55)
    except ValueError:
        pass


def decorate_axial_plane(ax, plane):
    if plane == "xy":
        ax.axhline(BODY_AXIS[0], color="white", linestyle="--", linewidth=0.75, alpha=0.9,
                   label="body axis")
        ax.axhline(NOZZLE_AXES[0][1], color="#00e5ff", linestyle=":", linewidth=0.9,
                   alpha=0.95, label="nozzle-0 axis")
        ax.set_ylabel("y [m]")
    elif plane == "xz":
        ax.axhline(BODY_AXIS[1], color="white", linestyle="--", linewidth=0.75, alpha=0.9,
                   label="body-axis z")
        ax.axhline(NOZZLE_AXES[1][2], color="#00e5ff", linestyle=":", linewidth=0.9,
                   alpha=0.95, label="nozzle-1 axis")
        ax.axhline(NOZZLE_AXES[2][2], color="#52ff65", linestyle=":", linewidth=0.9,
                   alpha=0.95, label="nozzle-2 axis")
        ax.set_ylabel("z [m]")
    else:
        raise ValueError("Unknown plane: {}".format(plane))
    ax.axvline(NOZZLE_X, color="white", linestyle=":", linewidth=0.65, alpha=0.8)
    ax.set_xlabel("x [m]")
    ax.set_aspect("equal")


def save_field_figure(path, x, transverse, data, solid, title, plane, cmap, norm=None,
                      vmin=None, vmax=None, colorbar_label=""):
    fig, ax = plt.subplots(figsize=(12.2, 7.2), constrained_layout=True)
    pcm = ax.pcolormesh(x, transverse, data, shading="auto", cmap=cmap, norm=norm,
                        vmin=vmin, vmax=vmax)
    add_solid_overlay(ax, x, transverse, solid)
    decorate_axial_plane(ax, plane)
    ax.set_title(title)
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label(colorbar_label)
    ax.legend(loc="lower right", framealpha=0.85, fontsize=8)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    reader = load_reader(args.reader)

    header = reader._parse_header_3d(args.plotfile)
    level = header["finest_level"]
    box = header["all_boxes"][level]
    nx = int(box[3]) - int(box[0]) + 1
    ny = int(box[4]) - int(box[1]) + 1
    nz = int(box[5]) - int(box[2]) + 1
    lo = np.asarray(header["lo"], dtype=float)
    hi = np.asarray(header["hi"], dtype=float)
    spacing = (hi - lo) / np.asarray([nx, ny, nz], dtype=float)
    x_coords = lo[0] + (np.arange(nx) + 0.5) * spacing[0]
    y_coords = lo[1] + (np.arange(ny) + 0.5) * spacing[1]
    z_coords = lo[2] + (np.arange(nz) + 0.5) * spacing[2]

    iy_needed = set()
    iy0, iy1, _ = bracket(y_coords, NOZZLE_12_PLANE_Y)
    iy_needed.update([iy0, iy1])
    z_targets = [BODY_AXIS[1]] + [axis[2] for axis in NOZZLE_AXES]
    iz_needed = set()
    for z_value in z_targets:
        iz0, iz1, _ = bracket(z_coords, z_value)
        iz_needed.update([iz0, iz1])

    requested = [
        "Density", "pressure", "x_velocity", "y_velocity", "z_velocity", "sld"
    ]
    print("Reading selected AMR slices from {}".format(args.plotfile), flush=True)
    slices, xc, yc, zc, dx, dy, dz, meta = reader.read_amrex_3d_multi_slices_multi_ix(
        args.plotfile,
        requested,
        ix_list=[],
        iy_list=sorted(iy_needed),
        iz_list=sorted(iz_needed),
        out_level=-1,
        prolong_method="linear",
    )
    print("Loaded grid {}x{}x{} at dx={:.8g} m".format(meta["Nx"], meta["Ny"], meta["Nz"], dx),
          flush=True)

    # z planes contain the body/nozzle-0 axial view and all requested x-lines;
    # y planes provide an axial view through nozzles 1 and 2.
    z_fields = {name: slices[name][2] for name in requested}
    y_fields = {name: slices[name][1] for name in requested}
    centre_xy = {
        name: interp_z_plane(z_fields[name], zc, BODY_AXIS[1]) for name in requested
    }
    nozzle12_xz = {
        name: interp_z_plane(y_fields[name], yc, NOZZLE_12_PLANE_Y) for name in requested
    }

    valid_xy = fluid_mask(centre_xy["sld"], centre_xy["Density"], centre_xy["pressure"])
    mach_xy = mach_from_primitives(
        centre_xy["Density"], centre_xy["pressure"], centre_xy["x_velocity"],
        centre_xy["y_velocity"], centre_xy["z_velocity"],
    )
    mach_xy[~valid_xy] = np.nan
    pressure_xy = centre_xy["pressure"].copy()
    pressure_xy[~valid_xy] = np.nan

    valid_xz = fluid_mask(
        nozzle12_xz["sld"], nozzle12_xz["Density"], nozzle12_xz["pressure"]
    )
    mach_xz = mach_from_primitives(
        nozzle12_xz["Density"], nozzle12_xz["pressure"], nozzle12_xz["x_velocity"],
        nozzle12_xz["y_velocity"], nozzle12_xz["z_velocity"],
    )
    mach_xz[~valid_xz] = np.nan
    pressure_xz = nozzle12_xz["pressure"].copy()
    pressure_xz[~valid_xz] = np.nan

    x_sel = (xc >= args.x_min) & (xc <= args.x_max)
    y_sel = (yc >= args.y_min) & (yc <= args.y_max)
    z_sel = (zc >= args.y_min) & (zc <= args.y_max)
    xp = xc[x_sel]
    yp = yc[y_sel]
    zp = zc[z_sel]
    crop_xy = np.ix_(x_sel, y_sel)
    crop_xz = np.ix_(x_sel, z_sel)
    mach_plot_xy = mach_xy[crop_xy].T
    pressure_plot_xy = pressure_xy[crop_xy].T
    valid_plot_xy = valid_xy[crop_xy].T
    solid_plot_xy = (centre_xy["sld"][crop_xy].T > 0.5)
    schlieren_plot_xy, schlieren_scale_xy = numerical_schlieren(
        centre_xy["Density"][crop_xy].T, valid_plot_xy, dy, dx, args.schlieren_k
    )
    mach_plot_xz = mach_xz[crop_xz].T
    pressure_plot_xz = pressure_xz[crop_xz].T
    valid_plot_xz = valid_xz[crop_xz].T
    solid_plot_xz = (nozzle12_xz["sld"][crop_xz].T > 0.5)
    schlieren_plot_xz, schlieren_scale_xz = numerical_schlieren(
        nozzle12_xz["Density"][crop_xz].T, valid_plot_xz, dz, dx, args.schlieren_k
    )

    sim_ms = float(meta["sim_time"]) * 1.0e3
    plane_z = BODY_AXIS[1]
    plot_name = os.path.basename(os.path.normpath(args.plotfile))
    prefix = "srp_tri_L2_{}".format(plot_name)
    title_xy = "{}; t={:.3f} ms; z={:.3f} m".format(prefix, sim_ms, plane_z)
    title_xz = "{}; t={:.3f} ms; y={:.3f} m".format(
        prefix, sim_ms, NOZZLE_12_PLANE_Y
    )

    field_paths = {
        "mach_xy": os.path.join(args.output_dir, prefix + "_mach_xy_z0p275.png"),
        "pressure_xy": os.path.join(args.output_dir, prefix + "_pressure_xy_z0p275.png"),
        "schlieren_xy": os.path.join(args.output_dir, prefix + "_schlieren_xy_z0p275.png"),
        "mach_xz": os.path.join(args.output_dir, prefix + "_mach_xz_y0p259.png"),
        "pressure_xz": os.path.join(args.output_dir, prefix + "_pressure_xz_y0p259.png"),
        "schlieren_xz": os.path.join(args.output_dir, prefix + "_schlieren_xz_y0p259.png"),
    }
    save_field_figure(
        field_paths["mach_xy"], xp, yp, mach_plot_xy, solid_plot_xy,
        "Mach number — " + title_xy, "xy", "turbo", vmin=0.0, vmax=args.mach_max,
        colorbar_label="Mach number",
    )
    save_field_figure(
        field_paths["pressure_xy"], xp, yp, pressure_plot_xy, solid_plot_xy,
        "Pressure — " + title_xy, "xy", "magma",
        norm=LogNorm(vmin=args.pressure_min, vmax=args.pressure_max),
        colorbar_label="Pressure [Pa] (log scale)",
    )
    save_field_figure(
        field_paths["schlieren_xy"], xp, yp, schlieren_plot_xy, solid_plot_xy,
        "Numerical schlieren — " + title_xy, "xy", "gray", vmin=0.0, vmax=1.0,
        colorbar_label="exp(-k |grad rho| / P99.5), k={:g}".format(args.schlieren_k),
    )
    save_field_figure(
        field_paths["mach_xz"], xp, zp, mach_plot_xz, solid_plot_xz,
        "Mach number — " + title_xz, "xz", "turbo", vmin=0.0, vmax=args.mach_max,
        colorbar_label="Mach number",
    )
    save_field_figure(
        field_paths["pressure_xz"], xp, zp, pressure_plot_xz, solid_plot_xz,
        "Pressure — " + title_xz, "xz", "magma",
        norm=LogNorm(vmin=args.pressure_min, vmax=args.pressure_max),
        colorbar_label="Pressure [Pa] (log scale)",
    )
    save_field_figure(
        field_paths["schlieren_xz"], xp, zp, schlieren_plot_xz, solid_plot_xz,
        "Numerical schlieren — " + title_xz, "xz", "gray", vmin=0.0, vmax=1.0,
        colorbar_label="exp(-k |grad rho| / P99.5), k={:g}".format(args.schlieren_k),
    )

    # One two-plane composite for quick inspection of all three nozzles.
    fig, axes = plt.subplots(3, 2, figsize=(18.0, 14.4), constrained_layout=True)
    plane_specs = [
        (yp, "xy", solid_plot_xy, mach_plot_xy, pressure_plot_xy, schlieren_plot_xy,
         "XY: z=0.275 m (body axis + nozzle 0)"),
        (zp, "xz", solid_plot_xz, mach_plot_xz, pressure_plot_xz, schlieren_plot_xz,
         "XZ: y=0.259 m (nozzles 1 + 2)"),
    ]
    field_specs = [
        ("Mach", "turbo", None, 0.0, args.mach_max, "Mach number"),
        ("Pressure", "magma", LogNorm(args.pressure_min, args.pressure_max), None, None,
         "Pressure [Pa] (log)"),
        ("Schlieren", "gray", None, 0.0, 1.0, "Numerical schlieren"),
    ]
    for col, (transverse, plane, solid, mach, pressure, schlieren, plane_title) in enumerate(plane_specs):
        data_by_row = [mach, pressure, schlieren]
        for row, (panel_name, cmap, norm, vmin, vmax, cblabel) in enumerate(field_specs):
            ax = axes[row, col]
            pcm = ax.pcolormesh(xp, transverse, data_by_row[row], shading="auto", cmap=cmap,
                                norm=norm, vmin=vmin, vmax=vmax)
            add_solid_overlay(ax, xp, transverse, solid)
            decorate_axial_plane(ax, plane)
            ax.set_title("{} — {}".format(panel_name, plane_title))
            fig.colorbar(pcm, ax=ax, pad=0.012).set_label(cblabel)
        axes[0, col].legend(loc="lower right", framealpha=0.85, fontsize=7)
    fig.suptitle("SRP tri-nozzle latest flow field — {}, t={:.3f} ms".format(
        plot_name, sim_ms), fontsize=14)
    composite_path = os.path.join(args.output_dir, prefix + "_flowfield_composite.png")
    fig.savefig(composite_path, dpi=210)
    plt.close(fig)

    profiles = {
        "body_axis": axis_profile(z_fields, yc, zc, BODY_AXIS[0], BODY_AXIS[1])
    }
    for name, y_value, z_value in NOZZLE_AXES:
        profiles[name] = axis_profile(z_fields, yc, zc, y_value, z_value)
    for profile in profiles.values():
        profile["density_gradient"] = segmented_abs_gradient(profile["density"], xc)

    colors = {
        "body_axis": "#111111", "nozzle_0": "#d62728",
        "nozzle_1": "#1f77b4", "nozzle_2": "#2ca02c",
    }
    labels = {
        "body_axis": "body axis (y=z=0.275)",
        "nozzle_0": "nozzle 0 (y=0.307, z=0.275)",
        "nozzle_1": "nozzle 1 (y=0.259, z=0.302)",
        "nozzle_2": "nozzle 2 (y=0.259, z=0.248)",
    }
    fig, axes = plt.subplots(3, 1, figsize=(12.5, 12.0), sharex=True, constrained_layout=True)
    for name, profile in profiles.items():
        style = "-" if name == "body_axis" else "--"
        width = 1.8 if name == "body_axis" else 1.25
        axes[0].plot(xc, profile["mach"], style, color=colors[name], linewidth=width, label=labels[name])
        axes[1].plot(xc, profile["pressure"], style, color=colors[name], linewidth=width)
        axes[2].plot(xc, profile["density_gradient"], style, color=colors[name], linewidth=width)
    for ax in axes:
        ax.axvline(NOZZLE_X, color="#666666", linestyle=":", linewidth=1.0,
                   label="nozzle plane" if ax is axes[0] else None)
        ax.set_xlim(args.x_min, args.x_max)
        ax.grid(True, which="both", alpha=0.22)
    axes[0].set_ylabel("Mach number")
    axes[0].set_ylim(bottom=0.0)
    axes[0].legend(loc="best", fontsize=8, ncol=2)
    axes[1].set_ylabel("Pressure [Pa]")
    axes[1].set_yscale("log")
    axes[1].set_ylim(args.pressure_min, args.pressure_max)
    axes[2].set_ylabel("|d rho / dx| [kg m^-4]\n(axis schlieren proxy)")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("x [m]")
    axes[2].set_ylim(bottom=1.0e-3)
    fig.suptitle("SRP tri-nozzle axis profiles — {}, t={:.3f} ms".format(
        plot_name, sim_ms), fontsize=14)
    axis_path = os.path.join(args.output_dir, prefix + "_axis_profiles.png")
    fig.savefig(axis_path, dpi=220)
    plt.close(fig)

    csv_path = os.path.join(args.output_dir, prefix + "_axis_profiles.csv")
    columns = ["x_m"]
    for name in profiles:
        columns.extend([
            name + "_mach", name + "_pressure_Pa", name + "_density_kgm3",
            name + "_abs_drhodx_kgm4", name + "_sld",
        ])
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for i, x_value in enumerate(xc):
            row = ["{:.12g}".format(x_value)]
            for name, profile in profiles.items():
                row.extend([
                    "{:.12g}".format(profile["mach"][i]),
                    "{:.12g}".format(profile["pressure"][i]),
                    "{:.12g}".format(profile["density"][i]),
                    "{:.12g}".format(profile["density_gradient"][i]),
                    "{:.12g}".format(profile["sld"][i]),
                ])
            writer.writerow(row)

    summary = {
        "plotfile": os.path.abspath(args.plotfile),
        "simulation_time_s": float(meta["sim_time"]),
        "simulation_time_ms": sim_ms,
        "grid": [int(meta["Nx"]), int(meta["Ny"]), int(meta["Nz"])],
        "spacing_m": [float(dx), float(dy), float(dz)],
        "centre_plane_z_m": plane_z,
        "nozzle_12_plane_y_m": NOZZLE_12_PLANE_Y,
        "body_axis_yz_m": list(BODY_AXIS),
        "nozzle_axes_yz_m": {name: [y, z] for name, y, z in NOZZLE_AXES},
        "nozzle_plane_x_m": NOZZLE_X,
        "field_ranges_fluid_crop": {
            "xy_z0p275": {
                "mach_min": float(np.nanmin(mach_plot_xy)),
                "mach_max": float(np.nanmax(mach_plot_xy)),
                "pressure_min_Pa": float(np.nanmin(pressure_plot_xy)),
                "pressure_max_Pa": float(np.nanmax(pressure_plot_xy)),
            },
            "xz_y0p259": {
                "mach_min": float(np.nanmin(mach_plot_xz)),
                "mach_max": float(np.nanmax(mach_plot_xz)),
                "pressure_min_Pa": float(np.nanmin(pressure_plot_xz)),
                "pressure_max_Pa": float(np.nanmax(pressure_plot_xz)),
            },
        },
        "schlieren": {
            "definition": "exp(-k*abs(grad(Density))/P99.5_fluid)",
            "k": float(args.schlieren_k),
            "xy_gradient_reference_kg_m4": schlieren_scale_xy,
            "xz_gradient_reference_kg_m4": schlieren_scale_xz,
        },
        "outputs": {
            "composite": composite_path,
            "axis_profiles_png": axis_path,
            "axis_profiles_csv": csv_path,
            "fields": field_paths,
        },
    }
    summary_path = os.path.join(args.output_dir, prefix + "_summary.json")
    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print("Saved outputs in {}".format(args.output_dir), flush=True)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
