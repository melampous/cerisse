#!/usr/bin/env python3
"""Plot Mach/pressure and axis-proximal profiles from a 2-D AMReX RZ plotfile."""

from __future__ import annotations

import argparse
import inspect
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import yt
from yt.frontends.boxlib import data_structures as boxlib_data


def patch_yt_cylindrical_readonly_edge() -> None:
    """Work around yt 4.3 + recent NumPy returning a read-only edge view."""
    source = textwrap.dedent(
        inspect.getsource(boxlib_data.BoxlibDataset._parse_header_file)
    )
    source = source.replace(
        "dre = self.domain_right_edge\n",
        "dre = self.domain_right_edge.copy()\n",
    )
    namespace = dict(boxlib_data.__dict__)
    exec(source, namespace)
    boxlib_data.BoxlibDataset._parse_header_file = namespace[
        "_parse_header_file"
    ]


def load_finest(plotfile: Path):
    patch_yt_cylindrical_readonly_edge()
    dataset = yt.load(str(plotfile))
    level = int(dataset.index.max_level)
    dims = np.asarray(dataset.domain_dimensions, dtype=int).copy()
    dims[:2] *= 2**level
    dims[2] = 1
    grid = dataset.covering_grid(level, dataset.domain_left_edge, dims)

    def field(name: str) -> np.ndarray:
        return grid[("boxlib", name)].to_ndarray()[:, :, 0]

    return dataset, level, field


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--failure-r-mm", type=float)
    parser.add_argument("--failure-x-mm", type=float)
    parser.add_argument("--x-min-mm", type=float, default=-80.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    dataset, level, field = load_finest(args.plotfile)
    time_us = float(dataset.current_time) * 1.0e6
    time_tag = f"{time_us:.2f}us".replace(".", "p")
    density = field("Density")
    pressure = field("pressure")
    radial_velocity = field("x_velocity")
    axial_velocity = field("y_velocity")
    sound_speed = np.sqrt(1.4 * pressure / density)
    mach = np.hypot(radial_velocity, axial_velocity) / sound_speed

    nr, nz = density.shape
    r_lo = float(dataset.domain_left_edge[0])
    r_hi = float(dataset.domain_right_edge[0])
    z_lo = float(dataset.domain_left_edge[1])
    z_hi = float(dataset.domain_right_edge[1])
    r_edges = np.linspace(r_lo, r_hi, nr + 1) * 1.0e3
    z_edges = np.linspace(z_lo, z_hi, nz + 1) * 1.0e3
    r_centres = 0.5 * (r_edges[:-1] + r_edges[1:])
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])

    p_inf = 534.5806702256042
    pressure_ratio = pressure / p_inf
    axis_ring = 0
    mach_axis = mach[axis_ring, :]
    pressure_axis = pressure_ratio[axis_ring, :]
    axial_axis = axial_velocity[axis_ring, :]

    mach_disk_mm = np.nan
    for j in range(nz - 2, 0, -1):
        if (
            axial_axis[j] < 0.0
            and mach_axis[j] >= 1.0
            and mach_axis[j - 1] < 1.0
            and axial_axis[j - 1] < 0.0
        ):
            fraction = (mach_axis[j] - 1.0) / (
                mach_axis[j] - mach_axis[j - 1]
            )
            mach_disk_mm = z_centres[j] + fraction * (
                z_centres[j - 1] - z_centres[j]
            )
            break

    zoom_extent = (args.x_min_mm, z_hi * 1.0e3, 0.0, 45.0)
    fig, axes = plt.subplots(2, 1, figsize=(11.2, 8.4), constrained_layout=True)
    mach_map = axes[0].pcolormesh(
        z_edges, r_edges, mach, shading="auto", cmap="turbo", vmin=0.0, vmax=6.0
    )
    axes[0].contour(
        z_centres, r_centres, mach, levels=[1.0], colors="white", linewidths=0.8
    )
    axes[0].set(xlim=zoom_extent[:2], ylim=zoom_extent[2:], ylabel="r [mm]",
                title=f"Mach number, t={time_us:.1f} us, uniform L{level}")
    fig.colorbar(mach_map, ax=axes[0], label="Mach")

    pressure_map = axes[1].pcolormesh(
        z_edges,
        r_edges,
        pressure_ratio,
        shading="auto",
        cmap="magma",
        norm=LogNorm(vmin=0.5, vmax=250.0),
    )
    axes[1].set(xlim=zoom_extent[:2], ylim=zoom_extent[2:], xlabel="axial x [mm]",
                ylabel="r [mm]", title="Static-pressure ratio p/p_inf")
    fig.colorbar(pressure_map, ax=axes[1], label="p/p_inf")
    if np.isfinite(mach_disk_mm):
        for axis in axes:
            axis.axvline(mach_disk_mm, color="cyan", linestyle="--", linewidth=1.0)
            axis.text(
                mach_disk_mm - 1.0,
                42.0,
                f"axis M=1 at {mach_disk_mm:.1f} mm",
                color="cyan",
                ha="right",
                va="top",
                fontsize=8,
            )
    if args.failure_r_mm is not None and args.failure_x_mm is not None:
        for axis in axes:
            axis.scatter(
                [args.failure_x_mm], [args.failure_r_mm], marker="x",
                s=70, linewidths=1.8, color="cyan", zorder=5,
                label="cell failing in next update",
            )
        axes[0].legend(loc="upper right", fontsize=8)
    field_figure = args.output_dir / f"mach_pressure_{time_tag}.png"
    fig.savefig(field_figure, dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.2), sharex=True,
                             constrained_layout=True)
    axes[0].plot(z_centres, mach_axis, color="black", linewidth=1.1)
    axes[0].axhline(1.0, color="0.5", linestyle=":")
    axes[0].set(ylabel="Mach", ylim=(0.0, 6.2))
    axes[1].semilogy(z_centres, pressure_axis, color="tab:red", linewidth=1.1)
    axes[1].set(ylabel="p/p_inf", ylim=(0.3, 300.0))
    axes[2].plot(z_centres, axial_axis, color="tab:blue", linewidth=1.1)
    axes[2].axhline(0.0, color="0.5", linestyle=":")
    axes[2].set(
        xlabel="axial x [mm]",
        ylabel="u_x [m/s]",
        xlim=(args.x_min_mm, z_hi * 1.0e3),
    )
    if np.isfinite(mach_disk_mm):
        for axis in axes:
            axis.axvline(mach_disk_mm, color="tab:green", linestyle="--")
    fig.suptitle(f"Axis-proximal ring r={r_centres[axis_ring]:.3f} mm")
    axis_figure = args.output_dir / f"axis_profiles_{time_tag}.png"
    fig.savefig(axis_figure, dpi=220)
    plt.close(fig)

    diagnostics = args.output_dir / f"diagnostics_{time_tag}.txt"
    diagnostics.write_text(
        "\n".join(
            [
                f"time_s={float(dataset.current_time):.17g}",
                f"level={level}",
                f"finest_shape_rz={nr}x{nz}",
                f"dr_mm={(r_edges[1]-r_edges[0]):.12g}",
                f"dz_mm={(z_edges[1]-z_edges[0]):.12g}",
                f"axis_ring_r_mm={r_centres[axis_ring]:.12g}",
                f"axis_mach1_location_mm={mach_disk_mm:.12g}",
                f"mach_min={mach.min():.12g}",
                f"mach_max={mach.max():.12g}",
                f"pressure_min_Pa={pressure.min():.12g}",
                f"pressure_max_Pa={pressure.max():.12g}",
            ]
        )
        + "\n"
    )
    print(field_figure)
    print(axis_figure)
    print(diagnostics)


if __name__ == "__main__":
    main()
