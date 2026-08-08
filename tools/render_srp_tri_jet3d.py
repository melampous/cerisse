#!/usr/bin/env python3
"""Render a true 3-D oblique view of the SRP tri-nozzle jet state.

The current plotfile has no passive jet tracer.  This script therefore uses
local isentropic total pressure as a source proxy, retains only components
connected to the three nozzle seeds, and overlays instantaneous streamlines
seeded separately from the upstream boundary and the three nozzles.

Outputs include a high-resolution two-view PNG, a rotating GIF, surface mesh
data, and a JSON file recording every threshold and classification caveat.
"""

import argparse
import glob
import json
import os
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from PIL import Image
from scipy import ndimage
import yt


GAMMA = 1.4
P_INF = 574.56
RHO_INF = 0.030799249530956845
U_INF = 743.3947269116185
P0_INF = 188222.90041860205
P0_JET = 4437901.4399999995
NOZZLE_X = 0.299
NOZZLE_RADIUS = 0.0028
NOZZLES = (
    ("jet 0", 0.307, 0.275, "#ff3b30"),
    ("jet 1", 0.259, 0.302, "#ff9500"),
    ("jet 2", 0.259, 0.248, "#ffd60a"),
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plotfile", required=True)
    parser.add_argument("--stl", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--level", type=int, default=2)
    parser.add_argument("--envelope-ratio", type=float, default=2.0)
    parser.add_argument("--core-ratio", type=float, default=5.0)
    parser.add_argument("--gif-frames", type=int, default=18)
    return parser.parse_args()


def read_ascii_stl(path):
    vertices = []
    with open(path, "r", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if text.startswith("vertex"):
                vertices.append([float(value) for value in text.split()[1:4]])
    array = np.asarray(vertices, dtype=np.float32)
    if array.size == 0 or array.shape[0] % 3:
        raise RuntimeError("Could not parse ASCII STL: {}".format(path))
    return array.reshape((-1, 3, 3))


def load_field(container, name):
    print("Loading {}".format(name), flush=True)
    array = np.asarray(container[("boxlib", name)].d, dtype=np.float32).copy()
    container.clear_data()
    return array


def seed_component_labels(labels, le, spacing):
    selected = []
    seed_report = {}
    nx, ny, nz = labels.shape
    for name, y0, z0, _ in NOZZLES:
        # Search a short cylinder immediately upstream of the nozzle face.
        x0 = NOZZLE_X - 0.003
        ic = int(round((x0 - le[0]) / spacing - 0.5))
        jc = int(round((y0 - le[1]) / spacing - 0.5))
        kc = int(round((z0 - le[2]) / spacing - 0.5))
        radius_cells = max(3, int(round(1.5 * NOZZLE_RADIUS / spacing)))
        ilo, ihi = max(0, ic - 4), min(nx, ic + 5)
        jlo, jhi = max(0, jc - radius_cells), min(ny, jc + radius_cells + 1)
        klo, khi = max(0, kc - radius_cells), min(nz, kc + radius_cells + 1)
        block = labels[ilo:ihi, jlo:jhi, klo:khi]
        values, counts = np.unique(block[block > 0], return_counts=True)
        if values.size:
            label = int(values[np.argmax(counts)])
            selected.append(label)
            seed_report[name] = label
        else:
            seed_report[name] = None
    return sorted(set(selected)), seed_report


def downsample_binary(mask):
    nx, ny, nz = (dimension // 2 for dimension in mask.shape)
    trimmed = mask[: 2 * nx, : 2 * ny, : 2 * nz]
    reduced = trimmed.reshape(nx, 2, ny, 2, nz, 2).mean(axis=(1, 3, 5))
    reduced = ndimage.gaussian_filter(reduced.astype(np.float32), sigma=0.70, mode="constant")
    reduced[0, :, :] = 0.0
    reduced[-1, :, :] = 0.0
    reduced[:, 0, :] = 0.0
    reduced[:, -1, :] = 0.0
    reduced[:, :, 0] = 0.0
    reduced[:, :, -1] = 0.0
    return reduced


def extract_surface(scalar, lo, hi, value=0.5):
    bbox = np.asarray([[lo[0], hi[0]], [lo[1], hi[1]], [lo[2], hi[2]]], dtype=float)
    uniform = yt.load_uniform_grid(
        {"indicator": scalar}, scalar.shape, bbox=bbox, nprocs=1,
        periodicity=(True, True, True),
    )
    surface = uniform.surface(uniform.all_data(), ("stream", "indicator"), value)
    triangles = np.asarray(surface.triangles, dtype=np.float32).copy()
    del surface, uniform
    return triangles


def sample_scalar(array, point, le, spacing):
    index = (np.asarray(point, dtype=float) - le) / spacing - 0.5
    return float(ndimage.map_coordinates(
        array, index[:, None], order=1, mode="constant", cval=np.nan,
        prefilter=False,
    )[0])


def sample_velocity(arrays, point, le, spacing):
    return np.asarray([
        sample_scalar(array, point, le, spacing) for array in arrays
    ], dtype=float)


def trace_streamline(seed, velocity, solid, ghost, le, hi, spacing,
                     step_length, max_steps):
    point = np.asarray(seed, dtype=float)
    points = [point.copy()]
    margin = 1.5 * spacing
    for _ in range(max_steps):
        if np.any(point <= le + margin) or np.any(point >= hi - margin):
            break
        if sample_scalar(solid, point, le, spacing) > 0.45:
            break
        if sample_scalar(ghost, point, le, spacing) > 0.45:
            break
        vector = sample_velocity(velocity, point, le, spacing)
        speed = float(np.linalg.norm(vector))
        if not np.all(np.isfinite(vector)) or speed < 5.0:
            break
        direction = vector / speed
        midpoint = point + 0.5 * step_length * direction
        vector_mid = sample_velocity(velocity, midpoint, le, spacing)
        speed_mid = float(np.linalg.norm(vector_mid))
        if not np.all(np.isfinite(vector_mid)) or speed_mid < 5.0:
            break
        point = point + step_length * vector_mid / speed_mid
        points.append(point.copy())
    return np.asarray(points)


def build_streamlines(velocity, solid, ghost, le, hi, spacing):
    jet_paths = []
    ring_angles = np.linspace(0.0, 2.0 * np.pi, 7, endpoint=False)
    disk_offsets = [(0.0, 0.0)] + [
        (0.00145 * np.cos(angle), 0.00145 * np.sin(angle)) for angle in ring_angles
    ]
    for nozzle_index, (name, y0, z0, color) in enumerate(NOZZLES):
        for dy, dz in disk_offsets:
            seed = (NOZZLE_X - 0.0032, y0 + dy, z0 + dz)
            path = trace_streamline(
                seed, velocity, solid, ghost, le, hi, spacing,
                step_length=0.00075, max_steps=360,
            )
            if len(path) >= 8:
                jet_paths.append((nozzle_index, color, path))

    free_paths = []
    radii = (0.090, 0.125)
    for radius in radii:
        for angle in np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False):
            y0 = 0.275 + radius * np.cos(angle)
            z0 = 0.275 + radius * np.sin(angle)
            seed = (le[0] + 0.0025, y0, z0)
            if not (le[1] < y0 < hi[1] and le[2] < z0 < hi[2]):
                continue
            path = trace_streamline(
                seed, velocity, solid, ghost, le, hi, spacing,
                step_length=0.0015, max_steps=220,
            )
            if len(path) >= 8:
                free_paths.append(path)
    return jet_paths, free_paths


def surface_collection(triangles, color, alpha, linewidth=0.0):
    collection = Poly3DCollection(
        triangles, facecolors=color, edgecolors="none", linewidths=linewidth,
        alpha=alpha, zsort="average",
    )
    return collection


def shaded_body_collection(triangles):
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    normals /= np.linalg.norm(normals, axis=1, keepdims=True) + 1.0e-12
    light = np.asarray([-0.55, 0.35, 0.75])
    light /= np.linalg.norm(light)
    shade = 0.28 + 0.58 * np.clip(normals @ light, 0.0, 1.0)
    colors = np.stack([shade, shade, shade, np.full_like(shade, 0.96)], axis=1)
    return Poly3DCollection(
        triangles, facecolors=colors, edgecolors="none", linewidths=0.0,
        alpha=0.96, zsort="average",
    )


def decorate_scene(ax, body_triangles, envelope_triangles, core_triangles,
                   jet_paths, free_paths, view):
    ax.add_collection3d(shaded_body_collection(body_triangles))
    ax.add_collection3d(surface_collection(envelope_triangles, "#ff9f0a", 0.20))
    ax.add_collection3d(surface_collection(core_triangles, "#ff2d20", 0.48))

    for _, color, path in jet_paths:
        ax.plot(path[:, 0], path[:, 1], path[:, 2], color=color,
                linewidth=1.55, alpha=0.92)
    for path in free_paths:
        ax.plot(path[:, 0], path[:, 1], path[:, 2], color="#1887ff",
                linewidth=0.85, alpha=0.52)

    theta = np.linspace(0.0, 2.0 * np.pi, 80)
    for index, (_, y0, z0, color) in enumerate(NOZZLES):
        ax.plot(
            np.full(theta.shape, NOZZLE_X + 0.00035),
            y0 + NOZZLE_RADIUS * np.cos(theta),
            z0 + NOZZLE_RADIUS * np.sin(theta),
            color=color, linewidth=1.8, alpha=0.95,
        )
        ax.text(NOZZLE_X + 0.002, y0, z0, str(index), color=color,
                fontsize=8, weight="bold")

    ax.quiver(0.075, 0.397, 0.397, 0.045, 0.0, 0.0, color="#1887ff",
              linewidth=1.7, arrow_length_ratio=0.22)
    ax.text(0.073, 0.405, 0.405, "freestream +x", color="#0969c8", fontsize=8)
    ax.quiver(0.294, 0.300, 0.368, -0.040, 0.0, 0.0, color="#ff3b30",
              linewidth=1.7, arrow_length_ratio=0.22)
    ax.text(0.250, 0.307, 0.374, "retro-jets -x", color="#d92519", fontsize=8)

    ax.set_xlim(0.055, 0.52)
    ax.set_ylim(0.14, 0.41)
    ax.set_zlim(0.14, 0.41)
    ax.set_box_aspect((0.465, 0.27, 0.27))
    ax.set_xlabel("x [m]", labelpad=0)
    ax.set_ylabel("y [m]", labelpad=0)
    ax.set_zlabel("z [m]", labelpad=0)
    ax.tick_params(labelsize=7, pad=-1)
    ax.grid(False)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.set_alpha(0.0)
    try:
        ax.set_proj_type("persp", focal_length=0.80)
    except TypeError:
        ax.set_proj_type("persp")
    ax.view_init(elev=view[0], azim=view[1])


def legend_handles(envelope_ratio, core_ratio):
    return [
        Patch(facecolor="#8a8a8a", label="aeroshell / nozzle STL"),
        Patch(facecolor="#ff9f0a", alpha=0.35,
              label="jet envelope: connected p0 > {:.0f} x p0_inf".format(envelope_ratio)),
        Patch(facecolor="#ff2d20", alpha=0.65,
              label="high-p0 jet core: p0 > {:.0f} x p0_inf".format(core_ratio)),
        Line2D([0], [0], color="#ff3b30", lw=2.0, label="nozzle-seeded jet streamlines"),
        Line2D([0], [0], color="#1887ff", lw=1.5, label="upstream-seeded freestream lines"),
    ]


def render_static(path, body, envelope, core, jet_paths, free_paths,
                  sim_ms, envelope_ratio, core_ratio):
    fig = plt.figure(figsize=(17.5, 8.2), constrained_layout=True)
    views = ((22, -145, "upstream oblique"), (20, -55, "side oblique"))
    for panel, (elevation, azimuth, title) in enumerate(views, start=1):
        ax = fig.add_subplot(1, 2, panel, projection="3d")
        decorate_scene(
            ax, body, envelope, core, jet_paths, free_paths,
            (elevation, azimuth),
        )
        ax.set_title(title, fontsize=11)
    fig.legend(handles=legend_handles(envelope_ratio, core_ratio), loc="lower center",
               ncol=3, framealpha=0.95, fontsize=9)
    fig.suptitle(
        "SRP tri-nozzle 3-D jet state — plt09000, t={:.3f} ms".format(sim_ms),
        fontsize=15,
    )
    fig.text(
        0.5, 0.018,
        "Jet surfaces use a nozzle-connected total-pressure proxy; strict source tracking requires a passive scalar.",
        ha="center", fontsize=8.5, color="#555555",
    )
    fig.savefig(path, dpi=210, facecolor="white")
    plt.close(fig)


def render_turntable(path, body, envelope, core, jet_paths, free_paths,
                     sim_ms, envelope_ratio, core_ratio, frame_count):
    with tempfile.TemporaryDirectory(prefix="srp_jet3d_") as directory:
        filenames = []
        fig = plt.figure(figsize=(11.5, 7.2))
        ax = fig.add_subplot(111, projection="3d")
        decorate_scene(ax, body, envelope, core, jet_paths, free_paths, (21, -155))
        ax.set_title("SRP tri-nozzle 3-D jet state — t={:.3f} ms".format(sim_ms), fontsize=12)
        fig.legend(handles=legend_handles(envelope_ratio, core_ratio), loc="lower center",
                   ncol=2, framealpha=0.93, fontsize=8)
        fig.subplots_adjust(left=0.0, right=1.0, top=0.93, bottom=0.13)
        for frame, azimuth in enumerate(np.linspace(-170.0, 190.0, frame_count, endpoint=False)):
            ax.view_init(elev=21.0, azim=float(azimuth))
            filename = os.path.join(directory, "frame_{:03d}.png".format(frame))
            fig.savefig(filename, dpi=105, facecolor="white")
            filenames.append(filename)
        plt.close(fig)
        images = [Image.open(filename).convert("P", palette=Image.ADAPTIVE) for filename in filenames]
        images[0].save(
            path, save_all=True, append_images=images[1:], duration=145,
            loop=0, optimize=False, disposal=2,
        )
        for image in images:
            image.close()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    yt.set_log_level(50)

    dataset = yt.load(args.plotfile)
    sim_ms = float(dataset.current_time) * 1.0e3
    level = min(args.level, int(dataset.max_level))
    domain_width = np.asarray(dataset.domain_right_edge.d - dataset.domain_left_edge.d)
    base_dims = np.asarray(dataset.domain_dimensions, dtype=int)
    refinement = 2 ** level
    spacing = float(domain_width[0] / (base_dims[0] * refinement))

    # Every bound is aligned with an L2 cell edge for the 192^3, r=2 case.
    le = np.asarray([0.060, 0.140, 0.140], dtype=float)
    requested_hi = np.asarray([0.305, 0.410, 0.410], dtype=float)
    dims = np.rint((requested_hi - le) / spacing).astype(int)
    hi = le + dims * spacing
    print("ROI level={} spacing={} dims={} lo={} hi={}".format(
        level, spacing, dims.tolist(), le.tolist(), hi.tolist()), flush=True)

    grid = dataset.covering_grid(
        level=level, left_edge=dataset.arr(le, "code_length"), dims=dims,
    )
    rho = load_field(grid, "Density")
    pressure = load_field(grid, "pressure")
    ux = load_field(grid, "x_velocity")
    uy = load_field(grid, "y_velocity")
    uz = load_field(grid, "z_velocity")
    solid = load_field(grid, "sld")
    ghost = load_field(grid, "ghs")
    del grid

    fluid = (
        (solid < 0.5) & (ghost < 0.5) & np.isfinite(rho) & np.isfinite(pressure)
        & (rho > 0.0) & (pressure > 0.0)
    )
    speed2 = ux * ux + uy * uy + uz * uz
    mach2 = speed2 * rho / np.maximum(GAMMA * pressure, 1.0e-30)
    total_pressure = pressure * np.power(1.0 + 0.2 * mach2, 3.5)
    total_pressure = total_pressure.astype(np.float32)

    pristine_free = (
        fluid & (ux > 0.95 * U_INF)
        & (np.abs(pressure / P_INF - 1.0) < 0.05)
        & (np.abs(rho / RHO_INF - 1.0) < 0.05)
    )
    recirculation = fluid & (ux < -50.0) & (total_pressure < 1.25 * P0_INF)

    envelope_threshold = args.envelope_ratio * P0_INF
    core_threshold = args.core_ratio * P0_INF
    candidate = fluid & (total_pressure > envelope_threshold)
    structure = np.ones((3, 3, 3), dtype=np.uint8)
    labels, component_count = ndimage.label(candidate, structure=structure)
    selected_labels, seed_report = seed_component_labels(labels, le, spacing)
    missing_seeds = [name for name, value in seed_report.items() if value is None]
    if missing_seeds:
        raise RuntimeError(
            "No connected high-p0 component found at nozzle seed(s): {}".format(
                ", ".join(missing_seeds)
            )
        )
    envelope = np.isin(labels, selected_labels)
    core = envelope & (total_pressure > core_threshold)
    retro_core = envelope & (ux < -150.0)
    del labels, candidate, mach2, speed2

    if not np.any(envelope) or not np.any(core):
        raise RuntimeError("Jet envelope/core mask is empty; check thresholds and seed locations")
    boundary_contacts = {
        "x_lo": int(np.count_nonzero(envelope[0, :, :])),
        "x_hi": int(np.count_nonzero(envelope[-1, :, :])),
        "y_lo": int(np.count_nonzero(envelope[:, 0, :])),
        "y_hi": int(np.count_nonzero(envelope[:, -1, :])),
        "z_lo": int(np.count_nonzero(envelope[:, :, 0])),
        "z_hi": int(np.count_nonzero(envelope[:, :, -1])),
    }
    touched = [name for name, count in boundary_contacts.items() if count]
    if touched:
        raise RuntimeError(
            "Jet envelope touches ROI boundary {}; expand ROI before rendering".format(
                ", ".join(touched)
            )
        )

    jet_paths, free_paths = build_streamlines(
        (ux, uy, uz), solid, ghost, le, hi, spacing
    )

    envelope_scalar = downsample_binary(envelope)
    core_scalar = downsample_binary(core)
    surface_hi = le + 2.0 * spacing * np.asarray(envelope_scalar.shape)
    print("Extracting jet-envelope surface", flush=True)
    envelope_triangles = extract_surface(envelope_scalar, le, surface_hi, 0.5)
    print("Extracting high-p0 core surface", flush=True)
    core_triangles = extract_surface(core_scalar, le, surface_hi, 0.5)
    body_triangles = read_ascii_stl(args.stl)

    prefix = "srp_tri_L2_plt09000_jet3d"
    static_path = os.path.join(args.output_dir, prefix + "_oblique.png")
    gif_path = os.path.join(args.output_dir, prefix + "_turntable.gif")
    mesh_path = os.path.join(args.output_dir, prefix + "_surfaces.npz")
    summary_path = os.path.join(args.output_dir, prefix + "_summary.json")

    print("Rendering oblique views", flush=True)
    render_static(
        static_path, body_triangles, envelope_triangles, core_triangles,
        jet_paths, free_paths, sim_ms, args.envelope_ratio, args.core_ratio,
    )
    print("Rendering rotating GIF", flush=True)
    gif_frames = args.gif_frames
    if len(envelope_triangles) + len(core_triangles) > 200000:
        gif_frames = min(gif_frames, 12)
    render_turntable(
        gif_path, body_triangles, envelope_triangles, core_triangles,
        jet_paths, free_paths, sim_ms, args.envelope_ratio, args.core_ratio,
        gif_frames,
    )
    np.savez_compressed(
        mesh_path,
        envelope_triangles_m=envelope_triangles,
        core_triangles_m=core_triangles,
    )

    cell_volume = spacing ** 3
    envelope_indices = np.argwhere(envelope)
    core_indices = np.argwhere(core)
    summary = {
        "plotfile": os.path.abspath(args.plotfile),
        "simulation_time_ms": sim_ms,
        "classification_is_exact_tracer": False,
        "classification_warning": (
            "The plotfile has no passive scalar. Jet surfaces are a nozzle-connected "
            "local-total-pressure proxy; seeded curves are instantaneous streamlines."
        ),
        "formulae": {
            "mach": "|u|/sqrt(1.4*pressure/Density)",
            "local_total_pressure_Pa": "pressure*(1+0.2*Mach^2)^3.5",
        },
        "reference_states": {
            "freestream_p0_Pa": P0_INF,
            "jet_p0_Pa": P0_JET,
            "freestream_u_m_s": U_INF,
        },
        "thresholds": {
            "jet_envelope_p0_Pa": envelope_threshold,
            "jet_core_p0_Pa": core_threshold,
            "retro_core_ux_m_s": -150.0,
            "solid_sld_max": 0.5,
            "ghost_ghs_max": 0.5,
        },
        "roi": {
            "level": level,
            "lo_m": le.tolist(),
            "hi_m": hi.tolist(),
            "dimensions": dims.tolist(),
            "spacing_m": spacing,
        },
        "connected_components": {
            "candidate_count": int(component_count),
            "selected_labels": selected_labels,
            "nozzle_seed_labels": seed_report,
            "jet_envelope_roi_boundary_contacts": boundary_contacts,
        },
        "classification_counts": {
            "fluid_cells": int(np.count_nonzero(fluid)),
            "jet_envelope_cells": int(np.count_nonzero(envelope)),
            "jet_core_cells": int(np.count_nonzero(core)),
            "retro_core_cells": int(np.count_nonzero(retro_core)),
            "pristine_freestream_cells": int(np.count_nonzero(pristine_free)),
            "low_p0_recirculation_cells": int(np.count_nonzero(recirculation)),
        },
        "jet_geometry": {
            "envelope_volume_m3": float(np.count_nonzero(envelope) * cell_volume),
            "core_volume_m3": float(np.count_nonzero(core) * cell_volume),
            "envelope_bbox_m": [
                (le + envelope_indices.min(axis=0) * spacing).tolist(),
                (le + (envelope_indices.max(axis=0) + 1) * spacing).tolist(),
            ],
            "core_bbox_m": [
                (le + core_indices.min(axis=0) * spacing).tolist(),
                (le + (core_indices.max(axis=0) + 1) * spacing).tolist(),
            ],
        },
        "rendering": {
            "envelope_triangles": int(len(envelope_triangles)),
            "core_triangles": int(len(core_triangles)),
            "jet_streamlines": int(len(jet_paths)),
            "freestream_streamlines": int(len(free_paths)),
            "surface_grid_spacing_m": 2.0 * spacing,
            "surface_geometry_warning": (
                "The displayed mesh is a Gaussian-smoothed 2x downsample of the L2 mask; "
                "it is qualitative and must not be used for exact area measurements."
            ),
            "gif_frames": int(gif_frames),
        },
        "outputs": {
            "oblique_png": static_path,
            "turntable_gif": gif_path,
            "surface_mesh_npz": mesh_path,
        },
    }
    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
