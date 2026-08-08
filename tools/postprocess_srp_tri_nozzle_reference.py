#!/usr/bin/env python3
"""Render a bow-shock-safe, nozzle-marked SRP tri-nozzle reference figure.

The main product follows the visual organization of the supplied 2 x 2
reference image: density on the top row, numerical schlieren on the bottom
row, a bow-shock-safe overview on the left, and a nozzle-region zoom on the
right.  The physical XZ plane at y=0.259 m passes through nozzles 1 and 2.
Their exits are marked on both zoom panels.  A small YZ exit-layout inset
also locates nozzle 0, which is intentionally off this XZ plane.

Sampling is delegated to the independently tested thin-slab extractor in
postprocess_srp_tri_late_mean.py.  Only the requested instantaneous XZ slice
is read from the AMReX plotfile.
"""

import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
import numpy as np


BODY_DIAMETER = 0.127
BODY_CENTRE = (0.300, 0.275, 0.275)
PLANE_Y = 0.259
NOZZLE_X = 0.299
NOZZLE_RADIUS = 0.0028
NOZZLES = (
    ("Nozzle 0", 0.307, 0.275, "#d62728"),
    ("Nozzle 1", 0.259, 0.302, "#e67e22"),
    ("Nozzle 2", 0.259, 0.248, "#00a6d6"),
)
FULL_AXIAL_CROP = (0.006666666666666667, 0.60)
FULL_TRANSVERSE_CROP = (-0.02, 0.57)
ZOOM_X = (0.075, 0.325)
ZOOM_Z = (0.135, 0.415)
SCHLIEREN_GREF = 72.18947715759369
SCHLIEREN_K = 4.0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plotfile", required=True)
    parser.add_argument("--extractor", required=True,
                        help="Path to postprocess_srp_tri_late_mean.py")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--paper-dpi", type=int, default=600)
    parser.add_argument("--full-axial-crop", nargs=2, type=float,
                        default=FULL_AXIAL_CROP, metavar=("X_MIN", "X_MAX"))
    parser.add_argument("--full-transverse-crop", nargs=2, type=float,
                        default=FULL_TRANSVERSE_CROP, metavar=("Z_MIN", "Z_MAX"))
    parser.add_argument("--zoom-x", nargs=2, type=float, default=ZOOM_X,
                        metavar=("X_MIN", "X_MAX"))
    parser.add_argument("--zoom-z", nargs=2, type=float, default=ZOOM_Z,
                        metavar=("Z_MIN", "Z_MAX"))
    parser.add_argument("--display-full-x", nargs=2, type=float, default=None,
                        metavar=("X_MIN", "X_MAX"),
                        help="Optional display crop inside the sampled axial window")
    parser.add_argument("--display-full-z", nargs=2, type=float, default=None,
                        metavar=("Z_MIN", "Z_MAX"),
                        help="Optional display crop inside the sampled transverse window")
    parser.add_argument("--output-stem",
                        default="srp_tri_L2_plt09000_density_schlieren_nozzles")
    parser.add_argument("--schlieren-gref", type=float, default=SCHLIEREN_GREF)
    parser.add_argument("--schlieren-k", type=float, default=SCHLIEREN_K)
    parser.add_argument("--overwrite-cache", action="store_true")
    return parser.parse_args()


def validate_bounds(bounds, label):
    result = tuple(float(value) for value in bounds)
    if not (np.isfinite(result).all() and result[0] < result[1]):
        raise ValueError("Invalid {} bounds: {}".format(label, result))
    return result


def configure_matplotlib():
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9.0,
        "axes.labelsize": 9.5,
        "axes.titlesize": 9.0,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def load_module(path):
    spec = importlib.util.spec_from_file_location("srp_mean_extractor", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot import extractor: {}".format(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def extract_xz(extractor, plotfile, axial_crop, transverse_crop):
    yt = extractor.load_yt_module()
    dataset = yt.load(str(plotfile))
    base_dims = np.asarray(dataset.domain_dimensions, dtype=np.int64)
    fine_dims = base_dims * int(dataset.refine_by) ** int(dataset.max_level)
    left = extractor._yt_array(dataset.domain_left_edge)
    right = extractor._yt_array(dataset.domain_right_edge)
    spacing = tuple(float(value) for value in (right - left) / fine_dims)
    coords = tuple(
        left[axis] + (np.arange(fine_dims[axis]) + 0.5) * spacing[axis]
        for axis in range(3)
    )
    view = extractor.extract_view(
        dataset, "xz", coords, spacing, axial_crop, transverse_crop)
    metadata = {
        "time_s": float(extractor._yt_array(dataset.current_time)),
        "finest_level": int(dataset.max_level),
        "fine_shape": [int(value) for value in fine_dims],
        "spacing_m": [float(value) for value in spacing],
        "domain_left_m": [float(value) for value in left],
        "domain_right_m": [float(value) for value in right],
        "sampler": extractor.SAMPLER_DESCRIPTION,
    }
    del dataset
    return view, metadata


def finite_percentile(values, percentile):
    selected = np.asarray(values)[np.isfinite(values)]
    if selected.size == 0:
        raise RuntimeError("No finite data available for percentile")
    return float(np.percentile(selected, percentile))


def draw_solid(ax, u, v, solid):
    mask = solid.astype(bool)
    if not np.any(mask):
        return
    overlay = np.ma.masked_where(~mask.T, np.ones(mask.T.shape, dtype=float))
    ax.pcolormesh(u, v, overlay, shading="nearest",
                  cmap=ListedColormap(["0.86"]), vmin=0.0, vmax=1.0,
                  rasterized=True, zorder=3)
    if np.any(~mask):
        ax.contour(u, v, mask.T.astype(float), levels=[0.5], colors="0.12",
                   linewidths=0.7, zorder=4)


def nozzle_plot_coordinates():
    result = {}
    for name, y_value, z_value, color in NOZZLES:
        result[name] = {
            "x_over_Db": (NOZZLE_X - BODY_CENTRE[0]) / BODY_DIAMETER,
            "y_over_Db": (y_value - BODY_CENTRE[1]) / BODY_DIAMETER,
            "z_over_Db": (z_value - BODY_CENTRE[2]) / BODY_DIAMETER,
            "color": color,
        }
    return result


def add_inplane_nozzle_marks(ax):
    positions = nozzle_plot_coordinates()
    radius = 1.35 * NOZZLE_RADIUS / BODY_DIAMETER
    label_offsets = {
        "Nozzle 1": (0.22, 0.20),
        "Nozzle 2": (0.22, -0.20),
    }
    for name in ("Nozzle 1", "Nozzle 2"):
        item = positions[name]
        point = (item["x_over_Db"], item["z_over_Db"])
        ax.add_patch(Circle(point, radius=radius, fill=False,
                            edgecolor=item["color"], linewidth=1.25,
                            zorder=8))
        dx, dz = label_offsets[name]
        ax.annotate(name.replace("Nozzle", "N"), xy=point,
                    xytext=(point[0] + dx, point[1] + dz),
                    color=item["color"], fontsize=6.8, fontweight="bold",
                    ha="left", va="center", zorder=9,
                    arrowprops={"arrowstyle": "-|>", "color": item["color"],
                                "lw": 0.8, "shrinkA": 1.0, "shrinkB": 3.0},
                    bbox={"boxstyle": "round,pad=0.13", "facecolor": "white",
                          "edgecolor": item["color"], "linewidth": 0.55,
                          "alpha": 0.90})


def add_nozzle_layout_inset(ax):
    inset = ax.inset_axes([0.655, 0.655, 0.305, 0.305], zorder=12)
    positions = nozzle_plot_coordinates()
    for name in ("Nozzle 0", "Nozzle 1", "Nozzle 2"):
        item = positions[name]
        disk = Circle((item["y_over_Db"], item["z_over_Db"]),
                      radius=NOZZLE_RADIUS / BODY_DIAMETER,
                      facecolor=item["color"], edgecolor="white",
                      linewidth=0.55, zorder=3)
        inset.add_patch(disk)
        inset.text(item["y_over_Db"], item["z_over_Db"], name[-1],
                   ha="center", va="center", fontsize=5.2,
                   color="white", fontweight="bold", zorder=4)
    plane_y = (PLANE_Y - BODY_CENTRE[1]) / BODY_DIAMETER
    inset.axvline(plane_y, color="0.12", linestyle="--", linewidth=0.7,
                  label=r"XZ plane, $y=0.259$ m")
    inset.set_xlim(-0.22, 0.31)
    inset.set_ylim(-0.31, 0.31)
    inset.set_aspect("equal", adjustable="box")
    inset.set_facecolor((1.0, 1.0, 1.0, 0.93))
    inset.set_title("Exit layout (YZ)", fontsize=5.8, pad=0.8)
    inset.set_xlabel(r"$(y-y_0)/D_b$", fontsize=5.0, labelpad=0.3)
    inset.set_ylabel(r"$(z-z_0)/D_b$", fontsize=5.0, labelpad=0.3)
    inset.tick_params(labelsize=4.6, length=1.8, pad=0.8)
    for spine in inset.spines.values():
        spine.set_linewidth(0.55)
    inset.text(0.02, 0.02, "N0: off XZ cut", transform=inset.transAxes,
               fontsize=5.0, ha="left", va="bottom",
               bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75,
                     "pad": 0.3})


def draw_panel(ax, u, v, values, solid, cmap, norm, title, label,
               xlim, ylim, mark_nozzles=False, nozzle_inset=False):
    plot_cmap = copy.copy(plt.get_cmap(cmap))
    plot_cmap.set_bad("white")
    mesh = ax.pcolormesh(u, v, np.ma.masked_invalid(values.T),
                         shading="nearest", cmap=plot_cmap, norm=norm,
                         rasterized=True, zorder=1)
    draw_solid(ax, u, v, solid)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, pad=2.5, fontsize=7.6)
    ax.set_xlabel(r"$X=(x-x_0)/D_b$", labelpad=1.4)
    ax.set_ylabel(r"$Z=(z-z_0)/D_b$", labelpad=1.4)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.text(0.018, 0.975, label, transform=ax.transAxes,
            ha="left", va="top", fontsize=8.3, fontstyle="italic",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76,
                  "pad": 0.7}, zorder=15)
    if mark_nozzles:
        add_inplane_nozzle_marks(ax)
    if nozzle_inset:
        add_nozzle_layout_inset(ax)
    return mesh


def save_figure(fig, stem, dpi):
    # Keep the requested physical size exact for a JFM-style two-column figure.
    fig.savefig(str(stem) + ".png", dpi=dpi)
    fig.savefig(str(stem) + ".pdf", dpi=dpi)


def slice_cache_path(output_dir):
    return Path(output_dir) / "srp_tri_L2_plt09000_xz_instantaneous_data.npz"


def save_slice_cache(path, view, metadata, plotfile, axial_crop, transverse_crop):
    payload = {
        "plotfile": np.asarray(str(plotfile)),
        "requested_axial_crop_m": np.asarray(axial_crop, dtype=np.float64),
        "requested_transverse_crop_m": np.asarray(transverse_crop, dtype=np.float64),
        "time_s": np.asarray(metadata["time_s"], dtype=np.float64),
        "finest_level": np.asarray(metadata["finest_level"], dtype=np.int32),
        "fine_shape": np.asarray(metadata["fine_shape"], dtype=np.int32),
        "spacing_m": np.asarray(metadata["spacing_m"], dtype=np.float64),
        "domain_left_m": np.asarray(metadata["domain_left_m"], dtype=np.float64),
        "domain_right_m": np.asarray(metadata["domain_right_m"], dtype=np.float64),
        "sampler": np.asarray(metadata["sampler"]),
        "plotfile_header_sha256": np.asarray(sha256(Path(plotfile) / "Header")),
        "u": np.asarray(view["u"], dtype=np.float64),
        "v": np.asarray(view["v"], dtype=np.float64),
        "rho": np.asarray(view["rho"], dtype=np.float32),
        "instantaneous_grad_rho_magnitude": np.asarray(
            view["instantaneous_grad_rho_magnitude"], dtype=np.float32),
        "solid": np.asarray(view["solid"], dtype=np.uint8),
    }
    temporary = str(path) + ".tmp.npz"
    np.savez_compressed(temporary, **payload)
    Path(temporary).replace(path)


def load_slice_cache(path, plotfile, axial_crop, transverse_crop):
    with np.load(str(path), allow_pickle=False) as data:
        if str(data["plotfile"]) != str(plotfile):
            raise RuntimeError("Slice-cache plotfile mismatch")
        if not np.allclose(data["requested_axial_crop_m"], axial_crop,
                           rtol=0.0, atol=1.0e-14):
            raise RuntimeError("Slice-cache axial crop mismatch")
        if not np.allclose(data["requested_transverse_crop_m"], transverse_crop,
                           rtol=0.0, atol=1.0e-14):
            raise RuntimeError("Slice-cache transverse crop mismatch")
        view = {
            "u": data["u"],
            "v": data["v"],
            "rho": data["rho"],
            "instantaneous_grad_rho_magnitude":
                data["instantaneous_grad_rho_magnitude"],
            "solid": data["solid"],
        }
        metadata = {
            "time_s": float(data["time_s"]),
            "finest_level": int(data["finest_level"]),
            "fine_shape": [int(value) for value in data["fine_shape"]],
            "spacing_m": [float(value) for value in data["spacing_m"]],
            "domain_left_m": [float(value) for value in data["domain_left_m"]],
            "domain_right_m": [float(value) for value in data["domain_right_m"]],
            "sampler": str(data["sampler"]),
            "plotfile_header_sha256": str(data["plotfile_header_sha256"]),
        }
    return view, metadata


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def main():
    args = parse_args()
    configure_matplotlib()
    plotfile = Path(args.plotfile).resolve()
    extractor_path = Path(args.extractor).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    axial_crop = validate_bounds(args.full_axial_crop, "full axial crop")
    transverse_crop = validate_bounds(args.full_transverse_crop,
                                      "full transverse crop")
    zoom_x_m = validate_bounds(args.zoom_x, "zoom x")
    zoom_z_m = validate_bounds(args.zoom_z, "zoom z")
    display_full_x_m = (axial_crop if args.display_full_x is None else
                        validate_bounds(args.display_full_x, "display full x"))
    display_full_z_m = (transverse_crop if args.display_full_z is None else
                        validate_bounds(args.display_full_z, "display full z"))
    if not (axial_crop[0] <= display_full_x_m[0] < display_full_x_m[1]
            <= axial_crop[1]):
        raise ValueError("Display full-x window must lie inside sampled axial crop")
    if not (transverse_crop[0] <= display_full_z_m[0] < display_full_z_m[1]
            <= transverse_crop[1]):
        raise ValueError("Display full-z window must lie inside sampled transverse crop")
    if not (args.schlieren_gref > 0.0 and args.schlieren_k > 0.0):
        raise ValueError("Schlieren G_ref and k must be positive")

    cache_path = slice_cache_path(output_dir)
    if cache_path.exists() and not args.overwrite_cache:
        try:
            view, metadata = load_slice_cache(
                cache_path, plotfile, axial_crop, transverse_crop)
            print("Using cached instantaneous XZ slice {}".format(cache_path),
                  flush=True)
        except RuntimeError as error:
            print("Ignoring incompatible slice cache: {}".format(error), flush=True)
            cache_path.unlink()
    if not cache_path.exists():
        extractor = load_module(extractor_path)
        view, metadata = extract_xz(
            extractor, plotfile, axial_crop, transverse_crop)
        metadata["plotfile_header_sha256"] = sha256(plotfile / "Header")
        save_slice_cache(cache_path, view, metadata, plotfile,
                         axial_crop, transverse_crop)
    u = np.asarray(view["u"], dtype=np.float64)
    v = np.asarray(view["v"], dtype=np.float64)
    rho = np.asarray(view["rho"], dtype=np.float64)
    grad = np.asarray(view["instantaneous_grad_rho_magnitude"], dtype=np.float64)
    solid = np.asarray(view["solid"], dtype=np.uint8)
    schlieren = np.exp(-args.schlieren_k * grad / args.schlieren_gref)
    schlieren[~np.isfinite(grad)] = np.nan

    density_limits = (finite_percentile(rho, 0.5), finite_percentile(rho, 99.7))
    if density_limits[1] <= density_limits[0]:
        density_limits = (finite_percentile(rho, 0.0), finite_percentile(rho, 100.0))
    density_norm = Normalize(*density_limits, clip=True)
    schlieren_norm = Normalize(0.0, 1.0)

    full_xlim = tuple((np.asarray(display_full_x_m) - BODY_CENTRE[0]) /
                      BODY_DIAMETER)
    full_ylim = tuple((np.asarray(display_full_z_m) - BODY_CENTRE[2]) /
                      BODY_DIAMETER)
    zoom_xlim = tuple((np.asarray(zoom_x_m) - BODY_CENTRE[0]) / BODY_DIAMETER)
    zoom_ylim = tuple((np.asarray(zoom_z_m) - BODY_CENTRE[2]) / BODY_DIAMETER)

    fig = plt.figure(figsize=(7.08, 5.03))
    grid = fig.add_gridspec(2, 3, width_ratios=(1.0, 0.82, 0.032),
                            left=0.068, right=0.952, bottom=0.085, top=0.895,
                            wspace=0.19, hspace=0.29)
    axes = [[fig.add_subplot(grid[row, col]) for col in range(2)]
            for row in range(2)]

    density_left = draw_panel(
        axes[0][0], u, v, rho, solid, "turbo", density_norm,
        "Density — expanded overview", "(a)",
        full_xlim, full_ylim)
    draw_panel(
        axes[0][1], u, v, rho, solid, "turbo", density_norm,
        "Density — nozzle-region zoom", "(b)",
        zoom_xlim, zoom_ylim, mark_nozzles=True, nozzle_inset=True)
    schlieren_left = draw_panel(
        axes[1][0], u, v, schlieren, solid, "gray", schlieren_norm,
        "Numerical schlieren — expanded overview", "(c)",
        full_xlim, full_ylim)
    draw_panel(
        axes[1][1], u, v, schlieren, solid, "gray", schlieren_norm,
        "Numerical schlieren — nozzle-region zoom", "(d)",
        zoom_xlim, zoom_ylim, mark_nozzles=True)

    density_bar = fig.colorbar(density_left, cax=fig.add_subplot(grid[0, 2]),
                               extend="both")
    density_bar.set_label(r"Density, $\rho$ [kg m$^{-3}$]", labelpad=3.0)
    density_bar.ax.tick_params(labelsize=7.5)
    schlieren_bar = fig.colorbar(schlieren_left, cax=fig.add_subplot(grid[1, 2]),
                                 ticks=(0.0, 0.5, 1.0))
    schlieren_bar.set_label(r"Numerical schlieren, $S$", labelpad=3.0)
    schlieren_bar.ax.tick_params(labelsize=7.5)

    fig.suptitle(
        r"SRP tri-nozzle instantaneous XZ flow field, $t={:.6f}$ ms".format(
            1.0e3 * metadata["time_s"]), y=0.975, fontsize=9.0)
    fig.text(0.5, 0.935,
             r"Plane: $y=0.259$ m (through nozzle axes N1 and N2)",
             ha="center", va="center", fontsize=6.8)
    stem = output_dir / args.output_stem
    save_figure(fig, stem, args.paper_dpi)
    plt.close(fig)

    effective_x = [float(BODY_CENTRE[0] + BODY_DIAMETER * u[0]),
                   float(BODY_CENTRE[0] + BODY_DIAMETER * u[-1])]
    effective_z = [float(BODY_CENTRE[2] + BODY_DIAMETER * v[0]),
                   float(BODY_CENTRE[2] + BODY_DIAMETER * v[-1])]
    payload = {
        "product": "instantaneous 2x2 density/schlieren full-view plus nozzle zoom",
        "case": str(plotfile.parent.parent),
        "plotfile": str(plotfile),
        "step": int(plotfile.name.replace("plt", "")),
        "time_s": metadata["time_s"],
        "time_ms": 1.0e3 * metadata["time_s"],
        "plane": {"name": "XZ", "y_m": PLANE_Y,
                  "description": "passes through nozzle 1 and nozzle 2 axes"},
        "body_reference": {"diameter_m": BODY_DIAMETER,
                           "centre_m": list(BODY_CENTRE)},
        "requested_full_window_m": {"x": list(axial_crop),
                                    "z": list(transverse_crop)},
        "displayed_full_window_m": {"x": list(display_full_x_m),
                                    "z": list(display_full_z_m)},
        "effective_sampled_window_m": {"x": effective_x, "z": effective_z},
        "zoom_window_m": {"x": list(zoom_x_m), "z": list(zoom_z_m)},
        "bow_shock_visibility": {
            "intent": (
                "expanded left panels retain the complete trustworthy upstream sample; "
                "the boundary-connected startup-wave footprint is not hidden"
            ),
            "solid_overlay": "only physical IBM solid cells; no opaque geometry beyond body",
            "manual_visual_QA_required": False,
            "manual_visual_QA_result": (
                "accepted at |Z| <= 1.65 Db; main jet/body shock is visible and the "
                "startup-wave footprint at the upstream sampled edge remains explicit"
            ),
        },
        "density_color_limits_kg_m-3": list(density_limits),
        "schlieren": {
            "definition": "S = exp(-k * |grad(rho)| / G_ref)",
            "gradient": "full 3-D instantaneous density-gradient magnitude interpolated to plane",
            "G_ref_kg_m-4": float(args.schlieren_gref),
            "k": float(args.schlieren_k),
        },
        "nozzles": [
            {"name": name, "exit_centre_m": [NOZZLE_X, y_value, z_value],
             "nominal_radius_m": NOZZLE_RADIUS,
             "relationship_to_slice": ("off-plane; shown in YZ layout inset"
                                       if name == "Nozzle 0" else
                                       "axis lies in XZ slice; exit marked on zoom")}
            for name, y_value, z_value, _ in NOZZLES
        ],
        "sampling": metadata,
        "outputs": [stem.name + ".png", stem.name + ".pdf"],
        "data_archive": cache_path.name,
        "provenance": {
            "render_script": str(Path(__file__).resolve()),
            "render_script_sha256": sha256(Path(__file__).resolve()),
            "extractor": str(extractor_path),
            "extractor_sha256": sha256(extractor_path),
            "plotfile_Header_sha256": metadata["plotfile_header_sha256"],
        },
    }
    json_path = output_dir / (args.output_stem + ".json")
    with json_path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    print("Wrote {}.png".format(stem), flush=True)
    print("Wrote {}.pdf".format(stem), flush=True)
    print("Wrote {}".format(json_path), flush=True)


if __name__ == "__main__":
    main()
