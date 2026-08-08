#!/usr/bin/env python3
"""Create thesis figures for the three selected inviscid discretisations.

The script reads the completed Y9000X result matrix.  It does not rerun the
solver or alter any result file.  AMR fields are evaluated at the true leaf
cell centres and are drawn over the physical extent of each leaf cell.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


os.environ.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")

SCRIPT_DIR = Path(__file__).resolve().parent
NUMERICS_DIR = SCRIPT_DIR.parents[1]
if str(NUMERICS_DIR) not in sys.path:
    sys.path.insert(0, str(NUMERICS_DIR))

from isentropic_vortex.analyze import (  # noqa: E402
    GAMMA,
    R_SPECIFIC,
    exact_state,
    latest_plotfile as latest_vortex_plotfile,
)
from shuosher.analyze import read_profile as read_shu_profile  # noqa: E402

try:
    import yt
except ImportError as exc:  # pragma: no cover - dependency diagnostic
    raise SystemExit("this script requires yt") from exc


SCHEMES = (
    ("llf-wenoz5", "LLF--WENO-Z5", "#0072B2", "o", "-"),
    ("llf-teno5", "LLF--TENO5", "#009E73", "s", "--"),
    (
        "afd-hllc-wenoz5",
        "AFD--HLLC--WENO-Z5",
        "#D55E00",
        "^",
        "-.",
    ),
)
SHU_RESOLUTIONS = (200, 400, 800)
VORTEX_RESOLUTIONS = (40, 80, 160)


def configure_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ("Latin Modern Roman",),
            "text.usetex": True,
            "text.latex.preamble": (
                r"\usepackage{lmodern}\usepackage{amsmath}"
            ),
            "font.size": 8.2,
            "axes.labelsize": 8.2,
            "axes.titlesize": 8.4,
            "axes.linewidth": 0.7,
            "legend.fontsize": 7.2,
            "xtick.labelsize": 7.3,
            "ytick.labelsize": 7.3,
            "lines.linewidth": 1.10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 400,
        }
    )
    return plt


def save_figure(figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    figure.savefig(png, dpi=400, bbox_inches="tight", pad_inches=0.02)
    figure.savefig(pdf, dpi=400, bbox_inches="tight", pad_inches=0.02)
    print(f"wrote {png}")
    print(f"wrote {pdf}")


def plot_shu_osher(
    root: Path,
    output_dir: Path,
    plt,
) -> None:
    profiles: dict[tuple[str, int], tuple[np.ndarray, np.ndarray]] = {}
    common_time: float | None = None
    extrema = []

    for scheme, _, _, _, _ in SCHEMES:
        for resolution in SHU_RESOLUTIONS:
            case = root / scheme / f"N{resolution}"
            plotfile, dataset, x, density = read_shu_profile(case)
            if x.size != resolution:
                raise ValueError(
                    f"{plotfile} contains {x.size} cells, expected {resolution}"
                )
            if int(dataset.index.max_level) != 0:
                raise ValueError(f"{plotfile} is not a uniform-grid result")
            time = float(dataset.current_time)
            if common_time is None:
                common_time = time
            elif not math.isclose(
                time,
                common_time,
                rel_tol=0.0,
                abs_tol=1.0e-12 * max(1.0, abs(common_time)),
            ):
                raise ValueError("Shu Osher results have different final times")
            profiles[(scheme, resolution)] = (x, density)
            visible = (x >= -1.5) & (x <= 2.5)
            extrema.extend(density[visible])

    y_min = float(np.min(extrema))
    y_max = float(np.max(extrema))
    y_margin = 0.035 * (y_max - y_min)

    figure, axes = plt.subplots(
        3,
        1,
        figsize=(5.9055, 5.15),
        sharex=True,
        sharey=True,
    )
    for panel, (axis, resolution) in enumerate(
        zip(axes, SHU_RESOLUTIONS)
    ):
        for scheme, label, colour, _, line_style in SCHEMES:
            x, density = profiles[(scheme, resolution)]
            axis.plot(
                x,
                density,
                color=colour,
                linestyle=line_style,
                label=label,
            )
        axis.set_xlim(-1.5, 2.5)
        axis.set_ylim(y_min - y_margin, y_max + y_margin)
        axis.set_ylabel(r"$\rho$")
        axis.tick_params(direction="in", which="both", top=True, right=True)
        axis.grid(
            axis="y",
            color="0.88",
            linestyle=":",
            linewidth=0.6,
        )
        axis.text(
            0.018,
            0.945,
            rf"$\mathbf{{({chr(ord('a') + panel)})}}$",
            transform=axis.transAxes,
            ha="left",
            va="top",
        )
        axis.text(
            0.985,
            0.945,
            rf"$N={resolution}$",
            transform=axis.transAxes,
            ha="right",
            va="top",
        )

    axes[-1].set_xlabel(r"$x$")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.995),
        ncol=3,
        frameon=False,
        columnspacing=1.2,
        handlelength=2.2,
    )
    figure.subplots_adjust(
        left=0.10,
        right=0.985,
        bottom=0.085,
        top=0.94,
        hspace=0.10,
    )
    save_figure(figure, output_dir, "shu_osher_three_scheme")
    plt.close(figure)
    print(f"Shu Osher common final time: {common_time:.16e}")


def read_vortex_rows(
    root: Path,
) -> dict[tuple[str, str, int], dict[str, float]]:
    indexed: dict[tuple[str, str, int], dict[str, float]] = {}
    for scheme, _, _, _, _ in SCHEMES:
        for suite in ("uniform", "amr"):
            source = root / scheme / suite / "errors.csv"
            if not source.is_file():
                raise FileNotFoundError(source)
            with source.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            for source_row in rows:
                resolution = int(source_row["base_n"])
                if resolution not in VORTEX_RESOLUTIONS:
                    continue
                expected_level = 0 if suite == "uniform" else 1
                if int(source_row["max_level"]) != expected_level:
                    raise ValueError(
                        f"{source} contains the wrong AMR level for {suite}"
                    )
                key = (scheme, suite, resolution)
                if key in indexed:
                    raise ValueError(f"duplicate result {key} in {source}")
                indexed[key] = {
                    "density_L2": float(source_row["density_L2"]),
                    "ke_deficit": abs(
                        1.0 - float(source_row["perturbation_ke_ratio"])
                    ),
                    "time": float(source_row["time"]),
                }

    expected = {
        (scheme, suite, resolution)
        for scheme, _, _, _, _ in SCHEMES
        for suite in ("uniform", "amr")
        for resolution in VORTEX_RESOLUTIONS
    }
    if set(indexed) != expected:
        missing = sorted(expected - set(indexed))
        raise ValueError(f"the vortex matrix is incomplete: {missing}")
    return indexed


def plot_vortex_accuracy(
    root: Path,
    output_dir: Path,
    plt,
) -> None:
    from matplotlib.lines import Line2D

    indexed = read_vortex_rows(root)
    suite_styles = {
        "uniform": ("-", True, "Uniform grid"),
        "amr": ("--", False, "AMR"),
    }
    figure, axes = plt.subplots(1, 2, figsize=(5.9055, 2.78))

    for scheme, _, colour, marker, _ in SCHEMES:
        for suite, (line_style, filled, _) in suite_styles.items():
            selected = [
                indexed[(scheme, suite, resolution)]
                for resolution in VORTEX_RESOLUTIONS
            ]
            marker_face = colour if filled else "white"
            common = {
                "color": colour,
                "linestyle": line_style,
                "marker": marker,
                "markersize": 4.5,
                "markerfacecolor": marker_face,
                "markeredgecolor": colour,
                "markeredgewidth": 0.8,
            }
            axes[0].plot(
                VORTEX_RESOLUTIONS,
                [row["density_L2"] for row in selected],
                **common,
            )
            axes[1].plot(
                VORTEX_RESOLUTIONS,
                [row["ke_deficit"] for row in selected],
                **common,
            )

    for panel, axis in enumerate(axes):
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xticks(
            VORTEX_RESOLUTIONS,
            [str(value) for value in VORTEX_RESOLUTIONS],
        )
        axis.set_xlabel(r"Base grid resolution, $N$")
        axis.tick_params(direction="in", which="both", top=True, right=True)
        axis.grid(
            axis="y",
            which="major",
            color="0.85",
            linestyle=":",
            linewidth=0.7,
        )
        axis.text(
            -0.10,
            1.035,
            rf"$\mathbf{{({chr(ord('a') + panel)})}}$",
            transform=axis.transAxes,
            ha="left",
            va="bottom",
        )

    axes[0].set_ylabel(
        r"Density error, $\|e_\rho\|_{2,h}$"
    )
    axes[1].set_ylabel(
        r"Energy loss, "
        r"$|1-K^\prime/K^{\prime\star}|$"
    )

    scheme_handles = [
        Line2D(
            [0],
            [0],
            color=colour,
            marker=marker,
            linestyle="-",
            markersize=4.5,
            label=label,
        )
        for _, label, colour, marker, _ in SCHEMES
    ]
    suite_handles = [
        Line2D(
            [0],
            [0],
            color="0.2",
            linestyle=line_style,
            marker="o",
            markerfacecolor="0.2" if filled else "white",
            markersize=4.5,
            label=label,
        )
        for line_style, filled, label in suite_styles.values()
    ]
    figure.legend(
        handles=scheme_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.015),
        ncol=3,
        frameon=False,
        columnspacing=1.4,
        handlelength=2.0,
    )
    axes[1].legend(
        handles=suite_handles,
        loc="lower left",
        frameon=False,
        handlelength=2.2,
    )
    figure.subplots_adjust(
        left=0.105,
        right=0.985,
        bottom=0.17,
        top=0.82,
        wspace=0.35,
    )
    save_figure(figure, output_dir, "vortex_accuracy_three_scheme")
    plt.close(figure)


def field_key(dataset, name: str):
    wanted = name.lower()
    for key in dataset.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"field {name!r} is unavailable in {dataset}")


def uniform_field(dataset, name: str) -> np.ndarray:
    dimensions = [int(value) for value in dataset.domain_dimensions]
    grid = dataset.covering_grid(
        level=0,
        left_edge=dataset.domain_left_edge,
        dims=dimensions,
    )
    values = np.asarray(grid[field_key(dataset, name)].d, dtype=float)
    if values.shape != tuple(dimensions):
        raise ValueError(
            f"unexpected {name} shape {values.shape}, expected {dimensions}"
        )
    return values[:, :, 0].T


def nondimensional_vorticity(
    dataset,
    resolution: int,
) -> np.ndarray:
    if int(dataset.index.max_level) != 0:
        raise ValueError("vorticity panels require uniform results")
    if int(dataset.domain_dimensions[0]) != resolution:
        raise ValueError("the vorticity result has the wrong resolution")
    velocity_x = uniform_field(dataset, "x_velocity")
    velocity_y = uniform_field(dataset, "y_velocity")
    lower = np.asarray(dataset.domain_left_edge.d, dtype=float)
    upper = np.asarray(dataset.domain_right_edge.d, dtype=float)
    dx = (upper[0] - lower[0]) / resolution
    dy = (upper[1] - lower[1]) / resolution
    dv_dx = (
        np.roll(velocity_y, -1, axis=1)
        - np.roll(velocity_y, 1, axis=1)
    ) / (2.0 * dx)
    du_dy = (
        np.roll(velocity_x, -1, axis=0)
        - np.roll(velocity_x, 1, axis=0)
    ) / (2.0 * dy)
    u_inf = 0.5 * math.sqrt(GAMMA * R_SPECIFIC * 300.0)
    return (dv_dx - du_dy) * 0.005 / u_inf


def exact_arguments() -> SimpleNamespace:
    return SimpleNamespace(
        mach=0.5,
        beta=0.2,
        p_inf=1.0e5,
        T_inf=300.0,
        radius=0.005,
        xc0=0.025,
        yc0=0.05,
        x_lo=0.0,
        x_hi=0.1,
        y_lo=0.0,
        y_hi=0.1,
    )


def leaf_density_error(dataset) -> dict[str, np.ndarray]:
    """Return relative density error at each true AMR leaf-cell centre."""

    if int(dataset.index.max_level) != 1:
        raise ValueError("the AMR field panel requires one refined level")
    data = dataset.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    y = np.asarray(data[("index", "y")].d, dtype=float)
    dx = np.asarray(data[("index", "dx")].d, dtype=float)
    dy = np.asarray(data[("index", "dy")].d, dtype=float)
    density = np.asarray(
        data[field_key(dataset, "Density")].d,
        dtype=float,
    )
    exact_density, *_ = exact_state(
        x,
        y,
        float(dataset.current_time),
        exact_arguments(),
    )
    rho_inf = 1.0e5 / (R_SPECIFIC * 300.0)
    error = 100.0 * (density - exact_density) / rho_inf

    lower = np.asarray(dataset.domain_left_edge.d, dtype=float)
    upper = np.asarray(dataset.domain_right_edge.d, dtype=float)
    domain_area = (upper[0] - lower[0]) * (upper[1] - lower[1])
    leaf_area = float(np.sum(dx * dy))
    if not math.isclose(
        leaf_area,
        domain_area,
        rel_tol=0.0,
        abs_tol=2.0e-12 * domain_area,
    ):
        raise ValueError(
            "AMR leaf-cell areas do not cover the physical domain exactly"
        )
    return {
        "x": x,
        "y": y,
        "dx": dx,
        "dy": dy,
        "error": error,
        "lower": lower,
        "upper": upper,
    }


def plot_field_matrix(
    root: Path,
    output_dir: Path,
    plt,
) -> None:
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.patches import Rectangle

    vorticity_fields = []
    error_fields = []
    for scheme, label, _, _, _ in SCHEMES:
        uniform_case = root / scheme / "uniform" / "N40"
        amr_case = root / scheme / "amr" / "N160"
        uniform_plotfile = latest_vortex_plotfile(uniform_case)
        amr_plotfile = latest_vortex_plotfile(amr_case)
        uniform_dataset = yt.load(str(uniform_plotfile))
        amr_dataset = yt.load(str(amr_plotfile))
        if int(amr_dataset.domain_dimensions[0]) != 160:
            raise ValueError(f"{amr_plotfile} is not the N=160 AMR result")
        vorticity_fields.append(
            (label, nondimensional_vorticity(uniform_dataset, 40))
        )
        error_fields.append((label, leaf_density_error(amr_dataset)))

    vorticity_min = min(float(np.min(field)) for _, field in vorticity_fields)
    vorticity_max = max(float(np.max(field)) for _, field in vorticity_fields)
    vorticity_limit_low = min(-0.07, vorticity_min)
    vorticity_limit_high = max(0.40, vorticity_max)
    error_limit = max(
        float(np.max(np.abs(field["error"])))
        for _, field in error_fields
    )
    if error_limit <= 0.0:
        raise ValueError("the AMR error range is degenerate")

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(5.9055, 4.18),
        gridspec_kw={"hspace": 0.28, "wspace": 0.14},
    )
    vorticity_norm = TwoSlopeNorm(
        vmin=vorticity_limit_low,
        vcenter=0.0,
        vmax=vorticity_limit_high,
    )
    error_norm = TwoSlopeNorm(
        vmin=-error_limit,
        vcenter=0.0,
        vmax=error_limit,
    )
    vorticity_image = None
    error_collection = None

    for column, ((label, vorticity), (_, leaf_field)) in enumerate(
        zip(vorticity_fields, error_fields)
    ):
        top_axis = axes[0, column]
        bottom_axis = axes[1, column]
        vorticity_image = top_axis.imshow(
            vorticity,
            origin="lower",
            extent=(0.0, 1.0, 0.0, 1.0),
            cmap="RdBu_r",
            norm=vorticity_norm,
            interpolation="nearest",
            rasterized=True,
        )

        lower = leaf_field["lower"]
        upper = leaf_field["upper"]
        length_x = upper[0] - lower[0]
        length_y = upper[1] - lower[1]
        rectangles = [
            Rectangle(
                (
                    (x - 0.5 * dx - lower[0]) / length_x,
                    (y - 0.5 * dy - lower[1]) / length_y,
                ),
                dx / length_x,
                dy / length_y,
            )
            for x, y, dx, dy in zip(
                leaf_field["x"],
                leaf_field["y"],
                leaf_field["dx"],
                leaf_field["dy"],
            )
        ]
        error_collection = PatchCollection(
            rectangles,
            cmap="RdBu_r",
            norm=error_norm,
            linewidth=0.0,
            edgecolor="none",
            rasterized=True,
        )
        error_collection.set_array(leaf_field["error"])
        bottom_axis.add_collection(error_collection)
        bottom_axis.axvline(
            0.5,
            color="black",
            linestyle="--",
            linewidth=0.75,
        )
        top_axis.set_title(label, pad=5.0)
        top_axis.set_xlim(0.30, 0.70)
        top_axis.set_ylim(0.30, 0.70)
        top_axis.set_xticks((0.35, 0.50, 0.65))
        top_axis.set_yticks((0.35, 0.50, 0.65))
        bottom_axis.set_xlim(0.0, 1.0)
        bottom_axis.set_ylim(0.0, 1.0)
        bottom_axis.set_xticks((0.25, 0.50, 0.75))
        bottom_axis.set_yticks((0.25, 0.50, 0.75))

        for row, axis in enumerate((top_axis, bottom_axis)):
            axis.set_aspect("equal")
            axis.tick_params(
                direction="in",
                which="both",
                top=True,
                right=True,
                labelleft=(column == 0),
            )
            axis.text(
                0.025,
                0.975,
                rf"$\mathbf{{({chr(ord('a') + 3 * row + column)})}}$",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=8.2,
                color="black",
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.72,
                    "pad": 1.0,
                },
            )
        bottom_axis.set_xlabel(r"$x/L$")

    axes[0, 0].set_ylabel(r"$y/L$")
    axes[1, 0].set_ylabel(r"$y/L$")
    figure.text(
        0.020,
        0.715,
        r"Uniform $N=40$",
        rotation=90,
        ha="center",
        va="center",
    )
    figure.text(
        0.020,
        0.270,
        r"AMR $N_0=160$",
        rotation=90,
        ha="center",
        va="center",
    )
    figure.subplots_adjust(
        left=0.095,
        right=0.885,
        bottom=0.105,
        top=0.925,
    )
    top_cbar_axis = figure.add_axes((0.910, 0.565, 0.017, 0.305))
    bottom_cbar_axis = figure.add_axes((0.910, 0.145, 0.017, 0.305))
    assert vorticity_image is not None
    assert error_collection is not None
    top_cbar = figure.colorbar(vorticity_image, cax=top_cbar_axis)
    bottom_cbar = figure.colorbar(error_collection, cax=bottom_cbar_axis)
    top_cbar.set_ticks((-0.05, 0.0, 0.1, 0.2, 0.3, 0.4))
    top_cbar.set_label(r"$\omega_z R_{\mathrm{v}}/U_\infty$")
    bottom_cbar.set_label(
        r"$100(\rho-\rho^\star)/\rho_\infty$"
    )
    top_cbar.ax.tick_params(labelsize=7.0, direction="in")
    bottom_cbar.ax.tick_params(labelsize=7.0, direction="in")

    save_figure(figure, output_dir, "vortex_fields_three_scheme")
    plt.close(figure)
    print(
        "common vorticity range: "
        f"[{vorticity_limit_low:.6e}, {vorticity_limit_high:.6e}]"
    )
    print(f"common AMR density-error limit: {error_limit:.6e} percent")
    for (label, vorticity), (_, leaf_field) in zip(
        vorticity_fields,
        error_fields,
    ):
        print(
            f"{label}: vorticity "
            f"[{np.min(vorticity):.6e}, {np.max(vorticity):.6e}], "
            f"AMR leaf cells={leaf_field['error'].size}, "
            f"density error "
            f"[{np.min(leaf_field['error']):.6e}, "
            f"{np.max(leaf_field['error']):.6e}] percent"
        )


def main() -> None:
    default_shu_root = (
        NUMERICS_DIR
        / "shuosher"
        / "results"
        / "y9000x_four_scheme_20260727_r2"
    )
    default_vortex_root = (
        NUMERICS_DIR
        / "isentropic_vortex"
        / "results"
        / "y9000x_four_scheme_20260727_r2"
    )
    default_output = SCRIPT_DIR / "output"

    parser = argparse.ArgumentParser()
    parser.add_argument("--shu-root", type=Path, default=default_shu_root)
    parser.add_argument(
        "--vortex-root",
        type=Path,
        default=default_vortex_root,
    )
    parser.add_argument("--output-dir", type=Path, default=default_output)
    args = parser.parse_args()

    yt.set_log_level(40)
    plt = configure_matplotlib()
    plot_shu_osher(args.shu_root, args.output_dir, plt)
    plot_vortex_accuracy(args.vortex_root, args.output_dir, plt)
    plot_field_matrix(args.vortex_root, args.output_dir, plt)


if __name__ == "__main__":
    main()
