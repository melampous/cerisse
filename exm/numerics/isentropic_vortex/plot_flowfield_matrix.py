#!/usr/bin/env python3
"""Plot uniform-grid vortex fields and AMR-interface errors for four schemes."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import yt

from analyze import (
    GAMMA,
    R_SPECIFIC,
    exact_state,
    latest_plotfile,
)


SCHEMES = (
    ("llf-wenoz5", "LLF WENO-Z5"),
    ("llf-teno5", "LLF TENO5"),
    ("afd-hllc-wenoz5", "AFD HLLC WENO-Z5"),
    ("afd-hllc-teno5", "AFD HLLC TENO5"),
)


def field_key(dataset, name: str):
    wanted = name.lower()
    for key in dataset.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"field {name!r} is unavailable in {dataset}")


def covering_field(dataset, field: str, level: int) -> np.ndarray:
    refinement = 2**level
    dimensions = [
        int(dataset.domain_dimensions[0]) * refinement,
        int(dataset.domain_dimensions[1]) * refinement,
        1,
    ]
    grid = dataset.covering_grid(
        level=level,
        left_edge=dataset.domain_left_edge,
        dims=dimensions,
    )
    values = np.asarray(grid[field_key(dataset, field)].d, dtype=float)
    if values.shape != tuple(dimensions):
        raise ValueError(
            f"unexpected {field} shape {values.shape}; expected {tuple(dimensions)}"
        )
    return values[:, :, 0].T


def composite_error_raster(
    dataset,
    exact_args: SimpleNamespace,
    rho_inf: float,
) -> np.ndarray:
    """Rasterise leaf-cell errors without inventing fine-grid sample points."""

    level = int(dataset.index.max_level)
    refinement = 2**level
    nx = int(dataset.domain_dimensions[0]) * refinement
    ny = int(dataset.domain_dimensions[1]) * refinement
    lower = np.asarray(dataset.domain_left_edge.d, dtype=float)
    upper = np.asarray(dataset.domain_right_edge.d, dtype=float)
    fine_dx = (upper[0] - lower[0]) / nx
    fine_dy = (upper[1] - lower[1]) / ny

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
        exact_args,
    )
    cell_error = 100.0 * (density - exact_density) / rho_inf

    raster = np.full((ny, nx), np.nan, dtype=float)
    coverage = np.zeros((ny, nx), dtype=np.uint8)
    for x_cell, y_cell, dx_cell, dy_cell, error in zip(
        x,
        y,
        dx,
        dy,
        cell_error,
    ):
        width = int(round(dx_cell / fine_dx))
        height = int(round(dy_cell / fine_dy))
        x_start = int(
            round((x_cell - 0.5 * dx_cell - lower[0]) / fine_dx)
        )
        y_start = int(
            round((y_cell - 0.5 * dy_cell - lower[1]) / fine_dy)
        )
        x_stop = x_start + width
        y_stop = y_start + height
        if (
            width < 1
            or height < 1
            or x_start < 0
            or y_start < 0
            or x_stop > nx
            or y_stop > ny
        ):
            raise ValueError(
                "a leaf cell cannot be mapped to the finest display raster"
            )
        if np.any(coverage[y_start:y_stop, x_start:x_stop] != 0):
            raise ValueError("overlapping leaf cells found in AMR display raster")
        raster[y_start:y_stop, x_start:x_stop] = error
        coverage[y_start:y_stop, x_start:x_stop] = 1

    if not np.all(coverage == 1) or not np.all(np.isfinite(raster)):
        raise ValueError("AMR leaf cells do not cover the display raster exactly")
    return raster


def load_case(path: Path, expected_level: int):
    plotfile = latest_plotfile(path)
    dataset = yt.load(str(plotfile))
    if int(dataset.index.max_level) != expected_level:
        raise ValueError(
            f"{plotfile} has max_level={dataset.index.max_level}; "
            f"expected {expected_level}"
        )
    return plotfile, dataset


def exact_arguments(xc0: float) -> SimpleNamespace:
    return SimpleNamespace(
        mach=0.5,
        beta=0.2,
        p_inf=1.0e5,
        T_inf=300.0,
        radius=0.005,
        xc0=xc0,
        yc0=0.05,
        x_lo=0.0,
        x_hi=0.1,
        y_lo=0.0,
        y_hi=0.1,
    )


def plot_vorticity_comparison(
    root: Path,
    resolution: int,
    png_output: Path,
    pdf_output: Path,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    u_inf = 0.5 * np.sqrt(GAMMA * R_SPECIFIC * 300.0)
    vortex_radius = 0.005
    fields = []
    extrema = []
    for scheme, label in SCHEMES:
        case = root / scheme / "uniform" / f"N{resolution}"
        plotfile, dataset = load_case(case, 0)
        if int(dataset.domain_dimensions[0]) != resolution:
            raise ValueError(f"{plotfile} is not the requested N={resolution} case")
        velocity_x = covering_field(dataset, "x_velocity", 0)
        velocity_y = covering_field(dataset, "y_velocity", 0)
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
        dimensionless_vorticity = (dv_dx - du_dy) * vortex_radius / u_inf
        fields.append((label, dimensionless_vorticity))
        extrema.append(
            (
                label,
                float(np.min(dimensionless_vorticity)),
                float(np.max(dimensionless_vorticity)),
            )
        )

    norm = TwoSlopeNorm(vmin=-0.07, vcenter=0.0, vmax=0.40)
    figure, axes = plt.subplots(
        1,
        len(SCHEMES),
        figsize=(7.4, 2.15),
        sharex=True,
        sharey=True,
    )
    image = None
    for column, (axis, (label, vorticity)) in enumerate(zip(axes, fields)):
        image = axis.imshow(
            vorticity,
            origin="lower",
            extent=(0.0, 1.0, 0.0, 1.0),
            cmap="RdBu_r",
            norm=norm,
            interpolation="nearest",
            rasterized=True,
        )
        axis.set_title(label, pad=6.0)
        axis.set_xlim(0.30, 0.70)
        axis.set_ylim(0.30, 0.70)
        axis.set_xticks((0.35, 0.50, 0.65))
        axis.set_yticks((0.35, 0.50, 0.65))
        axis.set_aspect("equal")
        axis.set_xlabel(r"$x/L$")
        axis.tick_params(direction="in", top=True, right=True)
        axis.text(
            0.035,
            0.965,
            f"({chr(ord('a') + column)})",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.0,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.72,
                "pad": 1.0,
            },
        )
    axes[0].set_ylabel(r"$y/L$")
    figure.subplots_adjust(
        left=0.07,
        right=0.91,
        bottom=0.20,
        top=0.86,
        wspace=0.12,
    )
    colorbar_axis = figure.add_axes((0.93, 0.22, 0.014, 0.60))
    assert image is not None
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_ticks((-0.05, 0.0, 0.1, 0.2, 0.3, 0.4))
    colorbar.set_label(r"Vorticity, $\omega_z R_v/u_\infty$", fontsize=8.0)
    colorbar.ax.tick_params(labelsize=7.0, direction="in")
    png_output.parent.mkdir(parents=True, exist_ok=True)
    pdf_output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png_output, dpi=300, bbox_inches="tight")
    figure.savefig(pdf_output, dpi=300, bbox_inches="tight")
    plt.close(figure)

    for label, minimum, maximum in extrema:
        print(
            f"{label} N={resolution} dimensionless vorticity range: "
            f"[{minimum:.6e}, {maximum:.6e}]"
        )
    print(f"wrote {png_output}")
    print(f"wrote {pdf_output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        type=Path,
        help="four-scheme result root produced by run_method_matrix.sh",
    )
    parser.add_argument("--resolution", type=int, default=160)
    parser.add_argument(
        "--png",
        type=Path,
        help="PNG output; default: ROOT/vortex_flowfield_NRES.png",
    )
    parser.add_argument(
        "--pdf",
        type=Path,
        help="PDF output; default: ROOT/vortex_flowfield_NRES.pdf",
    )
    parser.add_argument("--vorticity-resolution", type=int, default=40)
    parser.add_argument(
        "--vorticity-png",
        type=Path,
        help="vorticity PNG; default: ROOT/vortex_vorticity_NRES.png",
    )
    parser.add_argument(
        "--vorticity-pdf",
        type=Path,
        help="vorticity PDF; default: ROOT/vortex_vorticity_NRES.pdf",
    )
    args = parser.parse_args()

    if args.resolution <= 0 or args.vorticity_resolution <= 0:
        raise ValueError("resolutions must be positive")

    rho_inf = 1.0e5 / (R_SPECIFIC * 300.0)
    u_inf = 0.5 * np.sqrt(GAMMA * R_SPECIFIC * 300.0)
    uniform_fields = []
    amr_errors = []
    metadata = []

    for scheme, label in SCHEMES:
        uniform_case = args.root / scheme / "uniform" / f"N{args.resolution}"
        amr_case = args.root / scheme / "amr" / f"N{args.resolution}"
        uniform_plot, uniform_dataset = load_case(uniform_case, 0)
        amr_plot, amr_dataset = load_case(amr_case, 1)

        if int(uniform_dataset.domain_dimensions[0]) != args.resolution:
            raise ValueError(
                f"{uniform_plot} is not the requested N={args.resolution} case"
            )
        if int(amr_dataset.domain_dimensions[0]) != args.resolution:
            raise ValueError(
                f"{amr_plot} is not the requested N={args.resolution} case"
            )

        density = covering_field(uniform_dataset, "Density", 0)
        velocity_x = covering_field(uniform_dataset, "x_velocity", 0)
        velocity_y = covering_field(uniform_dataset, "y_velocity", 0)
        uniform_fields.append(
            {
                "density_percent": 100.0 * (density - rho_inf) / rho_inf,
                "u_perturbation": velocity_x - u_inf,
                "v_perturbation": velocity_y,
            }
        )

        amr_errors.append(
            composite_error_raster(
                amr_dataset,
                exact_arguments(0.025),
                rho_inf,
            )
        )
        metadata.append(
            {
                "scheme": scheme,
                "label": label,
                "uniform_plotfile": str(uniform_plot),
                "amr_plotfile": str(amr_plot),
            }
        )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top_min = min(float(np.min(field["density_percent"])) for field in uniform_fields)
    top_max = max(float(np.max(field["density_percent"])) for field in uniform_fields)
    error_limit = max(float(np.max(np.abs(error))) for error in amr_errors)
    if not top_min < top_max:
        raise ValueError("uniform density range is degenerate")
    if not error_limit > 0.0:
        raise ValueError("AMR density-error range is degenerate")

    rc_parameters = {
        "font.family": "serif",
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 8.5,
        "axes.linewidth": 0.7,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
    }
    with plt.rc_context(rc_parameters):
        figure, axes = plt.subplots(
            2,
            len(SCHEMES),
            figsize=(7.4, 4.25),
            sharex=True,
            sharey=True,
        )
        top_image = None
        error_image = None
        extent = (0.0, 1.0, 0.0, 1.0)

        for column, ((_, label), field, error) in enumerate(
            zip(SCHEMES, uniform_fields, amr_errors)
        ):
            top_axis = axes[0, column]
            error_axis = axes[1, column]
            top_image = top_axis.imshow(
                field["density_percent"],
                origin="lower",
                extent=extent,
                cmap="cividis",
                vmin=top_min,
                vmax=top_max,
                interpolation="nearest",
                rasterized=True,
            )

            ny, nx = field["density_percent"].shape
            sample_x = np.linspace(4, nx - 5, 15, dtype=int)
            sample_y = np.linspace(4, ny - 5, 15, dtype=int)
            xx = (sample_x + 0.5) / nx
            yy = (sample_y + 0.5) / ny
            x_quiver, y_quiver = np.meshgrid(xx, yy, indexing="xy")
            u_quiver = field["u_perturbation"][np.ix_(sample_y, sample_x)]
            v_quiver = field["v_perturbation"][np.ix_(sample_y, sample_x)]
            perturbation_speed = np.sqrt(u_quiver * u_quiver + v_quiver * v_quiver)
            weak_velocity = perturbation_speed < 0.01 * u_inf
            u_quiver = np.ma.masked_where(weak_velocity, u_quiver)
            v_quiver = np.ma.masked_where(weak_velocity, v_quiver)
            arrows = top_axis.quiver(
                x_quiver,
                y_quiver,
                u_quiver,
                v_quiver,
                color="black",
                alpha=0.62,
                angles="xy",
                scale_units="xy",
                scale=430.0,
                width=0.0045,
                headwidth=3.2,
                headlength=4.2,
                headaxislength=3.8,
            )
            if column == 0:
                top_axis.quiverkey(
                    arrows,
                    0.69,
                    0.88,
                    0.1 * u_inf,
                    r"$0.1u_\infty$",
                    labelpos="N",
                    coordinates="axes",
                    fontproperties={"size": 7.0},
                )

            error_image = error_axis.imshow(
                error,
                origin="lower",
                extent=extent,
                cmap="RdBu_r",
                vmin=-error_limit,
                vmax=error_limit,
                interpolation="nearest",
                rasterized=True,
            )
            error_axis.axvline(
                0.5,
                color="black",
                linestyle="--",
                linewidth=0.8,
                alpha=0.85,
            )
            error_axis.text(
                0.25,
                0.96,
                "level 0",
                transform=error_axis.transAxes,
                ha="center",
                va="top",
                fontsize=7.0,
                color="black",
            )
            error_axis.text(
                0.75,
                0.96,
                "level 1",
                transform=error_axis.transAxes,
                ha="center",
                va="top",
                fontsize=7.0,
                color="black",
            )

            top_axis.set_title(label, pad=8.0)
            error_axis.set_xlabel(r"$x/L$")
            for row, axis in enumerate((top_axis, error_axis)):
                axis.set_aspect("equal")
                axis.set_xlim(0.0, 1.0)
                axis.set_ylim(0.0, 1.0)
                axis.set_xticks((0.25, 0.50, 0.75))
                axis.tick_params(direction="in", top=True, right=True)
                axis.text(
                    0.025,
                    0.975 if row == 0 else 0.025,
                    f"({chr(ord('a') + row * len(SCHEMES) + column)})",
                    transform=axis.transAxes,
                    ha="left",
                    va="top" if row == 0 else "bottom",
                    color="white" if row == 0 else "black",
                    fontsize=8.0,
                )

        axes[0, 0].set_ylabel(r"$y/L$")
        axes[1, 0].set_ylabel(r"$y/L$")
        figure.text(
            0.012,
            0.705,
            "Uniform field",
            rotation=90,
            ha="center",
            va="center",
            fontsize=8.5,
        )
        figure.text(
            0.012,
            0.275,
            "AMR error",
            rotation=90,
            ha="center",
            va="center",
            fontsize=8.5,
        )

        figure.subplots_adjust(
            left=0.075,
            right=0.905,
            bottom=0.11,
            top=0.91,
            wspace=0.14,
            hspace=0.19,
        )
        colorbar_top_axis = figure.add_axes((0.925, 0.545, 0.014, 0.315))
        colorbar_error_axis = figure.add_axes((0.925, 0.135, 0.014, 0.315))
        assert top_image is not None
        assert error_image is not None
        top_colorbar = figure.colorbar(top_image, cax=colorbar_top_axis)
        top_colorbar.set_label(
            "Density perturbation (%)",
            fontsize=8.0,
        )
        error_colorbar = figure.colorbar(error_image, cax=colorbar_error_axis)
        error_colorbar.set_label(
            "Density error (%)",
            fontsize=8.0,
        )
        top_colorbar.ax.tick_params(labelsize=7.0, direction="in")
        error_colorbar.ax.tick_params(labelsize=7.0, direction="in")

        png_output = args.png or args.root / (
            f"vortex_flowfield_N{args.resolution}.png"
        )
        pdf_output = args.pdf or args.root / (
            f"vortex_flowfield_N{args.resolution}.pdf"
        )
        png_output.parent.mkdir(parents=True, exist_ok=True)
        pdf_output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(png_output, dpi=300, bbox_inches="tight")
        figure.savefig(pdf_output, dpi=300, bbox_inches="tight")
        plt.close(figure)

        vorticity_png = args.vorticity_png or args.root / (
            f"vortex_vorticity_N{args.vorticity_resolution}.png"
        )
        vorticity_pdf = args.vorticity_pdf or args.root / (
            f"vortex_vorticity_N{args.vorticity_resolution}.pdf"
        )
        plot_vorticity_comparison(
            args.root,
            args.vorticity_resolution,
            vorticity_png,
            vorticity_pdf,
        )

    print(f"uniform density range: [{top_min:.6e}, {top_max:.6e}] percent")
    print(f"common AMR density-error limit: {error_limit:.6e} percent")
    for item in metadata:
        print(
            f"{item['label']}: {item['uniform_plotfile']} ; "
            f"{item['amr_plotfile']}"
        )
    print(f"wrote {png_output}")
    print(f"wrote {pdf_output}")

if __name__ == "__main__":
    yt.set_log_level(40)
    main()
