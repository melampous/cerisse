#!/usr/bin/env python3
"""Plot the Cartesian and axisymmetric RZ MMS convergence results.

The script reads archived CSV files only.  It does not run CERISSE and does
not alter the archived data.  The default inputs are the final Y9000X archives
used for the Chapter 2 verification study.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from typing import Iterable


os.environ.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")

SCRIPT_DIR = Path(__file__).resolve().parent
NUMERICS_DIR = SCRIPT_DIR.parents[1]

DEFAULT_CARTESIAN_ROOT = (
    NUMERICS_DIR
    / "mms_cartesian"
    / "results"
    / "y9000x_cartesian_mms_afd_hllc_wenoz5_20260727_r2"
)
DEFAULT_RZ_ROOT = (
    NUMERICS_DIR
    / "mms_rz_navier_stokes"
    / "results"
    / "y9000x_point_rz_ns_mms_20260727_r1"
)
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output"

RESOLUTIONS = (16, 32, 64, 128)
VARIABLES = (
    ("density", r"$\rho$", "#0072B2", "o"),
    ("x_momentum", r"$\rho u_x$", "#D55E00", "s"),
    ("y_momentum", r"$\rho u_y$", "#009E73", "^"),
    ("energy", r"$\rho E$", "#CC79A7", "D"),
)


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
            "axes.titlesize": 8.5,
            "axes.linewidth": 0.7,
            "legend.fontsize": 7.3,
            "xtick.labelsize": 7.3,
            "ytick.labelsize": 7.3,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "lines.linewidth": 1.10,
            "lines.markersize": 4.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
        }
    )
    return plt


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    by_resolution: dict[int, dict[str, str]] = {}
    for row in rows:
        resolution = int(row["base_n"])
        if resolution in by_resolution:
            raise ValueError(f"duplicate N={resolution} row in {path}")
        by_resolution[resolution] = row
    missing = set(RESOLUTIONS) - set(by_resolution)
    extra = set(by_resolution) - set(RESOLUTIONS)
    if missing or extra:
        raise ValueError(
            f"unexpected resolution set in {path}: "
            f"missing={sorted(missing)}, extra={sorted(extra)}"
        )
    return [by_resolution[n] for n in RESOLUTIONS]


def finite_positive(values: Iterable[float], label: str) -> list[float]:
    checked = [float(value) for value in values]
    if not all(math.isfinite(value) and value > 0.0 for value in checked):
        raise ValueError(f"{label} contains a non-positive or non-finite value")
    return checked


def extract(rows: list[dict[str, str]], column: str) -> list[float]:
    try:
        values = [float(row[column]) for row in rows]
    except KeyError as exc:
        raise ValueError(f"missing column {column}") from exc
    return finite_positive(values, column)


def style_axis(axis, show_xlabel: bool = True) -> None:
    axis.grid(
        True,
        which="major",
        color="0.86",
        linestyle=":",
        linewidth=0.6,
        zorder=0,
    )
    axis.grid(
        True,
        which="minor",
        color="0.93",
        linestyle=":",
        linewidth=0.45,
        zorder=0,
    )
    axis.set_xlim(1.0 / 14.0, 1.0 / 142.0)
    axis.set_xticks([1.0 / n for n in RESOLUTIONS])
    axis.set_xticklabels([rf"$1/{n}$" for n in RESOLUTIONS])
    if show_xlabel:
        axis.set_xlabel(r"Grid spacing, $h$")


def plot_global_l2(
    axis,
    rows: list[dict[str, str]],
    *,
    cylindrical: bool = False,
) -> None:
    h = [1.0 / n for n in RESOLUTIONS]
    cylindrical_columns = {
        "density": "density",
        "x_momentum": "axial_momentum",
        "y_momentum": "radial_momentum",
        "energy": "energy",
    }
    for key, _, colour, marker in VARIABLES:
        column_key = cylindrical_columns[key] if cylindrical else key
        values = extract(rows, f"{column_key}_L2")
        axis.loglog(
            h,
            values,
            color=colour,
            marker=marker,
            markerfacecolor="white",
            markeredgewidth=0.9,
            zorder=3,
        )


def add_order_guide(
    axis,
    exponent: int,
    anchor_h: float,
    anchor_error: float,
    *,
    label_x: float,
    label_factor: float = 1.0,
) -> None:
    h = [1.0 / n for n in RESOLUTIONS]
    guide = [
        anchor_error * (grid_spacing / anchor_h) ** exponent
        for grid_spacing in h
    ]
    axis.loglog(
        h,
        guide,
        color="0.35",
        linestyle=(0, (1.5, 1.8)),
        linewidth=0.9,
        zorder=2,
    )
    label_y = (
        anchor_error
        * (label_x / anchor_h) ** exponent
        * label_factor
    )
    axis.text(
        label_x,
        label_y,
        rf"$\mathcal{{O}}(h^{exponent})$",
        color="0.25",
        ha="center",
        va="bottom",
        fontsize=7.3,
    )


def make_figure(
    cartesian_root: Path,
    rz_root: Path,
    output_dir: Path,
) -> tuple[Path, Path]:
    plt = configure_matplotlib()
    from matplotlib.lines import Line2D

    euler_rows = read_rows(cartesian_root / "euler_errors.csv")
    ns_rows = read_rows(cartesian_root / "navier_stokes_errors.csv")
    rz_evolution_rows = read_rows(rz_root / "evolution_errors.csv")
    rz_rhs_rows = read_rows(rz_root / "rhs_errors.csv")

    # The thesis text block is 150 mm wide.  Generate at the final printed
    # width so that the nominal font sizes are not reduced by LaTeX.
    figure = plt.figure(figsize=(5.9055, 4.82))
    grid = figure.add_gridspec(
        2,
        2,
        left=0.105,
        right=0.985,
        bottom=0.105,
        top=0.875,
        wspace=0.30,
        hspace=0.43,
    )
    euler_axis = figure.add_subplot(grid[0, 0])
    ns_axis = figure.add_subplot(grid[0, 1])
    rz_axis = figure.add_subplot(grid[1, 0])
    axis_ring_axis = figure.add_subplot(grid[1, 1])

    for axis, rows, cylindrical in (
        (euler_axis, euler_rows, False),
        (ns_axis, ns_rows, False),
        (rz_axis, rz_evolution_rows, True),
    ):
        plot_global_l2(axis, rows, cylindrical=cylindrical)
        style_axis(axis)

    euler_axis.set_title("(a) Cartesian Euler")
    ns_axis.set_title("(b) Cartesian Navier--Stokes")
    rz_axis.set_title("(c) Axisymmetric RZ Navier--Stokes")
    axis_ring_axis.set_title("(d) RZ first-ring residual")
    euler_axis.set_ylabel(r"Global $L_2$ error")
    ns_axis.set_ylabel(r"Global $L_2$ error")
    rz_axis.set_ylabel(r"Global $L_2$ error")
    axis_ring_axis.set_ylabel(
        r"$\|e_{\mathcal{L},\rho u_r}\|_{\infty,i=0}$"
    )

    add_order_guide(
        euler_axis,
        5,
        1.0 / 16.0,
        2.2e-7,
        label_x=1.0 / 24.0,
        label_factor=1.05,
    )
    add_order_guide(
        ns_axis,
        2,
        1.0 / 16.0,
        9.0e-7,
        label_x=1.0 / 24.0,
        label_factor=0.85,
    )

    h = [1.0 / n for n in RESOLUTIONS]
    first_ring_rhs = extract(
        rz_rhs_rows,
        "radial_momentum_first_ring_Linf",
    )
    axis_ring_axis.loglog(
        h,
        first_ring_rhs,
        color="0.10",
        linestyle=(0, (5.0, 2.1)),
        marker="X",
        markerfacecolor="0.10",
        markeredgewidth=0.8,
        linewidth=1.25,
        zorder=4,
    )
    add_order_guide(
        rz_axis,
        2,
        1.0 / 16.0,
        5.0e-7,
        label_x=1.0 / 25.0,
        label_factor=0.52,
    )
    add_order_guide(
        axis_ring_axis,
        1,
        1.0 / 16.0,
        3.2e-4,
        label_x=1.0 / 43.0,
        label_factor=1.18,
    )

    # Plotting additional log-log curves can replace Matplotlib's locator.
    # Reapply the fixed grid-spacing ticks after all curves are present.
    for axis in (euler_axis, ns_axis, rz_axis, axis_ring_axis):
        style_axis(axis)

    variable_handles = [
        Line2D(
            [],
            [],
            color=colour,
            marker=marker,
            markerfacecolor="white",
            markeredgewidth=0.9,
            label=label,
        )
        for _, label, colour, marker in VARIABLES
    ]
    figure.legend(
        handles=variable_handles,
        loc="upper center",
        bbox_to_anchor=(0.52, 0.995),
        ncol=4,
        frameon=False,
        columnspacing=1.10,
        handletextpad=0.45,
        handlelength=2.0,
    )
    rz_handles = [
        Line2D(
            [],
            [],
            color=colour,
            marker=marker,
            markerfacecolor="white",
            markeredgewidth=0.9,
            label=label,
        )
        for label, (_, _, colour, marker) in zip(
            (r"$\rho$", r"$\rho u_x$", r"$\rho u_r$", r"$\rho E$"),
            VARIABLES,
        )
    ]
    rz_axis.legend(
        handles=rz_handles,
        loc="lower left",
        ncol=2,
        frameon=False,
        columnspacing=0.9,
        handletextpad=0.35,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "mms_convergence.png"
    pdf_path = output_dir / "mms_convergence.pdf"
    figure.savefig(png_path, dpi=300)
    figure.savefig(pdf_path)
    plt.close(figure)

    return pdf_path, png_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cartesian-root",
        type=Path,
        default=DEFAULT_CARTESIAN_ROOT,
        help="Cartesian MMS result archive",
    )
    parser.add_argument(
        "--rz-root",
        type=Path,
        default=DEFAULT_RZ_ROOT,
        help="axisymmetric RZ MMS result archive",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="directory for the PDF and PNG",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pdf_path, png_path = make_figure(
        args.cartesian_root.resolve(),
        args.rz_root.resolve(),
        args.output_dir.resolve(),
    )
    print(f"wrote {pdf_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
