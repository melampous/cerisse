#!/usr/bin/env python3
"""Summarise and plot the low-temporal-error vortex convergence study.

The figure deliberately treats the 45-degree transport case as the primary
spatial-order test.  Axis-aligned density data are retained in a separate
panel because their finest-grid contact-mode anomaly is not representative of
the velocity and pressure convergence.

Example
-------
python3 plot_low_temporal_error_convergence.py \
    results/y9000x_afd_hllc_wenoz5_spatial_20260727 \
    results/y9000x_afd_hllc_wenoz5_spatial_angle45_20260727 \
    --output-dir results/y9000x_afd_hllc_wenoz5_spatial_summary_20260727
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from typing import Any


os.environ.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")


MAIN_RESOLUTIONS = (40, 80, 160, 320)
METRICS = (
    (
        "density_L2",
        r"Density, $\rho$",
        "#0072B2",
        "o",
    ),
    (
        "velocity_vector_L2",
        r"Velocity, $\mathbf{u}$",
        "#D55E00",
        "s",
    ),
    (
        "pressure_L2",
        r"Pressure, $p$",
        "#009E73",
        "^",
    ),
)


def read_csv_rows(source: Path) -> list[dict[str, str]]:
    """Read a non-empty CSV and preserve its original string fields."""

    if not source.is_file():
        raise FileNotFoundError(f"missing required CSV: {source}")
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{source} has no CSV header")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{source} has no data rows")
    return rows


def finite_float(row: dict[str, str], field: str, source: Path) -> float:
    """Return a finite floating-point CSV field with a useful error message."""

    try:
        value = float(row[field])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"{source} has an invalid {field!r} value") from exc
    if not math.isfinite(value):
        raise ValueError(f"{source} has a non-finite {field!r} value")
    return value


def integer_value(row: dict[str, str], field: str, source: Path) -> int:
    value = finite_float(row, field, source)
    integer = int(value)
    if value != integer:
        raise ValueError(f"{source} has a non-integer {field!r} value")
    return integer


def load_main_series(root: Path) -> dict[int, dict[str, Any]]:
    """Load and validate the four-grid main convergence sequence."""

    source = root / "errors.csv"
    source_rows = read_csv_rows(source)
    indexed: dict[int, dict[str, Any]] = {}
    for source_row in source_rows:
        resolution = integer_value(source_row, "base_n", source)
        base_ny = integer_value(source_row, "base_ny", source)
        max_level = integer_value(source_row, "max_level", source)
        if resolution != base_ny:
            raise ValueError(f"{source} contains a non-square grid at N={resolution}")
        if max_level != 0:
            raise ValueError(f"{source} contains max_level={max_level}; expected 0")
        if resolution in indexed:
            raise ValueError(f"{source} contains duplicate N={resolution} rows")

        row: dict[str, Any] = {
            "base_n": resolution,
            "base_ny": base_ny,
            "max_level": max_level,
            "source_csv": str(source),
            "time": finite_float(source_row, "time", source),
            "density_L1": finite_float(source_row, "density_L1", source),
            "density_L2": finite_float(source_row, "density_L2", source),
            "density_Linf": finite_float(source_row, "density_Linf", source),
            "velocity_vector_L2": finite_float(
                source_row, "velocity_vector_L2", source
            ),
            "pressure_L2": finite_float(source_row, "pressure_L2", source),
            "mass_relative_change_from_initial": finite_float(
                source_row, "mass_relative_change_from_initial", source
            ),
            "energy_relative_change_from_initial": finite_float(
                source_row, "energy_relative_change_from_initial", source
            ),
            "perturbation_ke_ratio": finite_float(
                source_row, "perturbation_ke_ratio", source
            ),
        }
        indexed[resolution] = row

    found = tuple(sorted(indexed))
    if found != MAIN_RESOLUTIONS:
        raise ValueError(
            f"{source} must contain N={MAIN_RESOLUTIONS}; found N={found}"
        )
    times = [float(indexed[resolution]["time"]) for resolution in MAIN_RESOLUTIONS]
    if not all(
        math.isclose(value, times[0], rel_tol=0.0, abs_tol=2.0e-16)
        for value in times[1:]
    ):
        raise ValueError(f"{source} does not use one common final time")
    return indexed


def observed_order(
    coarse_error: float,
    fine_error: float,
    coarse_n: int,
    fine_n: int,
) -> float:
    return math.log(coarse_error / fine_error) / math.log(fine_n / coarse_n)


def enrich_main_series(
    direction: str,
    indexed: dict[int, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Add normalised errors and observed orders to one main sequence."""

    output: list[dict[str, Any]] = []
    for index, resolution in enumerate(MAIN_RESOLUTIONS):
        row = dict(indexed[resolution])
        row["direction"] = direction
        row["previous_n"] = "" if index == 0 else MAIN_RESOLUTIONS[index - 1]
        row["flag"] = (
            "axis_density_contact_mode_anomaly"
            if direction == "axis" and resolution == 320
            else ""
        )
        for metric, _, _, _ in METRICS:
            row[f"{metric}_normalised_to_N40"] = float(row[metric]) / float(
                indexed[40][metric]
            )
            if index == 0:
                row[f"{metric}_observed_order"] = float("nan")
            else:
                coarse_n = MAIN_RESOLUTIONS[index - 1]
                row[f"{metric}_observed_order"] = observed_order(
                    float(indexed[coarse_n][metric]),
                    float(row[metric]),
                    coarse_n,
                    resolution,
                )
        output.append(row)
    return output


def write_main_summary(output: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "direction",
        "base_n",
        "previous_n",
        "time",
        "density_L1",
        "density_L2",
        "density_L2_normalised_to_N40",
        "density_L2_observed_order",
        "density_Linf",
        "velocity_vector_L2",
        "velocity_vector_L2_normalised_to_N40",
        "velocity_vector_L2_observed_order",
        "pressure_L2",
        "pressure_L2_normalised_to_N40",
        "pressure_L2_observed_order",
        "mass_relative_change_from_initial",
        "energy_relative_change_from_initial",
        "perturbation_ke_ratio",
        "flag",
        "source_csv",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)


def read_single_case(source: Path) -> dict[str, Any]:
    rows = read_csv_rows(source)
    if len(rows) != 1:
        raise ValueError(f"{source} must contain exactly one data row")
    source_row = rows[0]
    row: dict[str, Any] = {
        "base_n": integer_value(source_row, "base_n", source),
        "source_csv": str(source),
    }
    for metric, _, _, _ in METRICS:
        row[metric] = finite_float(source_row, metric, source)
    return row


def load_supplemental_axis_point(axis_root: Path) -> dict[str, Any] | None:
    source = axis_root / "N240" / "errors.csv"
    if not source.is_file():
        return None
    point = read_single_case(source)
    if int(point["base_n"]) != 240:
        raise ValueError(f"{source} is not the expected N=240 supplemental case")
    return point


def load_audits(
    direction: str,
    root: Path,
    main: dict[int, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Load optional one-row audit CSVs and compare them with their main runs."""

    audit_root = root / "audits"
    if not audit_root.is_dir():
        return []
    output: list[dict[str, Any]] = []
    for source in sorted(audit_root.glob("*/errors.csv")):
        row = read_single_case(source)
        resolution = int(row["base_n"])
        if resolution not in main:
            raise ValueError(f"{source} has no N={resolution} main-run reference")
        reference = main[resolution]
        output_row: dict[str, Any] = {
            "direction": direction,
            "run": source.parent.name,
            "base_n": resolution,
            "reference_run": f"main_N{resolution}",
            "source_csv": str(source),
        }
        for metric, _, _, _ in METRICS:
            error = float(row[metric])
            reference_error = float(reference[metric])
            output_row[metric] = error
            output_row[f"{metric}_relative_difference_from_main"] = (
                error / reference_error - 1.0
            )
        output.append(output_row)
    return output


def write_audit_summary(
    output: Path,
    audits: list[dict[str, Any]],
) -> None:
    fields = [
        "direction",
        "run",
        "base_n",
        "reference_run",
        "density_L2",
        "density_L2_relative_difference_from_main",
        "velocity_vector_L2",
        "velocity_vector_L2_relative_difference_from_main",
        "pressure_L2",
        "pressure_L2_relative_difference_from_main",
        "source_csv",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in audits)


def write_axis_anomaly_summary(
    output: Path,
    axis_main: dict[int, dict[str, Any]],
    supplemental: dict[str, Any] | None,
    axis_audits: list[dict[str, Any]],
) -> None:
    """Document the extra N=240 point and the finest-grid audit envelope."""

    n320_audits = [row for row in axis_audits if int(row["base_n"]) == 320]
    density_values = [float(axis_main[320]["density_L2"])] + [
        float(row["density_L2"]) for row in n320_audits
    ]
    density_spread = (
        (max(density_values) - min(density_values)) / float(axis_main[320]["density_L2"])
        if len(density_values) > 1
        else float("nan")
    )

    rows: list[dict[str, Any]] = []
    if supplemental is not None:
        rows.append(
            {
                "run": "supplemental_N240",
                "base_n": 240,
                "density_L2": supplemental["density_L2"],
                "comparison_coarse_n": 160,
                "density_L2_observed_order_from_comparison": observed_order(
                    float(axis_main[160]["density_L2"]),
                    float(supplemental["density_L2"]),
                    160,
                    240,
                ),
                "density_L2_ratio_to_N240": 1.0,
                "N320_audit_count_excluding_main": len(n320_audits),
                "N320_density_audit_relative_spread": density_spread,
                "source_csv": supplemental["source_csv"],
                "interpretation": "supplemental_resolution_before_anomaly",
            }
        )
        comparison_n = 240
        comparison_error = float(supplemental["density_L2"])
    else:
        comparison_n = 160
        comparison_error = float(axis_main[160]["density_L2"])

    rows.append(
        {
            "run": "main_N320",
            "base_n": 320,
            "density_L2": axis_main[320]["density_L2"],
            "comparison_coarse_n": comparison_n,
            "density_L2_observed_order_from_comparison": observed_order(
                comparison_error,
                float(axis_main[320]["density_L2"]),
                comparison_n,
                320,
            ),
            "density_L2_ratio_to_N240": (
                float(axis_main[320]["density_L2"])
                / float(supplemental["density_L2"])
                if supplemental is not None
                else float("nan")
            ),
            "N320_audit_count_excluding_main": len(n320_audits),
            "N320_density_audit_relative_spread": density_spread,
            "source_csv": axis_main[320]["source_csv"],
            "interpretation": "axis_density_contact_mode_anomaly",
        }
    )
    fields = [
        "run",
        "base_n",
        "density_L2",
        "comparison_coarse_n",
        "density_L2_observed_order_from_comparison",
        "density_L2_ratio_to_N240",
        "N320_audit_count_excluding_main",
        "N320_density_audit_relative_spread",
        "interpretation",
        "source_csv",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_figure(
    png_output: Path,
    pdf_output: Path,
    axis_main: dict[int, dict[str, Any]],
    diagonal_main: dict[int, dict[str, Any]],
    supplemental_axis: dict[str, Any] | None,
    axis_audits: list[dict[str, Any]],
    dpi: int,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    resolutions = list(MAIN_RESOLUTIONS)
    rc_parameters = {
        "font.family": "serif",
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9.0,
        "axes.linewidth": 0.8,
        "legend.fontsize": 7.4,
        "xtick.labelsize": 7.8,
        "ytick.labelsize": 7.8,
        "lines.linewidth": 1.35,
        "lines.markersize": 5.0,
        "mathtext.fontset": "dejavuserif",
        "savefig.facecolor": "white",
    }
    with plt.rc_context(rc_parameters):
        figure, axes = plt.subplots(1, 3, figsize=(9.4, 3.05))

        # (a) Primary diagonal-transport spatial convergence.
        metric_handles: list[Line2D] = []
        for metric, label, color, marker in METRICS:
            normalised = [
                float(diagonal_main[resolution][metric])
                / float(diagonal_main[40][metric])
                for resolution in resolutions
            ]
            axes[0].plot(
                resolutions,
                normalised,
                color=color,
                marker=marker,
                markerfacecolor="white",
                markeredgewidth=0.9,
                zorder=3,
            )
            metric_handles.append(
                Line2D(
                    [0],
                    [0],
                    color=color,
                    marker=marker,
                    markerfacecolor="white",
                    markeredgewidth=0.9,
                    label=label,
                )
            )
        fifth_order = [(resolution / 40.0) ** -5 for resolution in resolutions]
        axes[0].plot(
            resolutions,
            fifth_order,
            color="0.25",
            linestyle=(0, (4, 2)),
            linewidth=1.15,
            zorder=2,
        )
        metric_handles.append(
            Line2D(
                [0],
                [0],
                color="0.25",
                linestyle=(0, (4, 2)),
                linewidth=1.15,
                label=r"Fifth order, $N^{-5}$",
            )
        )
        axes[0].set_title(r"Diagonal transport ($45^\circ$)")
        axes[0].set_ylabel(r"Normalised $L_2$ error, $E_N/E_{40}$")

        # (b) Orders from the same primary diagonal sequence.
        fine_resolutions = resolutions[1:]
        for metric, _, color, marker in METRICS:
            rates = [
                observed_order(
                    float(diagonal_main[coarse_n][metric]),
                    float(diagonal_main[fine_n][metric]),
                    coarse_n,
                    fine_n,
                )
                for coarse_n, fine_n in zip(resolutions[:-1], resolutions[1:])
            ]
            axes[1].plot(
                range(len(rates)),
                rates,
                color=color,
                marker=marker,
                markerfacecolor="white",
                markeredgewidth=0.9,
                zorder=3,
            )
        axes[1].axhline(
            5.0,
            color="0.25",
            linestyle=(0, (4, 2)),
            linewidth=1.15,
            zorder=2,
        )
        axes[1].set_xticks(
            range(len(fine_resolutions)),
            [
                f"{coarse}\N{RIGHTWARDS ARROW}{fine}"
                for coarse, fine in zip(resolutions[:-1], fine_resolutions)
            ],
        )
        axes[1].set_ylim(4.2, 5.95)
        axes[1].set_ylabel(r"Observed order, $p$")
        axes[1].set_xlabel("Refinement interval")
        axes[1].set_title("Diagonal observed orders")

        # (c) Density-only comparison that exposes, but does not generalise,
        # the axis-aligned contact-mode anomaly.
        diagonal_density = [
            float(diagonal_main[resolution]["density_L2"])
            / float(diagonal_main[40]["density_L2"])
            for resolution in resolutions
        ]
        axis_resolutions = [40, 80, 160]
        if supplemental_axis is not None:
            axis_resolutions.append(240)
        axis_resolutions.append(320)
        axis_density = []
        for resolution in axis_resolutions:
            error = (
                float(supplemental_axis["density_L2"])
                if resolution == 240 and supplemental_axis is not None
                else float(axis_main[resolution]["density_L2"])
            )
            axis_density.append(error / float(axis_main[40]["density_L2"]))

        axes[2].plot(
            resolutions,
            diagonal_density,
            color="#0072B2",
            linestyle="-",
            marker="o",
            markerfacecolor="#0072B2",
            markeredgewidth=0.8,
            label=r"Diagonal ($45^\circ$)",
            zorder=3,
        )
        axes[2].plot(
            axis_resolutions,
            axis_density,
            color="#D55E00",
            linestyle="--",
            marker="s",
            markerfacecolor="white",
            markeredgewidth=0.9,
            label="Axis aligned",
            zorder=3,
        )
        axes[2].scatter(
            [320],
            [axis_density[-1]],
            s=68,
            facecolors="none",
            edgecolors="#D55E00",
            linewidths=1.15,
            zorder=4,
        )
        axes[2].set_title("Axis aligned density check")
        axes[2].set_ylabel(r"Normalised density $L_2$ error")
        axes[2].legend(
            loc="lower left",
            frameon=False,
            handlelength=2.2,
            borderaxespad=0.15,
        )

        axis_velocity_order = observed_order(
            float(axis_main[160]["velocity_vector_L2"]),
            float(axis_main[320]["velocity_vector_L2"]),
            160,
            320,
        )
        axis_pressure_order = observed_order(
            float(axis_main[160]["pressure_L2"]),
            float(axis_main[320]["pressure_L2"]),
            160,
            320,
        )
        n320_audits = [row for row in axis_audits if int(row["base_n"]) == 320]
        audit_values = [float(axis_main[320]["density_L2"])] + [
            float(row["density_L2"]) for row in n320_audits
        ]
        audit_spread_percent = (
            100.0
            * (max(audit_values) - min(audit_values))
            / float(axis_main[320]["density_L2"])
            if len(audit_values) > 1
            else float("nan")
        )
        audit_line = (
            f"{len(n320_audits)} audits span {audit_spread_percent:.3f}%"
            if math.isfinite(audit_spread_percent)
            else "No finest-grid audit available"
        )
        if supplemental_axis is not None:
            density_ratio = float(axis_main[320]["density_L2"]) / float(
                supplemental_axis["density_L2"]
            )
            density_ratio_line = (
                rf"$E_{{\rho,320}}/E_{{\rho,240}}={density_ratio:.2f}$"
            )
        else:
            density_ratio_line = r"No $N=240$ comparison"
        axes[2].text(
            0.97,
            0.95,
            "Axis $N=320$ result\n"
            + density_ratio_line
            + "\n"
            + rf"$160\to320:\ p_{{\mathbf{{u}},2}}={axis_velocity_order:.2f}$, "
            + rf"$p_p={axis_pressure_order:.2f}$"
            + "\n"
            + audit_line,
            transform=axes[2].transAxes,
            ha="right",
            va="top",
            fontsize=6.9,
            linespacing=1.25,
            bbox={
                "facecolor": "white",
                "edgecolor": "0.8",
                "linewidth": 0.55,
                "alpha": 0.94,
                "pad": 2.0,
            },
            zorder=5,
        )

        for index, axis in enumerate((axes[0], axes[2])):
            axis.set_xscale("log", base=2)
            axis.set_yscale("log")
            ticks = resolutions if index == 0 else sorted(set(axis_resolutions + resolutions))
            axis.set_xticks(ticks, [str(value) for value in ticks])
            axis.set_xlabel(r"Resolution, $N$")
        axes[0].set_ylim(8.0e-6, 1.8)
        axes[2].set_ylim(8.0e-6, 1.8)

        for index, axis in enumerate(axes):
            axis.grid(
                axis="y",
                which="major",
                color="0.84",
                linestyle=":",
                linewidth=0.65,
                zorder=0,
            )
            axis.tick_params(direction="in", which="both", top=True, right=True)
            axis.text(
                -0.12,
                1.03,
                f"({chr(ord('a') + index)})",
                transform=axis.transAxes,
                ha="left",
                va="bottom",
                fontsize=8.5,
            )

        figure.legend(
            handles=metric_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.995),
            ncol=4,
            frameon=False,
            handlelength=2.1,
            columnspacing=1.35,
        )
        figure.subplots_adjust(
            left=0.07,
            right=0.995,
            bottom=0.18,
            top=0.78,
            wspace=0.36,
        )
        png_output.parent.mkdir(parents=True, exist_ok=True)
        pdf_output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(
            png_output,
            dpi=dpi,
            bbox_inches="tight",
            metadata={
                "Title": "Low-temporal-error isentropic-vortex convergence",
                "Author": "Cerisse verification workflow",
            },
        )
        figure.savefig(
            pdf_output,
            dpi=dpi,
            bbox_inches="tight",
            metadata={
                "Title": "Low-temporal-error isentropic-vortex convergence",
                "Author": "Cerisse verification workflow",
                "Subject": "AFD-HLLC-WENO-Z5 spatial convergence",
            },
        )
        plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create reusable CSV summaries and a three-panel convergence figure "
            "from axis-aligned and 45-degree isentropic-vortex studies."
        )
    )
    parser.add_argument(
        "axis_root",
        type=Path,
        help="axis-aligned study root containing errors.csv",
    )
    parser.add_argument(
        "diagonal_root",
        type=Path,
        help="45-degree study root containing errors.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="output directory; default: DIAGONAL_ROOT.parent/spatial_summary",
    )
    parser.add_argument(
        "--stem",
        default="afd_hllc_wenoz5_low_temporal_error_convergence",
        help="figure filename stem",
    )
    parser.add_argument("--dpi", type=int, default=400, help="PNG resolution")
    args = parser.parse_args()
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")

    output_dir = args.output_dir or args.diagonal_root.parent / "spatial_summary"
    axis_main = load_main_series(args.axis_root)
    diagonal_main = load_main_series(args.diagonal_root)
    if not math.isclose(
        float(axis_main[40]["time"]),
        float(diagonal_main[40]["time"]),
        rel_tol=0.0,
        abs_tol=2.0e-16,
    ):
        raise ValueError("axis and diagonal studies do not share the same final time")

    main_rows = enrich_main_series("axis", axis_main)
    main_rows.extend(enrich_main_series("diagonal_45deg", diagonal_main))
    supplemental_axis = load_supplemental_axis_point(args.axis_root)
    axis_audits = load_audits("axis", args.axis_root, axis_main)
    diagonal_audits = load_audits(
        "diagonal_45deg", args.diagonal_root, diagonal_main
    )

    convergence_csv = output_dir / "convergence_summary.csv"
    audit_csv = output_dir / "audit_summary.csv"
    anomaly_csv = output_dir / "axis_density_contact_anomaly_summary.csv"
    png_output = output_dir / f"{args.stem}.png"
    pdf_output = output_dir / f"{args.stem}.pdf"
    write_main_summary(convergence_csv, main_rows)
    write_audit_summary(audit_csv, axis_audits + diagonal_audits)
    write_axis_anomaly_summary(
        anomaly_csv,
        axis_main,
        supplemental_axis,
        axis_audits,
    )
    plot_figure(
        png_output,
        pdf_output,
        axis_main,
        diagonal_main,
        supplemental_axis,
        axis_audits,
        args.dpi,
    )

    print("Diagonal 45-degree observed orders:")
    for metric, label, _, _ in METRICS:
        rates = [
            observed_order(
                float(diagonal_main[coarse_n][metric]),
                float(diagonal_main[fine_n][metric]),
                coarse_n,
                fine_n,
            )
            for coarse_n, fine_n in zip(
                MAIN_RESOLUTIONS[:-1], MAIN_RESOLUTIONS[1:]
            )
        ]
        print(f"  {label:24s}: " + ", ".join(f"{rate:.4f}" for rate in rates))
    print("Axis-aligned 160->320 observed orders:")
    for metric, label, _, _ in METRICS:
        rate = observed_order(
            float(axis_main[160][metric]),
            float(axis_main[320][metric]),
            160,
            320,
        )
        print(f"  {label:24s}: {rate:.4f}")
    for output in (
        convergence_csv,
        audit_csv,
        anomaly_csv,
        png_output,
        pdf_output,
    ):
        print(f"wrote {output}")


if __name__ == "__main__":
    main()
