#!/usr/bin/env python3
"""Combine and plot the four-scheme isentropic-vortex method matrix.

The input files are the ``errors.csv`` files produced by ``analyze.py`` and
``run_method_matrix.sh``.  The AMR sequence is reported against the base-grid
resolution.  Its observed rates assess the complete AMR algorithm and should
not be interpreted as the formal order of the interior reconstruction.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/tmp/cerisse-matplotlib")


SCHEMES = (
    "llf-wenoz5",
    "llf-teno5",
    "afd-hllc-wenoz5",
    "afd-hllc-teno5",
)
SUITES = ("uniform", "amr")
RESOLUTIONS = (40, 80, 160)

SCHEME_STYLES = {
    "llf-wenoz5": ("#0072B2", "o", "LLF WENO-Z5"),
    "llf-teno5": ("#009E73", "s", "LLF TENO5"),
    "afd-hllc-wenoz5": ("#D55E00", "^", "AFD HLLC WENO-Z5"),
    "afd-hllc-teno5": ("#CC79A7", "D", "AFD HLLC TENO5"),
}

DERIVED_FIELDS = (
    "density_L2_observed_rate",
    "perturbation_ke_deficit",
    "density_L2_ratio_to_same_scheme_uniform",
    "velocity_vector_L2_ratio_to_same_scheme_uniform",
    "pressure_L2_ratio_to_same_scheme_uniform",
    "perturbation_ke_deficit_ratio_to_same_scheme_uniform",
    "active_cells_ratio_to_same_scheme_uniform",
    "advance_seconds_ratio_to_same_scheme_uniform",
)


def finite_float(row: dict[str, str], field: str, source: Path) -> float:
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


def read_source_csv(source: Path) -> tuple[list[str], list[dict[str, str]]]:
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{source} has no CSV header")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{source} contains no data rows")
    return list(reader.fieldnames), rows


def matching_value(
    rows: list[dict[str, object]], field: str, description: str
) -> None:
    baseline = float(rows[0][field])
    for row in rows[1:]:
        value = float(row[field])
        scale = max(1.0, abs(baseline), abs(value))
        if not math.isclose(value, baseline, rel_tol=0.0, abs_tol=1.0e-12 * scale):
            raise ValueError(
                f"{description} must have matching {field}: "
                f"{baseline:.17g} and {value:.17g} differ"
            )


def load_matrix(
    root: Path,
) -> tuple[list[str], list[dict[str, object]]]:
    sources = {
        (scheme, suite): root / scheme / suite / "errors.csv"
        for scheme in SCHEMES
        for suite in SUITES
    }
    missing = [source for source in sources.values() if not source.is_file()]
    if missing:
        paths = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(f"the method matrix is incomplete:\n{paths}")

    common_header: list[str] | None = None
    matrix_rows: list[dict[str, object]] = []
    for scheme in SCHEMES:
        for suite in SUITES:
            source = sources[(scheme, suite)]
            header, source_rows = read_source_csv(source)
            if common_header is None:
                common_header = header
            elif header != common_header:
                raise ValueError(f"{source} does not use the common CSV schema")

            rows_by_resolution: dict[int, dict[str, str]] = {}
            for source_row in source_rows:
                base_n = integer_value(source_row, "base_n", source)
                base_ny = integer_value(source_row, "base_ny", source)
                max_level = integer_value(source_row, "max_level", source)
                if base_n != base_ny:
                    raise ValueError(
                        f"{source} contains a non-square base grid "
                        f"{base_n} by {base_ny}"
                    )
                if base_n in rows_by_resolution:
                    raise ValueError(f"{source} contains duplicate N={base_n} rows")
                rows_by_resolution[base_n] = source_row

                expected_level = 0 if suite == "uniform" else 1
                if max_level != expected_level:
                    raise ValueError(
                        f"{source} reports max_level={max_level}; "
                        f"{suite} requires {expected_level}"
                    )

            if tuple(sorted(rows_by_resolution)) != RESOLUTIONS:
                raise ValueError(
                    f"{source} must contain N={RESOLUTIONS}; found "
                    f"{tuple(sorted(rows_by_resolution))}"
                )

            previous_row: dict[str, object] | None = None
            for base_n in RESOLUTIONS:
                source_row = rows_by_resolution[base_n]
                row: dict[str, object] = {
                    "scheme": scheme,
                    "suite": suite,
                    "source_csv": str(source),
                }
                for field in header:
                    output_field = (
                        "source_density_L2_rate"
                        if field == "density_L2_rate"
                        else field
                    )
                    row[output_field] = source_row[field]

                density_l2 = finite_float(source_row, "density_L2", source)
                if previous_row is None:
                    observed_rate = float("nan")
                else:
                    coarse_n = int(previous_row["base_n"])
                    coarse_error = float(previous_row["density_L2"])
                    observed_rate = math.log(coarse_error / density_l2) / math.log(
                        base_n / coarse_n
                    )
                row["base_n"] = base_n
                row["base_ny"] = integer_value(source_row, "base_ny", source)
                row["max_level"] = integer_value(
                    source_row, "max_level", source
                )
                row["active_cells"] = integer_value(
                    source_row, "active_cells", source
                )
                for field in (
                    "time",
                    "x_lo",
                    "x_hi",
                    "y_lo",
                    "y_hi",
                    "advance_seconds",
                    "density_L1",
                    "density_L2",
                    "density_Linf",
                    "u_L2",
                    "v_L2",
                    "velocity_vector_L2",
                    "pressure_L2",
                    "mass_relative_change_from_initial",
                    "energy_relative_change_from_initial",
                    "perturbation_ke_ratio",
                ):
                    row[field] = finite_float(source_row, field, source)
                row["density_L2_observed_rate"] = observed_rate
                row["perturbation_ke_deficit"] = abs(
                    1.0 - float(row["perturbation_ke_ratio"])
                )
                matrix_rows.append(row)
                previous_row = row

    assert common_header is not None

    for field in ("time", "x_lo", "x_hi", "y_lo", "y_hi"):
        matching_value(matrix_rows, field, "all method-matrix rows")

    indexed = {
        (str(row["scheme"]), str(row["suite"]), int(row["base_n"])): row
        for row in matrix_rows
    }
    for row in matrix_rows:
        uniform = indexed[(str(row["scheme"]), "uniform", int(row["base_n"]))]
        row["density_L2_ratio_to_same_scheme_uniform"] = float(
            row["density_L2"]
        ) / float(uniform["density_L2"])
        row["velocity_vector_L2_ratio_to_same_scheme_uniform"] = float(
            row["velocity_vector_L2"]
        ) / float(uniform["velocity_vector_L2"])
        row["pressure_L2_ratio_to_same_scheme_uniform"] = float(
            row["pressure_L2"]
        ) / float(uniform["pressure_L2"])
        row["perturbation_ke_deficit_ratio_to_same_scheme_uniform"] = float(
            row["perturbation_ke_deficit"]
        ) / float(uniform["perturbation_ke_deficit"])
        row["active_cells_ratio_to_same_scheme_uniform"] = int(
            row["active_cells"]
        ) / int(uniform["active_cells"])
        row["advance_seconds_ratio_to_same_scheme_uniform"] = float(
            row["advance_seconds"]
        ) / float(uniform["advance_seconds"])

    output_header = [
        "scheme",
        "suite",
        "source_csv",
        *(
            "source_density_L2_rate" if field == "density_L2_rate" else field
            for field in common_header
        ),
        *DERIVED_FIELDS,
    ]
    return output_header, matrix_rows


def write_summary(
    output: Path, header: list[str], rows: list[dict[str, object]]
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


def plot_summary(output: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    indexed = {
        (str(row["scheme"]), str(row["suite"]), int(row["base_n"])): row
        for row in rows
    }
    suite_styles = {
        "uniform": ("-", "Uniform grid", True),
        "amr": ("--", "AMR", False),
    }
    rc_parameters = {
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.2,
    }
    with plt.rc_context(rc_parameters):
        figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))
        for scheme in SCHEMES:
            color, marker, _ = SCHEME_STYLES[scheme]
            for suite in SUITES:
                line_style, _, filled = suite_styles[suite]
                selected = [
                    indexed[(scheme, suite, resolution)]
                    for resolution in RESOLUTIONS
                ]
                marker_face = color if filled else "white"
                axes[0].plot(
                    RESOLUTIONS,
                    [float(row["density_L2"]) for row in selected],
                    color=color,
                    linestyle=line_style,
                    marker=marker,
                    markersize=4.5,
                    markerfacecolor=marker_face,
                    markeredgecolor=color,
                    markeredgewidth=0.8,
                )
                axes[1].plot(
                    RESOLUTIONS,
                    [
                        float(row["perturbation_ke_deficit"])
                        for row in selected
                    ],
                    color=color,
                    linestyle=line_style,
                    marker=marker,
                    markersize=4.5,
                    markerfacecolor=marker_face,
                    markeredgecolor=color,
                    markeredgewidth=0.8,
                )

        for index, axis in enumerate(axes):
            axis.set_xscale("log", base=2)
            axis.set_yscale("log")
            axis.set_xticks(RESOLUTIONS, [str(value) for value in RESOLUTIONS])
            axis.set_xlabel(r"Base-grid resolution, $N$")
            axis.grid(axis="y", which="major", color="0.85", linestyle=":", linewidth=0.7)
            axis.tick_params(direction="in", which="both", top=True, right=True)
            axis.text(
                -0.11,
                1.04,
                f"({chr(ord('a') + index)})",
                transform=axis.transAxes,
                ha="left",
                va="bottom",
            )
        axes[0].set_ylabel(r"Density error, $L_2(\rho-\rho_{\mathrm{exact}})$")
        axes[1].set_ylabel(
            r"$|1-K^\prime/K^\prime_{\mathrm{exact}}|$"
        )

        scheme_handles = [
            Line2D(
                [0],
                [0],
                color=SCHEME_STYLES[scheme][0],
                marker=SCHEME_STYLES[scheme][1],
                linestyle="-",
                markersize=4.5,
                label=SCHEME_STYLES[scheme][2],
            )
            for scheme in SCHEMES
        ]
        suite_handles = [
            Line2D(
                [0],
                [0],
                color="0.2",
                linestyle=suite_styles[suite][0],
                marker="o",
                markerfacecolor="0.2" if suite_styles[suite][2] else "white",
                markersize=4.5,
                label=suite_styles[suite][1],
            )
            for suite in SUITES
        ]
        figure.legend(
            handles=scheme_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.02),
            ncol=4,
            frameon=False,
            columnspacing=1.2,
            handlelength=2.0,
        )
        axes[1].legend(
            handles=suite_handles,
            loc="lower left",
            frameon=False,
            handlelength=2.2,
        )
        figure.subplots_adjust(
            left=0.10, right=0.985, bottom=0.17, top=0.83, wspace=0.34
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=300, bbox_inches="tight")
        plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        type=Path,
        help="matrix root produced by run_method_matrix.sh",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        help="combined CSV; default: ROOT/method_matrix_summary.csv",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        help="comparison figure; default: ROOT/method_matrix_accuracy.png",
    )
    args = parser.parse_args()

    output_csv = args.csv or args.root / "method_matrix_summary.csv"
    output_figure = args.figure or args.root / "method_matrix_accuracy.png"
    header, rows = load_matrix(args.root)
    write_summary(output_csv, header, rows)
    plot_summary(output_figure, rows)

    for row in rows:
        rate = float(row["density_L2_observed_rate"])
        rate_text = "   -" if math.isnan(rate) else f"{rate:4.2f}"
        print(
            f"{row['scheme']:20s} {row['suite']:7s} "
            f"N={int(row['base_n']):3d}  "
            f"rho_L2={float(row['density_L2']):.6e}  "
            f"rate={rate_text}  "
            f"Kloss={float(row['perturbation_ke_deficit']):.6e}"
        )
    print(f"wrote {output_csv}")
    print(f"wrote {output_figure}")


if __name__ == "__main__":
    main()
