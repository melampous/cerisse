#!/usr/bin/env python3
"""Summarise the four-scheme Shu--Osher resolution matrix.

The N=800 profile for each scheme is used only to quantify the change with
resolution.  These differences are not formal discretisation errors and N=800
is not presented as an exact or converged solution.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from analyze import (
    advance_time,
    read_profile,
    require_matching_uniform_domains,
)


SCHEMES = (
    "llf-wenoz5",
    "llf-teno5",
    "afd-hllc-wenoz5",
    "afd-hllc-teno5",
)
CELLS = (200, 400, 800)


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
        default=Path("method_matrix_metrics.csv"),
    )
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()

    profiles: dict[tuple[str, int], tuple] = {}
    for scheme in SCHEMES:
        for cells in CELLS:
            case = args.root / scheme / f"N{cells}"
            plotfile, dataset, x, density = read_profile(case)
            if x.size != cells:
                raise ValueError(
                    f"{plotfile} contains {x.size} cells; expected {cells}"
                )
            profiles[(scheme, cells)] = (plotfile, dataset, x, density, case)

    rows = []
    for scheme in SCHEMES:
        (
            reference_plotfile,
            reference_dataset,
            x_reference,
            density_reference,
            _,
        ) = profiles[(scheme, 800)]
        reference_time = float(reference_dataset.current_time)

        for cells in CELLS:
            plotfile, dataset, x, density, case = profiles[(scheme, cells)]
            require_matching_uniform_domains(
                dataset, reference_dataset, plotfile
            )
            time = float(dataset.current_time)
            if abs(time - reference_time) > 1.0e-10 * max(
                1.0, abs(reference_time)
            ):
                raise ValueError(
                    f"time mismatch: {plotfile} has t={time}, whereas "
                    f"{reference_plotfile} has t={reference_time}"
                )

            if cells == 800:
                difference = np.zeros_like(density)
                is_resolution_reference = 1
            else:
                difference = density - np.interp(
                    x, x_reference, density_reference
                )
                is_resolution_reference = 0

            rows.append(
                {
                    "scheme": scheme,
                    "cells": cells,
                    "time": time,
                    "advance_seconds": advance_time(case),
                    "density_min": float(np.min(density)),
                    "density_max": float(np.max(density)),
                    "density_total_variation": float(
                        np.sum(np.abs(np.diff(density)))
                    ),
                    "density_L1_difference_to_same_scheme_N800": float(
                        np.mean(np.abs(difference))
                    ),
                    "density_L2_difference_to_same_scheme_N800": float(
                        np.sqrt(np.mean(difference * difference))
                    ),
                    "density_Linf_difference_to_same_scheme_N800": float(
                        np.max(np.abs(difference))
                    ),
                    "same_scheme_resolution_reference": (
                        is_resolution_reference
                    ),
                    "plotfile": str(plotfile),
                }
            )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['scheme']:20s} N={row['cells']:4d} "
            f"TV={row['density_total_variation']:.6e} "
            f"L1(N800)={row['density_L1_difference_to_same_scheme_N800']:.6e} "
            f"advance={row['advance_seconds']:.6g} s"
        )
    print(f"wrote {args.csv}")

    if args.figure:
        import matplotlib.pyplot as plt

        styles = {
            "llf-wenoz5": ("#0072B2", "-", "LLF WENO-Z5"),
            "llf-teno5": ("#009E73", "--", "LLF TENO5"),
            "afd-hllc-wenoz5": ("#D55E00", "-.", "AFD HLLC WENO-Z5"),
            "afd-hllc-teno5": ("#CC79A7", ":", "AFD HLLC TENO5"),
        }
        fig, axes = plt.subplots(
            len(CELLS),
            1,
            figsize=(7.0, 7.2),
            sharex=True,
            constrained_layout=True,
        )
        for axis, cells in zip(axes, CELLS):
            for scheme in SCHEMES:
                _, _, x, density, _ = profiles[(scheme, cells)]
                color, line_style, label = styles[scheme]
                axis.plot(
                    x,
                    density,
                    color=color,
                    linestyle=line_style,
                    linewidth=1.0,
                    label=label,
                )
            axis.text(
                0.02,
                0.93,
                f"N = {cells}",
                transform=axis.transAxes,
                ha="left",
                va="top",
            )
            axis.set_ylabel(r"$\rho$")
            axis.set_xlim(-1.5, 2.5)
        axes[0].legend(
            frameon=False,
            ncol=2,
            loc="lower left",
            fontsize=8,
        )
        axes[-1].set_xlabel(r"$x$")
        args.figure.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.figure, dpi=300)
        print(f"wrote {args.figure}")


if __name__ == "__main__":
    main()
