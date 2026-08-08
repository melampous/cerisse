#!/usr/bin/env python3
"""Analyse the same-source four-scheme Shu--Osher matrix.

Each scheme's N=800 profile is used only as a numerical resolution reference
for that scheme.  No formal convergence order is inferred.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analyze import advance_time, read_profile, require_matching_uniform_domains


SCHEMES = (
    "llf-wenoz5",
    "llf-teno5",
    "afd-hllc-wenoz5",
    "skew",
)
CELLS = (200, 400, 800)
LABELS = {
    "llf-wenoz5": "LLF WENO-Z5",
    "llf-teno5": "LLF TENO5",
    "afd-hllc-wenoz5": "AFD HLLC WENO-Z5",
    "skew": "Skew-JST",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--window-lo", type=float, default=-1.5)
    parser.add_argument("--window-hi", type=float, default=1.5)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--png", type=Path)
    parser.add_argument("--pdf", type=Path)
    parser.add_argument("--metadata", type=Path)
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
            profiles[(scheme, cells)] = (
                plotfile,
                dataset,
                x,
                density,
                case,
            )

    rows = []
    for scheme in SCHEMES:
        ref_plotfile, ref_ds, x_ref, rho_ref, _ = profiles[(scheme, 800)]
        ref_time = float(ref_ds.current_time)
        ref_window = (x_ref >= args.window_lo) & (x_ref <= args.window_hi)
        ref_tv = float(np.sum(np.abs(np.diff(rho_ref[ref_window]))))
        ref_range = float(
            np.max(rho_ref[ref_window]) - np.min(rho_ref[ref_window])
        )

        for cells in CELLS:
            plotfile, dataset, x, density, case = profiles[(scheme, cells)]
            require_matching_uniform_domains(dataset, ref_ds, plotfile)
            time = float(dataset.current_time)
            if abs(time - ref_time) > 1.0e-10 * max(1.0, abs(ref_time)):
                raise ValueError(
                    f"time mismatch: {plotfile} has t={time}, whereas "
                    f"{ref_plotfile} has t={ref_time}"
                )
            window = (x >= args.window_lo) & (x <= args.window_hi)
            x_window = x[window]
            rho_window = density[window]
            rho_ref_on_grid = np.interp(x_window, x_ref, rho_ref)
            difference = rho_window - rho_ref_on_grid
            tv = float(np.sum(np.abs(np.diff(rho_window))))
            density_range = float(np.max(rho_window) - np.min(rho_window))

            rows.append(
                {
                    "scheme": scheme,
                    "scheme_label": LABELS[scheme],
                    "cells": cells,
                    "time": time,
                    "advance_seconds": advance_time(case),
                    "postshock_window_lo": args.window_lo,
                    "postshock_window_hi": args.window_hi,
                    "resolution_reference": f"{LABELS[scheme]} N=800",
                    "reference_is_exact": False,
                    "formal_order_inferred": False,
                    "postshock_density_L1_to_same_scheme_N800": float(
                        np.mean(np.abs(difference))
                    ),
                    "postshock_density_L2_to_same_scheme_N800": float(
                        np.sqrt(np.mean(difference * difference))
                    ),
                    "postshock_density_total_variation": tv,
                    "postshock_total_variation_ratio_to_same_scheme_N800": (
                        tv / ref_tv
                    ),
                    "postshock_density_peak_to_trough_range": density_range,
                    "postshock_range_ratio_to_same_scheme_N800": (
                        density_range / ref_range
                    ),
                    "postshock_density_min": float(np.min(rho_window)),
                    "postshock_density_max": float(np.max(rho_window)),
                    "full_domain_density_min": float(np.min(density)),
                    "full_domain_density_max": float(np.max(density)),
                    "plotfile": str(plotfile),
                }
            )

    csv_path = args.csv or args.root / "same_source_four_scheme_metrics.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['scheme']:18s} N={row['cells']:4d} "
            f"L1(N800)={row['postshock_density_L1_to_same_scheme_N800']:.6e} "
            f"TV={row['postshock_density_total_variation']:.6e} "
            f"range={row['postshock_density_peak_to_trough_range']:.6e} "
            f"advance={row['advance_seconds']:.6g} s"
        )
    print(f"wrote {csv_path}")

    metadata_path = args.metadata or args.root / "analysis_metadata.json"
    metadata = {
        "case": "canonical Shu-Osher shock-entropy-wave interaction",
        "schemes": [LABELS[scheme] for scheme in SCHEMES],
        "grids": list(CELLS),
        "postshock_window": [args.window_lo, args.window_hi],
        "resolution_reference": (
            "Each scheme's N=800 profile is used only for its own "
            "resolution-difference metrics."
        ),
        "skew_N6400_usage": (
            "The separate Skew-JST N=6400 profile is not used as a reference "
            "for LLF or AFD-HLLC."
        ),
        "reference_is_exact": False,
        "formal_order_inferred": False,
        "indicator_note": (
            "Post-shock total variation and peak-to-trough range are simple "
            "small-scale retention indicators. They are not accuracy orders."
        ),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"wrote {metadata_path}")

    if args.png or args.pdf:
        import matplotlib.pyplot as plt

        styles = {
            "llf-wenoz5": ("#0072B2", "-", 1.05),
            "llf-teno5": ("#009E73", "--", 1.05),
            "afd-hllc-wenoz5": ("#D55E00", "-.", 1.05),
            "skew": ("#000000", ":", 1.2),
        }
        with plt.rc_context(
            {
                "font.size": 9,
                "axes.labelsize": 10,
                "axes.titlesize": 9,
                "legend.fontsize": 8,
                "xtick.labelsize": 8.5,
                "ytick.labelsize": 8.5,
                "lines.solid_capstyle": "round",
            }
        ):
            fig, axes = plt.subplots(
                len(CELLS),
                1,
                figsize=(7.0, 7.1),
                sharex=True,
                constrained_layout=True,
            )
            for axis, cells in zip(axes, CELLS):
                for scheme in SCHEMES:
                    _, _, x, density, _ = profiles[(scheme, cells)]
                    color, line_style, width = styles[scheme]
                    axis.plot(
                        x,
                        density,
                        color=color,
                        linestyle=line_style,
                        linewidth=width,
                        label=LABELS[scheme],
                    )
                axis.text(
                    0.018,
                    0.93,
                    f"N = {cells}",
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                )
                axis.set_ylabel(r"$\rho$")
                axis.set_xlim(-1.5, 2.5)
                axis.grid(False)
            axes[0].legend(
                frameon=False,
                ncol=2,
                loc="lower left",
            )
            axes[-1].set_xlabel(r"$x$")
            if args.png:
                args.png.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(args.png, dpi=400)
                print(f"wrote {args.png}")
            if args.pdf:
                args.pdf.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(args.pdf)
                print(f"wrote {args.pdf}")


if __name__ == "__main__":
    try:
        import yt
    except ImportError as exc:
        raise SystemExit(
            "analyze_same_source_four_scheme.py requires yt"
        ) from exc
    yt.set_log_level(40)
    main()
