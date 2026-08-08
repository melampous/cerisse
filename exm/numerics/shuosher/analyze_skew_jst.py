#!/usr/bin/env python3
"""Quantify the Skew--JST Shu--Osher resolution study.

The selected reference is a numerical comparison profile.  It is not an exact
solution and is not used to infer a formal convergence order.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analyze import advance_time, read_profile, require_matching_uniform_domains


CELLS = (200, 400, 800)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument(
        "--reference",
        type=Path,
        help="reference directory or plotfile; default is ROOT/N800",
    )
    parser.add_argument(
        "--reference-label",
        default="Skew-JST N=800 numerical resolution reference",
    )
    parser.add_argument("--window-lo", type=float, default=-1.5)
    parser.add_argument("--window-hi", type=float, default=1.5)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()

    reference_path = args.reference or args.root / "N800"
    ref_plotfile, ref_ds, x_ref, rho_ref = read_profile(reference_path)
    if int(ref_ds.index.max_level) != 0:
        raise ValueError(f"reference is not uniform level 0: {ref_plotfile}")
    ref_time = float(ref_ds.current_time)
    ref_window = (x_ref >= args.window_lo) & (x_ref <= args.window_hi)
    if np.count_nonzero(ref_window) < 2:
        raise ValueError("the post-shock reference window contains fewer than two cells")
    ref_min = float(np.min(rho_ref[ref_window]))
    ref_max = float(np.max(rho_ref[ref_window]))
    ref_tv = float(np.sum(np.abs(np.diff(rho_ref[ref_window]))))
    ref_range = ref_max - ref_min

    profiles = []
    rows = []
    for cells in CELLS:
        case = args.root / f"N{cells}"
        plotfile, ds, x, rho = read_profile(case)
        if x.size != cells:
            raise ValueError(f"{plotfile} contains {x.size} cells; expected {cells}")
        require_matching_uniform_domains(ds, ref_ds, plotfile)
        time = float(ds.current_time)
        if abs(time - ref_time) > 1.0e-10 * max(1.0, abs(ref_time)):
            raise ValueError(
                f"time mismatch: {plotfile} has t={time}, reference has t={ref_time}"
            )

        window = (x >= args.window_lo) & (x <= args.window_hi)
        if np.count_nonzero(window) < 2:
            raise ValueError(f"the post-shock window contains fewer than two N={cells} cells")
        x_window = x[window]
        rho_window = rho[window]
        rho_reference = np.interp(x_window, x_ref, rho_ref)
        difference = rho_window - rho_reference
        rows.append(
            {
                "scheme": "skew-jst",
                "cells": cells,
                "time": time,
                "advance_seconds": advance_time(case),
                "postshock_window_lo": args.window_lo,
                "postshock_window_hi": args.window_hi,
                "reference_label": args.reference_label,
                "reference_plotfile": str(ref_plotfile),
                "postshock_density_L1_to_reference": float(
                    np.mean(np.abs(difference))
                ),
                "postshock_density_L2_to_reference": float(
                    np.sqrt(np.mean(difference * difference))
                ),
                "postshock_density_Linf_to_reference": float(
                    np.max(np.abs(difference))
                ),
                "postshock_density_total_variation": float(
                    np.sum(np.abs(np.diff(rho_window)))
                ),
                "postshock_density_total_variation_ratio_to_reference": float(
                    np.sum(np.abs(np.diff(rho_window))) / ref_tv
                ),
                "postshock_density_min": float(np.min(rho_window)),
                "postshock_density_max": float(np.max(rho_window)),
                "postshock_density_range_ratio_to_reference": float(
                    (float(np.max(rho_window)) - float(np.min(rho_window)))
                    / ref_range
                ),
                "postshock_peak_error_to_reference": float(
                    float(np.max(rho_window)) - ref_max
                ),
                "postshock_trough_error_to_reference": float(
                    float(np.min(rho_window)) - ref_min
                ),
                "postshock_overshoot_above_reference_max": float(
                    max(0.0, float(np.max(rho_window)) - ref_max)
                ),
                "postshock_undershoot_below_reference_min": float(
                    max(0.0, ref_min - float(np.min(rho_window)))
                ),
                "full_domain_density_min": float(np.min(rho)),
                "full_domain_density_max": float(np.max(rho)),
                "plotfile": str(plotfile),
            }
        )
        profiles.append((cells, x, rho))

    csv_path = args.csv or args.root / "skew_jst_metrics.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"N={row['cells']:4d} "
            f"L1={row['postshock_density_L1_to_reference']:.6e} "
            f"L2={row['postshock_density_L2_to_reference']:.6e} "
            f"TV={row['postshock_density_total_variation']:.6e} "
            f"TV/ref={row['postshock_density_total_variation_ratio_to_reference']:.4f} "
            f"range/ref={row['postshock_density_range_ratio_to_reference']:.4f} "
            f"over={row['postshock_overshoot_above_reference_max']:.6e} "
            f"under={row['postshock_undershoot_below_reference_min']:.6e}"
        )
    print(f"reference: {args.reference_label}")
    print(f"reference plotfile: {ref_plotfile}")
    print(f"wrote {csv_path}")

    metadata_path = args.metadata or args.root / "analysis_metadata.json"
    metadata = {
        "case": "canonical Shu-Osher shock-entropy-wave interaction",
        "scheme": "order-4 conservative skew-symmetric central flux with JST dissipation",
        "reference_label": args.reference_label,
        "reference_plotfile": str(ref_plotfile),
        "reference_is_exact": False,
        "formal_order_inferred": False,
        "postshock_window": [args.window_lo, args.window_hi],
        "metrics_note": (
            "L1 and L2 are discrete density differences after interpolation "
            "of the stated numerical reference. Overshoot and undershoot use "
            "the extrema of that reference within the same window. The total "
            "variation ratio and peak-to-trough range ratio are simple "
            "small-scale retention indicators, not formal accuracy orders."
        ),
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"wrote {metadata_path}")

    if args.figure:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.0, 3.8))
        ax.plot(
            x_ref,
            rho_ref,
            color="black",
            linewidth=1.0,
            label=args.reference_label,
            zorder=4,
        )
        styles = {
            200: ("#0072B2", "--"),
            400: ("#D55E00", "-."),
            800: ("#009E73", ":"),
        }
        for cells, x, rho in profiles:
            color, line_style = styles[cells]
            ax.plot(
                x,
                rho,
                color=color,
                linestyle=line_style,
                linewidth=0.9,
                label=f"Skew-JST, N = {cells}",
            )
        ax.axvspan(
            args.window_lo,
            args.window_hi,
            color="0.92",
            linewidth=0.0,
            zorder=0,
            label="metric window",
        )
        ax.set(xlabel=r"$x$", ylabel=r"$\rho$", xlim=(-2.0, 2.8))
        ax.legend(frameon=False, fontsize=7.5, ncol=2)
        fig.tight_layout()
        args.figure.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.figure, dpi=300)
        print(f"wrote {args.figure}")


if __name__ == "__main__":
    try:
        import yt
    except ImportError as exc:
        raise SystemExit("analyze_skew_jst.py requires yt") from exc
    yt.set_log_level(40)
    main()
