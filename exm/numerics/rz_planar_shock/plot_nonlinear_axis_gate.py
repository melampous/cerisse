#!/usr/bin/env python3
"""Plot the short R-Z nonlinear-axis matrix from case-local JSON metrics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "llf_wenoz5": "#1f77b4",
    "afd_hllc_wenoz5": "#ff7f0e",
    "skew4_central": "#2ca02c",
    "skew4_jst": "#d62728",
}


def split_case(name: str) -> tuple[str, int]:
    scheme, fallback = name.rsplit("_fallback", 1)
    return scheme, int(fallback)


def signed_arrays(records: list[dict], key: str) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.asarray([row[key]["min"] for row in records]),
        np.asarray([row[key]["max"] for row in records]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix_directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.matrix_directory / "nonlinear_axis_history.png"

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), sharex=True)
    ax_spread, ax_ur, ax_even, ax_odd = axes.flat
    for metrics_path in sorted(args.matrix_directory.glob("*/metrics.json")):
        case = metrics_path.parent.name
        scheme, fallback = split_case(case)
        data = json.loads(metrics_path.read_text())
        rows = data["plotfiles"]
        time = np.asarray([row["time"] for row in rows])
        color = COLORS[scheme]
        linestyle = "-" if fallback else "--"
        label = f"{scheme}, fallback {'on' if fallback else 'off'}"

        ax_spread.plot(
            time,
            [row["max_radial_pressure_spread_rel_local"] for row in rows],
            color=color,
            linestyle=linestyle,
            marker="o",
            markersize=2.5,
            label=label,
        )
        ax_ur.plot(
            time,
            [row["max_axis_abs_radial_velocity_over_sound"] for row in rows],
            color=color,
            linestyle=linestyle,
            marker="o",
            markersize=2.5,
        )
        even_min, even_max = signed_arrays(
            rows, "signed_axis_to_second_ring_pressure_rel"
        )
        odd_min, odd_max = signed_arrays(
            rows, "signed_odd_axis_extrapolation_ur_over_sound"
        )
        ax_even.plot(time, even_min, color=color, linestyle=linestyle)
        ax_even.plot(time, even_max, color=color, linestyle=linestyle)
        ax_odd.plot(time, odd_min, color=color, linestyle=linestyle)
        ax_odd.plot(time, odd_max, color=color, linestyle=linestyle)

    ax_spread.set_ylabel(r"$\max_z\,\Delta_r p/\max_r|p|$")
    ax_ur.set_ylabel(r"$\max_z |u_r|/a$ on ring 0")
    ax_even.set_ylabel(r"signed $(p_0-p_1)/\max_r|p|$")
    ax_odd.set_ylabel(r"signed $(3u_{r,0}-u_{r,1})/(2a_0)$")
    ax_even.set_xlabel("time")
    ax_odd.set_xlabel("time")
    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.25)
    ax_spread.set_yscale("symlog", linthresh=1.0e-14)
    ax_ur.set_yscale("symlog", linthresh=1.0e-14)
    ax_even.set_yscale("symlog", linthresh=1.0e-14)
    ax_odd.set_yscale("symlog", linthresh=1.0e-14)
    ax_spread.legend(fontsize=7, ncol=2)
    fig.suptitle(
        "M=10 stationary planar shock in R-Z, 16x64, CFL=0.25\n"
        "ring 0 is centred at r=dr/2 (not a sample at r=0)"
    )
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
