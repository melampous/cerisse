#!/usr/bin/env python3
"""Plot first-ring profiles and Mach-disk loci for a real A/B pair."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

import analyze_real_ab as analysis


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("control", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--level", type=int, default=4)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    cases = [
        (analysis.load_near_axis(args.control, args.level), "Control", "#1f77b4"),
        (
            analysis.load_near_axis(args.candidate, args.level),
            "Candidate v2",
            "#d62728",
        ),
    ]
    if args.baseline:
        cases.insert(
            0,
            (
                analysis.load_near_axis(args.baseline, args.level),
                "Baseline",
                "#555555",
            ),
        )
    for data, _, _ in cases[1:]:
        analysis.require_matching_grid(cases[0][0], data, ("level", "dr", "dz"))

    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    z_window = (cases[0][0]["z"] >= -0.125) & (cases[0][0]["z"] <= -0.108)
    for data, label, color in cases:
        z_mm = data["z"][z_window] * 1.0e3
        shock_axis_mm = data["shock"][0] * 1.0e3
        axes[0, 0].semilogy(
            z_mm, data["pressure"][0, z_window], color=color, label=label
        )
        axes[0, 1].plot(
            z_mm, data["radial_velocity"][0, z_window], color=color, label=label
        )
        axes[1, 0].plot(
            z_mm, data["vorticity"][0, z_window], color=color, label=label
        )
        for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
            axis.axvline(shock_axis_mm, color=color, linestyle=":", linewidth=1)

        radial_window = data["r"] <= 0.012
        axes[1, 1].plot(
            data["r"][radial_window] * 1.0e3,
            data["shock"][radial_window] * 1.0e3,
            color=color,
            label=label,
        )

    axes[0, 0].set_title("First-ring pressure")
    axes[0, 0].set_ylabel("pressure [Pa]")
    axes[0, 1].set_title("First-ring radial velocity")
    axes[0, 1].set_ylabel("$u_r$ [m/s]")
    axes[1, 0].set_title("First-ring azimuthal vorticity")
    axes[1, 0].set_ylabel("$\\omega_\\theta$ [1/s]")
    axes[1, 1].set_title("Mach-disk pressure-gradient locus")
    axes[1, 1].set_xlabel("r [mm]")
    axes[1, 1].set_ylabel("z [mm]")
    for axis in axes.flat[:3]:
        axis.set_xlabel("z [mm]")
    for axis in axes.flat:
        axis.grid(alpha=0.22)
        axis.legend(frameon=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)
    print(args.output)


if __name__ == "__main__":
    main()
