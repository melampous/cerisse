#!/usr/bin/env python3
"""Resolve competing axial compression peaks in the NPR15 center defect."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_npr15_defect_ab import load_variant


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--z-min", type=float, default=2.30)
    parser.add_argument("--z-max", type=float, default=2.55)
    parser.add_argument("--ring-count", type=int, default=16)
    return parser.parse_args()


def local_maxima(signal: np.ndarray, mask: np.ndarray) -> np.ndarray:
    interior = np.flatnonzero(
        mask[1:-1]
        & (signal[1:-1] > 0.0)
        & (signal[1:-1] >= signal[:-2])
        & (signal[1:-1] > signal[2:])
    )
    return interior + 1


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    variant = load_variant("baseline", args.plotfile)
    z_d = variant["z_d"]
    r_d = variant["r_d"]
    density_gradient = variant["rho_gradient"]
    pressure = variant["arrays"]["pressure"]
    density = variant["arrays"]["Density"]
    mach = variant["mach"]
    mask = (z_d >= args.z_min) & (z_d <= args.z_max)
    ring_count = min(args.ring_count, density.shape[0])

    rows = []
    peak_sets: dict[str, list[dict]] = {}
    for radial_index in range(ring_count):
        maxima = local_maxima(density_gradient[radial_index], mask)
        ordered = maxima[
            np.argsort(density_gradient[radial_index, maxima])[::-1]
        ]
        row_max = max(
            float(density_gradient[radial_index, ordered[0]])
            if ordered.size else 0.0,
            np.finfo(float).tiny,
        )
        peaks = []
        for rank, axial_index in enumerate(ordered[:8], start=1):
            lower = max(0, axial_index - 2)
            upper = min(z_d.size, axial_index + 3)
            before = slice(lower, axial_index)
            after = slice(axial_index + 1, upper)
            peak = {
                "radial_index": radial_index,
                "r_over_D": float(r_d[radial_index]),
                "rank": rank,
                "z_over_D": float(z_d[axial_index]),
                "density_gradient": float(
                    density_gradient[radial_index, axial_index]
                ),
                "relative_density_gradient": float(
                    density_gradient[radial_index, axial_index] / row_max
                ),
                "density_before": float(
                    np.mean(density[radial_index, before])
                ),
                "density_after": float(
                    np.mean(density[radial_index, after])
                ),
                "pressure_before_pa": float(
                    np.mean(pressure[radial_index, before])
                ),
                "pressure_after_pa": float(
                    np.mean(pressure[radial_index, after])
                ),
                "mach_before": float(np.mean(mach[radial_index, before])),
                "mach_after": float(np.mean(mach[radial_index, after])),
            }
            peaks.append(peak)
            rows.append(peak)
        peak_sets[str(radial_index)] = peaks

    with (args.output_dir / "npr15_peak_identity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "plotfile": str(args.plotfile),
        "time_s": variant["metrics"]["time_s"],
        "z_window_over_D": [args.z_min, args.z_max],
        "ring_count": ring_count,
        "peaks": peak_sets,
    }
    with (args.output_dir / "npr15_peak_identity.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 5.0), constrained_layout=True)
    selected = min(8, ring_count)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, selected))
    for radial_index, color in enumerate(colors):
        gradient = density_gradient[radial_index]
        scale = max(np.max(gradient[mask]), np.finfo(float).tiny)
        axes[0].plot(
            z_d[mask],
            gradient[mask] / scale + radial_index * 1.15,
            color=color,
            label=f"i={radial_index}, r/D={r_d[radial_index]:.4f}",
        )
        axes[1].plot(
            z_d[mask], mach[radial_index, mask], color=color,
            label=f"i={radial_index}",
        )
    axes[0].set(
        xlabel="z/D",
        ylabel="normalized density gradient + ring offset",
        title="Competing compression peaks",
    )
    axes[0].legend(fontsize=7)
    axes[1].set(
        xlabel="z/D", ylabel="Mach", title="Mach profiles"
    )
    axes[1].axhline(1.0, color="gray", lw=0.8)
    axes[1].legend(fontsize=7)

    radial_mask = r_d <= 0.08
    gradient_roi = np.maximum(
        density_gradient[np.ix_(radial_mask, mask)], 0.0
    )
    row_scale = np.maximum(
        np.max(gradient_roi, axis=1, keepdims=True),
        np.finfo(float).tiny,
    )
    image = axes[2].pcolormesh(
        z_d[mask],
        r_d[radial_mask],
        gradient_roi / row_scale,
        shading="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
    )
    figure.colorbar(image, ax=axes[2], label="row-normalized density gradient")
    axes[2].set(
        xlabel="z/D", ylabel="r/D", title="Near-axis compression ridges"
    )
    for axis in axes:
        axis.grid(alpha=0.2)
    figure.savefig(args.output_dir / "npr15_peak_identity.png", dpi=220)
    plt.close(figure)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
