#!/usr/bin/env python3
"""Compare frozen-checkpoint NPR=15 Mach-disk root-cause variants."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_VARIANTS = (
    "baseline",
    "viscous_annular",
    "transverse",
    "transverse_viscous",
    "directed",
)

LABELS = {
    "baseline": "baseline",
    "viscous_annular": "viscous annular",
    "transverse": "transverse LLF",
    "transverse_viscous": "transverse LLF + viscous annular",
    "directed": "normal + transverse LLF",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--target-time-ms", type=float, default=2.25)
    parser.add_argument("--variants", nargs="+", default=DEFAULT_VARIANTS)
    return parser.parse_args()


def read_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty metrics file: {path}")
    result: dict[str, np.ndarray] = {}
    for key in rows[0]:
        if key == "plotfile":
            result[key] = np.asarray([row[key] for row in rows], dtype=object)
        else:
            result[key] = np.asarray([float(row[key]) for row in rows])
    return result


def summarize(fields: dict[str, np.ndarray]) -> dict[str, float]:
    indent = np.abs(fields["density_center_minus_annulus_D"])
    return {
        "start_time_ms": float(fields["time_ms"][0]),
        "end_time_ms": float(fields["time_ms"][-1]),
        "frame_count": int(fields["time_ms"].size),
        "max_shape_range_D": float(np.nanmax(fields["density_shape_range_D"])),
        "median_shape_range_D": float(np.nanmedian(fields["density_shape_range_D"])),
        "max_abs_center_minus_annulus_D": float(np.nanmax(indent)),
        "median_shock_fwhm_D": float(
            np.nanmedian(fields["density_gradient_fwhm_D"])
        ),
        "max_disk_band_abs_ur_m_per_s": float(
            np.nanmax(fields["disk_band_radial_velocity_absmax_m_per_s"])
        ),
        "max_disk_band_abs_vorticity_per_s": float(
            np.nanmax(fields["disk_band_vorticity_absmax_per_s"])
        ),
        "mean_axis_front_x_over_D": float(
            np.nanmean(fields["density_center_x_over_D"])
        ),
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = {}
    summaries = {}
    compact = {}
    for variant in args.variants:
        root = args.audit_root / variant
        data[variant] = read_csv(root / "frame_metrics.csv")
        compact[variant] = np.load(root / "compact_fields.npz")
        summaries[variant] = summarize(data[variant])

    baseline = summaries[args.variants[0]]
    ratio_keys = (
        "max_shape_range_D",
        "max_abs_center_minus_annulus_D",
        "median_shock_fwhm_D",
        "max_disk_band_abs_ur_m_per_s",
        "max_disk_band_abs_vorticity_per_s",
    )
    for variant in args.variants:
        summaries[variant]["relative_to_baseline"] = {
            key: (
                summaries[variant][key] / baseline[key]
                if baseline[key] != 0.0
                else None
            )
            for key in ratio_keys
        }

    with (args.output_dir / "matrix_summary.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(summaries, stream, indent=2)
        stream.write("\n")

    scalar_keys = [
        key for key, value in baseline.items() if not isinstance(value, dict)
    ]
    with (args.output_dir / "matrix_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["variant", *scalar_keys])
        for variant in args.variants:
            writer.writerow(
                [variant, *(summaries[variant][key] for key in scalar_keys)]
            )

    figure, axes = plt.subplots(
        3, 2, figsize=(12.0, 12.0), constrained_layout=True
    )
    colors = plt.get_cmap("tab10").colors
    for index, variant in enumerate(args.variants):
        fields = data[variant]
        archive = compact[variant]
        color = colors[index % len(colors)]
        label = LABELS.get(variant, variant)
        time_index = int(
            np.argmin(np.abs(fields["time_ms"] - args.target_time_ms))
        )
        front = archive["density_fronts"][time_index]
        radial = archive["r_over_D"][: front.size]

        axes[0, 0].plot(front, radial, color=color, label=label)
        axes[0, 1].plot(
            fields["time_ms"], fields["density_shape_range_D"],
            "o-", ms=3, color=color, label=label,
        )
        axes[1, 0].plot(
            fields["time_ms"],
            fields["density_center_minus_annulus_D"],
            "o-", ms=3, color=color, label=label,
        )
        axes[1, 1].plot(
            fields["time_ms"], fields["density_gradient_fwhm_D"],
            "o-", ms=3, color=color, label=label,
        )
        axes[2, 0].plot(
            fields["time_ms"],
            fields["disk_band_radial_velocity_absmax_m_per_s"],
            "o-", ms=3, color=color, label=label,
        )
        axes[2, 1].semilogy(
            fields["time_ms"],
            fields["disk_band_vorticity_absmax_per_s"],
            "o-", ms=3, color=color, label=label,
        )

    axes[0, 0].set(
        xlabel="front x/D", ylabel="r/D",
        title=f"Instantaneous density front near {args.target_time_ms:.2f} ms",
    )
    axes[0, 1].set(
        xlabel="time [ms]", ylabel="radial front range [D]",
        title="Mach-disk transverse deformation",
    )
    axes[1, 0].axhline(0.0, color="black", lw=0.8)
    axes[1, 0].set(
        xlabel="time [ms]", ylabel="center - annulus [D]",
        title="Axis-local displacement",
    )
    axes[1, 1].set(
        xlabel="time [ms]", ylabel="axis shock FWHM [D]",
        title="Normal shock thickness",
    )
    axes[2, 0].set(
        xlabel="time [ms]", ylabel="max disk-band |u_r| [m/s]",
        title="Radial-velocity response",
    )
    axes[2, 1].set(
        xlabel="time [ms]", ylabel="max disk-band |vorticity| [1/s]",
        title="Vorticity response",
    )
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    figure.suptitle("NPR=15 frozen-checkpoint root-cause matrix", fontsize=15)
    figure.savefig(args.output_dir / "npr15_root_matrix.png", dpi=220)
    plt.close(figure)

    print(json.dumps(summaries, indent=2), flush=True)


if __name__ == "__main__":
    main()
