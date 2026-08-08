#!/usr/bin/env python3
"""Track the physical NPR15 Mach disk through time and in the mean field."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

from analyze_npr15_defect_ab import load_variant
from audit_npr15_disk import composite_region, positive_gradient
from npr15_shock_tracking import track_primary_mach_disk
from npr15_shock_tracking import subcell_peak


CORE_R_MAX_D = 0.08
DIAMETER = 0.0254


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def time_weights(times: np.ndarray) -> np.ndarray:
    if times.size == 1:
        return np.ones(1)
    weights = np.empty_like(times)
    weights[0] = 0.5 * (times[1] - times[0])
    weights[-1] = 0.5 * (times[-1] - times[-2])
    weights[1:-1] = 0.5 * (times[2:] - times[:-2])
    return weights / weights.sum()


def load_online_mean_track(
    plotfile: Path,
    reference_r_over_d: np.ndarray,
    reference_position: np.ndarray,
) -> dict:
    dataset = yt.load(str(plotfile))
    fields = ["DensityMEAN"]
    arrays, r, z, _, spacing, level = composite_region(
        dataset, fields, DIAMETER, 0.15, 1.8, 3.4
    )
    r_d = r / DIAMETER
    z_d = z / DIAMETER
    radial_count = int(np.count_nonzero(r_d <= CORE_R_MAX_D))
    gradient = positive_gradient(arrays["DensityMEAN"], z)
    local_reference = np.interp(
        r_d[:radial_count], reference_r_over_d, reference_position
    )
    position = np.full(radial_count, np.nan)
    for radial_index, expected in enumerate(local_reference):
        window = np.abs(z_d - expected) <= 0.08
        candidates = np.flatnonzero(window)
        if not candidates.size:
            continue
        axial_index = int(
            candidates[np.argmax(gradient[radial_index, candidates])]
        )
        position[radial_index] = subcell_peak(
            gradient[radial_index], z_d, axial_index
        )
    return {
        "position": position,
        "r_over_D": r_d[:radial_count],
        "spacing_over_D": float(spacing[1] / DIAMETER),
        "level": int(level),
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    yt.funcs.mylog.setLevel(50)

    plotfiles = []
    for segment in ("flush", "stat"):
        plotfiles.extend(sorted((args.case_dir / segment).glob("plt[0-9]*")))
    if not plotfiles:
        raise FileNotFoundError(f"no plotfiles below {args.case_dir}")

    records = []
    tracks: list[tuple[np.ndarray, np.ndarray]] = []
    for plotfile in plotfiles:
        variant = load_variant(plotfile.parent.name, plotfile)
        r_d = variant["r_d"]
        radial_count = int(np.count_nonzero(r_d <= CORE_R_MAX_D))
        primary = track_primary_mach_disk(
            variant["rho_gradient"],
            variant["arrays"]["Density"],
            variant["mach"],
            variant["z_d"],
            radial_count,
        )
        position = primary["position"]
        valid = np.isfinite(position)
        spacing_d = variant["metrics"]["axial_spacing_over_D"]
        record = {
            "segment": plotfile.parent.name,
            "plotfile": plotfile.name,
            "time_ms": variant["metrics"]["time_s"] * 1.0e3,
            "finest_level": variant["metrics"]["finest_level"],
            "spacing_over_D": spacing_d,
            "valid_fraction": float(np.mean(valid)),
            "axis_z_over_D": float(position[0]),
            "core_mean_z_over_D": float(
                np.nanmean(position) if np.any(valid) else np.nan
            ),
            "core_range_cells": float(
                np.ptp(position[valid]) / spacing_d
                if np.any(valid) else np.nan
            ),
            "axis_mach_before": float(primary["mach_before"][0]),
            "axis_mach_after": float(primary["mach_after"][0]),
        }
        records.append(record)
        tracks.append((r_d[:radial_count].copy(), position.copy()))
        print(
            f"{record['segment']}/{record['plotfile']} "
            f"t={record['time_ms']:.6f} ms "
            f"valid={record['valid_fraction']:.3f} "
            f"range={record['core_range_cells']:.3f} cells",
            flush=True,
        )

    usable = np.asarray(
        [
            record["valid_fraction"] == 1.0
            and record["finest_level"] == 4
            for record in records
        ]
    )
    if not np.any(usable):
        raise RuntimeError("no fully tracked L4 Mach-disk frames")
    target_index = int(
        np.argmax(
            [
                radial.size if keep else 0
                for keep, (radial, _) in zip(usable, tracks)
            ]
        )
    )
    target_r = tracks[target_index][0]
    usable_records = [
        record for keep, record in zip(usable, records) if keep
    ]
    usable_fronts = np.stack(
        [
            np.interp(target_r, radial, front)
            for keep, (radial, front) in zip(usable, tracks)
            if keep
        ]
    )
    usable_times = np.asarray(
        [record["time_ms"] for record in usable_records]
    )
    order = np.argsort(usable_times)
    usable_times = usable_times[order]
    usable_fronts = usable_fronts[order]
    usable_records = [usable_records[index] for index in order]

    stat_mask = np.asarray(
        [record["segment"] == "stat" for record in usable_records]
    )
    stat_times = usable_times[stat_mask]
    stat_fronts = usable_fronts[stat_mask]
    weights = time_weights(stat_times)
    mean_locus = np.tensordot(weights, stat_fronts, axes=(0, 0))
    variance_locus = np.tensordot(
        weights, (stat_fronts - mean_locus) ** 2, axes=(0, 0)
    )
    std_locus = np.sqrt(variance_locus)

    final_plot = sorted((args.case_dir / "stat").glob("plt[0-9]*"))[-1]
    online_mean = load_online_mean_track(
        final_plot, target_r, mean_locus
    )
    online_position = online_mean["position"]
    if not np.all(np.isfinite(online_position)):
        raise RuntimeError("online mean Mach disk was not fully tracked")

    with (args.output_dir / "npr15_primary_history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    np.savez_compressed(
        args.output_dir / "npr15_primary_history.npz",
        time_ms=usable_times,
        r_over_D=target_r,
        primary_front_z_over_D=usable_fronts,
        stat_mean_locus_z_over_D=mean_locus,
        stat_std_locus_over_D=std_locus,
        online_mean_r_over_D=online_mean["r_over_D"],
        online_mean_locus_z_over_D=online_position,
    )

    finest_spacing_d = online_mean["spacing_over_D"]
    stat_ranges = np.ptp(stat_fronts, axis=1) / finest_spacing_d
    summary = {
        "case_dir": str(args.case_dir),
        "physical_criterion": (
            "first compression with upstream Mach >= 1.5 and "
            "downstream Mach <= 1.5"
        ),
        "usable_L4_frame_count": int(usable_fronts.shape[0]),
        "stat_frame_count": int(stat_fronts.shape[0]),
        "stat_time_start_ms": float(stat_times[0]),
        "stat_time_end_ms": float(stat_times[-1]),
        "finest_spacing_over_D": finest_spacing_d,
        "instantaneous_stat_core_range_cells": {
            "median": float(np.median(stat_ranges)),
            "maximum": float(np.max(stat_ranges)),
        },
        "time_averaged_instantaneous_locus": {
            "core_range_cells": float(
                np.ptp(mean_locus) / finest_spacing_d
            ),
            "axis_z_over_D": float(mean_locus[0]),
            "outer_core_z_over_D": float(mean_locus[-1]),
            "maximum_temporal_std_over_D": float(np.max(std_locus)),
        },
        "online_mean_field_locus": {
            "source_plotfile": str(final_plot),
            "core_range_cells": float(
                np.ptp(online_position) / finest_spacing_d
            ),
            "axis_z_over_D": float(online_position[0]),
            "outer_core_z_over_D": float(online_position[-1]),
        },
    }
    with (args.output_dir / "npr15_primary_history_summary.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    figure, axes = plt.subplots(2, 2, figsize=(13.0, 8.5), constrained_layout=True)
    axes[0, 0].plot(
        usable_times, usable_fronts[:, 0], label="axis"
    )
    axes[0, 0].plot(
        usable_times, usable_fronts.mean(axis=1), label="core mean"
    )
    axes[0, 0].set(
        xlabel="time [ms]", ylabel="Mach-disk z/D", title="Physical disk breathing"
    )
    axes[0, 0].legend()

    ranges = np.ptp(usable_fronts, axis=1) / finest_spacing_d
    axes[0, 1].plot(usable_times, ranges)
    axes[0, 1].set(
        xlabel="time [ms]",
        ylabel="core radial range [L4 cells]",
        title="Instantaneous disk flatness",
    )

    offsets = (
        usable_fronts - usable_fronts.mean(axis=1, keepdims=True)
    ) / finest_spacing_d
    limit = max(np.percentile(np.abs(offsets), 99.0), 0.1)
    image = axes[1, 0].pcolormesh(
        usable_times,
        target_r,
        offsets.T,
        shading="nearest",
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit,
    )
    figure.colorbar(image, ax=axes[1, 0], label="front offset [L4 cells]")
    axes[1, 0].set(
        xlabel="time [ms]", ylabel="r/D", title="Core shape history"
    )

    axes[1, 1].plot(
        mean_locus, target_r, label="time mean of instantaneous loci"
    )
    axes[1, 1].fill_betweenx(
        target_r,
        mean_locus - std_locus,
        mean_locus + std_locus,
        alpha=0.2,
        label="temporal +/-1 std",
    )
    axes[1, 1].plot(
        online_position,
        online_mean["r_over_D"],
        "--",
        label="track on online mean field",
    )
    axes[1, 1].set(
        xlabel="Mach-disk z/D", ylabel="r/D", title="Mean disk locus"
    )
    axes[1, 1].legend(fontsize=8)

    for axis in axes.flat:
        axis.grid(alpha=0.22)
    figure.savefig(args.output_dir / "npr15_primary_history.png", dpi=220)
    plt.close(figure)

    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
