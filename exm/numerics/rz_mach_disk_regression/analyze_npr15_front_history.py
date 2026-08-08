#!/usr/bin/env python3
"""Track NPR15 Mach-disk radial modes from startup through the stat window."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def contiguous_axis_width(front: np.ndarray, threshold_cells: float,
                          spacing_d: float) -> int:
    difference_cells = np.abs(front - front[0]) / spacing_d
    exceed = np.flatnonzero(difference_cells > threshold_cells)
    return int(exceed[0]) if exceed.size else int(front.size)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    yt.funcs.mylog.setLevel(50)

    plotfiles = []
    for segment in ("flush", "stat"):
        plotfiles.extend(sorted((args.case_dir / segment).glob("plt[0-9]*")))
    if not plotfiles:
        raise FileNotFoundError(f"no plotfiles below {args.case_dir}")

    records: list[dict] = []
    fronts: list[tuple[np.ndarray, np.ndarray]] = []
    for index, plotfile in enumerate(plotfiles):
        variant = load_variant(plotfile.parent.name, plotfile)
        metrics = variant["metrics"]
        front = variant["density_front_subcell"]
        spacing_d = metrics["axial_spacing_over_D"]
        r_d = variant["r_d"][:front.size]

        annulus = (r_d >= 0.3) & (r_d <= 0.5)
        annulus_location = float(np.median(front[annulus]))
        record = {
            "sequence": index,
            "segment": plotfile.parent.name,
            "plotfile": plotfile.name,
            "time_ms": metrics["time_s"] * 1.0e3,
            "finest_level": metrics["finest_level"],
            "spacing_over_D": spacing_d,
            "axis_z_over_D": float(front[0]),
            "annulus_z_over_D": annulus_location,
            "center_minus_annulus_cells": float(
                (front[0] - annulus_location) / spacing_d
            ),
            "full_range_cells": float(np.ptp(front) / spacing_d),
            "first8_range_cells": float(np.ptp(front[:8]) / spacing_d),
            "first16_range_cells": float(np.ptp(front[:16]) / spacing_d),
            "axis_width_1cell_rings": contiguous_axis_width(
                front, 1.0, spacing_d
            ),
            "axis_width_4cell_rings": contiguous_axis_width(
                front, 4.0, spacing_d
            ),
            "axis_width_8cell_rings": contiguous_axis_width(
                front, 8.0, spacing_d
            ),
            "axis_width_4cell_over_D": contiguous_axis_width(
                front, 4.0, spacing_d
            ) * spacing_d,
            "axis_band_abs_ur_max_m_per_s":
                metrics["axis_band_abs_ur_max_m_per_s"],
            "first8_abs_ur_max_m_per_s":
                metrics["first8_band_abs_ur_max_m_per_s"],
            "first8_abs_vorticity_max_per_s":
                metrics["first8_band_abs_vorticity_max_per_s"],
        }
        records.append(record)
        fronts.append((r_d.copy(), front.copy()))
        print(
            f"{record['segment']}/{record['plotfile']} "
            f"t={record['time_ms']:.6f} ms "
            f"axis-annulus={record['center_minus_annulus_cells']:.3f} cells "
            f"range={record['full_range_cells']:.3f} cells",
            flush=True,
        )

    target_index = int(np.argmax([radial.size for radial, _ in fronts]))
    radial_coordinate = fronts[target_index][0]
    fronts_array = np.stack(
        [
            np.interp(radial_coordinate, radial, front)
            for radial, front in fronts
        ]
    )
    times_ms = np.asarray([record["time_ms"] for record in records])
    order = np.argsort(times_ms)
    times_ms = times_ms[order]
    fronts_array = fronts_array[order]
    records = [records[position] for position in order]
    finest_spacing_d = min(
        record["spacing_over_D"] for record in records
    )

    centered = fronts_array - fronts_array.mean(axis=0, keepdims=True)
    _, singular_values, spatial_modes = np.linalg.svd(
        centered, full_matrices=False
    )
    modal_energy = singular_values * singular_values
    modal_fraction = modal_energy / modal_energy.sum()

    with (args.output_dir / "npr15_front_history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    np.savez_compressed(
        args.output_dir / "npr15_front_history.npz",
        time_ms=times_ms,
        r_over_D=radial_coordinate,
        density_front_z_over_D=fronts_array,
        singular_values=singular_values,
        modal_fraction=modal_fraction,
        spatial_modes=spatial_modes,
    )

    signed_indent = np.asarray(
        [record["center_minus_annulus_cells"] for record in records]
    )
    full_range = np.asarray(
        [record["full_range_cells"] for record in records]
    )
    first8_range = np.asarray(
        [record["first8_range_cells"] for record in records]
    )
    axis_width = np.asarray(
        [record["axis_width_4cell_rings"] for record in records]
    )
    ur_axis = np.asarray(
        [record["axis_band_abs_ur_max_m_per_s"] for record in records]
    )
    vorticity = np.asarray(
        [record["first8_abs_vorticity_max_per_s"] for record in records]
    )
    stat_start = next(
        (
            record["time_ms"] for record in records
            if record["segment"] == "stat"
        ),
        float("nan"),
    )

    summary = {
        "case_dir": str(args.case_dir),
        "frame_count": len(records),
        "time_start_ms": float(times_ms[0]),
        "time_end_ms": float(times_ms[-1]),
        "stat_start_ms": stat_start,
        "finest_spacing_over_D": float(finest_spacing_d),
        "signed_center_minus_annulus_cells": {
            "minimum": float(signed_indent.min()),
            "maximum": float(signed_indent.max()),
            "median": float(np.median(signed_indent)),
            "sign_change_count": int(
                np.count_nonzero(signed_indent[1:] * signed_indent[:-1] < 0.0)
            ),
        },
        "front_range_cells": {
            "median": float(np.median(full_range)),
            "maximum": float(full_range.max()),
            "first8_median": float(np.median(first8_range)),
            "first8_maximum": float(first8_range.max()),
        },
        "axis_core_width_within_4cell_difference_rings": {
            "minimum": int(axis_width.min()),
            "maximum": int(axis_width.max()),
            "median": float(np.median(axis_width)),
        },
        "amr": {
            "minimum_finest_level": int(
                min(record["finest_level"] for record in records)
            ),
            "maximum_finest_level": int(
                max(record["finest_level"] for record in records)
            ),
            "first_L4_time_ms": float(
                next(
                    record["time_ms"] for record in records
                    if record["finest_level"] == 4
                )
            ),
        },
        "correlations": {
            "abs_indent_with_axis_abs_ur": float(
                np.corrcoef(np.abs(signed_indent), ur_axis)[0, 1]
            ),
            "abs_indent_with_first8_abs_vorticity": float(
                np.corrcoef(np.abs(signed_indent), vorticity)[0, 1]
            ),
        },
        "svd_modal_fraction_first8": modal_fraction[:8].tolist(),
    }
    with (args.output_dir / "npr15_front_history_summary.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    figure, axes = plt.subplots(2, 2, figsize=(13.0, 8.5), constrained_layout=True)
    axes[0, 0].plot(times_ms, fronts_array[:, 0], label="axis")
    annulus = (
        (radial_coordinate >= 0.3) & (radial_coordinate <= 0.5)
    )
    axes[0, 0].plot(
        times_ms, np.median(fronts_array[:, annulus], axis=1),
        label="r/D=0.3-0.5 median",
    )
    axes[0, 0].set(
        xlabel="time [ms]", ylabel="shock z/D", title="Disk breathing"
    )
    axes[0, 0].legend()

    axes[0, 1].plot(times_ms, signed_indent, label="center-annulus")
    axes[0, 1].plot(times_ms, full_range, label="full radial range")
    axes[0, 1].plot(times_ms, first8_range, label="first-8 range")
    axes[0, 1].axhline(0.0, color="black", lw=0.7)
    axes[0, 1].set(
        xlabel="time [ms]", ylabel="L4 cells", title="Signed shape metrics"
    )
    axes[0, 1].legend()

    front_offset = (
        fronts_array - np.median(fronts_array[:, annulus], axis=1)[:, None]
    ) / summary["finest_spacing_over_D"]
    image = axes[1, 0].pcolormesh(
        times_ms,
        radial_coordinate,
        front_offset.T,
        shading="nearest",
        cmap="coolwarm",
        vmin=-np.percentile(np.abs(front_offset), 99.0),
        vmax=np.percentile(np.abs(front_offset), 99.0),
    )
    figure.colorbar(image, ax=axes[1, 0], label="front offset [L4 cells]")
    axes[1, 0].set(
        xlabel="time [ms]", ylabel="r/D", title="Radial front mode"
    )

    for mode in range(min(4, spatial_modes.shape[0])):
        axes[1, 1].plot(
            radial_coordinate,
            spatial_modes[mode],
            label=f"mode {mode + 1}: {modal_fraction[mode]:.1%}",
        )
    axes[1, 1].set(
        xlabel="r/D", ylabel="normalized amplitude", title="Temporal SVD modes"
    )
    axes[1, 1].legend()

    for axis in axes.flat:
        axis.grid(alpha=0.22)
    if np.isfinite(stat_start):
        for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
            axis.axvline(stat_start, color="gray", ls=":", lw=0.8)
    figure.savefig(args.output_dir / "npr15_front_history.png", dpi=220)
    plt.close(figure)

    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
