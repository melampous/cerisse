#!/usr/bin/env python3
"""Compare the historical NPR15 pit with exact-time baseline and H restarts."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt


REPO_ROOT = Path(__file__).resolve().parents[3]
NPR_ANALYSIS_DIR = REPO_ROOT / "analysis" / "panda_benchmark" / "npr_sweep_case"
sys.path.insert(0, str(NPR_ANALYSIS_DIR))

from audit_npr15_disk import (  # noqa: E402
    GAMMA,
    RGAS,
    composite_region,
    gradient_fwhm,
    peak_index,
    positive_gradient,
    track_radial_front,
)


DIAMETER = 0.0254
TRACK_Z_MIN_D = 2.0
TRACK_Z_MAX_D = 2.65
SHAPE_R_MAX_D = 0.5


def subcell_peak(signal: np.ndarray, coordinate: np.ndarray, index: int) -> float:
    if index <= 0 or index + 1 >= signal.size:
        return float(coordinate[index])
    left, center, right = signal[index - 1:index + 2]
    denominator = left - 2.0 * center + right
    offset = 0.0
    if denominator != 0.0:
        offset = 0.5 * (left - right) / denominator
    offset = float(np.clip(offset, -0.5, 0.5))
    return float(
        coordinate[index] + offset * (coordinate[index + 1] - coordinate[index])
    )


def refine_front(
    gradient: np.ndarray, coordinate_d: np.ndarray, cell_front: np.ndarray
) -> np.ndarray:
    refined = np.empty_like(cell_front)
    for radial_index, position in enumerate(cell_front):
        index = int(np.argmin(np.abs(coordinate_d - position)))
        refined[radial_index] = subcell_peak(
            gradient[radial_index], coordinate_d, index
        )
    return refined


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("historical_plot", type=Path)
    parser.add_argument("baseline_plot", type=Path)
    parser.add_argument("h_plot", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def load_variant(label: str, plotfile: Path) -> dict:
    dataset = yt.load(str(plotfile))
    arrays, r, z, level_map, spacing, level = composite_region(
        dataset,
        ["Density", "pressure", "temperature", "x_velocity", "y_velocity"],
        DIAMETER,
        0.8,
        1.8,
        3.4,
    )
    r_d = r / DIAMETER
    z_d = z / DIAMETER
    shape_count = int(np.count_nonzero(r_d <= SHAPE_R_MAX_D))
    rho_gradient = positive_gradient(arrays["Density"], z)
    pressure_gradient = positive_gradient(arrays["pressure"], z)
    center_index = peak_index(
        rho_gradient[0], z_d, TRACK_Z_MIN_D, TRACK_Z_MAX_D
    )
    pressure_center_index = peak_index(
        pressure_gradient[0], z_d, TRACK_Z_MIN_D, TRACK_Z_MAX_D
    )
    density_front, density_peak = track_radial_front(
        rho_gradient, z_d, center_index, shape_count
    )
    pressure_front, pressure_peak = track_radial_front(
        pressure_gradient, z_d, pressure_center_index, shape_count
    )
    density_front_subcell = refine_front(
        rho_gradient[:shape_count], z_d, density_front
    )
    pressure_front_subcell = refine_front(
        pressure_gradient[:shape_count], z_d, pressure_front
    )
    spacing_d = float(spacing[1] / DIAMETER)
    annulus = (r_d[:shape_count] >= 0.3) & (r_d[:shape_count] <= 0.5)
    disk_band = np.abs(z_d - density_front[0]) <= 0.18
    radial_velocity = arrays["x_velocity"]
    axial_velocity = arrays["y_velocity"]
    vorticity = (
        np.gradient(radial_velocity, z, axis=1, edge_order=2)
        - np.gradient(axial_velocity, r, axis=0, edge_order=2)
    )
    mach = np.hypot(radial_velocity, axial_velocity) / np.sqrt(
        GAMMA * RGAS * arrays["temperature"]
    )
    level_roi = (
        (r_d[:, None] <= SHAPE_R_MAX_D)
        & (z_d[None, :] >= TRACK_Z_MIN_D)
        & (z_d[None, :] <= TRACK_Z_MAX_D)
    )

    def radial_range_cells(count: int) -> float:
        count = min(count, density_front.size)
        return float(np.ptp(density_front[:count]) / spacing_d)

    def radial_subcell_range_cells(count: int) -> float:
        count = min(count, density_front_subcell.size)
        return float(np.ptp(density_front_subcell[:count]) / spacing_d)

    metrics = {
        "label": label,
        "plotfile": str(plotfile),
        "time_s": float(dataset.current_time),
        "finest_level": int(level),
        "axial_spacing_over_D": spacing_d,
        "shape_radial_count": shape_count,
        "finest_fraction": float(np.mean(level_map[level_roi] == level)),
        "density_center_z_over_D": float(density_front[0]),
        "pressure_center_z_over_D": float(pressure_front[0]),
        "shape_range_over_D": float(np.ptp(density_front)),
        "shape_range_cells": float(np.ptp(density_front) / spacing_d),
        "center_minus_annulus_over_D": float(
            density_front[0] - np.median(density_front[annulus])
        ),
        "center_minus_annulus_cells": float(
            (density_front[0] - np.median(density_front[annulus])) / spacing_d
        ),
        "first8_range_cells": radial_range_cells(8),
        "first16_range_cells": radial_range_cells(16),
        "density_center_subcell_z_over_D": float(density_front_subcell[0]),
        "shape_range_subcell_cells": float(
            np.ptp(density_front_subcell) / spacing_d
        ),
        "center_minus_annulus_subcell_cells": float(
            (
                density_front_subcell[0]
                - np.median(density_front_subcell[annulus])
            )
            / spacing_d
        ),
        "first8_range_subcell_cells": radial_subcell_range_cells(8),
        "first16_range_subcell_cells": radial_subcell_range_cells(16),
        "axis_density_fwhm_cells": float(
            gradient_fwhm(rho_gradient[0], center_index, spacing_d) / spacing_d
        ),
        "axis_density_gradient_peak": float(density_peak[0]),
        "axis_pressure_gradient_peak": float(pressure_peak[0]),
        "axis_band_abs_ur_max_m_per_s": float(
            np.max(np.abs(radial_velocity[0, disk_band]))
        ),
        "first8_band_abs_ur_max_m_per_s": float(
            np.max(np.abs(radial_velocity[:8, :][:, disk_band]))
        ),
        "first8_band_abs_vorticity_max_per_s": float(
            np.max(np.abs(vorticity[:8, :][:, disk_band]))
        ),
        "density_min": float(np.min(arrays["Density"])),
        "pressure_min_pa": float(np.min(arrays["pressure"])),
    }
    return {
        "label": label,
        "metrics": metrics,
        "arrays": arrays,
        "r": r,
        "z": z,
        "r_d": r_d,
        "z_d": z_d,
        "density_front": density_front,
        "pressure_front": pressure_front,
        "density_front_subcell": density_front_subcell,
        "pressure_front_subcell": pressure_front_subcell,
        "rho_gradient": rho_gradient,
        "mach": mach,
        "vorticity": vorticity,
    }


def require_common_grid(variants: list[dict]) -> None:
    reference = variants[0]
    for variant in variants[1:]:
        if not (
            np.array_equal(reference["r"], variant["r"])
            and np.array_equal(reference["z"], variant["z"])
        ):
            raise RuntimeError("variant extraction grids differ")


def write_metrics(variants: list[dict], output_dir: Path) -> None:
    metrics = {
        "test": "NPR15 chk05000 to historical plt05098 exact-time A/B",
        "diameter_m": DIAMETER,
        "variants": {
            variant["label"]: variant["metrics"] for variant in variants
        },
    }
    times = [variant["metrics"]["time_s"] for variant in variants]
    metrics["max_time_difference_s"] = float(np.ptp(times))

    baseline = metrics["variants"]["baseline"]
    candidate = metrics["variants"]["h0125"]
    metrics["h_minus_baseline"] = {
        key: float(candidate[key] - baseline[key])
        for key in (
            "shape_range_cells",
            "center_minus_annulus_cells",
            "first8_range_cells",
            "first16_range_cells",
            "shape_range_subcell_cells",
            "center_minus_annulus_subcell_cells",
            "first8_range_subcell_cells",
            "first16_range_subcell_cells",
            "axis_density_fwhm_cells",
            "axis_band_abs_ur_max_m_per_s",
            "first8_band_abs_ur_max_m_per_s",
            "first8_band_abs_vorticity_max_per_s",
        )
    }
    historical_variant, baseline_variant, h_variant = variants
    metrics["field_difference_relative_l2"] = {}
    for field in ("Density", "pressure", "x_velocity", "y_velocity"):
        historical_field = historical_variant["arrays"][field]
        baseline_field = baseline_variant["arrays"][field]
        h_field = h_variant["arrays"][field]
        scale = max(np.linalg.norm(baseline_field), np.finfo(float).tiny)
        metrics["field_difference_relative_l2"][field] = {
            "baseline_minus_historical": float(
                np.linalg.norm(baseline_field - historical_field) / scale
            ),
            "h_minus_baseline": float(
                np.linalg.norm(h_field - baseline_field) / scale
            ),
        }
    with (output_dir / "npr15_defect_metrics.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(metrics, stream, indent=2)
        stream.write("\n")

    keys = list(variants[0]["metrics"])
    with (output_dir / "npr15_defect_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        for variant in variants:
            writer.writerow(variant["metrics"])


def plot_comparison(variants: list[dict], output_dir: Path) -> None:
    r_d = variants[0]["r_d"]
    z_d = variants[0]["z_d"]
    roi_z = (z_d >= 1.95) & (z_d <= 2.75)
    roi_r = r_d <= 0.65
    density_values = np.concatenate(
        [variant["arrays"]["Density"][np.ix_(roi_r, roi_z)].ravel()
         for variant in variants]
    )
    ur_limit = max(
        np.percentile(
            np.abs(variant["arrays"]["x_velocity"][np.ix_(roi_r, roi_z)]),
            99.5,
        )
        for variant in variants
    )
    gradient_values = np.concatenate(
        [np.abs(variant["rho_gradient"][np.ix_(roi_r, roi_z)]).ravel()
         for variant in variants]
    )
    gradient_floor = max(np.percentile(gradient_values, 5.0), 1.0e-30)
    gradient_ceiling = max(np.percentile(gradient_values, 99.8), gradient_floor)

    figure, axes = plt.subplots(
        len(variants), 4, figsize=(15.5, 10.0), constrained_layout=True
    )
    for row, variant in enumerate(variants):
        density = variant["arrays"]["Density"]
        radial_velocity = variant["arrays"]["x_velocity"]
        gradient = np.clip(
            np.abs(variant["rho_gradient"]), gradient_floor, gradient_ceiling
        )

        image = axes[row, 0].pcolormesh(
            z_d, r_d, density, shading="nearest", cmap="viridis",
            vmin=float(density_values.min()), vmax=float(density_values.max()),
        )
        axes[row, 0].plot(
            variant["density_front_subcell"],
            r_d[:variant["density_front_subcell"].size],
            color="white",
            lw=1.1,
        )
        figure.colorbar(image, ax=axes[row, 0], label="density")

        image = axes[row, 1].pcolormesh(
            z_d,
            r_d,
            np.log10(gradient),
            shading="nearest",
            cmap="magma",
            vmin=np.log10(gradient_floor),
            vmax=np.log10(gradient_ceiling),
        )
        figure.colorbar(image, ax=axes[row, 1], label="log10 abs(density_z)")

        image = axes[row, 2].pcolormesh(
            z_d,
            r_d,
            radial_velocity,
            shading="nearest",
            cmap="coolwarm",
            vmin=-ur_limit,
            vmax=ur_limit,
        )
        figure.colorbar(image, ax=axes[row, 2], label="radial velocity [m/s]")

        axis = axes[row, 3]
        for other in variants:
            axis.plot(
                other["density_front_subcell"],
                r_d[:other["density_front_subcell"].size],
                lw=1.0 if other is variant else 0.7,
                alpha=1.0 if other is variant else 0.45,
                label=other["label"],
            )
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
        axis.set_xlabel("shock z/D")
        axis.set_ylabel("r/D")

        for column in range(3):
            axes[row, column].set(
                xlim=(1.95, 2.75), ylim=(0.0, 0.65),
                xlabel="z/D", ylabel="r/D",
            )
        axes[row, 0].set_title(f"{variant['label']}: density + tracked front")
        axes[row, 1].set_title("shock image")
        axes[row, 2].set_title("radial velocity")
        axes[row, 3].set_title("front overlay")

    figure.savefig(output_dir / "npr15_defect_fields.png", dpi=220)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    for variant in variants:
        axes[0].plot(
            r_d[:variant["density_front_subcell"].size],
            variant["density_front_subcell"],
            label=variant["label"],
        )
    axes[0].set(xlabel="r/D", ylabel="shock z/D", title="Mach-disk locus")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    baseline = next(item for item in variants if item["label"] == "baseline")
    h_variant = next(item for item in variants if item["label"] == "h0125")
    count = min(
        baseline["density_front_subcell"].size,
        h_variant["density_front_subcell"].size,
    )
    spacing_d = baseline["metrics"]["axial_spacing_over_D"]
    axes[1].plot(
        r_d[:count],
        (h_variant["density_front_subcell"][:count]
         - baseline["density_front_subcell"][:count]) / spacing_d,
        color="black",
    )
    axes[1].axhline(0.0, color="gray", lw=0.8)
    axes[1].set(
        xlabel="r/D",
        ylabel="H minus baseline shock z [L4 cells]",
        title="Candidate displacement",
    )
    axes[1].grid(alpha=0.25)
    figure.savefig(output_dir / "npr15_defect_profiles.png", dpi=220)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    yt.funcs.mylog.setLevel(50)
    variants = [
        load_variant("historical", args.historical_plot),
        load_variant("baseline", args.baseline_plot),
        load_variant("h0125", args.h_plot),
    ]
    require_common_grid(variants)
    write_metrics(variants, args.output_dir)
    plot_comparison(variants, args.output_dir)
    print(
        json.dumps(
            {variant["label"]: variant["metrics"] for variant in variants},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
