#!/usr/bin/env python3
"""Re-evaluate NPR15 variants with a physical Mach-disk tracker."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_npr15_defect_ab import load_variant, require_common_grid
from npr15_shock_tracking import (
    track_primary_mach_disk,
    track_secondary_subsonic_compression,
)


CORE_R_MAX_D = 0.08
DISPLAY_R_MAX_D = 0.12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("historical_plot", type=Path)
    parser.add_argument("baseline_plot", type=Path)
    parser.add_argument("h_plot", type=Path)
    parser.add_argument("wb1_plot", type=Path)
    parser.add_argument("wb2_plot", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def add_physical_tracks(variant: dict) -> None:
    r_d = variant["r_d"]
    radial_count = int(np.count_nonzero(r_d <= DISPLAY_R_MAX_D))
    primary = track_primary_mach_disk(
        variant["rho_gradient"],
        variant["arrays"]["Density"],
        variant["mach"],
        variant["z_d"],
        radial_count,
    )
    secondary = track_secondary_subsonic_compression(
        variant["rho_gradient"],
        variant["mach"],
        variant["z_d"],
        primary,
        radial_count,
    )
    core = r_d[:radial_count] <= CORE_R_MAX_D
    spacing_d = variant["metrics"]["axial_spacing_over_D"]
    primary_core = primary["position"][core]
    if not np.all(np.isfinite(primary_core)):
        raise RuntimeError(
            f"{variant['label']} has an untracked primary shock in the core"
        )

    secondary_valid = np.isfinite(secondary["position"])
    dominance = (
        secondary_valid
        & np.isfinite(primary["gradient_peak"])
        & (secondary["gradient_peak"] > primary["gradient_peak"])
    )

    def finite_range_cells(values: np.ndarray) -> float:
        finite = values[np.isfinite(values)]
        return float(np.ptp(finite) / spacing_d)

    physical_metrics = {
        "core_r_max_over_D": CORE_R_MAX_D,
        "core_ring_count": int(np.count_nonzero(core)),
        "primary_core_valid_fraction": float(
            np.mean(np.isfinite(primary["position"][core]))
        ),
        "primary_core_mean_z_over_D": float(np.mean(primary_core)),
        "primary_core_range_cells": finite_range_cells(primary_core),
        "primary_first8_range_cells": finite_range_cells(
            primary["position"][:8]
        ),
        "primary_first16_range_cells": finite_range_cells(
            primary["position"][:16]
        ),
        "primary_axis_z_over_D": float(primary["position"][0]),
        "primary_axis_mach_before": float(primary["mach_before"][0]),
        "primary_axis_mach_after": float(primary["mach_after"][0]),
        "primary_axis_density_before": float(primary["density_before"][0]),
        "primary_axis_density_after": float(primary["density_after"][0]),
        "primary_axis_density_gradient": float(
            primary["gradient_peak"][0]
        ),
        "secondary_axis_z_over_D": float(secondary["position"][0]),
        "secondary_axis_mach_before": float(secondary["mach_before"][0]),
        "secondary_axis_mach_after": float(secondary["mach_after"][0]),
        "secondary_axis_density_gradient": float(
            secondary["gradient_peak"][0]
        ),
        "secondary_to_primary_axis_gradient_ratio": float(
            secondary["gradient_peak"][0] / primary["gradient_peak"][0]
        ),
        "secondary_stronger_radial_ring_count": int(
            np.count_nonzero(dominance)
        ),
        "secondary_stronger_r_max_over_D": float(
            np.max(r_d[:radial_count][dominance])
            if np.any(dominance) else 0.0
        ),
    }
    variant["primary"] = primary
    variant["secondary"] = secondary
    variant["physical_metrics"] = physical_metrics


def write_metrics(variants: list[dict], output_dir: Path) -> dict:
    by_label = {variant["label"]: variant for variant in variants}
    baseline = by_label["baseline"]
    keys = (
        "primary_core_mean_z_over_D",
        "primary_core_range_cells",
        "primary_first8_range_cells",
        "primary_first16_range_cells",
        "primary_axis_density_gradient",
        "secondary_axis_z_over_D",
        "secondary_axis_density_gradient",
        "secondary_to_primary_axis_gradient_ratio",
    )
    result = {
        "test": "NPR15 physical primary Mach-disk identity A/B",
        "criterion": (
            "first positive density-gradient peak with upstream Mach >= 1.5 "
            "and downstream Mach <= 1.5"
        ),
        "variants": {
            variant["label"]: variant["physical_metrics"]
            for variant in variants
        },
        "candidate_minus_baseline": {},
        "baseline_reproduction_relative_l2": {},
    }
    for label in ("h0125", "wb1", "wb2"):
        candidate = by_label[label]
        result["candidate_minus_baseline"][label] = {
            key: float(
                candidate["physical_metrics"][key]
                - baseline["physical_metrics"][key]
            )
            for key in keys
        }
    historical = by_label["historical"]
    for field in ("Density", "pressure", "x_velocity", "y_velocity"):
        reference = baseline["arrays"][field]
        scale = max(np.linalg.norm(reference), np.finfo(float).tiny)
        result["baseline_reproduction_relative_l2"][field] = float(
            np.linalg.norm(reference - historical["arrays"][field]) / scale
        )
    with (output_dir / "npr15_physical_disk_metrics.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    return result


def plot_tracks(variants: list[dict], output_dir: Path) -> None:
    by_label = {variant["label"]: variant for variant in variants}
    baseline = by_label["baseline"]
    r_d = baseline["r_d"]
    z_d = baseline["z_d"]
    display_count = baseline["primary"]["position"].size
    display_r = r_d[:display_count]
    spacing_d = baseline["metrics"]["axial_spacing_over_D"]

    figure, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), constrained_layout=True)
    for variant in variants:
        axes[0, 0].plot(
            variant["primary"]["position"],
            display_r,
            label=variant["label"],
        )
    axes[0, 0].set(
        xlabel="primary shock z/D",
        ylabel="r/D",
        title="Supersonic-to-subsonic Mach disk",
    )
    axes[0, 0].legend(fontsize=8)

    primary = baseline["primary"]
    secondary = baseline["secondary"]
    radial_mask = r_d <= DISPLAY_R_MAX_D
    axial_mask = (z_d >= 2.32) & (z_d <= 2.55)
    gradient = np.maximum(
        baseline["rho_gradient"][np.ix_(radial_mask, axial_mask)], 0.0
    )
    row_scale = np.maximum(
        np.max(gradient, axis=1, keepdims=True), np.finfo(float).tiny
    )
    image = axes[0, 1].pcolormesh(
        z_d[axial_mask],
        r_d[radial_mask],
        gradient / row_scale,
        shading="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
    )
    axes[0, 1].plot(
        primary["position"], display_r, color="cyan", lw=1.4,
        label="primary Mach disk",
    )
    axes[0, 1].plot(
        secondary["position"], display_r, color="white", lw=1.0,
        label="post-disk compression",
    )
    figure.colorbar(image, ax=axes[0, 1], label="row-normalized density gradient")
    axes[0, 1].set(
        xlabel="z/D", ylabel="r/D", title="Two distinct compression ridges"
    )
    axes[0, 1].legend(fontsize=8)

    for label in ("h0125", "wb1", "wb2"):
        candidate = by_label[label]
        displacement = (
            candidate["primary"]["position"] - primary["position"]
        ) / spacing_d
        axes[1, 0].plot(display_r, displacement, label=label)
    axes[1, 0].axhline(0.0, color="black", lw=0.8)
    axes[1, 0].set(
        xlabel="r/D",
        ylabel="candidate minus baseline [L4 cells]",
        title="Primary-disk displacement",
    )
    axes[1, 0].legend()

    for variant in variants[1:]:
        ratio = (
            variant["secondary"]["gradient_peak"]
            / variant["primary"]["gradient_peak"]
        )
        axes[1, 1].plot(display_r, ratio, label=variant["label"])
    axes[1, 1].axhline(1.0, color="black", lw=0.8)
    axes[1, 1].set(
        xlabel="r/D",
        ylabel="secondary/primary density-gradient peak",
        title="Why max-gradient tracking switches identity",
    )
    axes[1, 1].legend(fontsize=8)

    for axis in axes.flat:
        axis.grid(alpha=0.22)
    figure.savefig(output_dir / "npr15_physical_disk_tracks.png", dpi=220)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    variants = [
        load_variant("historical", args.historical_plot),
        load_variant("baseline", args.baseline_plot),
        load_variant("h0125", args.h_plot),
        load_variant("wb1", args.wb1_plot),
        load_variant("wb2", args.wb2_plot),
    ]
    require_common_grid(variants)
    for variant in variants:
        add_physical_tracks(variant)
    result = write_metrics(variants, args.output_dir)
    plot_tracks(variants, args.output_dir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
