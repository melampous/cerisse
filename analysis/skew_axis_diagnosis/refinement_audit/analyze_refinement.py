#!/usr/bin/env python3
"""Resolution audit of the residual near-axis line in the Run165 Skew RZ run.

This is deliberately a post-processing-only script.  It compares synchronous
uniform-L1 and uniform-L2 fields in physical coordinates and never changes a
plotfile or production source.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[3]
AXIS_READER = REPO / "analysis/skew_axis_diagnosis/analyze_axis_rings.py"
CASES = {
    "uniform_L1": HERE / "uniform_l1_200us/plt00500",
    "uniform_L2": HERE / "uniform_l2_200us/plt01000",
}
TIME_CONTROL = HERE / "uniform_l1_matched_fine_dt_200us/plt02000"
P_INF = 534.58067022560419
FIELDS = (
    "pressure",
    "density",
    "temperature",
    "radial_velocity",
    "axial_velocity",
    "mach",
)
FLOORS = {
    "pressure": 1.0,
    "density": 1.0e-12,
    "temperature": 1.0e-12,
    "radial_velocity": 1.0,
    "axial_velocity": 1.0,
    "mach": 1.0e-12,
}
SCALAR_FIELDS = ("pressure", "density", "temperature", "axial_velocity")


def load_reader():
    spec = importlib.util.spec_from_file_location("axis_reader", AXIS_READER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {AXIS_READER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.patch_yt_cylindrical_readonly_edge()
    return module


def load_case(reader, path: Path) -> dict:
    dataset, level, fields = reader.load_fields(path)
    nr, nz = fields["pressure"].shape
    rlo = float(dataset.domain_left_edge[0])
    rhi = float(dataset.domain_right_edge[0])
    zlo = float(dataset.domain_left_edge[1])
    zhi = float(dataset.domain_right_edge[1])
    dr = (rhi - rlo) / nr
    dz = (zhi - zlo) / nz
    return {
        "path": path,
        "dataset": dataset,
        "level": level,
        "r": rlo + (np.arange(nr) + 0.5) * dr,
        "z": zlo + (np.arange(nz) + 0.5) * dz,
        "dr": dr,
        "dz": dz,
        "fields": fields,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def radial_residual(values: np.ndarray, floor: float) -> np.ndarray:
    left = values[:-2]
    centre = values[1:-1]
    right = values[2:]
    scale = np.maximum.reduce(
        (np.abs(left), np.abs(centre), np.abs(right), np.full_like(centre, floor))
    )
    return (centre - 0.5 * (left + right)) / scale


def longest_true_run(mask: np.ndarray, spacing: float) -> float:
    padded = np.concatenate(([False], mask.astype(bool), [False]))
    jumps = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(jumps == 1)
    ends = np.flatnonzero(jumps == -1)
    return 0.0 if not starts.size else float(np.max(ends - starts) * spacing)


def contiguous_width(mask: np.ndarray, centre: int, spacing: float) -> float:
    lo = centre
    hi = centre
    while lo > 0 and mask[lo - 1]:
        lo -= 1
    while hi + 1 < mask.size and mask[hi + 1]:
        hi += 1
    return float((hi - lo + 1) * spacing)


def export_filament_metrics(data: dict[str, dict]) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    ring_rows: list[dict] = []
    for name, item in data.items():
        active = (item["z"] >= -0.060) & (item["z"] <= -0.020)
        nr_keep = 8 if item["level"] == 1 else 16
        for field in ("pressure", "density", "temperature", "radial_velocity", "mach"):
            residual = np.abs(
                radial_residual(item["fields"][field][:nr_keep], FLOORS[field])
            )
            p95 = np.quantile(residual[:, active], 0.95, axis=1)
            peak = np.max(residual[:, active], axis=1)
            mean = np.mean(residual[:, active], axis=1)
            for local_ring in range(residual.shape[0]):
                ring_rows.append(
                    {
                        "case": name,
                        "level": item["level"],
                        "dr_mm": f"{item['dr'] * 1.0e3:.12g}",
                        "field": field,
                        "centre_ring": local_ring + 1,
                        "r_mm": f"{item['r'][local_ring + 1] * 1.0e3:.12g}",
                        "mean_abs_residual": f"{mean[local_ring]:.17g}",
                        "p95_abs_residual": f"{p95[local_ring]:.17g}",
                        "peak_abs_residual": f"{peak[local_ring]:.17g}",
                        "longest_axial_over_2pct_mm": f"{longest_true_run(residual[local_ring, active] > 0.02, item['dz']) * 1.0e3:.17g}",
                    }
                )
            dominant = int(np.argmax(p95))
            mask = p95 > 0.02
            rows.append(
                {
                    "case": name,
                    "level": item["level"],
                    "dr_mm": f"{item['dr'] * 1.0e3:.12g}",
                    "field": field,
                    "dominant_centre_ring": dominant + 1,
                    "dominant_r_mm": f"{item['r'][dominant + 1] * 1.0e3:.12g}",
                    "dominant_p95_residual": f"{p95[dominant]:.17g}",
                    "dominant_peak_residual": f"{peak[dominant]:.17g}",
                    "dominant_longest_axial_over_2pct_mm": f"{longest_true_run(residual[dominant, active] > 0.02, item['dz']) * 1.0e3:.17g}",
                    "radial_width_p95_over_2pct_mm": f"{contiguous_width(mask, dominant, item['dr']) * 1.0e3:.17g}",
                    "first_r_mm_p95_over_2pct": (
                        f"{item['r'][np.flatnonzero(mask)[0] + 1] * 1.0e3:.12g}"
                        if np.any(mask) else "nan"
                    ),
                    "last_r_mm_p95_over_2pct": (
                        f"{item['r'][np.flatnonzero(mask)[-1] + 1] * 1.0e3:.12g}"
                        if np.any(mask) else "nan"
                    ),
                }
            )
    write_csv(HERE / "filament_resolution_metrics.csv", rows)
    write_csv(HERE / "first8_16_ring_metrics.csv", ring_rows)
    return rows, ring_rows


def export_parity_metrics(reader, data: dict[str, dict]) -> list[dict]:
    rows: list[dict] = []
    for name, item in data.items():
        active = (item["z"] >= -0.080) & (item["z"] <= -0.005)
        for field in SCALAR_FIELDS:
            _, residual, relative = reader.scalar_metrics(
                item["fields"][field], FLOORS[field]
            )
            local = np.abs(relative[active])
            rows.append(
                {
                    "case": name,
                    "level": item["level"],
                    "dr_mm": f"{item['dr'] * 1.0e3:.12g}",
                    "field": field,
                    "parity": "even",
                    "max_defect": f"{np.max(local):.17g}",
                    "p95_defect": f"{np.quantile(local, 0.95):.17g}",
                    "rms_defect": f"{np.sqrt(np.mean(local * local)):.17g}",
                    "max_signed_residual": f"{residual[active][np.argmax(local)]:.17g}",
                }
            )
        prediction = reader.odd_axis_prediction(item["fields"]["radial_velocity"])
        residual = item["fields"]["radial_velocity"][0] - prediction
        relative = np.abs(residual) / np.maximum(item["fields"]["sound_speed"][0], 1.0)
        local = relative[active]
        rows.append(
            {
                "case": name,
                "level": item["level"],
                "dr_mm": f"{item['dr'] * 1.0e3:.12g}",
                "field": "radial_velocity_over_sound",
                "parity": "odd",
                "max_defect": f"{np.max(local):.17g}",
                "p95_defect": f"{np.quantile(local, 0.95):.17g}",
                "rms_defect": f"{np.sqrt(np.mean(local * local)):.17g}",
                "max_signed_residual": f"{residual[active][np.argmax(local)]:.17g}",
            }
        )
    write_csv(HERE / "axis_parity_resolution_metrics.csv", rows)
    return rows


def export_time_step_control(reference: dict, control: dict) -> None:
    rows: list[dict] = []
    active = (reference["z"] >= -0.060) & (reference["z"] <= -0.020)
    for field in ("pressure", "density", "temperature", "radial_velocity", "mach"):
        summaries = []
        for item in (reference, control):
            residual = np.abs(radial_residual(item["fields"][field][:8], FLOORS[field]))
            p95 = np.quantile(residual[:, active], 0.95, axis=1)
            summaries.append((float(np.max(p95)), int(np.argmax(p95)) + 1))
        scale = np.maximum.reduce(
            (
                np.abs(reference["fields"][field][:8, active]),
                np.abs(control["fields"][field][:8, active]),
                np.full_like(reference["fields"][field][:8, active], FLOORS[field]),
            )
        )
        difference = np.abs(
            control["fields"][field][:8, active]
            - reference["fields"][field][:8, active]
        ) / scale
        rows.append(
            {
                "field": field,
                "reference_finest_dt_s": "2e-7",
                "control_finest_dt_s": "5e-8",
                "reference_dominant_ring": summaries[0][1],
                "control_dominant_ring": summaries[1][1],
                "reference_dominant_p95": f"{summaries[0][0]:.17g}",
                "control_dominant_p95": f"{summaries[1][0]:.17g}",
                "relative_change_in_dominant_p95": f"{summaries[1][0] / summaries[0][0] - 1.0:.17g}",
                "field_difference_p95": f"{np.quantile(difference, 0.95):.17g}",
                "field_difference_max": f"{np.max(difference):.17g}",
            }
        )
    write_csv(HERE / "l1_time_step_control.csv", rows)


def plot_common_fields(data: dict[str, dict]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.4), sharex=True, sharey=True,
                               constrained_layout=True)
    last_mach = None
    last_pressure = None
    for column, (name, item) in enumerate(data.items()):
        r_edges = (np.arange(item["r"].size + 1) * item["dr"]) * 1.0e3
        z_edges = (
            float(item["dataset"].domain_left_edge[1])
            + np.arange(item["z"].size + 1) * item["dz"]
        ) * 1.0e3
        last_mach = axes[0, column].pcolormesh(
            z_edges, r_edges, item["fields"]["mach"], shading="auto",
            cmap="turbo", vmin=0.0, vmax=6.0,
        )
        last_pressure = axes[1, column].pcolormesh(
            z_edges, r_edges, item["fields"]["pressure"] / P_INF,
            shading="auto", cmap="magma", norm=LogNorm(vmin=0.5, vmax=250.0),
        )
        axes[0, column].set_title(
            f"{name.replace('_', ' ')}, dr={item['dr'] * 1e3:.3f} mm"
        )
        axes[1, column].set_xlabel("axial z [mm]")
        for row in range(2):
            axes[row, column].set_xlim(-80.0, 0.5)
            axes[row, column].set_ylim(0.0, 12.0)
    axes[0, 0].set_ylabel("Mach\nr [mm]")
    axes[1, 0].set_ylabel("pressure\nr [mm]")
    figure.colorbar(last_mach, ax=axes[0, :], label="Mach")
    figure.colorbar(last_pressure, ax=axes[1, :], label=r"$p/p_\infty$")
    figure.suptitle("Run165 Skew K=0, synchronous t=200 us, common scales")
    figure.savefig(HERE / "mach_pressure_L1_L2_200us_common_scale.png", dpi=220)
    plt.close(figure)


def plot_residual_heatmaps(data: dict[str, dict]) -> None:
    fields = ("pressure", "density", "temperature", "mach")
    figure, axes = plt.subplots(len(fields), 2, figsize=(10.8, 10.0), sharex=True,
                               sharey=True, constrained_layout=True)
    image = None
    for column, (name, item) in enumerate(data.items()):
        keep = 10 if item["level"] == 1 else 20
        active = (item["z"] >= -0.080) & (item["z"] <= -0.005)
        z = item["z"][active] * 1.0e3
        z_edges = np.concatenate(
            ([z[0] - 0.5 * item["dz"] * 1.0e3],
             z + 0.5 * item["dz"] * 1.0e3)
        )
        for row, field in enumerate(fields):
            residual = np.abs(radial_residual(
                item["fields"][field][:keep], FLOORS[field]
            ))[:, active]
            r_edges = np.arange(1, keep) * item["dr"] * 1.0e3
            image = axes[row, column].pcolormesh(
                z_edges, r_edges, residual, shading="auto", cmap="magma",
                vmin=0.0, vmax=0.15,
            )
            if column == 0:
                axes[row, column].set_ylabel(f"{field}\nr [mm]")
            if row == 0:
                axes[row, column].set_title(name.replace("_", " "))
            if row == len(fields) - 1:
                axes[row, column].set_xlabel("axial z [mm]")
            axes[row, column].set_ylim(0.0, 9.0)
    figure.colorbar(image, ax=axes, label="absolute radial local-extremum residual")
    figure.suptitle("Near-axis line diagnostic at fixed physical radius, t=200 us")
    figure.savefig(HERE / "near_axis_residual_L1_L2_200us.png", dpi=220)
    plt.close(figure)


def plot_radial_residual_curves(ring_rows: list[dict]) -> None:
    fields = ("pressure", "density", "temperature", "radial_velocity", "mach")
    figure, axes = plt.subplots(2, 3, figsize=(11.2, 6.6), sharex=True,
                               constrained_layout=True)
    axes = axes.ravel()
    styles = {"uniform_L1": ("o-", "tab:blue"), "uniform_L2": (".-", "tab:orange")}
    for field, axis in zip(fields, axes):
        for case, (style, color) in styles.items():
            selected = [row for row in ring_rows if row["case"] == case and row["field"] == field]
            axis.plot(
                [float(row["r_mm"]) for row in selected],
                [float(row["p95_abs_residual"]) for row in selected],
                style, color=color, linewidth=1.3, markersize=4, label=case,
            )
        axis.axhline(0.02, color="0.45", linestyle=":", linewidth=0.9)
        axis.set_title(field)
        axis.set_xlim(0.0, 8.5)
        axis.set_ylim(bottom=0.0)
        axis.set_xlabel("physical radius r [mm]")
        axis.set_ylabel("p95 radial residual")
        axis.grid(alpha=0.22)
    axes[0].legend(fontsize=8)
    axes[-1].axis("off")
    figure.suptitle("First 8/16 rings compared at physical radius, -60 <= z <= -20 mm")
    figure.savefig(HERE / "radial_residual_vs_physical_radius_200us.png", dpi=220)
    plt.close(figure)


def main() -> None:
    reader = load_reader()
    data = {name: load_case(reader, path) for name, path in CASES.items()}
    time_control = load_case(reader, TIME_CONTROL)
    times = {name: float(item["dataset"].current_time) for name, item in data.items()}
    if max(times.values()) - min(times.values()) > 1.0e-12:
        raise RuntimeError(f"snapshots are not synchronous: {times}")
    filament, ring_rows = export_filament_metrics(data)
    parity = export_parity_metrics(reader, data)
    export_time_step_control(data["uniform_L1"], time_control)
    plot_common_fields(data)
    plot_residual_heatmaps(data)
    plot_radial_residual_curves(ring_rows)
    print("Filament dominant radii and widths (-60 <= z <= -20 mm):")
    for row in filament:
        print(
            f"{row['case']:10s} {row['field']:16s} "
            f"r={float(row['dominant_r_mm']):7.4f} mm "
            f"width={float(row['radial_width_p95_over_2pct_mm']):7.4f} mm "
            f"p95={float(row['dominant_p95_residual']):.5f}"
        )
    print("Axis parity max/p95 (-80 <= z <= -5 mm):")
    for row in parity:
        print(
            f"{row['case']:10s} {row['field']:28s} "
            f"max={float(row['max_defect']):.5f} p95={float(row['p95_defect']):.5f}"
        )


if __name__ == "__main__":
    main()
