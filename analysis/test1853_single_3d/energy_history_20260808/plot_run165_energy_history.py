#!/usr/bin/env python3
"""Create compact Run165 energy-history tables and plots from the AMR audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


NEAREST_WALL_NORMAL = np.array([-0.93887673, 0.19655647, -0.28262349])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit", type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--figure", required=True, type=Path)
    parser.add_argument("--zoom-figure", type=Path)
    return parser.parse_args()


def step_from_label(label: str) -> int:
    match = re.match(r"step(\d+)", label)
    if match is None:
        raise ValueError(f"cannot parse step from {label!r}")
    return int(match.group(1))


def threshold_fraction(record: dict, family: str, threshold: str) -> float:
    volume = record[f"{family}_threshold_volumes_m3"][threshold]
    return volume / record["fluid_volume_m3"]


def row_from_result(result: dict) -> dict:
    global_record = result["global_leaf_fluid"]
    probe = result["nearest_failure_probe_leaf_fluid_cell"]
    probe_region = result["regions_leaf_fluid"]["failure_probe_5mm"]
    low_temperature = result["lowest_temperature_cells"][0]
    low_position = low_temperature["position_m"]
    probe_position = probe["position_m"]
    probe_radius = math.hypot(probe_position[1], probe_position[2])
    probe_radial_velocity = (
        probe_position[1] * probe["y_velocity_m_s"]
        + probe_position[2] * probe["z_velocity_m_s"]
    ) / probe_radius
    probe_velocity = np.array(
        [
            probe["x_velocity_m_s"],
            probe["y_velocity_m_s"],
            probe["z_velocity_m_s"],
        ]
    )
    return {
        "step": step_from_label(result["label"]),
        "label": result["label"],
        "time_s": result["time_s"],
        "global_Tmin_K": global_record["minimum_temperature_K"],
        "global_T0min_K": global_record["minimum_total_temperature_K"],
        "global_Tmin_x_m": low_position[0],
        "global_Tmin_r_m": math.hypot(low_position[1], low_position[2]),
        "global_Tmin_cell_T0_K": low_temperature["total_temperature_K"],
        "global_Tmin_cell_speed_m_s": low_temperature["speed_m_s"],
        "probe_T_K": probe["temperature_K"],
        "probe_T0_K": probe["total_temperature_K"],
        "probe_speed_m_s": probe["speed_m_s"],
        "probe_x_velocity_m_s": probe["x_velocity_m_s"],
        "probe_y_velocity_m_s": probe["y_velocity_m_s"],
        "probe_z_velocity_m_s": probe["z_velocity_m_s"],
        "probe_radial_velocity_m_s": probe_radial_velocity,
        "probe_wall_normal_velocity_m_s": float(
            np.dot(probe_velocity, NEAREST_WALL_NORMAL)
        ),
        "probe_mach": probe["mach"],
        "probe_density_kg_m3": probe["density_kg_m3"],
        "probe_rhoE_J_m3": probe["total_energy_density_J_m3"],
        "probe_internal_energy_J_m3": probe["internal_energy_density_J_m3"],
        "probe_specific_total_energy_J_kg": (
            probe["total_energy_density_J_m3"] / probe["density_kg_m3"]
        ),
        "probe_specific_internal_energy_J_kg": (
            probe["internal_energy_density_J_m3"] / probe["density_kg_m3"]
        ),
        "probe5mm_Tmin_K": probe_region["minimum_temperature_K"],
        "probe5mm_T0min_K": probe_region["minimum_total_temperature_K"],
        "global_mass_kg": global_record["mass_kg"],
        "global_total_energy_J": global_record["total_energy_J"],
        "global_kinetic_energy_J": global_record["kinetic_energy_J"],
        "global_internal_energy_J": global_record["internal_energy_J"],
        "global_T0_lt_100_volume_fraction": threshold_fraction(
            global_record, "total_temperature", "100.0"
        ),
        "global_T0_lt_150_volume_fraction": threshold_fraction(
            global_record, "total_temperature", "150.0"
        ),
        "global_T_lt_30_volume_fraction": threshold_fraction(
            global_record, "temperature", "30.0"
        ),
        "probe5mm_T0_lt_150_volume_fraction": threshold_fraction(
            probe_region, "total_temperature", "150.0"
        ),
        "global_low_T0_inventory_proxy_J": global_record[
            "reference_total_enthalpy_deficit_proxy_J"
        ],
    }


def annotate_events(axis) -> None:
    axis.axvline(1339, color="0.45", linestyle=":", linewidth=1.0)
    axis.axvline(3367, color="#b2182b", linestyle="--", linewidth=1.0)


def main() -> None:
    args = parse_args()
    audit = json.loads(args.audit.read_text(encoding="ascii"))
    rows = sorted(
        (row_from_result(result) for result in audit["results"]),
        key=lambda row: row["step"],
    )

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    step = np.asarray([row["step"] for row in rows])
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), constrained_layout=True)

    axis = axes[0, 0]
    axis.plot(step, [row["global_Tmin_K"] for row in rows], "o-", label="global min T")
    axis.plot(step, [row["probe_T_K"] for row in rows], "s-", label="fixed near-wall probe T")
    axis.plot(step, [row["probe5mm_Tmin_K"] for row in rows], "^-", label="5 mm probe-region min T")
    axis.axhline(10.0, color="#b2182b", linewidth=1.0, label="10 K closure floor")
    annotate_events(axis)
    axis.set_ylabel("Static temperature [K]")
    axis.set_ylim(bottom=0.0)
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[0, 1]
    axis.plot(step, [row["probe_T0_K"] for row in rows], "s-", label="fixed probe T0")
    axis.plot(step, [row["probe5mm_T0min_K"] for row in rows], "^-", label="5 mm region min T0")
    axis.axhline(64.7298675, color="#2166ac", linestyle=":", label="initial quiescent T0")
    axis.axhline(338.6666667, color="#1b7837", linestyle=":", label="freestream T0")
    annotate_events(axis)
    axis.set_ylabel("Total-temperature diagnostic [K]")
    axis.set_ylim(bottom=0.0)
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[1, 0]
    axis.plot(step, [row["global_total_energy_J"] for row in rows], "o-", label="total")
    axis.plot(step, [row["global_kinetic_energy_J"] for row in rows], "s-", label="kinetic")
    axis.plot(step, [row["global_internal_energy_J"] for row in rows], "^-", label="internal")
    annotate_events(axis)
    axis.set_xlabel("Coarse step")
    axis.set_ylabel("Open-domain energy inventory [J]")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[1, 1]
    axis.plot(
        step,
        100.0 * np.asarray([row["global_T0_lt_150_volume_fraction"] for row in rows]),
        "o-",
        label="global volume T0 < 150 K",
    )
    axis.plot(
        step,
        100.0 * np.asarray([row["probe5mm_T0_lt_150_volume_fraction"] for row in rows]),
        "s-",
        label="5 mm region volume T0 < 150 K",
    )
    axis.plot(
        step,
        100.0 * np.asarray([row["global_T_lt_30_volume_fraction"] for row in rows]),
        "^-",
        label="global volume T < 30 K",
    )
    annotate_events(axis)
    axis.set_xlabel("Coarse step")
    axis.set_ylabel("Leaf-fluid volume fraction [%]")
    axis.set_ylim(bottom=0.0)
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    figure.suptitle(
        "Run165 full-field and near-wall energy history\n"
        "dotted: jet enabled at 1339; dashed: first strict 10 K crossing at 3367",
        fontsize=12,
    )
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180)
    plt.close(figure)

    if args.zoom_figure is not None:
        zoom_rows = [row for row in rows if 3200 <= row["step"] <= 3400]
        zoom_step = np.asarray([row["step"] for row in zoom_rows])
        zoom, zoom_axes = plt.subplots(
            2, 2, figsize=(12.0, 8.0), constrained_layout=True
        )

        axis = zoom_axes[0, 0]
        axis.plot(
            zoom_step,
            [row["probe_T_K"] for row in zoom_rows],
            "s-",
            label="fixed probe T",
        )
        axis.plot(
            zoom_step,
            [row["global_Tmin_K"] for row in zoom_rows],
            "o-",
            label="global min T",
        )
        axis.axhline(10.0, color="#b2182b", linewidth=1.0)
        annotate_events(axis)
        axis.set_ylabel("Static temperature [K]")
        axis.legend(fontsize=8)
        axis.grid(alpha=0.25)

        axis = zoom_axes[0, 1]
        axis.plot(
            zoom_step,
            [row["probe_T0_K"] for row in zoom_rows],
            "s-",
            label="fixed probe T0",
        )
        axis.plot(
            zoom_step,
            [row["probe5mm_T0min_K"] for row in zoom_rows],
            "^-",
            label="5 mm region min T0",
        )
        annotate_events(axis)
        axis.set_ylabel("Total-temperature diagnostic [K]")
        axis.legend(fontsize=8)
        axis.grid(alpha=0.25)

        axis = zoom_axes[1, 0]
        for key, label, marker in (
            ("probe_speed_m_s", "speed", "o-"),
            ("probe_x_velocity_m_s", "x velocity", "s-"),
            ("probe_radial_velocity_m_s", "radial velocity", "^-"),
            ("probe_wall_normal_velocity_m_s", "wall-normal velocity", "d-"),
        ):
            axis.plot(zoom_step, [row[key] for row in zoom_rows], marker, label=label)
        axis.axhline(0.0, color="0.35", linewidth=0.8)
        annotate_events(axis)
        axis.set_xlabel("Coarse step")
        axis.set_ylabel("Fixed-probe velocity [m/s]")
        axis.legend(fontsize=8)
        axis.grid(alpha=0.25)

        axis = zoom_axes[1, 1]
        total = np.asarray([row["probe_rhoE_J_m3"] for row in zoom_rows])
        internal = np.asarray(
            [row["probe_internal_energy_J_m3"] for row in zoom_rows]
        )
        axis.plot(zoom_step, total, "o-", label="rho E")
        axis.plot(zoom_step, total - internal, "s-", label="kinetic density")
        axis.plot(zoom_step, internal, "^-", label="internal-energy density")
        annotate_events(axis)
        axis.set_xlabel("Coarse step")
        axis.set_ylabel("Fixed-probe energy density [J/m3]")
        axis.legend(fontsize=8)
        axis.grid(alpha=0.25)

        for axis in zoom_axes.flat:
            axis.set_xlim(3180, 3420)

        zoom.suptitle(
            "Run165 near-wall failure window\n"
            "wall normal from the nearest STL triangle; negative velocity points into solid",
            fontsize=12,
        )
        args.zoom_figure.parent.mkdir(parents=True, exist_ok=True)
        zoom.savefig(args.zoom_figure, dpi=180)
        plt.close(zoom)


if __name__ == "__main__":
    main()
