#!/usr/bin/env python3
"""Audit Run165 thermodynamic energy over AMR leaf fluid cells."""

from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path

import numpy as np
import yt


GAMMA = 1.4
MOLECULAR_WEIGHT = 0.02896
UNIVERSAL_GAS_CONSTANT = 8.31446261815324
R_GAS = UNIVERSAL_GAS_CONSTANT / MOLECULAR_WEIGHT
CV = R_GAS / (GAMMA - 1.0)
CP = GAMMA * CV
REFERENCE_TOTAL_TEMPERATURE_K = 338.6666666666667
FAILURE_PROBE_M = np.array(
    [0.00290625, 0.005181818182, -0.007363636364], dtype=np.float64
)

REGIONS = {
    "jet_interaction": {"x": (-0.12, 0.02), "r_max": 0.08},
    "near_nozzle": {"x": (-0.005, 0.02), "r_max": 0.03},
    "forebody_nearfield": {"x": (-0.02, 0.10), "r_max": 0.08},
    "failure_probe_5mm": {"center": FAILURE_PROBE_M, "radius": 0.005},
    "failure_probe_10mm": {"center": FAILURE_PROBE_M, "radius": 0.010},
}

T_THRESHOLDS_K = (0.0, 5.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0, 100.0)
T0_THRESHOLDS_K = (0.0, 100.0, 150.0, 200.0, 250.0, 300.0, 330.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfiles", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--label", action="append", default=[])
    parser.add_argument("--low-records", type=int, default=32)
    return parser.parse_args()


def field(dataset, name: str):
    matches = [candidate for candidate in dataset.field_list if candidate[1] == name]
    if len(matches) != 1:
        raise RuntimeError(f"expected one field named {name!r}, found {matches}")
    return matches[0]


def empty_accumulator() -> dict:
    return {
        "leaf_cells": 0,
        "fluid_cells": 0,
        "fluid_volume_m3": 0.0,
        "mass_kg": 0.0,
        "total_energy_J": 0.0,
        "kinetic_energy_J": 0.0,
        "internal_energy_J": 0.0,
        "nonfinite_cells": 0,
        "rho_nonpositive_cells": 0,
        "internal_energy_nonpositive_cells": 0,
        "minimum_density_kg_m3": math.inf,
        "minimum_internal_energy_density_J_m3": math.inf,
        "minimum_temperature_K": math.inf,
        "minimum_total_temperature_K": math.inf,
        "maximum_mach": -math.inf,
        "temperature_threshold_counts": {str(value): 0 for value in T_THRESHOLDS_K},
        "temperature_threshold_volumes_m3": {
            str(value): 0.0 for value in T_THRESHOLDS_K
        },
        "total_temperature_threshold_counts": {
            str(value): 0 for value in T0_THRESHOLDS_K
        },
        "total_temperature_threshold_volumes_m3": {
            str(value): 0.0 for value in T0_THRESHOLDS_K
        },
        "reference_total_enthalpy_deficit_proxy_J": 0.0,
    }


def add_array_statistics(acc: dict, values: dict[str, np.ndarray], volume: float) -> None:
    rho = values["rho"]
    rho_energy = values["rho_energy"]
    kinetic = values["kinetic"]
    internal = values["internal"]
    temperature = values["temperature"]
    total_temperature = values["total_temperature"]
    mach = values["mach"]

    count = int(rho.size)
    acc["fluid_cells"] += count
    acc["fluid_volume_m3"] += count * volume
    acc["mass_kg"] += float(np.sum(rho, dtype=np.float64)) * volume
    acc["total_energy_J"] += float(np.sum(rho_energy, dtype=np.float64)) * volume
    acc["kinetic_energy_J"] += float(np.sum(kinetic, dtype=np.float64)) * volume
    acc["internal_energy_J"] += float(np.sum(internal, dtype=np.float64)) * volume

    finite = (
        np.isfinite(rho)
        & np.isfinite(rho_energy)
        & np.isfinite(internal)
        & np.isfinite(temperature)
        & np.isfinite(total_temperature)
        & np.isfinite(mach)
    )
    acc["nonfinite_cells"] += int(count - np.count_nonzero(finite))
    acc["rho_nonpositive_cells"] += int(np.count_nonzero(rho <= 0.0))
    acc["internal_energy_nonpositive_cells"] += int(np.count_nonzero(internal <= 0.0))

    if np.any(finite):
        acc["minimum_density_kg_m3"] = min(
            acc["minimum_density_kg_m3"], float(np.min(rho[finite]))
        )
        acc["minimum_internal_energy_density_J_m3"] = min(
            acc["minimum_internal_energy_density_J_m3"],
            float(np.min(internal[finite])),
        )
        acc["minimum_temperature_K"] = min(
            acc["minimum_temperature_K"], float(np.min(temperature[finite]))
        )
        acc["minimum_total_temperature_K"] = min(
            acc["minimum_total_temperature_K"],
            float(np.min(total_temperature[finite])),
        )
        acc["maximum_mach"] = max(acc["maximum_mach"], float(np.max(mach[finite])))

    for threshold in T_THRESHOLDS_K:
        selected = finite & (temperature < threshold)
        selected_count = int(np.count_nonzero(selected))
        acc["temperature_threshold_counts"][str(threshold)] += selected_count
        acc["temperature_threshold_volumes_m3"][str(threshold)] += (
            selected_count * volume
        )

    for threshold in T0_THRESHOLDS_K:
        selected = finite & (total_temperature < threshold)
        selected_count = int(np.count_nonzero(selected))
        acc["total_temperature_threshold_counts"][str(threshold)] += selected_count
        acc["total_temperature_threshold_volumes_m3"][str(threshold)] += (
            selected_count * volume
        )

    deficit = rho * CP * np.maximum(
        REFERENCE_TOTAL_TEMPERATURE_K - total_temperature, 0.0
    )
    acc["reference_total_enthalpy_deficit_proxy_J"] += (
        float(np.sum(deficit[finite], dtype=np.float64)) * volume
    )


def finalize_accumulator(acc: dict) -> dict:
    result = dict(acc)
    for key in (
        "minimum_density_kg_m3",
        "minimum_internal_energy_density_J_m3",
        "minimum_temperature_K",
        "minimum_total_temperature_K",
        "maximum_mach",
    ):
        if not math.isfinite(result[key]):
            result[key] = None
    if result["total_energy_J"] != 0.0:
        result["internal_energy_fraction"] = (
            result["internal_energy_J"] / result["total_energy_J"]
        )
        result["kinetic_energy_fraction"] = (
            result["kinetic_energy_J"] / result["total_energy_J"]
        )
    else:
        result["internal_energy_fraction"] = None
        result["kinetic_energy_fraction"] = None
    return result


def low_record(
    grid,
    local_index: tuple[int, int, int],
    arrays: dict[str, np.ndarray],
) -> dict:
    index = tuple(int(value) for value in local_index)
    xyz = [
        float(grid.LeftEdge[axis] + (index[axis] + 0.5) * grid.dds[axis])
        for axis in range(3)
    ]
    return {
        "level": int(grid.Level),
        "grid_id": int(grid.id),
        "local_index": list(index),
        "position_m": xyz,
        "density_kg_m3": float(arrays["rho"][index]),
        "total_energy_density_J_m3": float(arrays["rho_energy"][index]),
        "kinetic_energy_density_J_m3": float(arrays["kinetic"][index]),
        "internal_energy_density_J_m3": float(arrays["internal"][index]),
        "temperature_K": float(arrays["temperature"][index]),
        "total_temperature_K": float(arrays["total_temperature"][index]),
        "x_velocity_m_s": float(arrays["x_velocity"][index]),
        "y_velocity_m_s": float(arrays["y_velocity"][index]),
        "z_velocity_m_s": float(arrays["z_velocity"][index]),
        "speed_m_s": float(arrays["speed"][index]),
        "mach": float(arrays["mach"][index]),
        "sld": float(arrays["sld"][index]),
        "ghs": float(arrays["ghs"][index]),
    }


def region_mask(grid, specification: dict, fluid: np.ndarray) -> np.ndarray:
    x = float(grid.LeftEdge[0]) + (
        np.arange(fluid.shape[0], dtype=np.float64) + 0.5
    ) * float(grid.dds[0])
    y = float(grid.LeftEdge[1]) + (
        np.arange(fluid.shape[1], dtype=np.float64) + 0.5
    ) * float(grid.dds[1])
    z = float(grid.LeftEdge[2]) + (
        np.arange(fluid.shape[2], dtype=np.float64) + 0.5
    ) * float(grid.dds[2])
    if "center" in specification:
        center = specification["center"]
        distance_squared = (
            (x[:, None, None] - center[0]) ** 2
            + (y[None, :, None] - center[1]) ** 2
            + (z[None, None, :] - center[2]) ** 2
        )
        return fluid & (distance_squared <= specification["radius"] ** 2)
    radial_squared = y[None, :, None] ** 2 + z[None, None, :] ** 2
    return (
        fluid
        & (x[:, None, None] >= specification["x"][0])
        & (x[:, None, None] <= specification["x"][1])
        & (radial_squared <= specification["r_max"] ** 2)
    )


def nearest_probe_record(grid, fluid: np.ndarray, arrays: dict[str, np.ndarray]):
    lower = np.asarray(grid.LeftEdge, dtype=np.float64)
    upper = np.asarray(grid.RightEdge, dtype=np.float64)
    spacing = np.asarray(grid.dds, dtype=np.float64)
    if np.any(FAILURE_PROBE_M < lower - spacing) or np.any(
        FAILURE_PROBE_M > upper + spacing
    ):
        return None
    nominal = np.rint((FAILURE_PROBE_M - lower) / spacing - 0.5).astype(int)
    best = None
    for i in range(max(0, nominal[0] - 2), min(fluid.shape[0], nominal[0] + 3)):
        for j in range(max(0, nominal[1] - 2), min(fluid.shape[1], nominal[1] + 3)):
            for k in range(max(0, nominal[2] - 2), min(fluid.shape[2], nominal[2] + 3)):
                if not fluid[i, j, k]:
                    continue
                record = low_record(grid, (i, j, k), arrays)
                distance = float(
                    np.linalg.norm(np.asarray(record["position_m"]) - FAILURE_PROBE_M)
                )
                if best is None or distance < best[0]:
                    best = (distance, record)
    return best


def audit_plotfile(path: Path, label: str, low_record_count: int) -> dict:
    dataset = yt.load(str(path))
    fields = {
        name: field(dataset, name)
        for name in ("Xmom", "Ymom", "Zmom", "Energy", "Density", "sld", "ghs")
    }
    global_acc = empty_accumulator()
    level_acc = {
        str(level): empty_accumulator() for level in range(dataset.max_level + 1)
    }
    region_acc = {name: empty_accumulator() for name in REGIONS}
    low_temperature_heap: list[tuple[float, int, dict]] = []
    low_total_temperature_heap: list[tuple[float, int, dict]] = []
    probe_record = None
    serial = 0

    for grid in sorted(dataset.index.grids, key=lambda item: (item.Level, item.id)):
        leaf = np.asarray(grid.child_mask, dtype=bool)
        sld = np.asarray(grid[fields["sld"]], dtype=np.float64)
        ghs = np.asarray(grid[fields["ghs"]], dtype=np.float64)
        fluid = leaf & (sld < 0.5) & (ghs < 0.5)
        global_acc["leaf_cells"] += int(np.count_nonzero(leaf))
        level_acc[str(grid.Level)]["leaf_cells"] += int(np.count_nonzero(leaf))
        if not np.any(fluid):
            continue

        rho = np.asarray(grid[fields["Density"]], dtype=np.float64)
        rho_energy = np.asarray(grid[fields["Energy"]], dtype=np.float64)
        mx = np.asarray(grid[fields["Xmom"]], dtype=np.float64)
        my = np.asarray(grid[fields["Ymom"]], dtype=np.float64)
        mz = np.asarray(grid[fields["Zmom"]], dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_velocity = mx / rho
            y_velocity = my / rho
            z_velocity = mz / rho
            speed_squared = (mx * mx + my * my + mz * mz) / (rho * rho)
            speed = np.sqrt(speed_squared)
            kinetic = 0.5 * rho * speed_squared
            internal = rho_energy - kinetic
            temperature = internal / (rho * CV)
            total_temperature = temperature + speed_squared / (2.0 * CP)
            sound_speed_squared = GAMMA * R_GAS * temperature
            mach = speed / np.sqrt(sound_speed_squared)

        arrays = {
            "rho": rho,
            "rho_energy": rho_energy,
            "kinetic": kinetic,
            "internal": internal,
            "temperature": temperature,
            "total_temperature": total_temperature,
            "x_velocity": x_velocity,
            "y_velocity": y_velocity,
            "z_velocity": z_velocity,
            "speed": speed,
            "mach": mach,
            "sld": sld,
            "ghs": ghs,
        }
        selected_arrays = {name: values[fluid] for name, values in arrays.items()}
        cell_volume = float(np.prod(grid.dds.to_value()))
        add_array_statistics(global_acc, selected_arrays, cell_volume)
        add_array_statistics(level_acc[str(grid.Level)], selected_arrays, cell_volume)

        for name, specification in REGIONS.items():
            selected = region_mask(grid, specification, fluid)
            if not np.any(selected):
                continue
            region_acc[name]["leaf_cells"] += int(np.count_nonzero(selected))
            add_array_statistics(
                region_acc[name],
                {array_name: array[selected] for array_name, array in arrays.items()},
                cell_volume,
            )

        candidate_probe = nearest_probe_record(grid, fluid, arrays)
        if candidate_probe is not None and (
            probe_record is None or candidate_probe[0] < probe_record[0]
        ):
            probe_record = candidate_probe

        for quantity, heap in (
            ("temperature", low_temperature_heap),
            ("total_temperature", low_total_temperature_heap),
        ):
            data = arrays[quantity]
            flat_indices = np.flatnonzero(fluid & np.isfinite(data))
            if not flat_indices.size:
                continue
            candidates = min(low_record_count, flat_indices.size)
            local_order = np.argpartition(
                data.ravel()[flat_indices], candidates - 1
            )[:candidates]
            for flat_index in flat_indices[local_order]:
                index = np.unravel_index(flat_index, data.shape)
                record = low_record(grid, index, arrays)
                serial += 1
                item = (-record[f"{quantity}_K"], serial, record)
                if len(heap) < low_record_count:
                    heapq.heappush(heap, item)
                elif item > heap[0]:
                    heapq.heapreplace(heap, item)

    low_temperature_records = [item[2] for item in low_temperature_heap]
    low_temperature_records.sort(key=lambda item: item["temperature_K"])
    low_total_temperature_records = [
        item[2] for item in low_total_temperature_heap
    ]
    low_total_temperature_records.sort(
        key=lambda item: item["total_temperature_K"]
    )
    return {
        "label": label,
        "plotfile": str(path),
        "time_s": float(dataset.current_time),
        "max_level": int(dataset.max_level),
        "domain_dimensions": [int(value) for value in dataset.domain_dimensions],
        "global_leaf_fluid": finalize_accumulator(global_acc),
        "by_level_leaf_fluid": {
            level: finalize_accumulator(acc) for level, acc in level_acc.items()
        },
        "regions_leaf_fluid": {
            name: finalize_accumulator(acc) for name, acc in region_acc.items()
        },
        "failure_probe_m": FAILURE_PROBE_M.tolist(),
        "nearest_failure_probe_leaf_fluid_cell": (
            None
            if probe_record is None
            else {"distance_m": probe_record[0], **probe_record[1]}
        ),
        "lowest_temperature_cells": low_temperature_records,
        "lowest_total_temperature_cells": low_total_temperature_records,
    }


def main() -> None:
    args = parse_args()
    if args.label and len(args.label) != len(args.plotfiles):
        raise SystemExit("--label must be supplied once per plotfile or omitted")
    labels = args.label or [path.name for path in args.plotfiles]
    results = []
    for path, label in zip(args.plotfiles, labels):
        print(f"[energy-audit] loading {label}: {path}", flush=True)
        result = audit_plotfile(path, label, args.low_records)
        results.append(result)
        summary = result["global_leaf_fluid"]
        print(
            "[energy-audit] "
            f"{label} t={result['time_s']:.12g} "
            f"Tmin={summary['minimum_temperature_K']:.9g} "
            f"T0min={summary['minimum_total_temperature_K']:.9g} "
            f"T0<300K={summary['total_temperature_threshold_counts']['300.0']}",
            flush=True,
        )

    output = {
        "analysis": "run165_amr_leaf_fluid_energy_history",
        "gas_model": {
            "gamma": GAMMA,
            "molecular_weight_kg_mol": MOLECULAR_WEIGHT,
            "R_J_kg_K": R_GAS,
            "cv_J_kg_K": CV,
            "cp_J_kg_K": CP,
        },
        "selection": "AMR leaf cells with sld<0.5 and ghs<0.5",
        "reference_total_temperature_K": REFERENCE_TOTAL_TEMPERATURE_K,
        "reference_total_enthalpy_deficit_proxy_note": (
            "Integral rho*cp*max(T0_ref-T0,0)dV; diagnostic only, not a "
            "closed-domain conservation residual because Run165 has open boundaries."
        ),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="ascii")
    print(f"[energy-audit] wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
