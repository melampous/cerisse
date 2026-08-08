#!/usr/bin/env python3
"""Analyse the AFD HLLC shock-LLF switch for an advected contact."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np

try:
    import yt
except ImportError as exc:  # pragma: no cover
    raise SystemExit("analyze_ab.py requires yt") from exc


GAMMA = 1.4
RHO_LOW = 1.0
RHO_HIGH = 2.0
P0 = 1.0
UX0 = 0.5 * math.sqrt(GAMMA * P0 / RHO_LOW)
SMOOTHNESS_THRESHOLD = 0.08
SHOCK_PRESSURE_JUMP = 0.01
SHOCK_COMPRESSION = 0.001
GRIDS = (64, 128, 256, 512)
MODES = (0, 1)


def latest_plotfile(path: Path) -> Path:
    pattern = re.compile(r"^plt(\d+)$")
    candidates = [
        candidate
        for candidate in path.glob("plt*")
        if pattern.fullmatch(candidate.name)
        and (candidate / "Header").is_file()
    ]
    if not candidates:
        raise FileNotFoundError(f"no plotfile found in {path}")
    return max(
        candidates,
        key=lambda candidate: int(pattern.fullmatch(candidate.name)[1]),
    )


def field_key(ds, name: str):
    wanted = name.lower()
    for key in ds.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"{name!r} is absent from {ds.field_list}")


def ordered_state(plotfile: Path) -> dict[str, np.ndarray | float]:
    ds = yt.load(str(plotfile))
    if int(ds.index.max_level) != 0:
        raise ValueError(f"{plotfile} is not a Level-0-only calculation")
    data = ds.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    order = np.argsort(x)
    volume = np.asarray(data[("index", "cell_volume")].d, dtype=float)[order]
    rho = np.asarray(data[field_key(ds, "Density")].d, dtype=float)[order]
    xmom = np.asarray(data[field_key(ds, "Xmom")].d, dtype=float)[order]
    ymom = np.asarray(data[field_key(ds, "Ymom")].d, dtype=float)[order]
    zmom = np.asarray(data[field_key(ds, "Zmom")].d, dtype=float)[order]
    energy = np.asarray(data[field_key(ds, "Energy")].d, dtype=float)[order]
    ux = xmom / rho
    uy = ymom / rho
    uz = zmom / rho
    pressure = (GAMMA - 1.0) * (
        energy - 0.5 * rho * (ux * ux + uy * uy + uz * uz)
    )
    sound = np.sqrt(GAMMA * pressure / rho)
    return {
        "x": x[order],
        "volume": volume,
        "rho": rho,
        "ux": ux,
        "uy": uy,
        "uz": uz,
        "pressure": pressure,
        "sound": sound,
        "xmom": xmom,
        "ymom": ymom,
        "zmom": zmom,
        "energy": energy,
        "time": float(ds.current_time),
    }


def sensor_faces(state: dict[str, np.ndarray | float]) -> dict[str, object]:
    """Reproduce the periodic one-dimensional bulk AFD shock gate."""

    fields = (
        ("density", np.asarray(state["rho"])),
        ("pressure", np.asarray(state["pressure"])),
        ("x_velocity", np.asarray(state["ux"])),
        ("y_velocity", np.asarray(state["uy"])),
        ("z_velocity", np.asarray(state["uz"])),
    )
    sound = np.asarray(state["sound"])
    pressure = np.asarray(state["pressure"])
    normal_velocity = np.asarray(state["ux"])
    ncell = sound.size
    broad_faces = 0
    pressure_faces = 0
    compression_faces = 0
    pressure_compression_faces = 0
    fallback_faces = 0
    field_counts = {name: 0 for name, _ in fields}
    tiny = np.finfo(float).tiny

    for face in range(ncell):
        ids = np.array([(face + m - 3) % ncell for m in range(6)])
        sound_scale = float(np.max(np.abs(sound[ids])))
        broad = False
        for field_index, (name, values) in enumerate(fields):
            field_scale = float(np.max(np.abs(values[ids])))
            floor = sound_scale if field_index >= 2 else field_scale
            field_triggered = False
            for m in range(4):
                a, b, c = values[ids[m : m + 3]]
                denominator = abs(a) + 2.0 * abs(b) + abs(c) + floor + tiny
                if (
                    abs(a - 2.0 * b + c) / denominator
                    > SMOOTHNESS_THRESHOLD
                ):
                    field_triggered = True
                    broad = True
            if field_triggered:
                field_counts[name] += 1

        pressure_jump = any(
            abs(pressure[ids[m + 1]] - pressure[ids[m]])
            / max(pressure[ids[m]], pressure[ids[m + 1]], tiny)
            > SHOCK_PRESSURE_JUMP
            for m in range(5)
        )
        compression = any(
            max(
                normal_velocity[(ids[m] - 1) % ncell]
                - normal_velocity[(ids[m] + 1) % ncell],
                0.0,
            )
            / (
                2.0
                * max(
                    sound[(ids[m] - 1) % ncell],
                    sound[ids[m]],
                    sound[(ids[m] + 1) % ncell],
                    tiny,
                )
            )
            > SHOCK_COMPRESSION
            for m in range(1, 5)
        )
        broad_faces += int(broad)
        pressure_faces += int(pressure_jump)
        compression_faces += int(compression)
        pressure_compression_faces += int(pressure_jump and compression)
        fallback_faces += int(broad and pressure_jump and compression)

    return {
        "broad_faces": broad_faces,
        "pressure_faces": pressure_faces,
        "compression_faces": compression_faces,
        "pressure_compression_faces": pressure_compression_faces,
        "fallback_faces": fallback_faces,
        "field_counts": field_counts,
    }


def exact_density(x: np.ndarray) -> np.ndarray:
    return np.where((x >= 0.25) & (x < 0.75), RHO_HIGH, RHO_LOW)


def exact_initial_state(ncell: int) -> dict[str, np.ndarray]:
    x = (np.arange(ncell, dtype=float) + 0.5) / ncell
    rho = exact_density(x)
    pressure = np.full(ncell, P0)
    ux = np.full(ncell, UX0)
    zeros = np.zeros(ncell)
    return {
        "rho": rho,
        "pressure": pressure,
        "ux": ux,
        "uy": zeros,
        "uz": zeros,
        "sound": np.sqrt(GAMMA * pressure / rho),
    }


def analyse_case(case_dir: Path, mode: int, ncell: int) -> dict[str, object]:
    plotfile = latest_plotfile(case_dir)
    state = ordered_state(plotfile)
    x = np.asarray(state["x"])
    volume = np.asarray(state["volume"])
    rho = np.asarray(state["rho"])
    pressure = np.asarray(state["pressure"])
    xmom = np.asarray(state["xmom"])
    energy = np.asarray(state["energy"])
    total_volume = float(np.sum(volume))
    exact = exact_density(x)
    sensor = sensor_faces(state)
    initial_sensor = sensor_faces(exact_initial_state(ncell))
    field_counts = sensor["field_counts"]
    initial_field_counts = initial_sensor["field_counts"]

    density_error = rho - exact
    periodic_tv = float(np.sum(np.abs(rho - np.roll(rho, 1))))
    exact_tv = 2.0 * (RHO_HIGH - RHO_LOW)
    mixed = (rho > RHO_LOW + 0.01) & (rho < RHO_HIGH - 0.01)

    exact_mass = 1.5
    exact_momentum = exact_mass * UX0
    exact_energy = P0 / (GAMMA - 1.0) + 0.5 * exact_mass * UX0**2
    mass = float(np.sum(rho * volume))
    momentum = float(np.sum(xmom * volume))
    total_energy = float(np.sum(energy * volume))

    return {
        "afd_shock_llf": mode,
        "ncell": ncell,
        "cfl": 0.20,
        "afd_smoothness_threshold": SMOOTHNESS_THRESHOLD,
        "afd_shock_pressure_jump": SHOCK_PRESSURE_JUMP,
        "afd_shock_compression": SHOCK_COMPRESSION,
        "plotfile": str(plotfile),
        "time": float(state["time"]),
        "density_l1": float(np.sum(np.abs(density_error) * volume) / total_volume),
        "density_l2": float(
            math.sqrt(np.sum(density_error**2 * volume) / total_volume)
        ),
        "density_linf": float(np.max(np.abs(density_error))),
        "density_min": float(np.min(rho)),
        "density_max": float(np.max(rho)),
        "density_undershoot": float(max(0.0, RHO_LOW - np.min(rho))),
        "density_overshoot": float(max(0.0, np.max(rho) - RHO_HIGH)),
        "tv_ratio": periodic_tv / exact_tv,
        "mixed_cells_1pct": int(np.count_nonzero(mixed)),
        "mixed_width_fraction_1pct": float(np.count_nonzero(mixed) / ncell),
        "pressure_min": float(np.min(pressure)),
        "pressure_max": float(np.max(pressure)),
        "pressure_linf_error": float(np.max(np.abs(pressure - P0))),
        "offline_broad_nonsmooth_faces": sensor["broad_faces"],
        "offline_pressure_jump_faces": sensor["pressure_faces"],
        "offline_compression_faces": sensor["compression_faces"],
        "offline_pressure_compression_faces":
            sensor["pressure_compression_faces"],
        "offline_llf_fallback_faces": sensor["fallback_faces"],
        "offline_llf_fallback_fraction": sensor["fallback_faces"] / ncell,
        "initial_broad_nonsmooth_faces": initial_sensor["broad_faces"],
        "initial_pressure_jump_faces": initial_sensor["pressure_faces"],
        "initial_compression_faces": initial_sensor["compression_faces"],
        "initial_pressure_compression_faces":
            initial_sensor["pressure_compression_faces"],
        "initial_llf_fallback_faces": initial_sensor["fallback_faces"],
        "initial_llf_fallback_fraction":
            initial_sensor["fallback_faces"] / ncell,
        "initial_broad_density_trigger_faces":
            initial_field_counts["density"],
        "offline_broad_density_trigger_faces": field_counts["density"],
        "offline_broad_pressure_trigger_faces": field_counts["pressure"],
        "offline_broad_x_velocity_trigger_faces":
            field_counts["x_velocity"],
        "mass_relative_error": (mass - exact_mass) / exact_mass,
        "momentum_relative_error": (momentum - exact_momentum) / exact_momentum,
        "energy_relative_error": (total_energy - exact_energy) / exact_energy,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_comparison(path: Path, rows: list[dict[str, object]]) -> None:
    by_key = {(int(row["afd_shock_llf"]), int(row["ncell"])): row for row in rows}
    columns = (
        "density_l1",
        "mixed_cells_1pct",
        "tv_ratio",
        "pressure_linf_error",
        "offline_broad_nonsmooth_faces",
        "offline_pressure_jump_faces",
        "offline_compression_faces",
        "offline_pressure_compression_faces",
        "offline_llf_fallback_faces",
    )
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["ncell"]
            + [f"{name}_mode0" for name in columns]
            + [f"{name}_mode1" for name in columns]
            + ["density_l1_mode1_over_mode0"]
        )
        for ncell in GRIDS:
            off = by_key[(0, ncell)]
            on = by_key[(1, ncell)]
            writer.writerow(
                [ncell]
                + [off[name] for name in columns]
                + [on[name] for name in columns]
                + [float(on["density_l1"]) / float(off["density_l1"])]
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    rows = [
        analyse_case(
            args.root / "runs" / f"shock_llf_{mode}" / f"N{ncell}",
            mode,
            ncell,
        )
        for mode in MODES
        for ncell in GRIDS
    ]
    write_csv(args.root / "metrics.csv", rows)
    write_comparison(args.root / "comparison.csv", rows)


if __name__ == "__main__":
    main()
