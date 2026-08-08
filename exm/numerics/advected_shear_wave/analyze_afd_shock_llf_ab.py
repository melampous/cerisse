#!/usr/bin/env python3
"""Analyse the shock-LLF switch for the periodic advected shear wave."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np
import yt


GAMMA = 1.4
SMOOTHNESS_THRESHOLD = 0.08
SHOCK_PRESSURE_JUMP = 0.01
SHOCK_COMPRESSION = 0.001
PPW = (8, 12, 16, 24, 32)
MODES = (0, 1)
SUITES = {
    "linear_mode1": {"mode": 1, "amplitude_over_c0": 1.0e-5},
    "finite_mode4": {"mode": 4, "amplitude_over_c0": 0.5},
}


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


def ordered_primitives(plotfile: Path) -> dict[str, np.ndarray | float]:
    ds = yt.load(str(plotfile))
    data = ds.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    order = np.argsort(x)
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
        "volume": np.asarray(data[("index", "cell_volume")].d, dtype=float)[order],
        "rho": rho,
        "pressure": pressure,
        "ux": ux,
        "uy": uy,
        "uz": uz,
        "sound": sound,
        "xmom": xmom,
        "ymom": ymom,
        "zmom": zmom,
        "energy": energy,
        "time": float(ds.current_time),
    }


def sensor_faces(state: dict[str, np.ndarray]) -> dict[str, object]:
    fields = (
        ("density", state["rho"]),
        ("pressure", state["pressure"]),
        ("x_velocity", state["ux"]),
        ("y_velocity", state["uy"]),
        ("z_velocity", state["uz"]),
    )
    sound = state["sound"]
    pressure = state["pressure"]
    normal_velocity = state["ux"]
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


def wrapped_difference(angle: float, reference: float) -> float:
    return math.atan2(
        math.sin(angle - reference),
        math.cos(angle - reference),
    )


def exact_initial_state(
    wave_mode: int, amplitude_over_c0: float, ncell: int
) -> dict[str, np.ndarray]:
    x = (np.arange(ncell, dtype=float) + 0.5) / ncell
    c0 = math.sqrt(GAMMA)
    amplitude = amplitude_over_c0 * c0
    return {
        "rho": np.ones(ncell),
        "pressure": np.ones(ncell),
        "ux": np.full(ncell, 0.5 * c0),
        "uy": amplitude * np.sin(2.0 * math.pi * wave_mode * x),
        "uz": np.zeros(ncell),
        "sound": np.full(ncell, c0),
    }


def analyse_case(
    case_dir: Path,
    suite: str,
    mode_switch: int,
    ppw: int,
) -> dict[str, object]:
    setup = SUITES[suite]
    wave_mode = int(setup["mode"])
    amplitude_over_c0 = float(setup["amplitude_over_c0"])
    state = ordered_primitives(latest_plotfile(case_dir))
    x = np.asarray(state["x"])
    volume = np.asarray(state["volume"])
    rho = np.asarray(state["rho"])
    uy = np.asarray(state["uy"])
    pressure = np.asarray(state["pressure"])
    time = float(state["time"])
    total_volume = float(np.sum(volume))
    c0 = math.sqrt(GAMMA)
    ux0 = 0.5 * c0
    amplitude0 = amplitude_over_c0 * c0
    wave_number = 2.0 * math.pi * wave_mode

    sine_coefficient = float(
        2.0 * np.sum(uy * np.sin(wave_number * x) * volume) / total_volume
    )
    cosine_coefficient = float(
        2.0 * np.sum(uy * np.cos(wave_number * x) * volume) / total_volume
    )
    amplitude = math.hypot(sine_coefficient, cosine_coefficient)
    amplitude_ratio = amplitude / amplitude0
    measured_phase = math.atan2(-cosine_coefficient, sine_coefficient)
    exact_phase = wave_number * ux0 * time
    phase_error = wrapped_difference(measured_phase, exact_phase)
    phase_speed_error = phase_error / exact_phase
    uy_mean = float(np.sum(uy * volume) / total_volume)
    perturbation_ke = float(
        np.sum(0.5 * rho * (uy - uy_mean) ** 2 * volume)
    )
    initial_ke = 0.25 * amplitude0**2 * total_volume
    sensor = sensor_faces(state)
    initial_sensor = sensor_faces(
        exact_initial_state(wave_mode, amplitude_over_c0, ppw * wave_mode)
    )
    field_counts = sensor["field_counts"]
    initial_field_counts = initial_sensor["field_counts"]
    equivalent_viscosity = (
        -math.log(amplitude_ratio) / (wave_number**2 * time)
        if amplitude_ratio > 0.0
        else float("nan")
    )
    return {
        "suite": suite,
        "wave_mode": wave_mode,
        "amplitude_over_c0": amplitude_over_c0,
        "afd_shock_llf": mode_switch,
        "ppw": ppw,
        "ncell": ppw * wave_mode,
        "cfl": 0.05,
        "afd_smoothness_threshold": SMOOTHNESS_THRESHOLD,
        "afd_shock_pressure_jump": SHOCK_PRESSURE_JUMP,
        "afd_shock_compression": SHOCK_COMPRESSION,
        "plotfile": str(latest_plotfile(case_dir)),
        "time": time,
        "amplitude": amplitude,
        "amplitude_ratio": amplitude_ratio,
        "ke_ratio": perturbation_ke / initial_ke,
        "phase_error_rad": phase_error,
        "phase_speed_error_over_ux": phase_speed_error,
        "nu_num": equivalent_viscosity,
        "density_min": float(np.min(rho)),
        "density_max": float(np.max(rho)),
        "density_linf_error": float(np.max(np.abs(rho - 1.0))),
        "pressure_min": float(np.min(pressure)),
        "pressure_max": float(np.max(pressure)),
        "pressure_linf_error": float(np.max(np.abs(pressure - 1.0))),
        "offline_broad_nonsmooth_faces": sensor["broad_faces"],
        "offline_pressure_jump_faces": sensor["pressure_faces"],
        "offline_compression_faces": sensor["compression_faces"],
        "offline_pressure_compression_faces":
            sensor["pressure_compression_faces"],
        "offline_llf_fallback_faces": sensor["fallback_faces"],
        "offline_llf_fallback_fraction":
            sensor["fallback_faces"] / (ppw * wave_mode),
        "initial_broad_nonsmooth_faces": initial_sensor["broad_faces"],
        "initial_pressure_jump_faces": initial_sensor["pressure_faces"],
        "initial_compression_faces": initial_sensor["compression_faces"],
        "initial_pressure_compression_faces":
            initial_sensor["pressure_compression_faces"],
        "initial_llf_fallback_faces": initial_sensor["fallback_faces"],
        "initial_llf_fallback_fraction":
            initial_sensor["fallback_faces"] / (ppw * wave_mode),
        "initial_broad_y_velocity_trigger_faces":
            initial_field_counts["y_velocity"],
        "offline_broad_density_trigger_faces": field_counts["density"],
        "offline_broad_pressure_trigger_faces": field_counts["pressure"],
        "offline_broad_x_velocity_trigger_faces":
            field_counts["x_velocity"],
        "offline_broad_y_velocity_trigger_faces":
            field_counts["y_velocity"],
    }


def analyse(root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for suite in SUITES:
        for mode_switch in MODES:
            for ppw in PPW:
                case_dir = (
                    root
                    / "runs"
                    / suite
                    / f"shock_llf_{mode_switch}"
                    / f"PPW{ppw}"
                )
                rows.append(analyse_case(case_dir, suite, mode_switch, ppw))
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_comparison(
    path: Path, rows: list[dict[str, object]], root: Path
) -> None:
    by_key = {
        (str(row["suite"]), int(row["afd_shock_llf"]), int(row["ppw"])): row
        for row in rows
    }
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "suite",
                "wave_mode",
                "ppw",
                "ncell",
                "amplitude_ratio_mode0",
                "amplitude_ratio_mode1",
                "amplitude_ratio_difference",
                "ke_ratio_mode0",
                "ke_ratio_mode1",
                "ke_ratio_difference",
                "phase_error_difference_rad",
                "offline_broad_nonsmooth_faces_mode0",
                "offline_broad_nonsmooth_faces_mode1",
                "offline_pressure_jump_faces_mode0",
                "offline_pressure_jump_faces_mode1",
                "offline_compression_faces_mode0",
                "offline_compression_faces_mode1",
                "offline_pressure_compression_faces_mode0",
                "offline_pressure_compression_faces_mode1",
                "offline_llf_fallback_faces_mode0",
                "offline_llf_fallback_faces_mode1",
                "max_conservative_state_difference",
            ]
        )
        for suite, setup in SUITES.items():
            wave_mode = int(setup["mode"])
            for ppw in PPW:
                off = by_key[(suite, 0, ppw)]
                on = by_key[(suite, 1, ppw)]
                state_off = ordered_primitives(
                    latest_plotfile(
                        root
                        / "runs"
                        / suite
                        / "shock_llf_0"
                        / f"PPW{ppw}"
                    )
                )
                state_on = ordered_primitives(
                    latest_plotfile(
                        root
                        / "runs"
                        / suite
                        / "shock_llf_1"
                        / f"PPW{ppw}"
                    )
                )
                maximum = max(
                    float(np.max(np.abs(state_on[name] - state_off[name])))
                    for name in ("rho", "xmom", "ymom", "zmom", "energy")
                )
                writer.writerow(
                    [
                        suite,
                        wave_mode,
                        ppw,
                        ppw * wave_mode,
                        off["amplitude_ratio"],
                        on["amplitude_ratio"],
                        float(on["amplitude_ratio"])
                        - float(off["amplitude_ratio"]),
                        off["ke_ratio"],
                        on["ke_ratio"],
                        float(on["ke_ratio"]) - float(off["ke_ratio"]),
                        float(on["phase_error_rad"]) - float(off["phase_error_rad"]),
                        off["offline_broad_nonsmooth_faces"],
                        on["offline_broad_nonsmooth_faces"],
                        off["offline_pressure_jump_faces"],
                        on["offline_pressure_jump_faces"],
                        off["offline_compression_faces"],
                        on["offline_compression_faces"],
                        off["offline_pressure_compression_faces"],
                        on["offline_pressure_compression_faces"],
                        off["offline_llf_fallback_faces"],
                        on["offline_llf_fallback_faces"],
                        maximum,
                    ]
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    rows = analyse(args.root)
    write_csv(args.root / "metrics.csv", rows)
    write_comparison(args.root / "comparison.csv", rows, args.root)


if __name__ == "__main__":
    main()
