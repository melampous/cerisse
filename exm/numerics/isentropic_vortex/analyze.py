#!/usr/bin/env python3
"""Volume-weighted errors for the periodic convecting isentropic vortex."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np

try:
    import yt
except ImportError as exc:  # pragma: no cover - dependency diagnostic
    raise SystemExit("analyze.py requires yt (python3 -m pip install yt)") from exc


GAMMA = 1.4
R_UNIVERSAL = 6.02214076e23 * 1.380649e-23
MOLECULAR_WEIGHT = 28.96e-3
R_SPECIFIC = R_UNIVERSAL / MOLECULAR_WEIGHT


def latest_plotfile(path: Path) -> Path:
    if path.is_dir() and (path / "Header").is_file():
        return path
    pattern = re.compile(r"^plt(\d+)$")
    candidates = [
        p for p in path.glob("plt*") if pattern.fullmatch(p.name) and (p / "Header").is_file()
    ]
    if not candidates:
        candidates = [
            p
            for p in path.glob("**/plt*")
            if pattern.fullmatch(p.name) and (p / "Header").is_file()
        ]
    if not candidates:
        raise FileNotFoundError(f"no AMReX plotfile found below {path}")
    return max(candidates, key=lambda candidate: (int(pattern.fullmatch(candidate.name)[1]), str(candidate)))


def field_key(ds, name: str):
    wanted = name.lower()
    for key in ds.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"field {name!r} is unavailable; fields are {ds.field_list}")


def advance_time(path: Path) -> float:
    pattern = re.compile(r"Run Time advance\s*=\s*([0-9.eE+-]+)")
    for log in (path / "run.log", path.parent / "run.log"):
        if log.is_file():
            match = pattern.search(log.read_text(errors="replace"))
            if match:
                return float(match.group(1))
    return float("nan")


def periodic_delta(position, centre: float, lower: float, length: float):
    return (position - centre + 0.5 * length) % length - 0.5 * length


def exact_state(x, y, time: float, args):
    length_x = args.x_hi - args.x_lo
    length_y = args.y_hi - args.y_lo
    sound_speed = math.sqrt(GAMMA * R_SPECIFIC * args.T_inf)
    mean_speed = args.mach * sound_speed
    flow_angle = math.radians(getattr(args, "flow_angle_deg", 0.0))
    u_inf = mean_speed * math.cos(flow_angle)
    v_inf = mean_speed * math.sin(flow_angle)
    x_centre = args.x_lo + ((args.xc0 - args.x_lo + u_inf * time) % length_x)
    y_centre = args.y_lo + ((args.yc0 - args.y_lo + v_inf * time) % length_y)
    xr = periodic_delta(x, x_centre, args.x_lo, length_x) / args.radius
    yr = periodic_delta(y, y_centre, args.y_lo, length_y) / args.radius
    radius_sq = xr * xr + yr * yr
    envelope = np.exp(-0.5 * radius_sq)
    circulation_speed = mean_speed * args.beta
    u = u_inf - circulation_speed * yr * envelope
    v = v_inf + circulation_speed * xr * envelope
    temperature = args.T_inf - 0.5 * circulation_speed**2 / (
        GAMMA * R_SPECIFIC / (GAMMA - 1.0)
    ) * np.exp(-radius_sq)
    rho_inf = args.p_inf / (R_SPECIFIC * args.T_inf)
    rho = rho_inf * (temperature / args.T_inf) ** (1.0 / (GAMMA - 1.0))
    pressure = rho * R_SPECIFIC * temperature
    energy = pressure / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    return rho, u, v, pressure, energy, u_inf, v_inf


def weighted_lpnorm(error, volume, order: int) -> float:
    total_volume = np.sum(volume)
    if order == 1:
        return float(np.sum(np.abs(error) * volume) / total_volume)
    if order == 2:
        return float(np.sqrt(np.sum(error * error * volume) / total_volume))
    raise ValueError(order)


def analyse_case(path: Path, args):
    plotfile = latest_plotfile(path)
    ds = yt.load(str(plotfile))
    domain_lo = np.asarray(ds.domain_left_edge.d, dtype=float)
    domain_hi = np.asarray(ds.domain_right_edge.d, dtype=float)
    expected_lo = np.array([args.x_lo, args.y_lo])
    expected_hi = np.array([args.x_hi, args.y_hi])
    scale = max(1.0, float(np.max(np.abs(np.concatenate((expected_lo, expected_hi))))))
    if not (
        np.allclose(domain_lo[:2], expected_lo, rtol=0.0, atol=1.0e-12 * scale)
        and np.allclose(domain_hi[:2], expected_hi, rtol=0.0, atol=1.0e-12 * scale)
    ):
        raise ValueError(
            f"domain mismatch for {plotfile}: plotfile has "
            f"[{domain_lo[0]}, {domain_hi[0]}] x [{domain_lo[1]}, {domain_hi[1]}], "
            f"but the analysis expects [{args.x_lo}, {args.x_hi}] x "
            f"[{args.y_lo}, {args.y_hi}]"
        )
    data = ds.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    y = np.asarray(data[("index", "y")].d, dtype=float)
    volume = np.asarray(data[("index", "cell_volume")].d, dtype=float)
    rho = np.asarray(data[field_key(ds, "Density")].d, dtype=float)
    xmom = np.asarray(data[field_key(ds, "Xmom")].d, dtype=float)
    ymom = np.asarray(data[field_key(ds, "Ymom")].d, dtype=float)
    zmom = np.asarray(data[field_key(ds, "Zmom")].d, dtype=float)
    energy = np.asarray(data[field_key(ds, "Energy")].d, dtype=float)

    u = xmom / rho
    v = ymom / rho
    pressure = (GAMMA - 1.0) * (
        energy - 0.5 * (xmom * xmom + ymom * ymom + zmom * zmom) / rho
    )
    time = float(ds.current_time)
    rho_e, u_e, v_e, p_e, energy_e, u_inf, v_inf = exact_state(
        x, y, time, args
    )
    rho_initial, _, _, _, energy_initial, _, _ = exact_state(
        x, y, 0.0, args
    )

    velocity_error_sq = (u - u_e) ** 2 + (v - v_e) ** 2
    total_volume = np.sum(volume)
    mass = np.sum(rho * volume)
    mass_initial = np.sum(rho_initial * volume)
    total_energy = np.sum(energy * volume)
    total_energy_initial = np.sum(energy_initial * volume)
    perturbation_ke = np.sum(
        0.5 * rho * ((u - u_inf) ** 2 + (v - v_inf) ** 2) * volume
    )
    perturbation_ke_exact = np.sum(
        0.5
        * rho_e
        * ((u_e - u_inf) ** 2 + (v_e - v_inf) ** 2)
        * volume
    )

    return {
        "plotfile": str(plotfile),
        "base_n": int(ds.domain_dimensions[0]),
        "base_ny": int(ds.domain_dimensions[1]),
        "max_level": int(ds.index.max_level),
        "active_cells": int(rho.size),
        "time": time,
        "x_lo": float(domain_lo[0]),
        "x_hi": float(domain_hi[0]),
        "y_lo": float(domain_lo[1]),
        "y_hi": float(domain_hi[1]),
        "advance_seconds": advance_time(path),
        "density_L1": weighted_lpnorm(rho - rho_e, volume, 1),
        "density_L2": weighted_lpnorm(rho - rho_e, volume, 2),
        "density_Linf": float(np.max(np.abs(rho - rho_e))),
        "u_L2": weighted_lpnorm(u - u_e, volume, 2),
        "v_L2": weighted_lpnorm(v - v_e, volume, 2),
        "velocity_vector_L2": float(
            np.sqrt(np.sum(velocity_error_sq * volume) / total_volume)
        ),
        "pressure_L2": weighted_lpnorm(pressure - p_e, volume, 2),
        "mass_relative_change_from_initial": float((mass - mass_initial) / mass_initial),
        "energy_relative_change_from_initial": float(
            (total_energy - total_energy_initial) / total_energy_initial
        ),
        "perturbation_ke_ratio": float(perturbation_ke / perturbation_ke_exact),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument("--csv", type=Path, default=Path("vortex_errors.csv"))
    parser.add_argument("--mach", type=float, default=0.5)
    parser.add_argument("--beta", type=float, default=0.2)
    parser.add_argument("--p-inf", dest="p_inf", type=float, default=1.0e5)
    parser.add_argument("--T-inf", dest="T_inf", type=float, default=300.0)
    parser.add_argument("--radius", type=float, default=0.005)
    parser.add_argument("--flow-angle-deg", type=float, default=0.0)
    parser.add_argument("--xc0", type=float, default=0.05)
    parser.add_argument("--yc0", type=float, default=0.05)
    parser.add_argument("--x-lo", type=float, default=0.0)
    parser.add_argument("--x-hi", type=float, default=0.1)
    parser.add_argument("--y-lo", type=float, default=0.0)
    parser.add_argument("--y-hi", type=float, default=0.1)
    args = parser.parse_args()

    rows = [analyse_case(path, args) for path in args.cases]
    uniform_rows = sorted(
        (row for row in rows if row["max_level"] == 0), key=lambda row: row["base_n"]
    )
    for row in uniform_rows:
        if row["base_n"] != row["base_ny"]:
            raise ValueError(
                f"uniform convergence sequence requires square grids; got "
                f"{row['base_n']} x {row['base_ny']} in {row['plotfile']}"
            )
    if uniform_rows:
        baseline = uniform_rows[0]
        for row in uniform_rows[1:]:
            for key in ("time", "x_lo", "x_hi", "y_lo", "y_hi"):
                scale = max(1.0, abs(baseline[key]))
                if not math.isclose(row[key], baseline[key], rel_tol=0.0, abs_tol=1.0e-12 * scale):
                    raise ValueError(
                        f"uniform convergence cases must have matching {key}: "
                        f"{baseline['plotfile']} and {row['plotfile']} differ"
                    )
        resolutions = [row["base_n"] for row in uniform_rows]
        if len(set(resolutions)) != len(resolutions):
            raise ValueError(f"uniform convergence sequence has duplicate grids: {resolutions}")
    rates = {id(row): float("nan") for row in rows}
    for coarse, fine in zip(uniform_rows, uniform_rows[1:]):
        rates[id(fine)] = math.log(coarse["density_L2"] / fine["density_L2"]) / math.log(
            fine["base_n"] / coarse["base_n"]
        )
    for row in rows:
        row["density_L2_rate"] = rates[id(row)]

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"N0={row['base_n']:4d} L={row['max_level']}  "
            f"rho_L2={row['density_L2']:.6e}  "
            f"u_L2={row['u_L2']:.6e}  v_L2={row['v_L2']:.6e}  "
            f"Kprime/Kprime_exact={row['perturbation_ke_ratio']:.8f}  "
            f"advance={row['advance_seconds']:.6g} s"
        )
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    yt.set_log_level(40)
    main()
