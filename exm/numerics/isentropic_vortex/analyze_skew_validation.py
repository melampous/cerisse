#!/usr/bin/env python3
"""Uniform-grid error and smooth-vortex diagnostics for the Skew studies.

The error norms use the analytic convecting-vortex state at the plotfile time.
Vorticity is evaluated with the same periodic fourth-order centred derivative
for the numerical and analytic velocity fields.  The resulting ratios measure
loss relative to a grid-matched reference and do not include a derivative-grid
bias in the denominator.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

import analyze as base


def periodic_derivative(values: np.ndarray, spacing: float, axis: int) -> np.ndarray:
    """Return the fourth-order centred periodic derivative."""
    return (
        np.roll(values, 2, axis=axis)
        - 8.0 * np.roll(values, 1, axis=axis)
        + 8.0 * np.roll(values, -1, axis=axis)
        - np.roll(values, -2, axis=axis)
    ) / (12.0 * spacing)


def lpnorm(error: np.ndarray, order: int) -> float:
    if order == 1:
        return float(np.mean(np.abs(error)))
    if order == 2:
        return float(np.sqrt(np.mean(error * error)))
    if order == 0:
        return float(np.max(np.abs(error)))
    raise ValueError(order)


def structured_field(dataset, grid, name: str) -> np.ndarray:
    values = np.asarray(grid[base.field_key(dataset, name)].d, dtype=float)
    if values.ndim == 3 and values.shape[2] == 1:
        values = values[:, :, 0]
    if values.ndim != 2:
        raise ValueError(f"{name} is not a two-dimensional uniform field: {values.shape}")
    return values


def analyse_uniform_case(path: Path, args) -> dict[str, float | int | str]:
    plotfile = base.latest_plotfile(path)
    dataset = base.yt.load(str(plotfile))
    if int(dataset.index.max_level) != 0:
        raise ValueError(f"Skew validation requires a uniform grid: {plotfile}")

    dimensions = np.asarray(dataset.domain_dimensions, dtype=int)
    nx, ny = int(dimensions[0]), int(dimensions[1])
    if nx != ny:
        raise ValueError(f"Skew validation requires a square grid: {nx} x {ny}")

    lower = np.asarray(dataset.domain_left_edge.d, dtype=float)
    upper = np.asarray(dataset.domain_right_edge.d, dtype=float)
    grid = dataset.covering_grid(
        level=0, left_edge=dataset.domain_left_edge, dims=dataset.domain_dimensions
    )
    rho = structured_field(dataset, grid, "Density")
    xmom = structured_field(dataset, grid, "Xmom")
    ymom = structured_field(dataset, grid, "Ymom")
    zmom = structured_field(dataset, grid, "Zmom")
    energy = structured_field(dataset, grid, "Energy")
    x = np.asarray(grid[("index", "x")].d, dtype=float)[:, :, 0]
    y = np.asarray(grid[("index", "y")].d, dtype=float)[:, :, 0]

    u = xmom / rho
    v = ymom / rho
    pressure = (base.GAMMA - 1.0) * (
        energy - 0.5 * (xmom * xmom + ymom * ymom + zmom * zmom) / rho
    )
    time = float(dataset.current_time)
    rho_e, u_e, v_e, p_e, energy_e, u_inf, v_inf = base.exact_state(
        x, y, time, args
    )
    rho_initial, _, _, _, energy_initial, _, _ = base.exact_state(x, y, 0.0, args)

    dx = (upper[0] - lower[0]) / nx
    dy = (upper[1] - lower[1]) / ny
    omega = periodic_derivative(v, dx, axis=0) - periodic_derivative(
        u, dy, axis=1
    )
    omega_e = periodic_derivative(v_e, dx, axis=0) - periodic_derivative(
        u_e, dy, axis=1
    )

    length_x = upper[0] - lower[0]
    length_y = upper[1] - lower[1]
    x_centre = lower[0] + ((args.xc0 - lower[0] + u_inf * time) % length_x)
    y_centre = lower[1] + ((args.yc0 - lower[1] + v_inf * time) % length_y)
    radius = np.sqrt(
        base.periodic_delta(x, x_centre, lower[0], length_x) ** 2
        + base.periodic_delta(y, y_centre, lower[1], length_y) ** 2
    )
    core = radius <= math.sqrt(2.0) * args.radius

    perturbation_ke = np.sum(
        0.5 * rho * ((u - u_inf) ** 2 + (v - v_inf) ** 2)
    )
    perturbation_ke_exact = np.sum(
        0.5 * rho_e * ((u_e - u_inf) ** 2 + (v_e - v_inf) ** 2)
    )
    perturbation_ke_ratio = float(perturbation_ke / perturbation_ke_exact)
    ke_decay_rate = float(-math.log(perturbation_ke_ratio) / time)
    # For the Gaussian vortex used here, the continuous constant-density
    # relation is Z/K = 2/R^2.  Combining it with dK/dt = -2 nu Z gives this
    # scale-specific, time-integrated equivalent viscosity.  It is a compact
    # dissipation indicator, not a universal viscosity of the numerical method.
    equivalent_numerical_viscosity = float(
        -args.radius * args.radius
        * math.log(perturbation_ke_ratio)
        / (4.0 * time)
    )
    enstrophy = 0.5 * np.sum(omega * omega)
    enstrophy_exact = 0.5 * np.sum(omega_e * omega_e)
    core_circulation = np.sum(omega[core])
    core_circulation_exact = np.sum(omega_e[core])

    errors = {
        "density": rho - rho_e,
        "u": u - u_e,
        "v": v - v_e,
        "pressure": pressure - p_e,
    }
    row: dict[str, float | int | str] = {
        "plotfile": str(plotfile),
        "base_n": nx,
        "time": time,
        "advance_seconds": base.advance_time(path),
        "minimum_density": float(np.min(rho)),
        "minimum_pressure": float(np.min(pressure)),
        "mass_relative_change_from_initial": float(
            (np.sum(rho) - np.sum(rho_initial)) / np.sum(rho_initial)
        ),
        "energy_relative_change_from_initial": float(
            (np.sum(energy) - np.sum(energy_initial)) / np.sum(energy_initial)
        ),
        "perturbation_ke_ratio": perturbation_ke_ratio,
        "ke_decay_rate_per_second": ke_decay_rate,
        "equivalent_numerical_viscosity_from_ke": equivalent_numerical_viscosity,
        "peak_vorticity_ratio": float(np.max(omega) / np.max(omega_e)),
        "enstrophy_ratio": float(enstrophy / enstrophy_exact),
        "core_circulation_ratio": float(
            core_circulation / core_circulation_exact
        ),
    }
    for variable, error in errors.items():
        row[f"{variable}_L1"] = lpnorm(error, 1)
        row[f"{variable}_L2"] = lpnorm(error, 2)
        row[f"{variable}_Linf"] = lpnorm(error, 0)
    return row


def add_observed_orders(rows: list[dict[str, float | int | str]]) -> None:
    error_keys = [
        f"{variable}_{norm}"
        for variable in ("density", "u", "v", "pressure")
        for norm in ("L1", "L2", "Linf")
    ]
    for row in rows:
        for key in error_keys:
            row[f"{key}_order"] = float("nan")
    ordered = sorted(rows, key=lambda row: int(row["base_n"]))
    for coarse, fine in zip(ordered, ordered[1:]):
        refinement = float(fine["base_n"]) / float(coarse["base_n"])
        for key in error_keys:
            fine[f"{key}_order"] = math.log(
                float(coarse[key]) / float(fine[key])
            ) / math.log(refinement)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument("--csv", type=Path, required=True)
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

    rows = [analyse_uniform_case(path, args) for path in args.cases]
    reference = rows[0]
    for row in rows[1:]:
        if not math.isclose(
            float(row["time"]),
            float(reference["time"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise ValueError("all cases in one convergence analysis need the same time")
    add_observed_orders(rows)
    rows.sort(key=lambda row: int(row["base_n"]))

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"N={int(row['base_n']):4d} "
            f"rho_L2={float(row['density_L2']):.6e} "
            f"order={float(row['density_L2_order']):.4f} "
            f"K'/K'e={float(row['perturbation_ke_ratio']):.8f} "
            f"omega_max/omega_e={float(row['peak_vorticity_ratio']):.8f} "
            f"Z/Z_e={float(row['enstrophy_ratio']):.8f} "
            f"advance={float(row['advance_seconds']):.6g} s"
        )
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    base.yt.set_log_level(40)
    main()
