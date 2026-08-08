#!/usr/bin/env python3
"""Measure full R-Z MMS errors using cylindrical finite-volume weights."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yt


FIELD_NAMES = {
    "radial_momentum": "Xmom",
    "axial_momentum": "Ymom",
    "azimuthal_momentum": "Zmom",
    "energy": "Energy",
    "density": "Density",
}


@dataclass(frozen=True)
class Parameters:
    gamma: float = 1.4
    rho0: float = 1.0
    p0: float = 1.0
    uz0: float = 0.30
    rho_amp: float = 0.050
    ur_amp: float = 0.080
    uz_amp: float = 0.050
    p_amp: float = 0.080
    rho_r2: float = 0.20
    ur_r2: float = 0.15
    uz_r2: float = 0.10
    p_r2: float = 0.25
    omega: float = 1.0
    uz_phase: float = 0.35
    p_phase: float = 0.21


def load_boxlib_rz(plotfile: Path):
    from yt.data_objects.static_output import MutableAttribute

    original_set = MutableAttribute.__set__

    def writable_set(self, instance, value):
        self.data[instance] = value

    MutableAttribute.__set__ = writable_set
    try:
        return yt.load(str(plotfile))
    finally:
        MutableAttribute.__set__ = original_set


def raw_field(grid, name: str) -> np.ndarray:
    for namespace in ("boxlib", "gas", "stream"):
        try:
            return np.squeeze(np.asarray(grid[(namespace, name)]))
        except Exception:
            pass
    raise KeyError(name)


def point_values(
    r: np.ndarray, z: np.ndarray, time: float, p: Parameters
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    theta = 2.0 * np.pi * z - p.omega * time
    rho_shape = p.rho_amp * (1.0 + p.rho_r2 * r**2)
    ur_shape = p.ur_amp * r * (1.0 + p.ur_r2 * r**2)
    uz_shape = p.uz_amp * (1.0 + p.uz_r2 * r**2)
    pressure_shape = p.p_amp * (1.0 + p.p_r2 * r**2)

    rho = p.rho0 + rho_shape * np.sin(theta)
    ur = ur_shape * np.cos(theta)
    uz = p.uz0 + uz_shape * np.sin(theta + p.uz_phase)
    pressure = p.p0 + pressure_shape * np.cos(theta - p.p_phase)

    rho_t = -p.omega * rho_shape * np.cos(theta)
    ur_t = p.omega * ur_shape * np.sin(theta)
    uz_t = -p.omega * uz_shape * np.cos(theta + p.uz_phase)
    pressure_t = p.omega * pressure_shape * np.sin(theta - p.p_phase)

    velocity_squared = ur**2 + uz**2
    energy = pressure / (p.gamma - 1.0) + 0.5 * rho * velocity_squared
    energy_t = (
        pressure_t / (p.gamma - 1.0)
        + 0.5 * rho_t * velocity_squared
        + rho * (ur * ur_t + uz * uz_t)
    )

    state = {
        "density": rho,
        "radial_momentum": rho * ur,
        "axial_momentum": rho * uz,
        "azimuthal_momentum": np.zeros_like(rho),
        "energy": energy,
    }
    temporal_derivative = {
        "density": rho_t,
        "radial_momentum": rho_t * ur + rho * ur_t,
        "axial_momentum": rho_t * uz + rho * uz_t,
        "azimuthal_momentum": np.zeros_like(rho),
        "energy": energy_t,
    }
    return state, temporal_derivative


def cylindrical_cell_averages(
    r_center: np.ndarray,
    z_center: np.ndarray,
    dr: float,
    dz: float,
    time: float,
    quantity: str,
    parameters: Parameters,
) -> dict[str, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(4)
    result = {
        name: np.zeros((r_center.size, z_center.size), dtype=float)
        for name in FIELD_NAMES
    }
    normalization = np.zeros_like(result["density"])

    for radial_node, radial_weight in zip(nodes, weights):
        r = r_center[:, None] + 0.5 * dr * radial_node
        for axial_node, axial_weight in zip(nodes, weights):
            z = z_center[None, :] + 0.5 * dz * axial_node
            state, rhs = point_values(r, z, time, parameters)
            point = state if quantity == "state" else rhs
            weight = radial_weight * axial_weight * r
            normalization += weight
            for name in result:
                result[name] += weight * point[name]

    for name in result:
        result[name] /= normalization
    return result


def norm_metrics(error: np.ndarray, volume_weight: np.ndarray) -> dict[str, float]:
    denominator = float(np.sum(volume_weight))
    return {
        "l1": float(np.sum(np.abs(error) * volume_weight) / denominator),
        "l2": float(
            np.sqrt(np.sum(error**2 * volume_weight) / denominator)
        ),
        "linf": float(np.max(np.abs(error))),
        "axis_linf": float(np.max(np.abs(error[0, :]))),
        "first_three_rings_linf": float(
            np.max(np.abs(error[: min(3, error.shape[0]), :]))
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--quantity", choices=("state", "rhs"), default="state")
    parser.add_argument(
        "--exact-time",
        type=float,
        help="override plotfile time (required for order_rk=0 RHS output)",
    )
    parser.add_argument("--uz0", type=float, default=0.30)
    parser.add_argument("--rho-amp", type=float, default=0.050)
    parser.add_argument("--ur-amp", type=float, default=0.080)
    parser.add_argument("--uz-amp", type=float, default=0.050)
    parser.add_argument("--p-amp", type=float, default=0.080)
    parser.add_argument("--rho-r2", type=float, default=0.20)
    parser.add_argument("--ur-r2", type=float, default=0.15)
    parser.add_argument("--uz-r2", type=float, default=0.10)
    parser.add_argument("--p-r2", type=float, default=0.25)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_boxlib_rz(args.plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    numerical = {
        name: raw_field(grid, plot_name)
        for name, plot_name in FIELD_NAMES.items()
    }

    reference_shape = numerical["density"].shape
    if reference_shape != tuple(dims[:2]):
        transposed_shape = tuple(dims[:2][::-1])
        if reference_shape != transposed_shape:
            raise RuntimeError(
                f"unexpected field shape {reference_shape}; domain is {tuple(dims[:2])}"
            )
        numerical = {name: values.T for name, values in numerical.items()}

    lo = np.asarray(ds.domain_left_edge, dtype=float)
    hi = np.asarray(ds.domain_right_edge, dtype=float)
    nr, nz = numerical["density"].shape
    dr = (hi[0] - lo[0]) / nr
    dz = (hi[1] - lo[1]) / nz
    r_center = lo[0] + (np.arange(nr) + 0.5) * dr
    z_center = lo[1] + (np.arange(nz) + 0.5) * dz
    plot_time = float(ds.current_time)
    exact_time = plot_time if args.exact_time is None else args.exact_time

    parameters = Parameters(
        uz0=args.uz0,
        rho_amp=args.rho_amp,
        ur_amp=args.ur_amp,
        uz_amp=args.uz_amp,
        p_amp=args.p_amp,
        rho_r2=args.rho_r2,
        ur_r2=args.ur_r2,
        uz_r2=args.uz_r2,
        p_r2=args.p_r2,
    )
    exact = cylindrical_cell_averages(
        r_center, z_center, dr, dz, exact_time, args.quantity, parameters
    )
    r_lo = r_center - 0.5 * dr
    r_hi = r_center + 0.5 * dr
    volume_weight = np.broadcast_to(
        (r_hi**2 - r_lo**2)[:, None], (nr, nz)
    )

    errors = {
        name: norm_metrics(numerical[name] - exact[name], volume_weight)
        for name in FIELD_NAMES
    }
    report: dict[str, object] = {
        "schema_version": 1,
        "plotfile": str(args.plotfile),
        "quantity": args.quantity,
        "state_semantics": "cylindrical_volume_average",
        "quadrature": "tensor Gauss-Legendre order 4 with radial Jacobian",
        "nr": nr,
        "nz": nz,
        "dr": float(dr),
        "dz": float(dz),
        "plot_time": plot_time,
        "exact_time": exact_time,
        "errors": errors,
    }

    if args.quantity == "state":
        density = numerical["density"]
        kinetic = 0.5 * (
            numerical["radial_momentum"] ** 2
            + numerical["axial_momentum"] ** 2
            + numerical["azimuthal_momentum"] ** 2
        ) / density
        pressure = (parameters.gamma - 1.0) * (
            numerical["energy"] - kinetic
        )
        report["admissibility"] = {
            "minimum_density": float(np.min(density)),
            "minimum_cell_pressure": float(np.min(pressure)),
        }

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
