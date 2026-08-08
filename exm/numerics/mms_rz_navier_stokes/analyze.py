#!/usr/bin/env python3
"""Point-sample error norms for the full R-Z Navier--Stokes MMS."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np

try:
    import yt
except ImportError as exc:
    raise SystemExit("analyze.py requires yt") from exc


GAMMA = 1.4


def latest_plotfile(path: Path) -> Path:
    if path.is_dir() and (path / "Header").is_file():
        return path
    pattern = re.compile(r"^plt(\d+)$")
    candidates = [
        item
        for item in path.glob("**/plt*")
        if pattern.fullmatch(item.name) and (item / "Header").is_file()
    ]
    if not candidates:
        raise FileNotFoundError(f"no plotfile below {path}")
    return max(
        candidates,
        key=lambda item: (int(pattern.fullmatch(item.name)[1]), str(item)),
    )


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


def field_key(ds, name: str):
    wanted = name.lower()
    for key in ds.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"field {name!r} not found in {ds.field_list}")


def exact_state(r, z, time):
    phase = 2.0 * np.pi * z - time
    rho = 1.0 + 0.050 * (1.0 + 0.20 * r**2) * np.sin(phase)
    temperature = (
        1.0
        + 0.060
        * (1.0 + 0.25 * r**2)
        * np.cos(phase - 0.21)
    )
    ur = (
        0.080
        * r
        * (1.0 + 0.15 * r**2)
        * np.cos(phase + 0.13)
    )
    uz = (
        0.30
        + 0.050
        * (1.0 + 0.10 * r**2)
        * np.sin(phase + 0.35)
    )
    pressure = rho * temperature
    return {
        "density": rho,
        "radial_momentum": rho * ur,
        "axial_momentum": rho * uz,
        "energy": (
            pressure / (GAMMA - 1.0)
            + 0.5 * rho * (ur * ur + uz * uz)
        ),
    }


def exact_time_derivative(r, z, time):
    complex_step = 1.0e-30
    complex_state = exact_state(r, z, time + 1j * complex_step)
    return {
        name: np.imag(value) / complex_step
        for name, value in complex_state.items()
    }


def weighted_metrics(error, volume, radial_index, nr):
    total = np.sum(volume)
    first_ring = radial_index == 0
    first_three = radial_index < 3
    interior = radial_index < max(nr - 3, 1)
    interior_volume = volume[interior]
    interior_error = error[interior]
    return {
        "L1": float(np.sum(np.abs(error) * volume) / total),
        "L2": float(np.sqrt(np.sum(error * error * volume) / total)),
        "Linf": float(np.max(np.abs(error))),
        "interior_L1": float(
            np.sum(np.abs(interior_error) * interior_volume)
            / np.sum(interior_volume)
        ),
        "interior_L2": float(
            np.sqrt(
                np.sum(interior_error * interior_error * interior_volume)
                / np.sum(interior_volume)
            )
        ),
        "interior_Linf": float(np.max(np.abs(interior_error))),
        "first_ring_Linf": float(np.max(np.abs(error[first_ring]))),
        "first_three_rings_Linf": float(
            np.max(np.abs(error[first_three]))
        ),
    }


def analyze_case(path: Path, quantity: str, physics: str, exact_time_override):
    plotfile = latest_plotfile(path)
    ds = load_boxlib_rz(plotfile)
    data = ds.all_data()
    r = np.asarray(data[("index", "r")].d, dtype=float)
    z = np.asarray(data[("index", "z")].d, dtype=float)
    volume = np.asarray(data[("index", "cell_volume")].d, dtype=float)
    numerical = {
        "density": np.asarray(
            data[field_key(ds, "Density")].d, dtype=float
        ),
        "radial_momentum": np.asarray(
            data[field_key(ds, "Xmom")].d, dtype=float
        ),
        "axial_momentum": np.asarray(
            data[field_key(ds, "Ymom")].d, dtype=float
        ),
        "energy": np.asarray(
            data[field_key(ds, "Energy")].d, dtype=float
        ),
    }
    nr = int(ds.domain_dimensions[0])
    nz = int(ds.domain_dimensions[1])
    radial_lo = float(ds.domain_left_edge[0])
    radial_hi = float(ds.domain_right_edge[0])
    dr = (radial_hi - radial_lo) / nr
    radial_index = np.rint((r - radial_lo) / dr - 0.5).astype(int)
    plot_time = float(ds.current_time)
    exact_time_value = (
        plot_time if exact_time_override is None else exact_time_override
    )
    exact = (
        exact_state(r, z, exact_time_value)
        if quantity == "state"
        else exact_time_derivative(r, z, exact_time_value)
    )

    row = {
        "physics": physics,
        "quantity": quantity,
        "state_semantics": "point_sample",
        "plotfile": str(plotfile),
        "base_n": nr,
        "base_nz": nz,
        "plot_time": plot_time,
        "exact_time": exact_time_value,
    }
    for variable in (
        "density",
        "radial_momentum",
        "axial_momentum",
        "energy",
    ):
        metrics = weighted_metrics(
            numerical[variable] - exact[variable],
            volume,
            radial_index,
            nr,
        )
        for metric, value in metrics.items():
            row[f"{variable}_{metric}"] = value
    return row


def add_rates(rows):
    variables = (
        "density",
        "radial_momentum",
        "axial_momentum",
        "energy",
    )
    metrics = (
        "L1",
        "L2",
        "Linf",
        "interior_L1",
        "interior_L2",
        "interior_Linf",
        "first_ring_Linf",
        "first_three_rings_Linf",
    )
    for row in rows:
        for variable in variables:
            for metric in metrics:
                row[f"{variable}_{metric}_rate"] = float("nan")
    for coarse, fine in zip(rows, rows[1:]):
        ratio = fine["base_n"] / coarse["base_n"]
        for variable in variables:
            for metric in metrics:
                key = f"{variable}_{metric}"
                fine[f"{key}_rate"] = math.log(
                    coarse[key] / fine[key]
                ) / math.log(ratio)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument(
        "--quantity", choices=("state", "rhs"), default="state"
    )
    parser.add_argument(
        "--physics",
        choices=("euler", "navier-stokes"),
        default="navier-stokes",
    )
    parser.add_argument("--exact-time", type=float)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    rows = sorted(
        (
            analyze_case(
                case, args.quantity, args.physics, args.exact_time
            )
            for case in args.cases
        ),
        key=lambda row: row["base_n"],
    )
    if any(row["base_n"] != row["base_nz"] for row in rows):
        raise ValueError("MMS convergence requires square grids")
    add_rates(rows)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{args.physics:14s} {args.quantity:5s} "
            f"N={row['base_n']:4d} "
            f"rho_L2={row['density_L2']:.6e} "
            f"mr_L2={row['radial_momentum_L2']:.6e} "
            f"mz_L2={row['axial_momentum_L2']:.6e} "
            f"E_L2={row['energy_L2']:.6e} "
            f"mr_axis={row['radial_momentum_first_ring_Linf']:.6e}"
        )
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    yt.set_log_level(40)
    main()
