#!/usr/bin/env python3
"""Point-sample error norms for the Cartesian Euler and NS MMS."""

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


def field_key(ds, name: str):
    wanted = name.lower()
    for key in ds.field_list:
        if key[1].lower() == wanted:
            return key
    raise KeyError(f"field {name!r} not found in {ds.field_list}")


def exact_state(x, y, time):
    twopi = 2.0 * np.pi
    rho = (
        1.0
        + 0.08 * np.sin(twopi * x - 0.70 * time)
        + 0.05 * np.cos(twopi * y + 0.30 * time)
    )
    u = (
        0.45
        + 0.07 * np.sin(twopi * x + 0.20 * time)
        + 0.06 * np.cos(twopi * y - 0.40 * time)
    )
    v = (
        0.30
        + 0.05 * np.cos(twopi * x - 0.60 * time)
        + 0.04 * np.sin(twopi * y + 0.50 * time)
    )
    pressure = (
        1.0
        + 0.07 * np.cos(twopi * x + 0.30 * time)
        + 0.06 * np.sin(twopi * y - 0.20 * time)
    )
    xmom = rho * u
    ymom = rho * v
    energy = pressure / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    return {
        "density": rho,
        "x_momentum": xmom,
        "y_momentum": ymom,
        "energy": energy,
    }


def norm(error, volume, order):
    total = np.sum(volume)
    if order == 1:
        return float(np.sum(np.abs(error) * volume) / total)
    if order == 2:
        return float(np.sqrt(np.sum(error * error * volume) / total))
    if order == np.inf:
        return float(np.max(np.abs(error)))
    raise ValueError(order)


def analyze_case(path: Path, physics: str):
    plotfile = latest_plotfile(path)
    ds = yt.load(str(plotfile))
    data = ds.all_data()
    x = np.asarray(data[("index", "x")].d, dtype=float)
    y = np.asarray(data[("index", "y")].d, dtype=float)
    volume = np.asarray(data[("index", "cell_volume")].d, dtype=float)
    numerical = {
        "density": np.asarray(data[field_key(ds, "Density")].d, dtype=float),
        "x_momentum": np.asarray(data[field_key(ds, "Xmom")].d, dtype=float),
        "y_momentum": np.asarray(data[field_key(ds, "Ymom")].d, dtype=float),
        "energy": np.asarray(data[field_key(ds, "Energy")].d, dtype=float),
    }
    time = float(ds.current_time)
    exact = exact_state(x, y, time)
    row = {
        "physics": physics,
        "plotfile": str(plotfile),
        "base_n": int(ds.domain_dimensions[0]),
        "base_ny": int(ds.domain_dimensions[1]),
        "time": time,
    }
    for variable in ("density", "x_momentum", "y_momentum", "energy"):
        error = numerical[variable] - exact[variable]
        row[f"{variable}_L1"] = norm(error, volume, 1)
        row[f"{variable}_L2"] = norm(error, volume, 2)
        row[f"{variable}_Linf"] = norm(error, volume, np.inf)
    return row


def add_rates(rows):
    variables = ("density", "x_momentum", "y_momentum", "energy")
    norms = ("L1", "L2", "Linf")
    for row in rows:
        for variable in variables:
            for norm_name in norms:
                row[f"{variable}_{norm_name}_rate"] = float("nan")
    for coarse, fine in zip(rows, rows[1:]):
        ratio = fine["base_n"] / coarse["base_n"]
        for variable in variables:
            for norm_name in norms:
                key = f"{variable}_{norm_name}"
                fine[f"{key}_rate"] = math.log(coarse[key] / fine[key]) / math.log(
                    ratio
                )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="+", type=Path)
    parser.add_argument("--physics", choices=("euler", "navier-stokes"), required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    rows = sorted(
        (analyze_case(case, args.physics) for case in args.cases),
        key=lambda row: row["base_n"],
    )
    if any(row["base_n"] != row["base_ny"] for row in rows):
        raise ValueError("MMS convergence requires square grids")
    reference_time = rows[0]["time"]
    if any(
        not math.isclose(
            row["time"], reference_time, rel_tol=0.0, abs_tol=1.0e-12
        )
        for row in rows
    ):
        raise ValueError("MMS cases have different final times")
    add_rates(rows)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{args.physics:14s} N={row['base_n']:4d} "
            f"rho_L2={row['density_L2']:.6e} "
            f"mx_L2={row['x_momentum_L2']:.6e} "
            f"my_L2={row['y_momentum_L2']:.6e} "
            f"E_L2={row['energy_L2']:.6e}"
        )
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    yt.set_log_level(40)
    main()
