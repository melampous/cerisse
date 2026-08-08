#!/usr/bin/env python3
"""Isolate the viscous R-Z truncation error from paired RHS plotfiles."""

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
    raise SystemExit("analyze_operator_difference.py requires yt") from exc


FIELD_NAMES = {
    "density": "Density",
    "radial_momentum": "Xmom",
    "axial_momentum": "Ymom",
    "energy": "Energy",
}


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


def raw_field(grid, name: str) -> np.ndarray:
    for namespace in ("boxlib", "gas", "stream"):
        try:
            return np.squeeze(np.asarray(grid[(namespace, name)]))
        except Exception:
            pass
    raise KeyError(name)


def load_case(path: Path):
    plotfile = latest_plotfile(path)
    ds = load_boxlib_rz(plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    fields = {
        name: raw_field(grid, plot_name)
        for name, plot_name in FIELD_NAMES.items()
    }
    shape = fields["density"].shape
    if shape != tuple(dims[:2]):
        if shape != tuple(dims[:2][::-1]):
            raise RuntimeError(
                f"unexpected field shape {shape}; domain is {tuple(dims[:2])}"
            )
        fields = {name: value.T for name, value in fields.items()}
    return plotfile, ds, fields


def metrics(error: np.ndarray):
    nr, nz = error.shape
    radial_index = np.arange(nr)[:, None]
    radial_weight = 2.0 * radial_index + 1.0
    volume = np.broadcast_to(radial_weight, (nr, nz))
    interior = np.broadcast_to(radial_index < max(nr - 3, 1), (nr, nz))
    axis = np.broadcast_to(radial_index == 0, (nr, nz))
    first_three = np.broadcast_to(radial_index < 3, (nr, nz))

    def l1(mask):
        return float(
            np.sum(np.abs(error[mask]) * volume[mask])
            / np.sum(volume[mask])
        )

    def l2(mask):
        return float(
            np.sqrt(
                np.sum(error[mask] ** 2 * volume[mask])
                / np.sum(volume[mask])
            )
        )

    full = np.ones_like(error, dtype=bool)
    return {
        "L1": l1(full),
        "L2": l2(full),
        "Linf": float(np.max(np.abs(error))),
        "interior_L1": l1(interior),
        "interior_L2": l2(interior),
        "interior_Linf": float(np.max(np.abs(error[interior]))),
        "first_ring_Linf": float(np.max(np.abs(error[axis]))),
        "first_three_rings_Linf": float(
            np.max(np.abs(error[first_three]))
        ),
    }


def analyze_pair(euler_path: Path, ns_path: Path):
    euler_plot, euler_ds, euler = load_case(euler_path)
    ns_plot, ns_ds, navier_stokes = load_case(ns_path)
    euler_dims = tuple(np.asarray(euler_ds.domain_dimensions, dtype=int)[:2])
    ns_dims = tuple(np.asarray(ns_ds.domain_dimensions, dtype=int)[:2])
    if euler_dims != ns_dims:
        raise ValueError(f"grid mismatch: {euler_dims} versus {ns_dims}")
    if not math.isclose(
        float(euler_ds.current_time),
        float(ns_ds.current_time),
        rel_tol=0.0,
        abs_tol=1.0e-14,
    ):
        raise ValueError("paired RHS plotfiles have different times")

    row = {
        "component": "viscous_truncation_error",
        "base_n": euler_dims[0],
        "base_nz": euler_dims[1],
        "euler_plotfile": str(euler_plot),
        "navier_stokes_plotfile": str(ns_plot),
    }
    for variable in FIELD_NAMES:
        # Both forced operators have exact RHS dU/dt. Their difference is the
        # discrete viscous operator minus its exact continuous counterpart.
        difference = navier_stokes[variable] - euler[variable]
        for metric, value in metrics(difference).items():
            row[f"{variable}_{metric}"] = value
    return row


def add_rates(rows):
    metric_suffixes = (
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
        for variable in FIELD_NAMES:
            for metric in metric_suffixes:
                row[f"{variable}_{metric}_rate"] = float("nan")
    for coarse, fine in zip(rows, rows[1:]):
        ratio = fine["base_n"] / coarse["base_n"]
        for variable in FIELD_NAMES:
            for metric in metric_suffixes:
                key = f"{variable}_{metric}"
                if coarse[key] == 0.0 or fine[key] == 0.0:
                    continue
                fine[f"{key}_rate"] = math.log(
                    coarse[key] / fine[key]
                ) / math.log(ratio)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--euler", nargs="+", type=Path, required=True)
    parser.add_argument(
        "--navier-stokes", nargs="+", type=Path, required=True
    )
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()
    if len(args.euler) != len(args.navier_stokes):
        raise ValueError("Euler and Navier-Stokes case counts must match")

    rows = sorted(
        (
            analyze_pair(euler, ns)
            for euler, ns in zip(args.euler, args.navier_stokes)
        ),
        key=lambda row: row["base_n"],
    )
    add_rates(rows)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"viscous N={row['base_n']:4d} "
            f"rho_Linf={row['density_Linf']:.3e} "
            f"mr_L2={row['radial_momentum_L2']:.6e} "
            f"mz_L2={row['axial_momentum_L2']:.6e} "
            f"E_L2={row['energy_L2']:.6e} "
            f"mr_axis={row['radial_momentum_first_ring_Linf']:.6e}"
        )
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    yt.set_log_level(40)
    main()
