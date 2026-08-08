#!/usr/bin/env python3
"""Compare smooth radial-flow Euler RHS fields with analytic values."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yt


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


def cylindrical_moment(r: np.ndarray, dr: float, power: int) -> np.ndarray:
    rlo = r - 0.5 * dr
    rhi = r + 0.5 * dr
    return (
        2.0
        / (power + 2.0)
        * (rhi ** (power + 2) - rlo ** (power + 2))
        / (rhi**2 - rlo**2)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--mode", type=int, choices=(0, 1), required=True)
    parser.add_argument("--rho", type=float, default=1.0)
    parser.add_argument("--a", type=float, default=0.20)
    parser.add_argument("--b", type=float, default=0.10)
    parser.add_argument("--point-values", action="store_true")
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_boxlib_rz(args.plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    density_rhs = raw_field(grid, "Density")
    radial_momentum_rhs = raw_field(grid, "Xmom")
    if density_rhs.shape[0] != dims[0] and density_rhs.shape[1] == dims[0]:
        density_rhs = density_rhs.T
        radial_momentum_rhs = radial_momentum_rhs.T

    lo = np.asarray(ds.domain_left_edge, dtype=float)
    hi = np.asarray(ds.domain_right_edge, dtype=float)
    dr = (hi[0] - lo[0]) / density_rhs.shape[0]
    r = lo[0] + (np.arange(density_rhs.shape[0]) + 0.5) * dr
    rr = np.broadcast_to(r[:, None], density_rhs.shape)

    if args.point_values:
        m1 = rr
        m2 = rr**2
        m3 = rr**3
        m5 = rr**5
    else:
        m1 = cylindrical_moment(rr, dr, 1)
        m2 = cylindrical_moment(rr, dr, 2)
        m3 = cylindrical_moment(rr, dr, 3)
        m5 = cylindrical_moment(rr, dr, 5)

    density_exact = -args.rho * (2.0 * args.a)
    momentum_exact = -args.rho * 3.0 * args.a**2 * m1
    if args.mode == 1:
        density_exact = -args.rho * (2.0 * args.a + 4.0 * args.b * m2)
        momentum_exact = -args.rho * (
            3.0 * args.a**2 * m1
            + 10.0 * args.a * args.b * m3
            + 7.0 * args.b**2 * m5
        )

    density_error = density_rhs - density_exact
    momentum_error = radial_momentum_rhs - momentum_exact
    metrics = {
        "plotfile": str(args.plotfile),
        "state_semantics": (
            "point" if args.point_values else "cylindrical_volume_average"
        ),
        "mode": args.mode,
        "nr": int(density_rhs.shape[0]),
        "nz": int(density_rhs.shape[1]),
        "dr": float(dr),
        "density_error_linf": float(np.max(np.abs(density_error))),
        "density_axis_error_linf": float(np.max(np.abs(density_error[0]))),
        "momentum_error_linf": float(np.max(np.abs(momentum_error))),
        "momentum_axis_error_linf": float(
            np.max(np.abs(momentum_error[0]))
        ),
        "momentum_axis_rhs_over_r_linf": float(
            np.max(np.abs(radial_momentum_rhs[0] / r[0]))
        ),
    }
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
