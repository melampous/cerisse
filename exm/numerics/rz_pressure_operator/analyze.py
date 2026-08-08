#!/usr/bin/env python3
"""Compare the computed radial-momentum RHS with -dp/dr."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--profile", type=int, choices=(0, 1), required=True)
    parser.add_argument("--p0", type=float, default=2.0)
    parser.add_argument("--amplitude", type=float, default=0.25)
    parser.add_argument("--shock-z", type=float, default=0.5)
    parser.add_argument("--curvature", type=float, default=0.20)
    parser.add_argument("--width", type=float, default=0.04)
    parser.add_argument(
        "--point-values",
        action="store_true",
        help="compare against the analytic RHS at cell centres",
    )
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_boxlib_rz(args.plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    rhs = raw_field(grid, "Xmom")
    if rhs.shape[0] != dims[0] and rhs.shape[1] == dims[0]:
        rhs = rhs.T

    lo = np.asarray(ds.domain_left_edge, dtype=float)
    hi = np.asarray(ds.domain_right_edge, dtype=float)
    dr = (hi[0] - lo[0]) / rhs.shape[0]
    dz = (hi[1] - lo[1]) / rhs.shape[1]
    r = lo[0] + (np.arange(rhs.shape[0]) + 0.5) * dr
    z = lo[1] + (np.arange(rhs.shape[1]) + 0.5) * dz
    rr, zz = np.meshgrid(r, z, indexing="ij")

    if args.profile == 0:
        if args.point_values:
            exact = -2.0 * args.amplitude * rr
        else:
            rlo = rr - 0.5 * dr
            rhi = rr + 0.5 * dr
            exact = (
                -(4.0 / 3.0)
                * args.amplitude
                * (rhi**3 - rlo**3)
                / (rhi**2 - rlo**2)
            )
        mask = np.ones_like(exact, dtype=bool)
    else:
        if args.point_values:
            arg = (
                zz - args.shock_z - args.curvature * rr**2
            ) / args.width
            exact = (
                2.0
                * args.amplitude
                * args.curvature
                * rr
                / args.width
                / np.cosh(arg) ** 2
            )
        else:
            nodes, weights = np.polynomial.legendre.leggauss(12)
            numerator = np.zeros_like(rr)
            denominator = np.zeros_like(rr)
            for node_r, weight_r in zip(nodes, weights):
                rq = rr + 0.5 * dr * node_r
                for node_z, weight_z in zip(nodes, weights):
                    zq = zz + 0.5 * dz * node_z
                    arg = (
                        zq - args.shock_z - args.curvature * rq**2
                    ) / args.width
                    weight = weight_r * weight_z * rq
                    numerator += weight * (
                        2.0
                        * args.amplitude
                        * args.curvature
                        * rq
                        / args.width
                        / np.cosh(arg) ** 2
                    )
                    denominator += weight
            exact = numerator / denominator
        arg = (zz - args.shock_z - args.curvature * rr**2) / args.width
        mask = np.abs(arg) <= 3.0

    error = rhs - exact
    first = np.arange(rhs.shape[0])[:, None] < min(8, rhs.shape[0])
    first = np.broadcast_to(first, rhs.shape)
    weighted_l2 = np.sqrt(np.sum(rr[mask] * error[mask] ** 2) / np.sum(rr[mask]))
    metrics = {
        "plotfile": str(args.plotfile),
        "state_semantics": "point" if args.point_values else "cylindrical_volume_average",
        "profile": args.profile,
        "nr": int(rhs.shape[0]),
        "nz": int(rhs.shape[1]),
        "dr": float(dr),
        "rhs_absmax": float(np.max(np.abs(rhs))),
        "exact_absmax": float(np.max(np.abs(exact))),
        "error_linf": float(np.max(np.abs(error[mask]))),
        "error_weighted_l2": float(weighted_l2),
        "axis_error_linf": float(np.max(np.abs(error[0]))),
        "first8_error_linf": float(np.max(np.abs(error[first & mask]))),
        "axis_rhs_over_r_linf": float(np.max(np.abs(rhs[0] / r[0]))),
    }
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
