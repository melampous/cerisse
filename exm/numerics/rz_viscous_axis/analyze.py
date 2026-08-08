#!/usr/bin/env python3
"""Compare manufactured R-Z viscous RHS fields with analytic values."""

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
    parser.add_argument("--mode", type=int, choices=(0, 1, 2, 3), required=True)
    parser.add_argument("--mu", type=float, default=0.01)
    parser.add_argument("--b", type=float, default=0.10)
    parser.add_argument("--c", type=float, default=0.20)
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
    xrhs = raw_field(grid, "Xmom")
    yrhs = raw_field(grid, "Ymom")
    if xrhs.shape[0] != dims[0] and xrhs.shape[1] == dims[0]:
        xrhs = xrhs.T
        yrhs = yrhs.T

    lo = np.asarray(ds.domain_left_edge, dtype=float)
    hi = np.asarray(ds.domain_right_edge, dtype=float)
    dr = (hi[0] - lo[0]) / xrhs.shape[0]
    r = lo[0] + (np.arange(xrhs.shape[0]) + 0.5) * dr
    rr = np.broadcast_to(r[:, None], xrhs.shape)

    xexact = np.zeros_like(xrhs)
    yexact = np.zeros_like(yrhs)
    if args.mode == 1:
        if args.point_values:
            r_target = rr
        else:
            rlo = rr - 0.5 * dr
            rhi = rr + 0.5 * dr
            r_target = (
                (2.0 / 3.0)
                * (rhi**3 - rlo**3)
                / (rhi**2 - rlo**2)
            )
        xexact = (32.0 / 3.0) * args.mu * args.b * r_target
    if args.mode == 2:
        yexact.fill(4.0 * args.mu * args.c)

    xerror = xrhs - xexact
    yerror = yrhs - yexact
    metrics = {
        "plotfile": str(args.plotfile),
        "state_semantics": "point" if args.point_values else "cylindrical_volume_average",
        "mode": args.mode,
        "nr": int(xrhs.shape[0]),
        "nz": int(xrhs.shape[1]),
        "dr": float(dr),
        "x_rhs_absmax": float(np.max(np.abs(xrhs))),
        "y_rhs_absmax": float(np.max(np.abs(yrhs))),
        "x_error_linf": float(np.max(np.abs(xerror))),
        "y_error_linf": float(np.max(np.abs(yerror))),
        "x_axis_error_linf": float(np.max(np.abs(xerror[0]))),
        "y_axis_error_linf": float(np.max(np.abs(yerror[0]))),
    }
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
