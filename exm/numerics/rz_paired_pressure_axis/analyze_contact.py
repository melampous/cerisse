#!/usr/bin/env python3
"""Measure pressure/velocity pollution in the R-Z constant-p contact gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yt


def load_plotfile(path: Path):
    # Work around the MutableAttribute incompatibility in the yt version used
    # by the local Cerisse environment (same workaround as analyze.py).
    from yt.data_objects.static_output import MutableAttribute

    original_set = MutableAttribute.__set__

    def writable_set(self, instance, value):
        self.data[instance] = value

    MutableAttribute.__set__ = writable_set
    try:
        return yt.load(str(path))
    finally:
        MutableAttribute.__set__ = original_set


def field(grid, *names: str) -> np.ndarray:
    for name in names:
        for namespace in ("boxlib", "gas", "stream"):
            try:
                return np.squeeze(np.asarray(grid[(namespace, name)],
                                             dtype=float))
            except Exception:
                pass
    raise KeyError(f"none of {names!r} found")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--p0", type=float, default=2.0)
    parser.add_argument("--rho-inner", type=float, default=4.0)
    parser.add_argument("--rho-outer", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_plotfile(args.plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    pressure = field(grid, "pressure")
    density = field(grid, "density", "Density")
    ur = field(grid, "x_velocity", "velocity_x", "radial_velocity")
    uz = field(grid, "y_velocity", "velocity_y", "axial_velocity")

    sound_ref = np.sqrt(args.gamma * args.p0 /
                        min(args.rho_inner, args.rho_outer))
    dr = float(ds.domain_width[0]) / int(ds.domain_dimensions[0])
    near_axis = np.zeros_like(pressure, dtype=bool)
    near_axis[: min(4, pressure.shape[0]), :] = True
    dp = pressure - args.p0
    metrics = {
        "plotfile": str(args.plotfile.resolve()),
        "time": float(ds.current_time),
        "cells": int(pressure.size),
        "pressure_linf_abs": float(np.max(np.abs(dp))),
        "pressure_linf_rel": float(np.max(np.abs(dp)) / args.p0),
        "pressure_span_rel": float((np.max(pressure) - np.min(pressure)) /
                                   args.p0),
        "radial_velocity_linf": float(np.max(np.abs(ur))),
        "radial_mach_linf": float(np.max(np.abs(ur)) / sound_ref),
        "axial_velocity_linf": float(np.max(np.abs(uz))),
        "near_axis_pressure_linf_rel":
            float(np.max(np.abs(dp[near_axis])) / args.p0),
        "near_axis_radial_mach_linf":
            float(np.max(np.abs(ur[near_axis])) / sound_ref),
        "density_min": float(np.min(density)),
        "density_max": float(np.max(density)),
    }
    output = json.dumps(metrics, indent=2, sort_keys=True)
    print(output)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(output + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
