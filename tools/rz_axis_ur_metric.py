#!/usr/bin/env python3
"""Measure near-axis radial velocity in a 2D AMReX RZ plotfile.

The metric is intentionally explicit: choose a number of radial cell layers and
an axial window, then report Linf and RMS/L2 of x_velocity over that mask.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", help="AMReX plotfile directory")
    parser.add_argument("--ncells", type=int, default=2,
                        help="number of radial cell layers from r=0")
    parser.add_argument("--zmin", type=float, default=0.0)
    parser.add_argument("--zmax", type=float, default=None)
    parser.add_argument("--mask", choices=("all", "fluid", "non_ghost", "fluid_non_ghost"),
                        default="fluid_non_ghost")
    parser.add_argument("--out-level", type=int, default=-1)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo / "tools"))
    import schlieren as sch  # pylint: disable=import-error,import-outside-toplevel

    names = ["x_velocity"]
    if args.mask in ("fluid", "fluid_non_ghost"):
        names.append("sld")
    if args.mask in ("non_ghost", "fluid_non_ghost"):
        names.append("ghs")

    fields, r, z, dr, dz, meta = sch.read_amrex_plotfile_multi(
        args.plotfile, names, out_level=args.out_level)
    ur = fields["x_velocity"]

    mask = np.zeros_like(ur, dtype=bool)
    mask[:args.ncells, :] = True
    mask &= z[None, :] >= args.zmin
    if args.zmax is not None:
        mask &= z[None, :] <= args.zmax
    if args.mask in ("fluid", "fluid_non_ghost"):
        mask &= fields["sld"] < 0.5
    if args.mask in ("non_ghost", "fluid_non_ghost"):
        mask &= fields["ghs"] < 0.5

    if not np.any(mask):
        raise SystemExit("empty diagnostic mask")

    vals = ur[mask]
    linf = float(np.max(np.abs(vals)))
    l2 = float(np.sqrt(np.mean(vals * vals)))
    mean = float(np.mean(vals))
    print(f"plotfile: {args.plotfile}")
    print(f"time: {meta['sim_time']:.16e}")
    print(f"level: {meta['output_level']}  shape: {ur.shape}  dr: {dr:.8e}  dz: {dz:.8e}")
    print(f"mask: ncells={args.ncells} zmin={args.zmin:.8e} "
          f"zmax={args.zmax if args.zmax is not None else 'None'} mask={args.mask} n={vals.size}")
    print(f"ur_Linf: {linf:.10e}")
    print(f"ur_L2:   {l2:.10e}")
    print(f"ur_mean: {mean:.10e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
