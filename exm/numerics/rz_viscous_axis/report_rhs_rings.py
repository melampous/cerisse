#!/usr/bin/env python3
"""Print the first radial rings of every stored manufactured RHS field."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--rings", type=int, default=4)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_boxlib_rz(args.plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    lo = np.asarray(ds.domain_left_edge, dtype=float)
    hi = np.asarray(ds.domain_right_edge, dtype=float)
    dr = (hi[0] - lo[0]) / dims[0]
    radius = lo[0] + (np.arange(dims[0]) + 0.5) * dr

    result = {"plotfile": str(args.plotfile), "dr": float(dr), "rings": []}
    for ring in range(min(args.rings, dims[0])):
        row = {"ring": ring, "r": float(radius[ring]), "fields": {}}
        for name in ("Density", "Xmom", "Ymom", "Zmom", "Energy"):
            value = None
            for namespace in ("boxlib", "gas", "stream"):
                try:
                    value = np.squeeze(np.asarray(grid[(namespace, name)]))
                    break
                except Exception:
                    pass
            if value is None:
                continue
            if value.shape[0] != dims[0] and value.shape[1] == dims[0]:
                value = value.T
            radial_row = value[ring, :]
            row["fields"][name] = {
                "mean": float(np.mean(radial_row)),
                "min": float(np.min(radial_row)),
                "max": float(np.max(radial_row)),
                "absmax": float(np.max(np.abs(radial_row))),
            }
        result["rings"].append(row)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
