#!/usr/bin/env python3
"""Report first-ring R-Z velocity/vorticity values without treating r=dr/2 as r=0."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yt


def load_plotfile(path: Path):
    from yt.data_objects.static_output import MutableAttribute

    original_set = MutableAttribute.__set__

    def writable_set(self, instance, value):
        self.data[instance] = value

    MutableAttribute.__set__ = writable_set
    try:
        return yt.load(str(path))
    finally:
        MutableAttribute.__set__ = original_set


def raw_field(grid, name: str) -> np.ndarray:
    for namespace in ("boxlib", "gas", "stream"):
        try:
            return np.squeeze(np.asarray(grid[(namespace, name)]))
        except Exception:
            pass
    raise KeyError(name)


def radial_first(values: np.ndarray, nr: int) -> np.ndarray:
    if values.shape[0] == nr:
        return values
    if values.shape[1] == nr:
        return values.T
    raise ValueError(f"Cannot identify radial dimension in {values.shape}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--rings", type=int, default=4)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_plotfile(args.plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int)
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    dr = float((ds.domain_right_edge[0] - ds.domain_left_edge[0]) / dims[0])

    density = radial_first(raw_field(grid, "Density"), int(dims[0]))
    radial_momentum = radial_first(raw_field(grid, "Xmom"), int(dims[0]))
    axial_momentum = radial_first(raw_field(grid, "Ymom"), int(dims[0]))
    arrays = {
        "u_r": radial_momentum / density,
        "u_z": axial_momentum / density,
        "mag_vorticity": radial_first(
            raw_field(grid, "magvort"), int(dims[0])
        ),
    }
    rows = []
    for ring in range(min(args.rings, int(dims[0]))):
        row = {"ring_i": ring, "r": (ring + 0.5) * dr}
        for label, values in arrays.items():
            line = values[ring, :]
            row[f"{label}_min"] = float(np.min(line))
            row[f"{label}_max"] = float(np.max(line))
            row[f"{label}_absmax"] = float(np.max(np.abs(line)))
        rows.append(row)

    result = {
        "plotfile": str(args.plotfile),
        "note": "ring_i=0 is centred at r=dr/2; it is not an r=0 sample",
        "dr": dr,
        "rows": rows,
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
