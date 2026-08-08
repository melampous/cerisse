#!/usr/bin/env python3
"""Audit uniform R-Z dU/dt globally and in the first four radial rings."""

from __future__ import annotations

import argparse
import csv
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


def orient_radial_first(values: np.ndarray, nr: int) -> np.ndarray:
    if values.ndim != 2:
        raise ValueError(f"Expected a 2-D field, got shape {values.shape}")
    if values.shape[0] == nr:
        return values
    if values.shape[1] == nr:
        return values.T
    raise ValueError(f"Cannot identify radial dimension in {values.shape}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--rings", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_plotfile(args.plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)

    fields = {
        "R_rho": "Density",
        "R_rho_ur": "Xmom",
        "R_rho_uz": "Ymom",
        "R_E": "Energy",
    }
    arrays = {
        label: orient_radial_first(raw_field(grid, raw_name), int(dims[0]))
        for label, raw_name in fields.items()
    }

    rows: list[dict[str, float | int | str]] = []
    ring_count = min(args.rings, int(dims[0]))
    for ring in range(ring_count):
        for label, values in arrays.items():
            line = values[ring, :]
            rows.append(
                {
                    "radial_ring": ring,
                    "field": label,
                    "sample_at_axial_mid": float(line[line.size // 2]),
                    "min": float(np.min(line)),
                    "max": float(np.max(line)),
                    "absmax": float(np.max(np.abs(line))),
                    "l2": float(np.sqrt(np.mean(line * line))),
                }
            )

    global_absmax = {
        label: float(np.max(np.abs(values))) for label, values in arrays.items()
    }
    finite = all(np.all(np.isfinite(values)) for values in arrays.values())
    worst = max(global_absmax.values(), default=0.0)
    result = {
        "plotfile": str(args.plotfile),
        "nr": int(dims[0]),
        "nz": int(dims[1]),
        "rings_checked": ring_count,
        "tolerance": args.tolerance,
        "finite": bool(finite),
        "global_absmax": global_absmax,
        "first_rings": rows,
        "worst_absmax": worst,
        "pass": bool(finite and worst <= args.tolerance),
    }

    encoded = json.dumps(result, indent=2)
    print(encoded)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(encoded + "\n", encoding="utf-8")
    if args.csv is not None:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    raise SystemExit(0 if result["pass"] else 1)


if __name__ == "__main__":
    main()
