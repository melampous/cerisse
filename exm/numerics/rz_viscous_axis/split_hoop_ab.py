#!/usr/bin/env python3
"""Recover the R-Z radial viscous face/source split from a hoop on/off A/B.

The executable writes the complete spatial RHS.  With all other inputs held
fixed, ``cns.rz_visc_hoop=0`` removes only ``-tau_theta_theta/r``.  Therefore
the off result is the metric face-flux divergence and on-minus-off is the
geometric hoop source, without adding instrumentation to the production RHS.
"""

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


def radial_momentum(plotfile: Path) -> tuple[np.ndarray, np.ndarray]:
    ds = load_boxlib_rz(plotfile)
    dims = np.asarray(ds.domain_dimensions, dtype=int).copy()
    grid = ds.covering_grid(0, ds.domain_left_edge, dims)
    field = None
    for namespace in ("boxlib", "gas", "stream"):
        try:
            field = np.squeeze(np.asarray(grid[(namespace, "Xmom")]))
            break
        except Exception:
            pass
    if field is None:
        raise KeyError("Xmom")
    if field.shape[0] != dims[0] and field.shape[1] == dims[0]:
        field = field.T
    lo = np.asarray(ds.domain_left_edge, dtype=float)
    hi = np.asarray(ds.domain_right_edge, dtype=float)
    dr = (hi[0] - lo[0]) / field.shape[0]
    radius = lo[0] + (np.arange(field.shape[0]) + 0.5) * dr
    return field, radius


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("hoop_on", type=Path)
    parser.add_argument("hoop_off", type=Path)
    parser.add_argument("--rings", type=int, default=4)
    parser.add_argument("--mu", type=float, default=0.01)
    parser.add_argument("--a", type=float, default=0.20)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    total, radius = radial_momentum(args.hoop_on)
    face, radius_off = radial_momentum(args.hoop_off)
    if total.shape != face.shape or not np.array_equal(radius, radius_off):
        raise ValueError("hoop on/off grids do not match")
    source = total - face
    expected_face = (2.0 / 3.0) * args.mu * args.a / radius

    rows = []
    for ring in range(min(args.rings, total.shape[0])):
        rows.append(
            {
                "ring": ring,
                "r": float(radius[ring]),
                "face_mean": float(np.mean(face[ring, :])),
                "hoop_source_mean": float(np.mean(source[ring, :])),
                "total_mean": float(np.mean(total[ring, :])),
                "face_z_spread": float(np.ptp(face[ring, :])),
                "source_z_spread": float(np.ptp(source[ring, :])),
                "total_z_absmax": float(np.max(np.abs(total[ring, :]))),
                "continuous_face_target": float(expected_face[ring]),
                "continuous_source_target": float(-expected_face[ring]),
            }
        )

    print(
        json.dumps(
            {
                "hoop_on": str(args.hoop_on),
                "hoop_off": str(args.hoop_off),
                "interpretation": {
                    "face": "R_r,viscous-face = RHS(hoop off)",
                    "hoop_source": "R_tau,geometric = RHS(on) - RHS(off)",
                    "total": "face + hoop_source",
                },
                "rings": rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
