#!/usr/bin/env python3
"""Analyze the radial-momentum RHS from the paired R-Z pressure operator."""

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


def analytic_rhs(
    rr: np.ndarray,
    zz: np.ndarray,
    profile: int,
    amplitude: float,
    quartic: float,
    layer_z: float,
    curvature: float,
    width: float,
) -> np.ndarray:
    if profile == 0:
        return -2.0 * amplitude * rr
    if profile == 2:
        return -2.0 * amplitude * rr - 4.0 * quartic * rr**3
    arg = (zz - layer_z - curvature * rr**2) / width
    return (
        2.0
        * amplitude
        * curvature
        * rr
        / width
        / np.cosh(arg) ** 2
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--profile", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--p0", type=float, default=2.0)
    parser.add_argument("--amplitude", type=float, default=0.25)
    parser.add_argument("--quartic", type=float, default=0.10)
    parser.add_argument("--layer-z", type=float, default=0.5)
    parser.add_argument("--curvature", type=float, default=0.20)
    parser.add_argument("--width", type=float, default=0.04)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()

    yt.funcs.mylog.setLevel(50)
    ds = load_plotfile(args.plotfile)
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
    exact = analytic_rhs(
        rr,
        zz,
        args.profile,
        args.amplitude,
        args.quartic,
        args.layer_z,
        args.curvature,
        args.width,
    )
    error = rhs - exact
    first4 = np.zeros_like(error, dtype=bool)
    first4[: min(4, rhs.shape[0]), :] = True
    pressure_scale = max(abs(args.p0) / dr, 1.0)
    weighted_l2 = np.sqrt(np.sum(rr * error**2) / np.sum(rr))

    metrics = {
        "plotfile": str(args.plotfile),
        "profile": args.profile,
        "nr": int(rhs.shape[0]),
        "nz": int(rhs.shape[1]),
        "dr": float(dr),
        "rhs_absmax": float(np.max(np.abs(rhs))),
        "exact_absmax": float(np.max(np.abs(exact))),
        "error_linf": float(np.max(np.abs(error))),
        "error_weighted_l2": float(weighted_l2),
        "axis_rhs_absmax": float(np.max(np.abs(rhs[0]))),
        "axis_exact_absmax": float(np.max(np.abs(exact[0]))),
        "axis_error_linf": float(np.max(np.abs(error[0]))),
        "first4_error_linf": float(np.max(np.abs(error[first4]))),
        "rhs_over_p_by_dr": float(np.max(np.abs(rhs)) / pressure_scale),
        "axis_error_over_p_by_dr": float(
            np.max(np.abs(error[0])) / pressure_scale
        ),
        "finite": bool(np.all(np.isfinite(rhs))),
    }

    encoded = json.dumps(metrics, indent=2)
    print(encoded)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(encoded + "\n", encoding="utf-8")

    if args.figure is not None:
        import matplotlib.pyplot as plt

        jline = int(np.argmax(np.abs(exact[0])))
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
        axes[0].plot(r, rhs[:, jline], label="computed RHS")
        axes[0].plot(r, exact[:, jline], "--", label="exact -dp/dr")
        axes[0].set_xlabel("r")
        axes[0].set_ylabel("radial-momentum RHS")
        axes[0].set_title(f"z={z[jline]:.6g}")
        axes[0].grid(alpha=0.25)
        axes[0].legend()
        image = axes[1].imshow(
            error.T,
            origin="lower",
            aspect="auto",
            extent=(r[0], r[-1], z[0], z[-1]),
            cmap="coolwarm",
        )
        axes[1].set_xlabel("r")
        axes[1].set_ylabel("z")
        axes[1].set_title("RHS error")
        fig.colorbar(image, ax=axes[1])
        fig.tight_layout()
        args.figure.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.figure, dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
