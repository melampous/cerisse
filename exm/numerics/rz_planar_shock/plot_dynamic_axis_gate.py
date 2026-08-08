#!/usr/bin/env python3
"""Plot the dynamic R-Z axis-pressure gate for all saved test cases.

The figure is deliberately based only on plotfile fields.  It therefore gives
an independent visual check of the scalar metrics written by
``analyze_current_axis_gate.py``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yt
import yt.frontends.boxlib.data_structures as boxlib_data_structures


# yt 4.3 plus NumPy 2 can expose a read-only padded domain edge while loading
# a 2-D cylindrical plotfile.  Geometry is irrelevant to these direct reads of
# cell-centred arrays, so use the same safe workaround as the gate analyser.
_geometry = boxlib_data_structures.Geometry
boxlib_data_structures.Geometry = lambda name: _geometry(
    "cartesian" if name == "cylindrical" else name
)
yt.funcs.mylog.setLevel(40)


@dataclass
class Snapshot:
    step: int
    time: float
    z: np.ndarray
    pressure: np.ndarray
    radial_velocity: np.ndarray
    sound_speed: np.ndarray


def read_snapshot(plotfile: Path, gamma: float) -> Snapshot:
    dataset = yt.load(str(plotfile))
    grid = dataset.covering_grid(
        level=0,
        left_edge=dataset.domain_left_edge,
        dims=dataset.domain_dimensions,
    )

    def field(name: str) -> np.ndarray:
        values = np.asarray(grid[("boxlib", name)].d, dtype=np.float64)
        return np.squeeze(values)

    pressure = field("pressure")
    density = field("density")
    radial_velocity = field("x_velocity")
    if pressure.ndim != 2:
        raise ValueError(f"expected a 2-D pressure field in {plotfile}, got {pressure.shape}")

    dimensions = np.asarray(dataset.domain_dimensions, dtype=int)
    left = np.asarray(dataset.domain_left_edge.d, dtype=float)
    right = np.asarray(dataset.domain_right_edge.d, dtype=float)
    dz = (right[1] - left[1]) / dimensions[1]
    z = left[1] + (np.arange(dimensions[1]) + 0.5) * dz
    return Snapshot(
        step=int(plotfile.name.removeprefix("plt")),
        time=float(dataset.current_time),
        z=z,
        pressure=pressure,
        radial_velocity=radial_velocity,
        sound_speed=np.sqrt(gamma * pressure / density),
    )


def find_plotfiles(case_directory: Path) -> list[Path]:
    plots = sorted(
        path
        for path in case_directory.glob("plt[0-9]*")
        if (path / "Header").is_file()
    )
    if not plots:
        raise FileNotFoundError(f"no plotfiles found in {case_directory}")
    return plots


def relative_ring_error(pressure: np.ndarray) -> np.ndarray:
    axis = pressure[0, :]
    second = pressure[1, :]
    scale = np.maximum.reduce(
        (np.abs(axis), np.abs(second), np.full_like(axis, np.finfo(float).tiny))
    )
    return np.abs(axis - second) / scale


def history(snapshot_series: list[Snapshot]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    times = []
    radial_pressure_spreads = []
    radial_mach_errors = []
    for snapshot in snapshot_series:
        pressure = snapshot.pressure
        local_scale = np.maximum(
            np.max(np.abs(pressure), axis=0), np.finfo(float).tiny
        )
        times.append(snapshot.time)
        radial_pressure_spreads.append(
            np.max(np.ptp(pressure, axis=0) / local_scale)
        )
        radial_mach_errors.append(
            np.max(np.abs(snapshot.radial_velocity) / snapshot.sound_speed)
        )
    return (
        np.asarray(times),
        np.asarray(radial_pressure_spreads),
        np.asarray(radial_mach_errors),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "results_directory",
        nargs="?",
        type=Path,
        default=Path("results/dynamic_axis_gate_20260807_current"),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--scheme", default="llf-wenoz5")
    args = parser.parse_args()

    root = args.results_directory.resolve()
    output = (
        args.output.resolve()
        if args.output
        else root / "dynamic_axis_pressure_gate_summary.png"
    )
    cases = (
        ("m3_plus", r"Mach 3, $u_z>0$"),
        ("m10_plus", r"Mach 10, $u_z>0$"),
        ("m10_minus", r"Mach 10, $u_z<0$"),
    )

    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(
        3,
        3,
        figsize=(15.2, 10.2),
        sharex=False,
        constrained_layout=True,
    )
    pressure_colors = ("#1f77b4", "#d62728", "#2ca02c")
    pressure_labels = ("axis ring", "second ring", "outer ring")
    history_floor = 1.0e-18

    for column, (case_name, title) in enumerate(cases):
        plotfiles = find_plotfiles(root / case_name)
        snapshots = [read_snapshot(plotfile, args.gamma) for plotfile in plotfiles]
        final = snapshots[-1]
        ring_indices = (0, 1, final.pressure.shape[0] - 1)

        top = axes[0, column]
        for radial_index, color, label in zip(
            ring_indices, pressure_colors, pressure_labels
        ):
            top.plot(
                final.z,
                final.pressure[radial_index, :],
                color=color,
                linewidth=2.1 if radial_index == 0 else 1.35,
                linestyle="-" if radial_index == 0 else ("--" if radial_index == 1 else ":"),
                label=label,
            )
        top.set_title(f"{title}\nfinal: step {final.step}, t={final.time:.5g}")
        if column == 0:
            top.set_ylabel("pressure")
        top.legend(loc="best", frameon=True, fontsize=9)

        middle = axes[1, column]
        ring_error = relative_ring_error(final.pressure)
        middle.semilogy(
            final.z,
            np.maximum(ring_error, history_floor),
            color="#7f3c8d",
            linewidth=1.45,
        )
        middle.axhline(
            1.0e-10,
            color="#d62728",
            linewidth=1.1,
            linestyle="--",
            label=r"gate $10^{-10}$",
        )
        middle.set_ylim(5.0e-18, 3.0e-9)
        middle.text(
            0.03,
            0.10,
            rf"max $={np.max(ring_error):.2e}$",
            transform=middle.transAxes,
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "fc": "white", "alpha": 0.85},
        )
        if column == 0:
            middle.set_ylabel(r"$|p_0-p_1|/p_{local}$")
        middle.set_xlim(final.z[0], final.z[-1])
        middle.set_xlabel("z")
        middle.legend(loc="upper right", fontsize=8)

        times, pressure_history, radial_mach_history = history(snapshots)
        bottom = axes[2, column]
        bottom.semilogy(
            times,
            np.maximum(pressure_history, history_floor),
            "o-",
            color="#11a579",
            linewidth=1.5,
            markersize=4.0,
            label=r"max$_z$ $(p_{max,r}-p_{min,r})/p_{local}$",
        )
        bottom.semilogy(
            times,
            np.maximum(radial_mach_history, history_floor),
            "s--",
            color="#e73f74",
            linewidth=1.35,
            markersize=3.8,
            label=r"max $|u_r|/a$",
        )
        bottom.axhline(1.0e-10, color="#d62728", linewidth=1.1, linestyle="--")
        bottom.set_ylim(5.0e-18, 3.0e-9)
        if times[-1] > times[0]:
            margin = 0.025 * (times[-1] - times[0])
            bottom.set_xlim(times[0] - margin, times[-1] + margin)
        bottom.set_xlabel("time")
        if column == 0:
            bottom.set_ylabel("radial-uniformity metric")
        bottom.legend(loc="best", frameon=True, fontsize=8)

    figure.suptitle(
        f"{args.scheme} R-Z dynamic axis-pressure gate — "
        "all-fluid, uniform L0, no IBM",
        fontsize=14,
        fontweight="semibold",
        y=1.025,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    main()
