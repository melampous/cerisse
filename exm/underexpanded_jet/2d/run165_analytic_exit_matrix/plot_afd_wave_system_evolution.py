#!/usr/bin/env python3
"""Plot the developing counterflow-jet wave system before the AFD failure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_hllc_muscl_afd_weno_pointwise as pointwise


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_hllc_muscl_vs_afd_weno"
OUTPUT = RESULTS / "afd_wave_system_evolution_68_73p6us.png"
STEPS = (170, 176, 182, 184)
P_INF = 534.5806702256042
GAMMA = 1.4


def axis_negative_flow_mach_one(
    mach: np.ndarray, axial_velocity: np.ndarray, x_centres: np.ndarray
) -> float:
    for j in range(mach.shape[1] - 2, 0, -1):
        if (
            axial_velocity[0, j] < 0.0
            and axial_velocity[0, j - 1] < 0.0
            and mach[0, j] >= 1.0
            and mach[0, j - 1] < 1.0
        ):
            fraction = (mach[0, j] - 1.0) / (
                mach[0, j] - mach[0, j - 1]
            )
            return float(
                x_centres[j]
                + fraction * (x_centres[j - 1] - x_centres[j])
            )
    return np.nan


def main() -> None:
    pointwise.patch_yt_cylindrical_readonly_edge()
    pointwise.yt.funcs.mylog.setLevel(50)

    snapshots = []
    for step in STEPS:
        fields, time_us, r_edges, x_edges, r_centres, x_centres = (
            pointwise.load_fields(pointwise.plotfile_for("afd_off", step))
        )
        pressure = fields["pressure"]
        sound_speed = np.sqrt(GAMMA * pressure / fields["rho"])
        mach = np.hypot(fields["ur"], fields["uz"]) / sound_speed
        disk_x = axis_negative_flow_mach_one(mach, fields["uz"], x_centres)
        local_js = np.arange(168, 175)
        peak_j = int(local_js[np.argmax(mach[4, local_js])])
        snapshots.append(
            (
                time_us,
                pressure,
                mach,
                r_edges,
                x_edges,
                r_centres,
                x_centres,
                disk_x,
                peak_j,
            )
        )

    fig, axes = plt.subplots(
        2, len(STEPS), figsize=(15.2, 6.8), sharex=True, sharey=True,
        constrained_layout=True
    )
    mach_image = pressure_image = None
    for column, snapshot in enumerate(snapshots):
        (
            time_us,
            pressure,
            mach,
            r_edges,
            x_edges,
            r_centres,
            x_centres,
            disk_x,
            peak_j,
        ) = snapshot
        mach_image = axes[0, column].pcolormesh(
            x_edges, r_edges, mach, shading="auto", cmap="turbo",
            vmin=0.0, vmax=7.0
        )
        pressure_image = axes[1, column].pcolormesh(
            x_edges,
            r_edges,
            np.log10(np.maximum(pressure / P_INF, 1.0e-8)),
            shading="auto",
            cmap="magma",
            vmin=0.0,
            vmax=2.4,
        )
        for axis in axes[:, column]:
            axis.contour(
                x_centres, r_centres, mach, levels=[1.0], colors="white",
                linewidths=0.7
            )
            axis.scatter(
                [x_centres[peak_j]], [r_centres[4]], marker="x",
                color="cyan", s=45, linewidths=1.5
            )
            if np.isfinite(disk_x):
                axis.axvline(disk_x, color="cyan", linestyle="--", linewidth=0.8)
            axis.set(xlim=(-35.0, 2.5), ylim=(0.0, 20.0))
        axes[0, column].set_title(
            f"t={time_us:.1f} us\naxis jet M=1: x={disk_x:.1f} mm"
        )
        axes[1, column].set_xlabel("axial x [mm]")

    axes[0, 0].set_ylabel("r [mm]\nMach")
    axes[1, 0].set_ylabel("r [mm]\nlog10(p/p_inf)")
    axes[0, 0].annotate(
        "jet direction",
        xy=(-14.0, 18.6),
        xytext=(-3.0, 18.6),
        arrowprops={"arrowstyle": "->", "color": "white"},
        color="white",
        ha="center",
        va="center",
    )
    if mach_image is not None:
        fig.colorbar(mach_image, ax=axes[0, :], label="Mach", shrink=0.9)
    if pressure_image is not None:
        fig.colorbar(
            pressure_image, ax=axes[1, :], label="log10(p/p_inf)", shrink=0.9
        )
    fig.suptitle(
        "AFD fallback-off: counterflow-jet wave-system development\n"
        "white: M=1; cyan x: moving local Mach peak at i=4; dashed: axis jet M=1",
        fontsize=13,
    )
    fig.savefig(OUTPUT, dpi=220)
    print(OUTPUT)


if __name__ == "__main__":
    main()
