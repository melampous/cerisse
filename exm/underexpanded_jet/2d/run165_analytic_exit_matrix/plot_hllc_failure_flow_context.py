#!/usr/bin/env python3
"""Plot the AFD-HLLC-WENO failure face in its Run165 RZ flow context."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_hllc_muscl_afd_weno_pointwise as pointwise


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "compare_hllc_muscl_vs_afd_weno"
PLOTFILE = RESULTS / "afd_fallback_off" / "plt00184"
OUTPUT = RESULTS / "same_state_reconstruction" / "hllc_failure_flow_context.png"
TARGET = (4, 171)
GAMMA = 1.4


def main() -> None:
    pointwise.patch_yt_cylindrical_readonly_edge()
    fields, time_us, r_edges, x_edges, r_centres, x_centres = pointwise.load_fields(PLOTFILE)
    rho = fields["rho"]
    pressure = fields["pressure"]
    speed = np.hypot(fields["ur"], fields["uz"])
    sound = np.sqrt(GAMMA * pressure / rho)
    mach = speed / sound
    dr = (r_centres[1] - r_centres[0]) * 1.0e-3
    dx = (x_centres[1] - x_centres[0]) * 1.0e-3
    grad_rho = np.hypot(*np.gradient(rho, dr, dx))
    schlieren = np.log10(grad_rho / max(float(grad_rho.max()), 1.0e-300) + 1.0e-6)

    ti, tj = TARGET
    target_r = r_centres[ti]
    target_x = x_centres[tj]
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.8), constrained_layout=True, sharex=True, sharey=True)
    panels = [
        (mach, "turbo", 0.0, 7.0, "Mach number"),
        (np.log10(pressure), "viridis", None, None, r"$\log_{10}(p/\mathrm{Pa})$"),
        (schlieren, "gray", -6.0, 0.0, r"normalized $\log_{10}|\nabla\rho|$"),
    ]
    for axis, (field, cmap, vmin, vmax, title) in zip(axes, panels):
        image = axis.pcolormesh(x_edges, r_edges, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        axis.contour(x_centres, r_centres, mach, levels=[1.0], colors="white", linewidths=0.7, alpha=0.8)
        axis.plot(target_x, target_r, marker="x", markersize=11, markeredgewidth=2.5, color="red")
        axis.annotate(
            "failure cell\n(4,171)",
            (target_x, target_r),
            xytext=(target_x + 3.0, target_r + 5.0),
            arrowprops={"arrowstyle": "->", "color": "red", "lw": 1.3},
            color="red",
            fontsize=8,
        )
        axis.set_title(title)
        axis.set_xlabel("axial x [mm]")
        axis.set_xlim(-35.0, 2.5)
        axis.set_ylim(0.0, 20.0)
        axis.grid(alpha=0.12)
        fig.colorbar(image, ax=axis, shrink=0.86)
    axes[0].set_ylabel("radius r [mm]")
    axes[0].annotate(
        "retro-jet flow direction",
        xy=(-18.0, 17.5),
        xytext=(-5.0, 17.5),
        arrowprops={"arrowstyle": "->", "lw": 1.4, "color": "black"},
        ha="center",
        fontsize=8,
    )
    fig.suptitle(f"AFD-HLLC-WENO pre-failure flow context, t={time_us:.2f} us")
    OUTPUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)
    print(OUTPUT)
    print(f"target: r={target_r:.9g} mm, x={target_x:.9g} mm, M={mach[ti,tj]:.9g}, p={pressure[ti,tj]:.9g} Pa")


if __name__ == "__main__":
    main()
