#!/usr/bin/env python3
"""Plot the near-field mean density and pressure used in thesis Figure 3.11."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
FIELD_DIR = HERE / "bc_variants" / "fields"

RHO_J = 1.6418
P_INF = 99_780.0
ETA_0 = 3.2735
P_E = ETA_0 * P_INF / 1.892929158737854

STATIONS = np.array([0.60, 0.75, 0.90, 1.05, 1.25, 1.55, 2.00, 3.05])


def load_field(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.load(FIELD_DIR / f"SLM_llfweno_{name}.npz")
    return data["x"], data["y"], data["field"]


def main() -> None:
    x, y, rho = load_field("DensityMEAN")
    xp, yp, pressure = load_field("pressureMEAN")
    xu, yu, ux = load_field("x_velocityMEAN")
    xv, yv, uy = load_field("y_velocityMEAN")

    if not (
        np.array_equal(x, xp)
        and np.array_equal(x, xu)
        and np.array_equal(x, xv)
        and np.array_equal(y, yp)
        and np.array_equal(y, yu)
        and np.array_equal(y, yv)
    ):
        raise ValueError("The extracted mean fields do not share a common grid.")

    ix = (x >= 0.0) & (x <= 3.10)
    iy = (y >= -1.0) & (y <= 1.0)
    xc = x[ix]
    yc = y[iy]
    rho_c = rho[np.ix_(ix, iy)].T / RHO_J
    pressure_c = pressure[np.ix_(ix, iy)].T / P_E
    ux_c = ux[np.ix_(ix, iy)].T
    uy_c = uy[np.ix_(ix, iy)].T

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10.2, 4.2),
        sharey=True,
        constrained_layout=True,
    )

    density_image = axes[0].imshow(
        rho_c,
        origin="lower",
        extent=(xc[0], xc[-1], yc[0], yc[-1]),
        aspect="equal",
        interpolation="bilinear",
        cmap="inferno",
        vmin=0.35,
        vmax=1.55,
    )

    # A modestly downsampled grid preserves the mean streamline pattern while
    # keeping the rendering deterministic and inexpensive.
    ds = 4
    axes[0].streamplot(
        xc[::ds],
        yc[::ds],
        ux_c[::ds, ::ds],
        uy_c[::ds, ::ds],
        color="white",
        density=0.78,
        linewidth=0.45,
        arrowsize=0.65,
        minlength=0.08,
    )

    pressure_image = axes[1].imshow(
        pressure_c,
        origin="lower",
        extent=(xc[0], xc[-1], yc[0], yc[-1]),
        aspect="equal",
        interpolation="bilinear",
        cmap="viridis",
        vmin=0.30,
        vmax=1.05,
    )

    for station in STATIONS:
        axes[1].axvline(
            station,
            color="white",
            linestyle=(0, (4, 3)),
            linewidth=0.9,
            alpha=0.95,
        )
        axes[1].text(
            station,
            -0.965,
            f"{station:g}",
            rotation=90,
            ha="center",
            va="bottom",
            color="white",
            fontsize=9.5,
        )

    axes[0].set_title(r"$\langle\rho\rangle/\rho_j$")
    axes[1].set_title(r"$\langle p\rangle/p_e$")

    for ax in axes:
        ax.set_xlim(0.0, 3.10)
        ax.set_ylim(-1.0, 1.0)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"$x/D_e$")
        ax.set_xticks(np.arange(0.0, 3.01, 0.5))
        ax.set_yticks(np.arange(-1.0, 1.01, 0.25))
        ax.tick_params(length=3.5, width=0.8)

    axes[0].set_ylabel(r"$y/D_e$")
    axes[1].tick_params(labelleft=False)

    density_bar = fig.colorbar(
        density_image,
        ax=axes[0],
        fraction=0.040,
        pad=0.035,
    )
    density_bar.set_label(r"$\langle\rho\rangle/\rho_j$", fontsize=12)
    density_bar.ax.tick_params(labelsize=10.5, length=3)

    pressure_bar = fig.colorbar(
        pressure_image,
        ax=axes[1],
        fraction=0.040,
        pad=0.035,
    )
    pressure_bar.set_label(r"$\langle p\rangle/p_e$", fontsize=12)
    pressure_bar.ax.tick_params(labelsize=10.5, length=3)

    output = HERE / "l5near_rho_p"
    fig.savefig(
        output.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.03,
    )
    fig.savefig(
        output.with_suffix(".pdf"),
        bbox_inches="tight",
        pad_inches=0.03,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
