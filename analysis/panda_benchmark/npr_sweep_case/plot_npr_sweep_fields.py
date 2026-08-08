#!/usr/bin/env python3
"""Pressure-ratio sweep field figure in the thesis Fig. 3.14 layout
(vsweep_fields_reliable.png conventions), rebuilt from the m142 NPR sweep:
rows eta0 = 3.27, 4.2, 4.9, 7.5, 15; columns mean Mach, temporal rms Mach,
numerical schlieren.  Data: plateau-window online statistics (D_e/256).

M_rms from stored second moments by first-order propagation,
  Var(M) = Var(q)/c^2 + (M^2/4T^2) Var(T),  q = |u|,
with Var(q) = (u^2 Var(u) + v^2 Var(v) + 2 u v Cov(u,v))/q^2; the q-T
covariance is not stored and is neglected.

Runs locally or on the head node (LITBASE points at the npz directory).
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D_E = 0.0254
GAM, RGAS = 1.4, 287.0
BASE = os.environ.get(
    "LITBASE",
    "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark/lit_compare")
OUTDIR = Path(os.environ.get("LITOUT", BASE))

CASES_5 = [  # (mean npz stem, eta0 label, actual NPR)
    ("sweep_3.273", "3.27", 3.273446),
    ("sweep_4.200", "4.2", 4.200),
    ("sweep_4.900", "4.9", 4.900),
    ("nprfix_7.5", "7.5", 7.5),
    ("nprfix_15", "15", 15.0),
]
CASES_FULL = [
    ("sweep_1.893", "1.89", 1.89293),
    ("sweep_2.394", "2.39", 2.394),
    ("sweep_2.800", "2.80", 2.800),
    ("sweep_3.273", "3.27", 3.273446),
    ("sweep_3.671", "3.67", 3.671),
    ("sweep_4.200", "4.2", 4.200),
    ("sweep_4.900", "4.9", 4.900),
    ("sweep_5.746", "5.75", 5.746),
    ("nprfix_7.5", "7.5", 7.5),
    ("nprfix_10", "10", 10.0),
    ("nprfix_15", "15", 15.0),
    ("nprfix_20", "20", 20.0),
]
FULL = len(sys.argv) > 1 and sys.argv[1] == "full"
CASES = CASES_FULL if FULL else CASES_5

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 9,
})


def load_case(stem):
    # prefer window-differenced clean statistics (nprfix cases: the
    # [2.0,2.75] ms segment contains decaying start-up excursions)
    try:
        m = np.load(f"{BASE}/data/LIT2D_{stem}_clean.npz")
        q = m
        print(f"  ({stem}: clean tail window)")
    except FileNotFoundError:
        m = np.load(f"{BASE}/data/LIT2D_{stem}.npz")
        try:
            q = np.load(f"{BASE}/data/LIT2D_{stem}_mrms.npz")
        except FileNotFoundError:
            q = np.load(f"{BASE}/data/LIT2D_{stem}_sqr.npz")
    dx = float(m["dx"])
    nr, nz = m["DensityMEAN"].shape
    r = (np.arange(nr) + 0.5) * dx / D_E
    z = (np.arange(nz) + 0.5) * dx / D_E
    ur, uz = m["x_velocityMEAN"], m["y_velocityMEAN"]
    T = m["temperatureMEAN"]
    qbar2 = ur**2 + uz**2
    qbar = np.sqrt(qbar2)
    c2 = GAM * RGAS * T
    mean_M = qbar / np.sqrt(c2)
    var_u = np.clip(q["x_velocitySQR"] - ur**2, 0, None)
    var_v = np.clip(q["y_velocitySQR"] - uz**2, 0, None)
    cov_uv = q["xy_velocityMEAN"] - ur * uz
    var_T = np.clip(q["temperatureSQR"] - T**2, 0, None)
    var_q = (ur**2 * var_u + uz**2 * var_v + 2 * ur * uz * cov_uv) / \
        np.maximum(qbar2, 1.0)
    var_M = np.clip(var_q, 0, None) / c2 + (mean_M**2 / (4 * T**2)) * var_T
    std_M = np.sqrt(np.clip(var_M, 0, None))
    return r, z, mean_M, std_M, m["Density"], T, uz


def disk_location(z, uz_axis, T_axis):
    """Axis M=1 downcross behind the first-cell expansion (mean field);
    None when there is no genuine deceleration below M=0.95."""
    M = uz_axis / np.sqrt(GAM * RGAS * T_axis)
    sel = z < 6.0
    zs, Ms = z[sel], M[sel]
    ipk = int(np.argmax(Ms))
    if Ms[ipk] < 1.25:
        return None
    sub = np.where(Ms[ipk:] < 1.0)[0]
    if len(sub) == 0:
        return None
    i0 = ipk + sub[0]
    up = np.where(Ms[i0:] > 1.0)[0]
    i1 = len(Ms) - 1 if len(up) == 0 else i0 + up[0]
    if float(Ms[i0:i1 + 1].min()) > 0.95:
        return None
    return float(zs[i0])


def mirror(field, radius):
    return (np.vstack((field[:0:-1], field)),
            np.concatenate((-radius[:0:-1], radius)))


def schlieren(density, radius, axial):
    grad_r, grad_x = np.gradient(density, radius, axial, edge_order=2)
    magnitude = np.sqrt(grad_r**2 + grad_x**2)
    reference = np.nanpercentile(magnitude[:, axial <= 6.0], 99.5)
    return np.log1p(10.0 * magnitude / max(reference, 1.0e-30))


def main():
    height = 8.1 if len(CASES) == 5 else 1.52 * len(CASES) + 0.9
    fig, axes = plt.subplots(len(CASES), 3, figsize=(7.5, height))
    mach_image = rms_image = schlieren_image = None

    for row, (stem, label, npr) in enumerate(CASES):
        r, z, mean_M, std_M, rho_inst, T, uz = load_case(stem)
        disk = disk_location(z, uz[0], T[0])
        print(f"eta0={label}: max meanM={np.nanmax(mean_M):.2f} "
              f"p99.9 Mrms={np.nanpercentile(std_M, 99.9):.2f} "
              f"disk={disk}")

        mean_mach, mr = mirror(mean_M, r)
        rms_mach, _ = mirror(std_M, r)
        schl, _ = mirror(schlieren(rho_inst, r * D_E, z * D_E), r)
        extent = [z.min(), z.max(), mr.min(), mr.max()]

        mach_image = axes[row, 0].imshow(
            mean_mach, origin="lower", extent=extent, aspect="equal",
            cmap="turbo", vmin=0.0, vmax=4.8)
        axes[row, 0].contour(z, mr, mean_mach, levels=[1.0],
                             colors="white", linewidths=0.55)
        rms_image = axes[row, 1].imshow(
            rms_mach, origin="lower", extent=extent, aspect="equal",
            cmap="magma", vmin=0.0, vmax=1.1)
        axes[row, 1].contour(z, mr, mean_mach, levels=[1.0],
                             colors="white", linewidths=0.5)
        schlieren_image = axes[row, 2].imshow(
            schl, origin="lower", extent=extent, aspect="equal",
            cmap="gray_r", vmin=0.0, vmax=2.6)

        if disk is not None:
            axes[row, 0].axvline(disk, color="white", ls=":", lw=0.9)
            axes[row, 1].axvline(disk, color="white", ls=":", lw=0.9)
            axes[row, 2].axvline(disk, color="black", ls=":", lw=0.9)

        axes[row, 0].set_ylabel(rf"$\eta_0={label}$" + "\n" + r"$r/D_e$")
        for column in range(3):
            ax = axes[row, column]
            ax.set_xlim(0.0, 6.0)
            ax.set_ylim(-2.0, 2.0)
            ax.tick_params(direction="out", length=2.5, pad=1.5, labelsize=8)
            if row == len(CASES) - 1:
                ax.set_xlabel(r"$x/D_e$")
            else:
                ax.set_xticklabels([])
            if column > 0:
                ax.set_yticklabels([])

    axes[0, 0].set_title(r"$\overline{M}$")
    axes[0, 1].set_title(r"$M_{\mathrm{rms}}$")
    axes[0, 2].set_title("numerical schlieren")

    cb_frac = 8.1 / height
    fig.subplots_adjust(left=0.08, right=0.995, bottom=0.13 * cb_frac,
                        top=1 - 0.025 * cb_frac, wspace=0.15, hspace=0.06)
    colorbar_axes = [
        fig.add_axes([0.08, 0.045 * cb_frac, 0.27, 0.018 * cb_frac]),
        fig.add_axes([0.38, 0.045 * cb_frac, 0.27, 0.018 * cb_frac]),
        fig.add_axes([0.68, 0.045 * cb_frac, 0.27, 0.018 * cb_frac]),
    ]
    colorbars = [
        fig.colorbar(mach_image, cax=colorbar_axes[0], orientation="horizontal"),
        fig.colorbar(rms_image, cax=colorbar_axes[1], orientation="horizontal"),
        fig.colorbar(schlieren_image, cax=colorbar_axes[2],
                     orientation="horizontal"),
    ]
    for colorbar, label in zip(
            colorbars,
            [r"$\overline{M}$", r"$M_{\mathrm{rms}}$",
             r"$\log(1+k|\nabla\rho|)$"]):
        colorbar.set_label(label, labelpad=-1)
        colorbar.ax.tick_params(labelsize=8, length=2.5, pad=1.5)

    out = OUTDIR / ("npr_sweep_fields_m142_full" if FULL else
                    "npr_sweep_fields_m142")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"saved {out}.png/.pdf")


if __name__ == "__main__":
    main()
