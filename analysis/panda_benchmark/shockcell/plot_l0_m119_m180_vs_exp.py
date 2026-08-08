#!/usr/bin/env python3
"""m119 / m180 2D RZ L0 (dx=1584 um, quasi-steady) TRUE statistics vs experiment.
Radial density at the exact experimental stations (fig6 Mj=1.19, fig5c Mj=1.80).
Normalization: rho_j = rho_amb*(1+0.2*Mj^2), rho_amb=1.1590 (same convention
that yields the established 1.6413 for m142/Mj=1.4426).
Outputs: fig44_m119L0_radial_vs_exp.png, fig45_m180L0_radial_vs_exp.png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
JOBS = [
    ("m119", "RADT_m119_2d_L0x.npz", "panda_fig6_mj119_EPAPS.npz", 1.4873,
     "fig44_m119L0_radial_vs_exp.png",
     "m119 (Mj=1.19) 2D L0 quasi-steady TRUE-stat [8.0,9.5]ms vs Panda fig6"),
    ("m180", "RADT_m180_2d_L0x.npz", "panda_fig5c_mj180_EPAPS.npz", 1.9100,
     "fig45_m180L0_radial_vs_exp.png",
     "m180 (Mj=1.80) 2D L0 quasi-steady TRUE-stat [11.5,13.0]ms vs Panda fig5c"),
]
for tag, simf, expf, rhoj_sim, outf, title in JOBS:
    d = np.load(HERE / simf, allow_pickle=True)
    exp = np.load(HERE.parent / expf, allow_pickle=True)
    stations = d["stations"]
    rr = d["rr"]
    nst = len(stations)
    ncol = 4
    nrow = int(np.ceil(nst / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(16, 3.6 * nrow), sharey=True)
    for i, st in enumerate(stations):
        ax = axs.flat[i]
        prof = d["rad_DensityMEAN"][i] / rhoj_sim      # (noff, NR)
        mid = prof[15]                                  # offset +0.00
        lo = np.nanmin(prof, axis=0); hi = np.nanmax(prof, axis=0)
        m = rr <= 2.0
        ax.fill_between(rr[m], lo[m], hi[m], color="#4292c6", alpha=0.15, lw=0)
        ax.plot(rr[m], mid[m], color="#084594", lw=1.6,
                label="2D L0 (dx=1584 µm)" if i == 0 else None)
        key = next((k for k in exp.files
                    if k.startswith("x") and abs(float(k[1:]) - st) < 1e-6), None)
        if key is not None:
            e = exp[key]
            r = np.abs(e[:, 0]); rho = e[:, 1] / float(exp["rhoj"])
            o = np.argsort(r)
            ax.plot(r[o], rho[o], "ko", ms=3.5, mfc="none",
                    label="Panda & Seasholtz (1999)" if i == 0 else None)
        ax.set_title(f"$x/D_e$ = {st}", fontsize=10)
        ax.grid(alpha=0.25, lw=0.4)
        if i // ncol == nrow - 1:
            ax.set_xlabel("$r/D_e$")
        if i % ncol == 0:
            ax.set_ylabel("$\\langle\\rho\\rangle/\\rho_j$")
    for j in range(nst, nrow * ncol):
        axs.flat[j].axis("off")
    fig.legend(loc="lower right", fontsize=9, frameon=False)
    fig.suptitle(title + "  (band: station ±0.30$D_e$ axial offsets)", y=0.995)
    fig.tight_layout(rect=[0, 0.02, 1, 0.98])
    fig.savefig(HERE / outf, dpi=160)
    print("wrote", outf)
