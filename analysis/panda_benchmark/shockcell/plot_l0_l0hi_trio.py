#!/usr/bin/env python3
"""Three jets (m142/m119/m180), 2D RZ single-level quasi-steady TRUE statistics:
L0 (dx=1584 um) vs L0hi (base doubled, dx=792 um) vs Panda & Seasholtz (1999).
Radial mean-density profiles at the experimental stations of each jet.
rho_j convention: rho_amb*(1+0.2*Mj^2), rho_amb=1.1590.
Outputs: fig48 (m142), fig49 (m119), fig50 (m180)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).parent
JOBS = [
    ("m142", "panda_fig5a_mj144_EPAPS.npz", 1.6413,
     [("RADT_L0lo.npz", "L0 (1584 µm) [8.0,9.5]ms", "#4292c6"),
      ("RADT_m142_2d_L0hi.npz", "L0hi (792 µm, base×2) [8.0,9.5]ms", "#d94801")],
     "fig48_m142_l0_l0hi_vs_exp.png",
     "m142 (Mj=1.44) 2D quasi-steady L0 vs L0hi vs Panda fig5a"),
    ("m119", "panda_fig6_mj119_EPAPS.npz", 1.4873,
     [("RADT_m119_2d_L0x.npz", "L0 (1584 µm) [8.0,9.5]ms", "#4292c6"),
      ("RADT_m119_2d_L0hix.npz", "L0hi (792 µm, base×2) [8.0,9.5]ms", "#d94801")],
     "fig49_m119_l0_l0hi_vs_exp.png",
     "m119 (Mj=1.19) 2D quasi-steady L0 vs L0hi vs Panda fig6"),
    ("m180", "panda_fig5c_mj180_EPAPS.npz", 1.9100,
     [("RADT_m180_2d_L0x.npz", "L0 (1584 µm) [11.5,13.0]ms", "#4292c6"),
      ("RADT_m180_2d_L0hix.npz", "L0hi (792 µm, base×2) [8.0,9.5]ms", "#d94801")],
     "fig50_m180_l0_l0hi_vs_exp.png",
     "m180 (Mj=1.80) 2D quasi-steady L0 vs L0hi vs Panda fig5c"),
]
for tag, expf, rhoj, cases, outf, title in JOBS:
    avail = [(fn, lab, col) for fn, lab, col in cases if (HERE / fn).exists()]
    if not avail:
        print(f"{tag}: no sim npz yet, skipped"); continue
    exp = np.load(HERE.parent / expf, allow_pickle=True)
    d0 = np.load(HERE / avail[0][0], allow_pickle=True)
    stations = d0["stations"]
    nst = len(stations)
    ncol = 4
    nrow = int(np.ceil(nst / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(16, 3.6 * nrow), sharey=True)
    for i, st in enumerate(stations):
        ax = axs.flat[i]
        for fn, lab, col in avail:
            d = np.load(HERE / fn, allow_pickle=True)
            rr = d["rr"]
            prof = d["rad_DensityMEAN"][i] / rhoj
            mid = prof[15]
            lo = np.nanmin(prof, axis=0); hi = np.nanmax(prof, axis=0)
            m = rr <= 2.0
            ax.fill_between(rr[m], lo[m], hi[m], color=col, alpha=0.12, lw=0)
            ax.plot(rr[m], mid[m], color=col, lw=1.5, label=lab if i == 0 else None)
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
