#!/usr/bin/env python3
"""Publication figure: m142 2D RZ AMR box layouts for max_level = 1..6.
2x3 panels, no title, panel tags (a)-(f), shared level legend.
Filled translucent patches per level (finest on top) + thin outlines."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from pathlib import Path

plt.rcParams.update({
    "font.size": 10.5, "axes.linewidth": 0.65,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254
d = np.load(HERE / "GRID_boxes_2d_L1to6.npz", allow_pickle=True)
CMAP = plt.cm.Blues
FILL = [CMAP(0.12 + 0.13 * l) for l in range(7)]
EDGE = "#40506080"
TAGS = ["L1", "L2", "L3", "L4", "L5", "L6"]
PANEL = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
XMAX, RMAX = 16.0, 3.0

PW, ML, MR, MB_LEG = 3.05, 0.55, 0.10, 0.22
PH = PW * RMAX / XMAX
WSP, HSP = 0.12, 0.14
FW = ML + 3 * PW + 2 * WSP + MR
FH = MB_LEG + 0.38 + 2 * PH + HSP + 0.06
fig = plt.figure(figsize=(FW, FH))
axs = []
for row in range(2):
    for col in range(3):
        left = (ML + col * (PW + WSP)) / FW
        bottom = (MB_LEG + 0.38 + (1 - row) * (PH + HSP)) / FH
        axs.append(fig.add_axes([left, bottom, PW / FW, PH / FH]))
for k, (tag, ax) in enumerate(zip(TAGS, axs)):
    finest = int(d[f"{tag}_finest"])
    for lev in range(finest + 1):
        for (rlo, rhi, zlo, zhi) in d[f"{tag}_lev{lev}"]:
            if zlo / D > XMAX or rlo / D > RMAX:
                continue
            ax.add_patch(Rectangle((zlo / D, rlo / D),
                                   (zhi - zlo) / D, (rhi - rlo) / D,
                                   facecolor=FILL[lev], edgecolor=EDGE,
                                   linewidth=0.3, zorder=lev + 1))
    ax.set_xlim(0, XMAX); ax.set_ylim(0, RMAX)
    ax.text(0.015, 0.93, f"{PANEL[k]} L{tag[1]}",
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5, zorder=20,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.2))
    ax.tick_params(length=2.5)
    if k < 3:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel("$x/D_e$", labelpad=1.5)
    if k % 3 == 0:
        ax.set_ylabel("$r/D_e$", labelpad=1.5)
    else:
        ax.set_yticklabels([])
handles = [Patch(facecolor=FILL[l], edgecolor="#405060", linewidth=0.4,
                 label=f"$\\ell={l}$") for l in range(7)]
fig.legend(handles=handles, ncol=7, loc="lower center", frameon=False,
           fontsize=9.5, bbox_to_anchor=(0.5, 0.01), handlelength=1.2,
           handleheight=0.9, columnspacing=1.1)
fig.savefig(HERE / "paperfig_grid2d_L1to6.png", dpi=300)
fig.savefig(HERE / "paperfig_grid2d_L1to6.pdf")
print("wrote paperfig_grid2d_L1to6.png/.pdf")
