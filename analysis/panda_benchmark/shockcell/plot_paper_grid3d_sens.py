"""Paper figure: 3D m142 grid-sensitivity family — AMR grid distributions.

Five independent runs (lmax = 1..5), each panel drawn from that run's own
plotfile BoxArray (y = 0 meridional slice). 2x3 layout, panel (f) slot hosts
the shared level legend. Style unified with paperfig_grid2d_L1to6.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

D = 0.0254
BASE = "/home/qiaoj/testcerisse/cerisse/analysis/panda_benchmark"
d = np.load(f"{BASE}/shockcell/GRID_boxes_3d_sens.npz")

FILL = [plt.cm.Blues(0.12 + 0.13 * l) for l in range(7)]
EDGE, ELW = "0.25", 0.45
# L0: no dedicated run, but the base grid is identical in every case — reuse it
TAGS = ["L0", "L1", "L2", "L3", "L4", "L5"]
PANEL = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

plt.rcParams.update({"font.size": 10.5, "axes.linewidth": 0.65,
                     "xtick.direction": "in", "ytick.direction": "in"})

fig, axes = plt.subplots(2, 3, figsize=(10.0, 3.55),
                         gridspec_kw={"wspace": 0.06, "hspace": 0.22})

for k, tag in enumerate(TAGS):
    ax = axes.flat[k]
    finest = 0 if tag == "L0" else int(d[f"{tag}_finest"])
    for l in range(finest + 1):
        B = d["L1_lev0"] if tag == "L0" else d[f"{tag}_lev{l}"]
        # boxes crossing the y=0 plane (first cell layer above the axis plane)
        yc = 2e-5
        sel = B[(B[:, 2] <= yc) & (B[:, 3] >= yc)]
        for b in sel:
            ax.add_patch(Rectangle((b[0] / D, b[4] / D),
                                   (b[1] - b[0]) / D, (b[5] - b[4]) / D,
                                   facecolor=FILL[l], edgecolor=EDGE,
                                   linewidth=ELW))
    ax.set_xlim(0, 16); ax.set_ylim(-3, 3)
    ax.set_aspect("equal")
    ax.text(0.022, 0.88, f"{PANEL[k]} L{finest}",
            transform=ax.transAxes, fontsize=10.5,
            bbox=dict(facecolor="white", edgecolor="0.4", linewidth=0.6,
                      boxstyle="square,pad=0.25"))
    if k >= 3:
        ax.set_xlabel(r"$x/D_e$")
    else:
        ax.set_xticklabels([])
    if k % 3 == 0:
        ax.set_ylabel(r"$z/D_e$")
    else:
        ax.set_yticklabels([])

handles = [Patch(facecolor=FILL[l], edgecolor=EDGE, linewidth=ELW,
                 label=f"$\\ell={l}$") for l in range(6)]
fig.legend(handles=handles, ncol=6, loc="lower center", frameon=False,
           fontsize=9.5, bbox_to_anchor=(0.5, 0.005), columnspacing=1.2,
           handlelength=1.3, handletextpad=0.5)

fig.subplots_adjust(left=0.045, right=0.995, top=0.985, bottom=0.175)
for ext in ("png", "pdf"):
    fig.savefig(f"{BASE}/m142_3d_grid_levels.{ext}", dpi=300,
                bbox_inches="tight", pad_inches=0.05)
print("saved m142_3d_grid_levels.png/pdf")
