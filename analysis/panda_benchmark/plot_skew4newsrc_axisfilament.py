#!/usr/bin/env python3
"""Radial cuts of the [2,6] ms mean density and radial velocity for
b2d_skew4_newsrc_L4 (new-tree skew-4 + JST, G = rF port), showing the
grid-attached axis filament: the anomaly is confined to the first 1-2
cells at the local finest spacing, while the off-axis jet is regular.
Statistics corrected for the stats_start_time dilution (factor 1.5)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SCR = "/tmp/claude-1000/-home-qiaoj-testcerisse-cerisse/b1e704a6-bfa2-43bd-8eef-98cfa73815d3/scratchpad"
d = np.load(SCR + "/radial_cuts.npz")
names = [str(n) for n in d["names"]]
D, RHO_J = 0.0254, 1.6413
imean, iur = names.index("DensityMEAN"), names.index("x_velocityMEAN")

stations = [("z2.5mm", "$x/D_e = 0.1$", "#7a0e0e"),
            ("z25.4mm", "$x/D_e = 1$", "#c23616"),
            ("z101.6mm", "$x/D_e = 4$", "#e67e22"),
            ("z152.4mm", "$x/D_e = 6$", "#f0b429")]

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
fig.subplots_adjust(left=0.075, right=0.985, top=0.93, bottom=0.14, wspace=0.24)

for tag, lab, col in stations:
    r = d["r_" + tag]/D
    rho = d["dat_" + tag][imean]*1.5/RHO_J
    ur = d["dat_" + tag][iur]*1.5
    w = r <= 0.31
    a1.plot(r[w], rho[w], color=col, lw=1.6, label=lab,
            path_effects=None)
    a2.plot(r[w], ur[w], color=col, lw=1.6, label=lab)
dxf = float(d["dxf"])/D
for ax in (a1, a2):
    for k in (1, 2, 4):
        ax.axvline(k*dxf, color="0.75", lw=0.7, ls=":")
    ax.set_xlim(0, 0.31)
    ax.grid(True, ls="--", lw=0.5, color="0.93")
    ax.set_axisbelow(True)
    ax.set_xlabel("$r/D_e$")
a1.set_ylabel(r"$\overline{\rho}/\rho_j$")
a1.set_title("(a) mean density, [2, 6] ms", fontsize=10)
a1.legend(fontsize=8.6)
a2.set_ylabel(r"$\overline{u}_r$  [m/s]")
a2.set_title("(b) mean radial velocity (axis regularity requires $O(r)$)",
             fontsize=10)
a2.axhline(0, color="0.35", lw=0.9, ls=":")
fig.savefig(HERE/"skew4newsrc_axis_filament.png", dpi=250)
fig.savefig(HERE/"skew4newsrc_axis_filament.pdf")
print("wrote skew4newsrc_axis_filament.png/pdf")
