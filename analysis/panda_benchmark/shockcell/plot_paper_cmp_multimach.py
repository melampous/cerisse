#!/usr/bin/env python3
"""Multi-Mach first-cell metrics: 2D NPR-sweep (99 um, quasi-steady
staircase statistics) vs the five Panda & Seasholtz (1999) centreline
density traverses. x_1 and S_1 are density-peak based (same criteria as
the single-condition analysis); Mj = 1.0 leg excluded (no cell system)."""
import numpy as np, csv
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True})
HERE = Path(__file__).parent
D = 0.0254

def boxfilt(a, x, wD):
    w = max(3, int(round(wD/(x[1]-x[0])))) | 1
    o = np.convolve(a, np.ones(w)/w, mode="same")
    o[:w//2] = a[:w//2]; o[-(w//2):] = a[-(w//2):]
    return o

def refine(x, y, i):
    if i <= 0 or i >= len(y)-1: return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])

def firstcell(x, r):
    rs = boxfilt(r, x, 0.05)
    dist = max(1, int(round(0.5/(x[1]-x[0]))))
    pmin,_ = find_peaks(np.where(x > 0.35, -rs, -np.inf), prominence=0.02, distance=dist)
    pmax,_ = find_peaks(np.where(x > 0.35, rs, -np.inf), prominence=0.02, distance=dist)
    ev = sorted([(i, "min") for i in pmin] + [(i, "max") for i in pmax])
    seq = []; expect = "min"
    for i, t in ev:
        if t == expect:
            seq.append((t, i)); expect = "max" if t == "min" else "min"
    maxs = [i for t, i in seq if t == "max"]
    if len(maxs) < 2: return np.nan, np.nan, np.nan
    x1 = refine(x, rs, maxs[0]); x2 = refine(x, rs, maxs[1])
    return x1, x2, x2-x1

sw = np.load(HERE / "SWEEP_axis_rhoMEAN.npz")
sim = []
for tag in ["2p394", "2p800", "3p273", "3p671", "4p200", "4p900", "5p746"]:
    eta = float(tag.replace("p", "."))
    Mj = np.sqrt(5.0*(eta**(2.0/7.0)-1.0))
    rho = sw["rho_"+tag]; dx = float(sw["dx_"+tag])
    x = (np.arange(len(rho))+0.5)*dx/D
    rhoj = 1.1590*(1+0.2*Mj**2)
    sim.append((Mj, *firstcell(x, rho/rhoj)))
sim = np.array(sim)
exp = []
with open(HERE.parent.parent.parent/"U_jet/analysis/panda_benchmark/panda_multiMach_exp_params.csv") as f:
    for row in csv.DictReader(f):
        exp.append((float(row["Mj"]), float(row["peak1_x"]),
                    float(row["Ls12_De"]), float(row["Pack_Ls"])))
exp = np.array(exp)

fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7))
axs[0].plot(sim[:, 0], sim[:, 1], "o-", color="#1f4e9c", ms=4, lw=0.9,
            label="2D sweep (99 µm)")
axs[0].plot(exp[:, 0], exp[:, 1], "s", color="k", ms=4.5, mfc="none", mew=0.9,
            label="Panda & Seasholtz (1999)")
axs[0].set_xlabel("$M_j$", labelpad=2); axs[0].set_ylabel("$x_1/D_e$", labelpad=2)
axs[0].legend(fontsize=6.6, frameon=False, loc="lower right", handlelength=1.5)
axs[1].plot(sim[:, 0], sim[:, 3], "o-", color="#1f4e9c", ms=4, lw=0.9,
            label="2D sweep (99 µm)")
axs[1].plot(exp[:, 0], exp[:, 2], "s", color="k", ms=4.5, mfc="none", mew=0.9,
            label="Exp. $S_1$")
axs[1].plot(exp[:, 0], exp[:, 3], ":", color="0.4", lw=1.0, label="Pack (1950)")
axs[1].set_xlabel("$M_j$", labelpad=2); axs[1].set_ylabel("$S_1/D_e$", labelpad=2)
axs[1].legend(fontsize=6.6, frameon=False, loc="upper left",
              bbox_to_anchor=(0.12, 1.0), handlelength=1.5)
for k, ax in enumerate(axs):
    ax.tick_params(length=2.5)
    ax.text(0.03, 0.95, "(%s)" % chr(97+k), transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
fig.tight_layout(w_pad=1.4)
fig.savefig(HERE / "paperfig_cmp_multimach.png", dpi=300)
fig.savefig(HERE / "paperfig_cmp_multimach.pdf")
print("wrote paperfig_cmp_multimach")
