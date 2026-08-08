#!/usr/bin/env python3
"""Effective dissipation of the mean shock-cell structure per numerical
scheme: oscillation envelope of the centreline mean density (detrended by
a one-cell 1.2 D_e moving average, Hilbert envelope, exponential fit over
1-6 D_e) for the experiment, baseline WENO-Z5 NS, TENO5 NS and Euler.
Larger decay length L_d = cells persist longer = less effective
dissipation of the mean cell structure."""
import numpy as np
from scipy.signal import hilbert
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})
HERE = Path(__file__).parent
D = 0.0254
RHO_J = 1.6413
FITLO, FITHI = 1.0, 6.0
CB, CT, CE = "#b3251e", "#167c72", "#6a51a3"

def moving_average(v, x, w):
    n = max(3, int(round(w / (x[1] - x[0]))))
    if n % 2 == 0:
        n += 1
    out = np.convolve(v, np.ones(n) / n, mode="same")
    out[:n//2] = v[:n//2]; out[-(n//2):] = v[-(n//2):]
    return out

xg = np.arange(0.60, 7.0, 0.04)          # uniform evaluation grid

def envelope(x, r):
    rr = np.interp(xg, x, r)
    osc = rr - moving_average(rr, xg, 1.2)
    env = np.abs(hilbert(osc))
    m = (xg >= FITLO) & (xg <= FITHI)
    c = np.polyfit(xg[m], np.log(env[m]), 1)
    return env, c

CURVES = []
xd = np.load(HERE.parent / "panda_m142den_axis_EPAPS.npz")
CURVES.append(("experiment", "k", xd["x"], xd["rho"]))
a3 = np.load(HERE / "AX3T_L5.npz", allow_pickle=True)
CURVES.append(("WENO-Z5 NS", CB, a3["x"] / D, a3["ax_DensityMEAN"] / RHO_J))
for lab, col, fn in (("TENO5 NS", CT, "SIM_teno_rhoMEAN_z0.npz"),
                     ("Euler (WENO-Z5)", CE, "SIM_euler_rhoMEAN_z0.npz")):
    d = np.load(HERE / fn)
    jc = np.argsort(np.abs(d["y"]))[:2]
    CURVES.append((lab, col, d["x"], d["rho"][:, jc].mean(1) / RHO_J))

fig, ax = plt.subplots(figsize=(7.0, 2.9))
ax.set_axisbelow(True)
ax.grid(True, ls="--", lw=0.5, color="0.87")
print("scheme            L_d [D_e]")
for lab, col, x, r in CURVES:
    env, c = envelope(x, r)
    Ld = -1.0 / c[0]
    print("%-16s  %6.2f" % (lab, Ld))
    if lab == "experiment":
        ax.semilogy(xg, env, "o", color="k", ms=2.6, mew=0,
                    label="experiment ($L_d=%.1f\\,D_e$)" % Ld)
    else:
        ax.semilogy(xg, env, color=col, lw=1.2,
                    label="%s ($L_d=%.1f\\,D_e$)" % (lab, Ld))
    xx = np.linspace(FITLO, FITHI, 40)
    ax.semilogy(xx, np.exp(np.polyval(c, xx)), color=col, lw=0.8, ls=":")
ax.set_xlim(0.8, 6.5)
ax.set_ylim(2e-2, 6e-1)
ax.set_xlabel("$x/D_e$", labelpad=2)
ax.set_ylabel("oscillation envelope $/\\rho_j$", labelpad=2)
ax.tick_params(length=2.5)
ax.legend(fontsize=6.8, frameon=False, loc="lower left", handlelength=1.5)
fig.tight_layout()
fig.savefig(HERE / "paperfig_schemes_envelope.png", dpi=300)
fig.savefig(HERE / "paperfig_schemes_envelope.pdf")
print("wrote paperfig_schemes_envelope")
