"""F3: m180 axis virtual-Pitot pressure vs Zapryagaev et al. (2015) Fig. 7.
F4: m180 radial virtual-Pitot profiles at the Zapryagaev stations (Fig. 6).

Primary curve: 2D RZ NPR=5.746 plateau statistics (D_e/256, 1.25 ms window).
Secondary: m180 3D instantaneous slice (D_e/128, t=6.5 ms, record_stats off).
Experimental condition: NPR=5.0, M_j=1.71 -> quantitative offsets expected;
reported feature positions annotated, no digitised experimental curve drawn.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lit_common import (D, GAM, RGAS, PINF, BASE, load2d, axis_state,
                        pitot_ratio, probe_avg, load_m180_slice,
                        m180_axis_inst, mach_field)

plt.rcParams.update({"font.size": 11, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True,
                     "ytick.right": True})

NPR = 5.746
P0 = NPR * PINF
c = load2d("LIT2D_sweep_5.746.npz")
s = load_m180_slice()

# ---------------- F3: axis ----------------
z, p, T, rho, uz, M = axis_state(c)
pt2d = probe_avg(pitot_ratio(p, np.abs(uz) / np.sqrt(GAM * RGAS * T)) / P0,
                 z[1] - z[0])
xi, pi, Ti, ui, Mi = m180_axis_inst(s)
pt3d = probe_avg(pitot_ratio(pi, np.abs(Mi)) / P0, xi[1] - xi[0])

fig, ax = plt.subplots(figsize=(9.8, 4.6))
ax.plot(z, pt2d, color="#08306b", lw=1.9,
        label="RZ mean, NPR=5.746 ($D_e/256$, 1.25 ms stats)")
ax.plot(xi, pt3d, color="#74a9cf", lw=1.1, alpha=0.9,
        label="3D instantaneous, NPR=5.746 ($D_e/128$, $t$=6.5 ms)")
# Zapryagaev reported features (NPR=5.0): first cell ~1.5D, minima ~2.8, 4.5
feats = [(1.5, "first cell end\n(exp., NPR=5.0)"),
         (2.8, "2nd minimum (exp.)"), (4.5, "3rd minimum (exp.)")]
for xf, lab in feats:
    ax.axvline(xf, color="#e7298a", lw=1.0, ls="--", alpha=0.8)
    ax.text(xf + 0.05, 0.97, lab, rotation=90, va="top", ha="left",
            fontsize=7.5, color="#e7298a")
sc = np.sqrt(NPR / 5.0)
for xf, _ in feats:
    ax.axvline(xf * sc, color="#e7298a", lw=1.0, ls=":", alpha=0.5)
ax.plot([], [], color="#e7298a", ls=":", alpha=0.5,
        label=r"exp. features $\times\sqrt{5.746/5.0}$")
ax.set_xlim(0, 8); ax.set_ylim(0, 1.02)
ax.set_xlabel(r"$x/D_e$")
ax.set_ylabel(r"virtual Pitot pressure $p_{t2}/p_0$")
ax.legend(fontsize=8.5, loc="lower right")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{BASE}/litfig_m180_pitot_axis.{ext}", dpi=250)
plt.close(fig)

# ---------------- F4: radial stations ----------------
STATIONS = [0.015, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]
Mf2 = mach_field(c)
pt2 = pitot_ratio(c["pressureMEAN"], Mf2) / P0            # (nr, nz)
Mf3 = (np.sqrt(s["x_velocity"]**2 + s["y_velocity"]**2) /
       np.sqrt(GAM * RGAS * s["temperature"]))
pt3 = pitot_ratio(s["pressure"], Mf3) / P0                # (ny, nx)
ypos = s["y"] > 0

fig, axes = plt.subplots(3, 4, figsize=(12.6, 8.2), sharey=True)
for k, (st, ax) in enumerate(zip(STATIONS, axes.flat)):
    j2 = np.argmin(np.abs(c["z"] - st))
    prof2 = probe_avg(pt2[:, j2], c["r"][1] - c["r"][0])
    ax.plot(c["r"], prof2, color="#08306b", lw=1.7,
            label="RZ mean" if k == 0 else None)
    i3 = np.argmin(np.abs(s["x"] - st))
    yp = s["y"][ypos]; pp = pt3[ypos, i3]
    yn = -s["y"][~ypos][::-1]; pn = pt3[~ypos, i3][::-1]
    ax.plot(yp, probe_avg(pp, yp[1] - yp[0]), color="#74a9cf", lw=1.0,
            label="3D inst. ($y>0$)" if k == 0 else None)
    ax.plot(yn, probe_avg(pn, yn[1] - yn[0]), color="#a63603", lw=0.8,
            alpha=0.7, label="3D inst. ($y<0$)" if k == 0 else None)
    ax.text(0.95, 0.90, rf"$x/D_e={st:g}$", transform=ax.transAxes,
            ha="right", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="0.6", lw=0.5,
                      boxstyle="square,pad=0.2"))
    ax.set_xlim(0, 1.6); ax.set_ylim(0, 1.05)
    if k >= 8:
        ax.set_xlabel(r"$r/D_e$")
    if k % 4 == 0:
        ax.set_ylabel(r"$p_{t2}/p_0$")
fig.legend(loc="lower center", ncol=3, fontsize=9, frameon=False,
           bbox_to_anchor=(0.5, 0.0))
fig.suptitle("Radial virtual-Pitot profiles, NPR=5.746 — station set of "
             "Zapryagaev et al. (2015, exp. NPR=5.0)", fontsize=11)
fig.tight_layout(rect=[0, 0.035, 1, 0.97])
for ext in ("png", "pdf"):
    fig.savefig(f"{BASE}/litfig_m180_pitot_radial.{ext}", dpi=250)
plt.close(fig)
print("saved litfig_m180_pitot_axis + litfig_m180_pitot_radial")
