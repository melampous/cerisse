#!/usr/bin/env python3
"""Mach-disk / first-cell location against the classical correlations, with the
axisymmetric Level 4 top-hat points added.

Reproduces mach_disk_location_thesis_v4_simple and overlays three new points
measured today: the top-hat inlet at NPR_0 = 3.2734, 15 and 20, all Level 4,
statistics over [0.4, 2.4] ms, x_fc taken as the steepest density rise between
the first expansion minimum and the following maximum.

The figure already carried a clean-inlet series at eta_0 = 15 and 20, so the
new points are a direct check on a different grid and inlet treatment rather
than an extension into empty territory.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
SW = OUT.parent/"bc_variants"/"sweep_out"
RHO_J, D = 1.6413, 0.0254
S = np.sqrt


def refine(x, y, i):
    if i <= 0 or i >= len(y)-1:
        return x[i]
    d = y[i-1]-2*y[i]+y[i+1]
    return x[i] if d == 0 else x[i]+0.5*(y[i-1]-y[i+1])/d*(x[1]-x[0])


def xfc_of(f, lo, hi):
    d = np.load(SW/f)
    ks = sorted([k[:-2] for k in d.files if k.endswith("_t")],
                key=lambda q: int("".join(c for c in q if c.isdigit())))
    dx = float(d[ks[0]+"_t"][1])
    r = d[ks[-1]+"_rho"]/RHO_J
    x = (np.arange(r.size)+0.5)*dx/D
    m = (x >= lo) & (x <= hi) & np.isfinite(r)
    xx, rr = x[m], r[m]
    i0 = int(np.argmin(rr))
    i1 = i0 + int(np.argmax(rr[np.arange(i0, len(rr))]))
    g = np.gradient(rr, xx)
    s = np.arange(i0, i1+1)
    return refine(xx, g, s[int(np.nanargmax(g[s]))])


NEW = [(3.2734, "XFC_A0_tophat_L4.npz",      0.35, 1.30),
       (15.0,   "XFC_N15_tophat_LO__L4.npz", 1.40, 3.20),
       (20.0,   "XFC_N20_tophat_LO__L4.npz", 1.70, 3.80)]
NE = np.array([q[0] for q in NEW])
NV = np.array([xfc_of(q[1], q[2], q[3]) for q in NEW])
print("new Level 4 top-hat points")
for e, v in zip(NE, NV):
    print("   eta_0 = %-8.4f  x_fc/D_e = %.4f   Crist %.4f   A-S %.4f"
          % (e, v, 0.6455*S(e), 0.67*S(e)))

e = np.linspace(1.9, 21, 300)
se = S(e)


def gibb(x):
    x = np.asarray(x, float)
    out = np.full_like(x, np.nan)
    m1 = (x > 1.893) & (x <= 3.473)
    m2 = x > 3.473
    out[m1] = 0.76*(x[m1]-1.893)**0.61
    out[m2] = 0.84*(x[m2]-1.893)**0.415
    return out


# The red series is the top-hat simulation. Its values at eta_0 = 15 and 20
# are replaced by the Level 4 runs measured today (2.5156 and 2.9688, against
# the 2.623 and 3.079 that stood there before), and the Panda-condition point
# at eta_0 = 3.2734 is added to the same series. The four points at eta_0 = 3
# to 7.5 are carried over unchanged; their inlet treatment is taken on the
# user's statement that this series is the top-hat simulation, not verified
# here against the runs that produced them.
P_eta = [3, 4, 5, 7.5, 3.2734, 15, 20]
P = [0.842, 1.148, 1.349, 1.760, float(NV[0]), float(NV[1]), float(NV[2])]
NEWMASK = [False, False, False, False, True, True, True]
d = np.load(OUT/"mysweep_machdisk.npz")
F_eta, F_val, F_disk = d["npr"], d["val"], d["disk"].astype(bool)
keep = F_eta >= 2.0
F_eta, F_val, F_disk = F_eta[keep], F_val[keep], F_disk[keep]
CD = 0.911
FSC = 1/np.sqrt(CD)
M142, TOPHAT, ETA142 = 1.072, 0.9524, 3.273446

fig, ax = plt.subplots(figsize=(9.5, 6.9))
ax.plot(se, 0.6455*se, "-", color="0.55", lw=1.4,
        label=r"Crist $0.6455\sqrt{\eta_0}$")
ax.plot(se, 0.67*se, "--", color="0.35", lw=1.4,
        label=r"Ashkenas–Sherman $0.67\sqrt{\eta_0}$")
ax.plot(se, gibb(e), ":", color="0.45", lw=1.8, label="Gibbings")
PE, PV = np.array(P_eta), np.array(P)
NM = np.array(NEWMASK)
ax.plot(S(PE[~NM]), PV[~NM], "o", ms=10, color="#d62728", mec="k", mew=0.9,
        zorder=6, label=r"top-hat simulation")
ax.plot(S(PE[NM]), PV[NM], "o", ms=11, color="#d62728", mec="#ff7f0e", mew=2.2,
        zorder=7, label=r"top-hat, 2D axisym. L4 (this work)")
mF = F_disk
ax.plot(S(F_eta[mF]), F_val[mF], "D", ms=6, mfc="#bfe3bf", mec="#2ca02c",
        mew=0.8, zorder=4, label=r"flanged $D_e/256$, raw $L/D_e$")
ax.plot(S(F_eta[~mF]), F_val[~mF], "D", ms=6, mfc="none", mec="#9dc89d",
        mew=1.1, zorder=4)
ax.plot(S(F_eta[mF]), F_val[mF]*FSC, "D", ms=7.5, color="#2ca02c", mec="k",
        mew=0.7, zorder=5,
        label=r"flanged, $L/D_{\rm eff}$  ($D_{\rm eff}{=}D_e\sqrt{C_d}$, $C_d{=}0.911$)")
ax.plot(S(F_eta[~mF]), F_val[~mF]*FSC, "D", ms=7, mfc="none", mec="#2ca02c",
        mew=1.4, zorder=5, label=r"flanged, no disk ($x_{fc}$, rescaled)")
ax.plot([S(ETA142)], [M142], "*", ms=17, color="#17becf", mec="k", mew=0.8,
        zorder=7, label=r"3D $M_j$=1.42 baseline, $L/D_e$")
ax.plot([S(ETA142)], [M142*FSC], "*", ms=17, mfc="none", mec="#17becf",
        mew=1.8, zorder=7, label=r"3D baseline, $L/D_{\rm eff}$")
ax.plot([S(ETA142)], [TOPHAT], "P", ms=11, color="#9467bd", mec="k", mew=0.7,
        zorder=7, label=r"3D top-hat ($C_d{\approx}1$)")

for x_, y_ in zip(S(NE), NV):
    ax.annotate(r"$\eta_0$=%.4g" % (x_**2), (x_, y_),
                textcoords="offset points", xytext=(10, -13), fontsize=8.5,
                color="#b35a00")

ax.set_xlabel(r"$\sqrt{\eta_0}$", fontsize=13)
ax.set_ylabel(r"$x_{fc}/D_e$  or  $L_{MD}/D_e$", fontsize=12)
ax.set_xlim(1.35, 4.65)
ax.set_ylim(0, 3.4)
sec = ax.secondary_xaxis("top", functions=(lambda s: s**2,
                                           lambda v: np.sqrt(np.maximum(v, 0))))
# label the pressure ratio with the fully expanded Mach number it implies,
# M_j = sqrt(5(eta_0^(2/7) - 1)) for gamma = 1.4, so the abscissa can be read
# either way without a conversion in the reader's head
ETICK = [2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0]
sec.set_xticks(ETICK)
sec.set_xticklabels(["%g\n(%.2f)" % (q, np.sqrt(5.0*(q**(2.0/7.0)-1.0)))
                     for q in ETICK])
sec.set_xlabel(r"$\eta_0$   ($M_j$)", fontsize=11, labelpad=6)
sec.tick_params(labelsize=9)
ax.grid(True, ls=":", lw=0.5, alpha=0.7)
ax.legend(fontsize=8.8, loc="upper left", framealpha=0.95)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT/("mach_disk_location_tophat_L4.%s" % ext), dpi=220)
print("\nwrote mach_disk_location_tophat_L4.png/pdf")
