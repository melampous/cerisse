#!/usr/bin/env python3
"""Mach-disk / first-cell location against the classical correlations.

Panel (a) is the conventional plot. On it several points are unavoidably close:
at eta_0 = 3.27 four results lie within 0.20 D_e, and at eta_0 = 15 and 20 the
top-hat and the flanged-rescaled values differ by less than 0.07 D_e. At a
vertical scale that has to span 0 to 3.4 D_e those separations are smaller than
the markers.

Panel (b) removes that by dividing out the correlation itself. Plotting
x_fc/(0.6455 sqrt(eta_0)) leaves the same data on an axis whose full range is
0.45 instead of 3.4, so a 0.02 D_e difference becomes legible. This is a change
of ordinate, not of position: no point is displaced, nothing is jittered, and
the abscissa is untouched, which is the only way to decongest a quantitative
correlation plot without misreporting where the data sit.

Marker sizes are ordered largest to smallest with the drawing order to match,
so a large symbol never hides a small one, and the large ones carry a partly
transparent face so an overlap still reads as an overlap.
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
NV = np.array([xfc_of(q[1], q[2], q[3]) for q in NEW])

CRIST = lambda q: 0.6455*S(q)


def gibb(x):
    x = np.asarray(x, float)
    out = np.full_like(x, np.nan)
    m1 = (x > 1.893) & (x <= 3.473)
    m2 = x > 3.473
    out[m1] = 0.76*(x[m1]-1.893)**0.61
    out[m2] = 0.84*(x[m2]-1.893)**0.415
    return out


# --- series. eta_0 = 3..7.5 carried over; 3.2734, 15, 20 measured here.
TH_E = np.array([3, 4, 5, 7.5, 3.2734, 15, 20], float)
TH_V = np.array([0.842, 1.148, 1.349, 1.760, NV[0], NV[1], NV[2]])
TH_NEW = np.array([0, 0, 0, 0, 1, 1, 1], bool)   # bookkeeping only, not drawn

d = np.load(OUT/"mysweep_machdisk.npz")
F_eta, F_val, F_disk = d["npr"], d["val"], d["disk"].astype(bool)
keep = F_eta >= 2.0
F_eta, F_val, F_disk = F_eta[keep], F_val[keep], F_disk[keep]
CD = 0.911
FSC = 1/np.sqrt(CD)
M142, TOPHAT, ETA142 = 1.072, 0.9524, 3.273446

print("%-34s %8s %8s %8s" % ("series", "eta_0", "x_fc", "/Crist"))
for e_, v_, n_ in zip(TH_E, TH_V, TH_NEW):
    print("%-34s %8.4f %8.4f %8.4f"
          % ("top-hat" + ("  <- measured here" if n_ else ""), e_, v_,
             v_/CRIST(e_)))

# ------------------------------------------------------------------ drawing
# largest first so nothing large is drawn over something small
SER = [
    dict(k="star_raw", e=[ETA142], v=[M142], m="*", ms=17, fc="#17becf",
         ec="k", mew=0.8, a=0.85, lab=r"3D $M_j$=1.42 baseline, $L/D_e$"),
    dict(k="star_eff", e=[ETA142], v=[M142*FSC], m="*", ms=17, fc="none",
         ec="#17becf", mew=1.8, a=1.0, lab=r"3D baseline, $L/D_{\rm eff}$"),
    dict(k="plus3d", e=[ETA142], v=[TOPHAT], m="P", ms=11, fc="#9467bd",
         ec="k", mew=0.7, a=0.85, lab=r"3D top-hat ($C_d{\approx}1$)"),
    dict(k="tophat", e=TH_E, v=TH_V, m="o", ms=9.5,
         fc="#d62728", ec="k", mew=0.9, a=0.85, lab="top-hat simulation"),
    dict(k="f_eff", e=F_eta[F_disk], v=F_val[F_disk]*FSC, m="D", ms=6.6,
         fc="#2ca02c", ec="k", mew=0.7, a=0.9,
         lab=r"flanged, $L/D_{\rm eff}$ ($C_d{=}0.911$)"),
    dict(k="f_eff_nd", e=F_eta[~F_disk], v=F_val[~F_disk]*FSC, m="D", ms=6.2,
         fc="none", ec="#2ca02c", mew=1.4, a=1.0,
         lab=r"flanged, no disk ($x_{fc}$, rescaled)"),
    dict(k="f_raw", e=F_eta[F_disk], v=F_val[F_disk], m="D", ms=5.0,
         fc="#bfe3bf", ec="#2ca02c", mew=0.8, a=0.95,
         lab=r"flanged $D_e/256$, raw $L/D_e$"),
    dict(k="f_raw_nd", e=F_eta[~F_disk], v=F_val[~F_disk], m="D", ms=4.8,
         fc="none", ec="#9dc89d", mew=1.1, a=1.0, lab=None),
]

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.8,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": False, "ytick.right": True,
                     "legend.frameon": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.8, 5.5))
fig.subplots_adjust(left=0.058, right=0.988, top=0.845, bottom=0.245, wspace=0.185)

e = np.linspace(1.9, 21, 400)
se = S(e)
for ax, norm in ((a1, False), (a2, True)):
    n = CRIST(e) if norm else 1.0
    ax.plot(se, 0.6455*se/n, "-", color="0.55", lw=1.3,
            label=r"Crist $0.6455\sqrt{\eta_0}$")
    ax.plot(se, 0.67*se/n, "--", color="0.35", lw=1.3,
            label=r"Ashkenas–Sherman $0.67\sqrt{\eta_0}$")
    ax.plot(se, gibb(e)/n, ":", color="0.45", lw=1.7, label="Gibbings")
    for z, s_ in enumerate(SER):
        ee = np.asarray(s_["e"], float)
        vv = np.asarray(s_["v"], float)/(CRIST(np.asarray(s_["e"], float))
                                         if norm else 1.0)
        ax.plot(S(ee), vv, s_["m"], ms=s_["ms"], mfc=s_["fc"], mec=s_["ec"],
                mew=s_["mew"], alpha=s_["a"], ls="", zorder=5+z,
                label=s_["lab"] if (not norm and s_["lab"]) else None)
    ax.set_xlim(1.35, 4.65)
    ax.set_xlabel(r"$\sqrt{\eta_0}$", fontsize=11.5)
    ax.grid(True, ls=":", lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)
    sec = ax.secondary_xaxis("top", functions=(lambda s: s**2,
                                               lambda v: S(np.maximum(v, 0))))
    ET = [2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0]
    sec.set_xticks(ET)
    sec.set_xticklabels(["%g\n(%.2f)" % (q, S(5.0*(q**(2.0/7.0)-1.0)))
                         for q in ET])
    sec.set_xlabel(r"$\eta_0$   ($M_j$)", fontsize=10.5, labelpad=5)
    sec.tick_params(labelsize=8.5)

a1.set_ylim(0, 3.4)
a1.set_ylabel(r"$x_{fc}/D_e$  or  $L_{MD}/D_e$", fontsize=11.5)
a1.text(0.030, 0.955, "(a)", transform=a1.transAxes, fontsize=11, va="top")
a2.axhline(1.0, color="0.55", lw=1.3)
a2.set_ylim(0.60, 1.14)
a2.set_ylabel(r"$x_{fc}\,/\,0.6455\sqrt{\eta_0}$", fontsize=11.5)
a2.text(0.030, 0.955, "(b)", transform=a2.transAxes, fontsize=11, va="top")

h, l = a1.get_legend_handles_labels()
fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=3,
           fontsize=9.0, handlelength=1.9, columnspacing=2.0, labelspacing=0.45)
for ext in ("png", "pdf"):
    fig.savefig(OUT/("mach_disk_tophat_paper.%s" % ext), dpi=300)
print("\nwrote mach_disk_tophat_paper.png/pdf")
