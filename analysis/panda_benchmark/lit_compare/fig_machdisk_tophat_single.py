#!/usr/bin/env python3
"""Mach-disk / first-cell location against the classical correlations.

Single panel, in the layout of the original figure.

Points that would otherwise hide one another are handled inside the panel
rather than by adding a second one: markers are ordered largest to smallest
with the drawing order to match, so a large symbol is never laid over a small
one, and every filled symbol carries a partly transparent face so a genuine
overlap still reads as an overlap instead of as a single point. No marker is
displaced and no jitter is applied; the positions are the measured ones.

The upper abscissa carries the fully expanded Mach number that each pressure
ratio implies, M_j = sqrt(5(eta_0^(2/7) - 1)) for gamma = 1.4.

The four regime nodes quoted in the flow-regimes section of the chapter are
drawn as vertical markers: the
Muraoka-Hiejima first-Mach-disk band 3.03-3.12 and the Gibbings regular-to-Mach
reflection transition at 3.473. The first is drawn as a band because the
chapter states these boundaries are approximate and geometry dependent rather
than hard thresholds.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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


NEW = [(3.2734, "XFC_A0_tophat_L4.npz", 0.35, 1.30)]
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


# --- series. The six sweep locations use the same mean-field extraction as
# Figure 3.21; eta_0=3.2734 is the separate top-hat comparison at the Panda
# operating condition.
TH_E = np.array([3, 3.2734, 4, 5, 7.5, 15, 20], float)
TH_V = np.array([0.824, NV[0], 1.1385, 1.3321, 1.7323, 2.6054, 3.0464])
TH_NEW = np.array([1, 1, 1, 1, 1, 1, 1], bool)   # bookkeeping only, not drawn

d = np.load(OUT/"mysweep_machdisk.npz")
F_eta, F_val, F_disk = d["npr"], d["val"], d["disk"].astype(bool)
keep = F_eta >= 2.0
F_eta, F_val, F_disk = F_eta[keep], F_val[keep], F_disk[keep]
# Effective discharge coefficient of the baseline shifted-tanh profile
# (Table 3.2).  Use the same value here when converting D_e to D_eff.
CD = 0.8789
FSC = 1/np.sqrt(CD)
ETA142 = 3.273446
# The four Cartesian profiles, x_fc/D_e straight from Table 3.6 of the chapter.
# The values the figure carried before (1.072 for the baseline, 0.9524 for the
# top hat) are not in that table and 1.072 does not correspond to any defined
# feature of the baseline profile: re-extracting the centreline with the
# chapter's own definition, the steepest positive gradient of the first
# recompression, gives 0.8828, and the first mean-density maximum gives 1.1836,
# against 0.883 and 1.176 in the table. At x = 1.072 the profile is simply at
# rho/rho_j = 1.326 partway up the recompression.
# Only the baseline and the top hat are drawn; the wall-attached and thin
# shifted profiles of Table 3.6 (0.930 and 0.922) are left out on request.
TAB36 = [("baseline", 0.883), ("top-hat", 0.945)]

print("%-34s %8s %8s %8s" % ("series", "eta_0", "x_fc", "/Crist"))
for e_, v_, n_ in zip(TH_E, TH_V, TH_NEW):
    print("%-34s %8.4f %8.4f %8.4f"
          % ("top-hat" + ("  <- measured here" if n_ else ""), e_, v_,
             v_/CRIST(e_)))

# ------------------------------------------------------------------ drawing
# Largest first so nothing large is drawn over something small. One fill per
# series and nothing else: every marker is opaque, carries the same black edge
# of the same width, and is painted in exactly the colour of its own legend
# swatch. Nothing about a point's fill varies with the point.
SER = [
    dict(k="d3_base", e=[ETA142], v=[0.883], m="*", ms=16, fc="#17becf",
         ec="k", mew=0.8, lab=r"Cartesian baseline"),
    dict(k="d3_hat", e=[ETA142], v=[0.945], m="P", ms=11, fc="#9467bd",
         ec="k", mew=0.8, lab=r"Cartesian top-hat"),
    dict(k="tophat", e=TH_E, v=TH_V, m="o", ms=9.5,
         fc="#d62728", ec="k", mew=0.8, lab="RZ top-hat sweep"),
    dict(k="f_eff", e=F_eta, v=F_val*FSC, m="D", ms=6.6,
         fc="#2ca02c", ec="k", mew=0.8,
         lab=r"RZ baseline profile, $L/D_{\rm eff}$"),
    dict(k="f_raw", e=F_eta, v=F_val, m="D", ms=5.0,
         fc="#bfe3bf", ec="k", mew=0.8,
         lab=r"RZ baseline profile, $L/D_e$"),
]

plt.rcParams.update({"font.size": 13.2, "axes.linewidth": 0.9,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": False, "ytick.right": True,
                     "legend.frameon": True})
fig, ax = plt.subplots(figsize=(9.6, 7.4))
fig.subplots_adjust(left=0.098, right=0.980, top=0.850, bottom=0.110)

e = np.linspace(1.9, 21, 400)
se = S(e)
ax.plot(se, 0.6455*se, "-", color="0.55", lw=1.4,
        label=r"Crist $0.6455\sqrt{\eta_0}$")
ax.plot(se, 0.67*se, "--", color="0.35", lw=1.4,
        label=r"Ashkenas–Sherman $0.67\sqrt{\eta_0}$")
ax.plot(se, gibb(e), ":", color="0.45", lw=1.8, label="Gibbings")
for z, s_ in enumerate(SER):
    ax.plot(S(np.asarray(s_["e"], float)), np.asarray(s_["v"], float),
            s_["m"], ms=s_["ms"], mfc=s_["fc"], mec=s_["ec"], mew=s_["mew"],
            alpha=1.0, ls="", zorder=5+z, label=s_["lab"])

# ---- regime nodes, values as quoted in the chapter's flow-regimes text
ax.axvspan(S(3.03), S(3.12), color="#f4d35e", alpha=0.50, zorder=1, lw=0)
ax.axvline(S(3.473), color="#b05f00", ls=(0, (4, 2)), lw=1.1, alpha=0.85, zorder=2)
ax.set_xlim(1.35, 4.65)
ax.set_ylim(0, 3.4)
ax.set_xlabel(r"$\sqrt{\eta_0}$", fontsize=15)
ax.set_ylabel(r"$x_{fc}/D$  or  $L_{MD}/D$", fontsize=15)
ax.grid(True, ls=":", lw=0.5, alpha=0.7)
ax.set_axisbelow(True)

sec = ax.secondary_xaxis("top", functions=(lambda s: s**2,
                                           lambda v: S(np.maximum(v, 0))))
ET = [2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0]
sec.set_xticks(ET)
sec.set_xticklabels(["%g\n(%.2f)" % (q, S(5.0*(q**(2.0/7.0)-1.0))) for q in ET])
sec.set_xlabel(r"$\eta_0$   ($M_j$)", fontsize=14, labelpad=7)
sec.tick_params(labelsize=12.0)

NODEKEY = [
    Patch(facecolor="#f4d35e", alpha=0.50,
          label=r"first Mach disk, $\eta_0=3.03$–$3.12$ (Muraoka–Hiejima)"),
    Line2D([], [], color="#b05f00", ls=(0, (4, 2)), lw=1.1,
           label=r"regular$\to$Mach reflection, $\eta_0=3.47$ (Gibbings)")]
# Two columns, grouped by what the entry is rather than by drawing order:
# the left column is everything taken from the literature - the three
# correlations and the two regime nodes - and the right column is everything
# computed here. Matplotlib fills a legend column by column, so five then five
# with ncol=2 puts the split exactly on that boundary.
h, l = ax.get_legend_handles_labels()
CORR = [(hh, ll) for hh, ll in zip(h, l)
        if any(k in ll for k in ("Crist", "Ashkenas", "Gibbings"))]
SIM = [(hh, ll) for hh, ll in zip(h, l)
       if not any(k in ll for k in ("Crist", "Ashkenas", "Gibbings"))]
LEFT = [x[0] for x in CORR] + NODEKEY
RIGHT = [x[0] for x in SIM]
LL = [x[1] for x in CORR] + [x.get_label() for x in NODEKEY]
RL = [x[1] for x in SIM]
# Five entries per column, so the two columns line up row for row. Spacing is
# opened up: matplotlib sizes each legend row across both columns, so the rows
# stay level once the counts match.
ax.legend(handles=LEFT+RIGHT, labels=LL+RL, loc="lower right", ncol=2,
          fontsize=9.5, framealpha=0.95, borderpad=0.95, labelspacing=0.78,
          handlelength=1.9, handletextpad=0.85, columnspacing=2.6,
          borderaxespad=0.6)
for ext in ("png", "pdf"):
    fig.savefig(OUT/("mach_disk_tophat_single.%s" % ext), dpi=300,
                bbox_inches="tight", pad_inches=0.03)
print("\nwrote mach_disk_tophat_single.png/pdf")
